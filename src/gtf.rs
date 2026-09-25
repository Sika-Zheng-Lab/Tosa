//! GTF annotation file parser.
//!
//! Parses a GTF file, extracts exon records, groups them by transcript,
//! and derives intron coordinates (and thus exon-intron boundary coordinates).

use log::info;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::boundary::{intron_to_boundaries, BoundaryIndex};
use crate::types::Strand;

/// An exon record from a GTF file.
#[derive(Debug, Clone)]
struct ExonRecord {
    chrom: String,
    start: i64, // 0-based start
    end: i64,   // 0-based exclusive end
    strand: Strand,
}

/// Parse a GTF file and build a BoundaryIndex.
///
/// The `boundary_anchor_length` parameter controls how many bases on each side
/// of a splice site the boundary interval spans (default 1). A read's aligned
/// segment must fully contain this interval to be counted as a boundary read.
///
/// Algorithm:
/// 1. Read all "exon" records from the GTF
/// 2. Group exons by transcript_id
/// 3. For each transcript, sort exons by start position
/// 4. Derive introns from gaps between adjacent exons
/// 5. For each intron, create 5' and 3' boundary entries
/// 6. Deduplicate boundaries across all transcripts
pub fn parse_gtf(
    path: &str,
    boundary_anchor_length: i64,
) -> Result<BoundaryIndex, Box<dyn std::error::Error>> {
    info!("Parsing GTF file: {}", path);

    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Collect exons grouped by transcript_id
    let mut transcript_exons: HashMap<String, Vec<ExonRecord>> = HashMap::new();
    let mut line_count = 0;
    let mut exon_count = 0;

    for line in reader.lines() {
        let line = line?;
        line_count += 1;

        // Skip comment lines
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        // Only process exon features
        if fields[2] != "exon" {
            continue;
        }

        let chrom = fields[0].to_string();
        // GTF is 1-based inclusive, convert to 0-based half-open
        let start: i64 = fields[3].parse::<i64>()? - 1;
        let end: i64 = fields[4].parse::<i64>()?;
        let strand = match fields[6] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => Strand::Unknown,
        };

        // Extract transcript_id from attributes (field 8)
        let transcript_id = extract_attribute(fields[8], "transcript_id");
        if let Some(tid) = transcript_id {
            transcript_exons.entry(tid).or_default().push(ExonRecord {
                chrom,
                start,
                end,
                strand,
            });
            exon_count += 1;
        }
    }

    info!(
        "Parsed {} lines, {} exon records, {} transcripts",
        line_count,
        exon_count,
        transcript_exons.len()
    );

    // Derive introns and build boundary index
    let mut boundary_index = BoundaryIndex::new();
    let mut seen_boundaries: HashSet<String> = HashSet::new();
    let mut intron_count = 0;

    for (_tid, mut exons) in transcript_exons {
        // Sort exons by start position
        exons.sort_by_key(|e| e.start);

        // Derive introns from gaps between adjacent exons
        for i in 0..exons.len().saturating_sub(1) {
            let intron_start = exons[i].end;
            let intron_end = exons[i + 1].start;

            // Skip if exons overlap or are adjacent (no intron)
            if intron_start >= intron_end {
                continue;
            }

            let chrom = &exons[i].chrom;
            let strand = exons[i].strand;

            let (five_p, three_p) = intron_to_boundaries(
                chrom,
                intron_start,
                intron_end,
                strand,
                boundary_anchor_length,
            );

            // Deduplicate: same boundary can arise from multiple transcripts
            if seen_boundaries.insert(five_p.boundary_id.clone()) {
                boundary_index.add(chrom, five_p);
            }
            if seen_boundaries.insert(three_p.boundary_id.clone()) {
                boundary_index.add(chrom, three_p);
            }
            intron_count += 1;
        }
    }

    info!(
        "Derived {} introns, {} unique boundary coordinates",
        intron_count,
        seen_boundaries.len()
    );

    Ok(boundary_index)
}

/// Extract an attribute value from a GTF attributes string.
///
/// GTF attributes format: `key "value"; key "value"; ...`
fn extract_attribute(attributes: &str, key: &str) -> Option<String> {
    let search = format!("{} \"", key);
    if let Some(start_idx) = attributes.find(&search) {
        let value_start = start_idx + search.len();
        if let Some(end_idx) = attributes[value_start..].find('"') {
            return Some(attributes[value_start..value_start + end_idx].to_string());
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::BoundaryType;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_extract_attribute() {
        let attrs = r#"gene_id "GENE1"; transcript_id "TX1"; exon_number "1";"#;
        assert_eq!(
            extract_attribute(attrs, "transcript_id"),
            Some("TX1".to_string())
        );
        assert_eq!(
            extract_attribute(attrs, "gene_id"),
            Some("GENE1".to_string())
        );
        assert_eq!(extract_attribute(attrs, "missing_key"), None);
    }

    #[test]
    fn test_parse_gtf_simple() {
        let _ = env_logger::Builder::from_default_env()
            .filter_level(log::LevelFilter::Info)
            .try_init();
        // Create a minimal GTF file with one gene, one transcript, three exons
        let gtf_content = "\
chr1\tensembl\tgene\t1001\t5000\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";
chr1\tensembl\ttranscript\t1001\t5000\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";
chr1\tensembl\texon\t1001\t1200\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; exon_number \"1\";
chr1\tensembl\texon\t2001\t2200\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; exon_number \"2\";
chr1\tensembl\texon\t3001\t3500\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\"; exon_number \"3\";
";
        let mut tmpfile = NamedTempFile::new().unwrap();
        write!(tmpfile, "{}", gtf_content).unwrap();

        let boundary_index = parse_gtf(tmpfile.path().to_str().unwrap(), 1).unwrap();

        // Exon 1: [1000, 1200), Exon 2: [2000, 2200), Exon 3: [3000, 3500)
        // Intron 1: 1200-2000 → 5' boundary: 1199-1201, 3' boundary: 1999-2001
        // Intron 2: 2200-3000 → 5' boundary: 2199-2201, 3' boundary: 2999-3001

        // Check 5' boundary of intron 1
        let overlapping = boundary_index.find_overlapping("chr1", 1100, 1300);
        assert_eq!(overlapping.len(), 1);
        assert_eq!(overlapping[0].boundary_id, "chr1:1199-1201");
        assert_eq!(overlapping[0].boundary_type, BoundaryType::FivePrime);

        // Check 3' boundary of intron 1
        let overlapping = boundary_index.find_overlapping("chr1", 1990, 2010);
        assert_eq!(overlapping.len(), 1);
        assert_eq!(overlapping[0].boundary_id, "chr1:1999-2001");
        assert_eq!(overlapping[0].boundary_type, BoundaryType::ThreePrime);

        // Check 5' boundary of intron 2
        let overlapping = boundary_index.find_overlapping("chr1", 2100, 2300);
        assert_eq!(overlapping.len(), 1);
        assert_eq!(overlapping[0].boundary_id, "chr1:2199-2201");

        // Check 3' boundary of intron 2
        let overlapping = boundary_index.find_overlapping("chr1", 2990, 3010);
        assert_eq!(overlapping.len(), 1);
        assert_eq!(overlapping[0].boundary_id, "chr1:2999-3001");
    }

    #[test]
    fn test_parse_gtf_dedup_across_transcripts() {
        // Two transcripts sharing the same intron
        let gtf_content = "\
chr1\tensembl\texon\t1001\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
chr1\tensembl\texon\t2001\t2200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
chr1\tensembl\texon\t1001\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX2\";
chr1\tensembl\texon\t2001\t2500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX2\";
";
        let mut tmpfile = NamedTempFile::new().unwrap();
        write!(tmpfile, "{}", gtf_content).unwrap();

        let boundary_index = parse_gtf(tmpfile.path().to_str().unwrap(), 1).unwrap();

        // Both transcripts have intron 1200-2000, boundaries should be deduplicated
        let overlapping = boundary_index.find_overlapping("chr1", 1100, 1300);
        assert_eq!(overlapping.len(), 1); // Not 2
    }

    #[test]
    fn test_parse_gtf_comment_lines() {
        let gtf_content = "\
#!genome-build GRCh38
#!genome-version GRCh38
chr1\tensembl\texon\t1001\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
chr1\tensembl\texon\t2001\t2200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
";
        let mut tmpfile = NamedTempFile::new().unwrap();
        write!(tmpfile, "{}", gtf_content).unwrap();

        let boundary_index = parse_gtf(tmpfile.path().to_str().unwrap(), 1).unwrap();
        let overlapping = boundary_index.find_overlapping("chr1", 1100, 1300);
        assert_eq!(overlapping.len(), 1);
    }

    #[test]
    fn test_parse_gtf_short_lines() {
        // Lines with fewer than 9 tab-separated fields should be skipped
        let gtf_content = "\
short_line\twithout\tenough\tfields
chr1\tensembl\texon\t1001\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
chr1\tensembl\texon\t2001\t2200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
";
        let mut tmpfile = NamedTempFile::new().unwrap();
        write!(tmpfile, "{}", gtf_content).unwrap();

        let boundary_index = parse_gtf(tmpfile.path().to_str().unwrap(), 1).unwrap();
        let overlapping = boundary_index.find_overlapping("chr1", 1100, 1300);
        assert_eq!(overlapping.len(), 1);
    }

    #[test]
    fn test_parse_gtf_unknown_strand() {
        // Exons with strand "." should be parsed as Strand::Unknown
        let gtf_content = "\
chr1\tensembl\texon\t1001\t1200\t.\t.\t.\tgene_id \"G1\"; transcript_id \"TX1\";
chr1\tensembl\texon\t2001\t2200\t.\t.\t.\tgene_id \"G1\"; transcript_id \"TX1\";
";
        let mut tmpfile = NamedTempFile::new().unwrap();
        write!(tmpfile, "{}", gtf_content).unwrap();

        let boundary_index = parse_gtf(tmpfile.path().to_str().unwrap(), 1).unwrap();
        let overlapping = boundary_index.find_overlapping("chr1", 1098, 1302);
        assert_eq!(overlapping.len(), 1);
        // Verify the boundary has Unknown strand
        assert_eq!(overlapping[0].strand, Strand::Unknown);
    }

    #[test]
    fn test_parse_gtf_overlapping_exons() {
        // Adjacent/overlapping exons where intron_start >= intron_end should produce no intron
        let gtf_content = "\
chr1\tensembl\texon\t1001\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
chr1\tensembl\texon\t1200\t1400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
";
        let mut tmpfile = NamedTempFile::new().unwrap();
        write!(tmpfile, "{}", gtf_content).unwrap();

        let boundary_index = parse_gtf(tmpfile.path().to_str().unwrap(), 1).unwrap();
        assert!(
            !boundary_index.boundaries.contains_key("chr1")
                || boundary_index.find_overlapping("chr1", 0, 10000).is_empty()
        );
    }

    #[test]
    fn test_parse_gtf_single_exon_transcript() {
        // A transcript with only one exon should produce no introns/boundaries
        let gtf_content = "\
chr1\tensembl\texon\t1001\t1200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\";
";
        let mut tmpfile = NamedTempFile::new().unwrap();
        write!(tmpfile, "{}", gtf_content).unwrap();

        let boundary_index = parse_gtf(tmpfile.path().to_str().unwrap(), 1).unwrap();
        assert!(
            !boundary_index.boundaries.contains_key("chr1")
                || boundary_index.find_overlapping("chr1", 0, 10000).is_empty()
        );
    }
}
