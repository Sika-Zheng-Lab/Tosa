//! Integration tests for Tosa.
//!
//! These tests use a small synthetic example BAM (26 reads on chr1)
//! and a matching GTF annotation with 3 exons / 2 introns.
//! CRAM tests use the same data converted to CRAM format with a reference FASTA.
//!
//! Expected junctions (unstranded):
//!   chr1:1201-1499  .  13
//!   chr1:1701-1999  .   5
//!
//! Expected junctions (XS strand):
//!   chr1:1201-1499  +  10
//!   chr1:1201-1499  -   3
//!   chr1:1701-1999  +   5
//!
//! Expected boundaries (unstranded, with GTF):
//!   chr1:1199-1201  5p  .  2
//!   chr1:1498-1500  3p  .  2
//!   chr1:1699-1701  5p  .  1
//!   chr1:1998-2000  3p  .  1
//!
//! Single-cell barcodes: AAAA-1, BBBB-1, CCCC-1
//! Expected per-cell junction counts (unstranded, after UMI dedup):
//!   chr1:1201-1499  .  AAAA-1   6  (10 reads, 6 unique UMIs)
//!   chr1:1201-1499  .  BBBB-1   2  (3 reads, 2 unique UMIs)
//!   chr1:1701-1999  .  CCCC-1   3  (5 reads, 3 unique UMIs)

use std::collections::HashSet;
use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("examples");
    path
}

fn test_bam_path() -> String {
    test_data_dir()
        .join("example.bam")
        .to_str()
        .unwrap()
        .to_string()
}

fn test_cram_path() -> String {
    test_data_dir()
        .join("example.cram")
        .to_str()
        .unwrap()
        .to_string()
}

fn test_gtf_path() -> String {
    test_data_dir()
        .join("annotation.gtf")
        .to_str()
        .unwrap()
        .to_string()
}

fn test_barcodes_path() -> String {
    test_data_dir()
        .join("barcodes.tsv")
        .to_str()
        .unwrap()
        .to_string()
}

// ---------------------------------------------------------------------------
// Helper: build a RunConfig for bulk mode
// ---------------------------------------------------------------------------
fn bulk_config(strand: tosa::types::StrandMode, gtf: Option<String>) -> tosa::types::RunConfig {
    tosa::types::RunConfig {
        mode: tosa::types::Mode::Bulk,
        bam_file: test_bam_path(),
        output_prefix: "/dev/null".to_string(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: None,
        strand_mode: strand,
        gtf_file: gtf,
        verbose: false,
        threads: 1,
    }
}

// ===========================================================================
// 1. Total read count
// ===========================================================================
#[test]
fn test_count_total_reads() {
    let total = tosa::bam_reader::count_total_reads(&test_bam_path(), 1).unwrap();
    assert_eq!(total, 26, "Synthetic BAM has exactly 26 records");
}

// ===========================================================================
// 2. Unstranded junction counting – exact values
// ===========================================================================
#[test]
fn test_bulk_junction_unstranded() {
    let config = bulk_config(tosa::types::StrandMode::Unstranded, None);
    let result = tosa::bam_reader::process_bam_records(&config, &HashSet::new(), None).unwrap();

    assert_eq!(
        result.junction_totals.len(),
        2,
        "Should find exactly 2 junctions"
    );

    // The composite key includes the strand suffix, e.g. "chr1:1201-1499:."
    let j1 = result
        .junction_totals
        .iter()
        .find(|(k, _)| k.starts_with("chr1:1201-1499"))
        .map(|(_, &v)| v);
    let j2 = result
        .junction_totals
        .iter()
        .find(|(k, _)| k.starts_with("chr1:1701-1999"))
        .map(|(_, &v)| v);

    assert_eq!(j1, Some(13), "chr1:1201-1499 should have 13 reads");
    assert_eq!(j2, Some(5), "chr1:1701-1999 should have 5 reads");
}

// ===========================================================================
// 3. XS strand junction counting – strand-split values
// ===========================================================================
#[test]
fn test_bulk_junction_xs_strand() {
    let config = bulk_config(tosa::types::StrandMode::XS, None);
    let result = tosa::bam_reader::process_bam_records(&config, &HashSet::new(), None).unwrap();

    // XS splits the first junction into + and -, second only +
    assert_eq!(
        result.junction_totals.len(),
        3,
        "XS mode should yield 3 junction rows"
    );

    let get = |key: &str| result.junction_totals.get(key).copied().unwrap_or(0);

    assert_eq!(get("chr1:1201-1499:+"), 10, "chr1:1201-1499 + strand = 10");
    assert_eq!(get("chr1:1201-1499:-"), 3, "chr1:1201-1499 - strand = 3");
    assert_eq!(get("chr1:1701-1999:+"), 5, "chr1:1701-1999 + strand = 5");
}

// ===========================================================================
// 4. Boundary counting with GTF (unstranded)
// ===========================================================================
#[test]
fn test_bulk_boundary_with_gtf() {
    let config = bulk_config(tosa::types::StrandMode::Unstranded, Some(test_gtf_path()));
    let boundary_index = tosa::gtf::parse_gtf(&test_gtf_path(), 1).unwrap();
    let result =
        tosa::bam_reader::process_bam_records(&config, &HashSet::new(), Some(&boundary_index))
            .unwrap();

    // Junctions must still be present
    assert_eq!(result.junction_totals.len(), 2);

    // Boundary counts
    assert!(
        !result.boundary_totals.is_empty(),
        "Should have boundary counts"
    );

    let get_b = |key: &str| result.boundary_totals.get(key).copied().unwrap_or(0);

    assert_eq!(get_b("chr1:1199-1201"), 2, "5' boundary at intron 1 start");
    assert_eq!(get_b("chr1:1498-1500"), 2, "3' boundary at intron 1 end");
    assert_eq!(get_b("chr1:1699-1701"), 1, "5' boundary at intron 2 start");
    assert_eq!(get_b("chr1:1998-2000"), 1, "3' boundary at intron 2 end");
}

// ===========================================================================
// 5. GTF parsing – structure check
// ===========================================================================
#[test]
fn test_gtf_parsing() {
    let boundary_index = tosa::gtf::parse_gtf(&test_gtf_path(), 1).unwrap();

    // Should have boundaries on chr1
    assert!(
        boundary_index.boundaries.contains_key("chr1"),
        "GTF contains chr1 annotations"
    );

    // GTF: exons at 1000-1200, 1500-1700, 2000-2300 (1-based inclusive)
    // Intron 1: 1201-1499 → 5' boundary 1199-1201, 3' boundary 1498-1500
    // Intron 2: 1701-1999 → 5' boundary 1699-1701, 3' boundary 1998-2000
    let overlapping_5p = boundary_index.find_overlapping("chr1", 1199, 1202);
    assert!(
        !overlapping_5p.is_empty(),
        "Should find 5' boundary near exon1-intron1 junction"
    );

    let overlapping_3p = boundary_index.find_overlapping("chr1", 1497, 1500);
    assert!(
        !overlapping_3p.is_empty(),
        "Should find 3' boundary near intron1-exon2 junction"
    );
}

// ===========================================================================
// 6. Output round-trip (write + verify)
// ===========================================================================
#[test]
fn test_bulk_output_write() {
    use std::collections::HashMap;
    use tosa::types::Strand;

    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir
        .path()
        .join("test_output")
        .to_str()
        .unwrap()
        .to_string();

    let mut junction_totals = HashMap::new();
    junction_totals.insert("chr1:100-200:+".to_string(), 5u32);
    junction_totals.insert("chr1:300-400:-".to_string(), 3u32);

    let mut junction_strands = HashMap::new();
    junction_strands.insert("chr1:100-200:+".to_string(), Strand::Plus);
    junction_strands.insert("chr1:300-400:-".to_string(), Strand::Minus);

    tosa::output::write_junction_bulk(&prefix, &junction_totals, &junction_strands).unwrap();

    // Verify the output file was created
    let output_path = format!("{}_junction.tsv.gz", prefix);
    assert!(
        std::path::Path::new(&output_path).exists(),
        "Output file should exist"
    );

    // Read and verify content
    use flate2::read::GzDecoder;
    use std::io::Read;
    let file = std::fs::File::open(&output_path).unwrap();
    let mut decoder = GzDecoder::new(file);
    let mut content = String::new();
    decoder.read_to_string(&mut content).unwrap();

    assert!(
        content.contains("Junction\tStrand\tCount"),
        "Should have header"
    );
    assert!(
        content.contains("chr1:100-200\t+\t5"),
        "Should have junction entry"
    );
    assert!(
        content.contains("chr1:300-400\t-\t3"),
        "Should have junction entry"
    );
}

// ---------------------------------------------------------------------------
// Helper: build a RunConfig for single-cell mode
// ---------------------------------------------------------------------------
fn single_config(
    strand: tosa::types::StrandMode,
    gtf: Option<String>,
    barcode_file: Option<String>,
) -> tosa::types::RunConfig {
    tosa::types::RunConfig {
        mode: tosa::types::Mode::Single,
        bam_file: test_bam_path(),
        output_prefix: "/dev/null".to_string(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: barcode_file,
        strand_mode: strand,
        gtf_file: gtf,
        verbose: false,
        threads: 1,
    }
}

/// Load barcodes of interest from the example barcodes.tsv.
fn load_barcodes() -> HashSet<String> {
    let path = test_barcodes_path();
    tosa::data_loader::load_cell_barcodes(Some(&path)).unwrap()
}

// ===========================================================================
// 7. Single-cell junction counting (unstranded)
// ===========================================================================
#[test]
fn test_single_junction_unstranded() {
    let barcodes = load_barcodes();
    let config = single_config(
        tosa::types::StrandMode::Unstranded,
        None,
        Some(test_barcodes_path()),
    );
    let result = tosa::bam_reader::process_bam_records(&config, &barcodes, None).unwrap();

    // 3 barcodes observed
    assert_eq!(
        result.cell_barcodes.len(),
        3,
        "Should discover 3 cell barcodes"
    );
    assert!(result.cell_barcodes.contains("AAAA-1"));
    assert!(result.cell_barcodes.contains("BBBB-1"));
    assert!(result.cell_barcodes.contains("CCCC-1"));

    // 2 junctions
    assert_eq!(result.junction_counts.len(), 2, "Should find 2 junctions");

    // chr1:1201-1499  AAAA-1=6 (UMI dedup), BBBB-1=2 (UMI dedup)
    let j1_key = result
        .junction_counts
        .keys()
        .find(|k| k.starts_with("chr1:1201-1499"))
        .unwrap();
    let j1 = &result.junction_counts[j1_key];
    assert_eq!(j1.get("AAAA-1").copied().unwrap_or(0), 6);
    assert_eq!(j1.get("BBBB-1").copied().unwrap_or(0), 2);

    // chr1:1701-1999  CCCC-1=3 (UMI dedup)
    let j2_key = result
        .junction_counts
        .keys()
        .find(|k| k.starts_with("chr1:1701-1999"))
        .unwrap();
    let j2 = &result.junction_counts[j2_key];
    assert_eq!(j2.get("CCCC-1").copied().unwrap_or(0), 3);
}

// ===========================================================================
// 8. Single-cell junction counting (XS strand)
// ===========================================================================
#[test]
fn test_single_junction_xs_strand() {
    let barcodes = load_barcodes();
    let config = single_config(
        tosa::types::StrandMode::XS,
        None,
        Some(test_barcodes_path()),
    );
    let result = tosa::bam_reader::process_bam_records(&config, &barcodes, None).unwrap();

    // XS splits junction 1 into + and - → 3 junction-strand combos
    assert_eq!(
        result.junction_counts.len(),
        3,
        "XS mode should yield 3 junction rows"
    );

    assert_eq!(
        result
            .junction_counts
            .get("chr1:1201-1499:+")
            .and_then(|m| m.get("AAAA-1"))
            .copied()
            .unwrap_or(0),
        6
    );
    assert_eq!(
        result
            .junction_counts
            .get("chr1:1201-1499:-")
            .and_then(|m| m.get("BBBB-1"))
            .copied()
            .unwrap_or(0),
        2
    );
    assert_eq!(
        result
            .junction_counts
            .get("chr1:1701-1999:+")
            .and_then(|m| m.get("CCCC-1"))
            .copied()
            .unwrap_or(0),
        3
    );
}

// ===========================================================================
// 9. Single-cell boundary counting with GTF
// ===========================================================================
#[test]
fn test_single_boundary_with_gtf() {
    let barcodes = load_barcodes();
    let config = single_config(
        tosa::types::StrandMode::Unstranded,
        Some(test_gtf_path()),
        Some(test_barcodes_path()),
    );
    let boundary_index = tosa::gtf::parse_gtf(&test_gtf_path(), 1).unwrap();
    let result =
        tosa::bam_reader::process_bam_records(&config, &barcodes, Some(&boundary_index)).unwrap();

    // All boundary reads are on CCCC-1
    let get_bc = |key: &str| {
        result
            .boundary_counts
            .get(key)
            .and_then(|m| m.get("CCCC-1"))
            .copied()
            .unwrap_or(0)
    };

    assert_eq!(
        get_bc("chr1:1199-1201"),
        1,
        "5' boundary intron 1, CCCC-1 (UMI dedup)"
    );
    assert_eq!(
        get_bc("chr1:1498-1500"),
        1,
        "3' boundary intron 1, CCCC-1 (UMI dedup)"
    );
    assert_eq!(get_bc("chr1:1699-1701"), 1, "5' boundary intron 2, CCCC-1");
    assert_eq!(get_bc("chr1:1998-2000"), 1, "3' boundary intron 2, CCCC-1");
}

// ===========================================================================
// 10. Single-cell output round-trip (MatrixMarket + barcodes + features)
// ===========================================================================
#[test]
fn test_single_output_write() {
    use std::collections::HashMap;
    use tosa::types::Strand;

    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir.path().join("sc_out").to_str().unwrap().to_string();

    // Simulate 2 junctions × 2 barcodes
    let mut junction_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    junction_counts
        .entry("chr1:1201-1499:.".to_string())
        .or_default()
        .insert("AAAA-1".to_string(), 10);
    junction_counts
        .entry("chr1:1201-1499:.".to_string())
        .or_default()
        .insert("BBBB-1".to_string(), 3);
    junction_counts
        .entry("chr1:1701-1999:.".to_string())
        .or_default()
        .insert("CCCC-1".to_string(), 5);

    let mut cell_barcodes = HashSet::new();
    cell_barcodes.insert("AAAA-1".to_string());
    cell_barcodes.insert("BBBB-1".to_string());
    cell_barcodes.insert("CCCC-1".to_string());

    let mut junction_strands = HashMap::new();
    junction_strands.insert("chr1:1201-1499:.".to_string(), Strand::Unknown);
    junction_strands.insert("chr1:1701-1999:.".to_string(), Strand::Unknown);

    tosa::output::write_junction_single(
        &prefix,
        &junction_counts,
        &cell_barcodes,
        &junction_strands,
    )
    .unwrap();

    use flate2::read::GzDecoder;
    use std::io::Read;

    // Verify matrix.mtx.gz
    let mtx_path = format!("{}_matrix.mtx.gz", prefix);
    assert!(std::path::Path::new(&mtx_path).exists());
    let mut content = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(&mtx_path).unwrap())
        .read_to_string(&mut content)
        .unwrap();
    assert!(content.contains("%%MatrixMarket"));
    // 2 features × 3 barcodes × 3 non-zero entries
    assert!(
        content.contains("2 3 3"),
        "Matrix dimensions should be 2×3 with 3 entries"
    );

    // Verify barcodes.tsv.gz
    let bc_path = format!("{}_barcodes.tsv.gz", prefix);
    let mut bc_content = String::new();
    GzDecoder::new(std::fs::File::open(&bc_path).unwrap())
        .read_to_string(&mut bc_content)
        .unwrap();
    assert!(bc_content.contains("AAAA-1"));
    assert!(bc_content.contains("BBBB-1"));
    assert!(bc_content.contains("CCCC-1"));

    // Verify features.tsv.gz
    let feat_path = format!("{}_features.tsv.gz", prefix);
    let mut feat_content = String::new();
    GzDecoder::new(std::fs::File::open(&feat_path).unwrap())
        .read_to_string(&mut feat_content)
        .unwrap();
    assert!(feat_content.contains("chr1:1201-1499"));
    assert!(feat_content.contains("chr1:1701-1999"));

    // Verify junction_barcodes.tsv.gz
    let jb_path = format!("{}_junction_barcodes.tsv.gz", prefix);
    let mut jb_content = String::new();
    GzDecoder::new(std::fs::File::open(&jb_path).unwrap())
        .read_to_string(&mut jb_content)
        .unwrap();
    assert!(jb_content.contains("Feature\tStrand\tBarcode\tCount"));
    assert!(jb_content.contains("AAAA-1\t10"));
    assert!(jb_content.contains("BBBB-1\t3"));
    assert!(jb_content.contains("CCCC-1\t5"));
}

// ===========================================================================
// 11. Bulk boundary output round-trip (write + verify)
// ===========================================================================
#[test]
fn test_bulk_boundary_output_write() {
    use std::collections::HashMap;
    use tosa::types::{BoundaryType, Strand};

    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir
        .path()
        .join("test_boundary")
        .to_str()
        .unwrap()
        .to_string();

    let mut boundary_totals = HashMap::new();
    boundary_totals.insert("chr1:1200-1201".to_string(), 5u32);
    boundary_totals.insert("chr1:1498-1499".to_string(), 3u32);

    let mut boundary_types = HashMap::new();
    boundary_types.insert("chr1:1200-1201".to_string(), BoundaryType::FivePrime);
    boundary_types.insert("chr1:1498-1499".to_string(), BoundaryType::ThreePrime);

    let mut boundary_strands = HashMap::new();
    boundary_strands.insert("chr1:1200-1201".to_string(), Strand::Plus);
    boundary_strands.insert("chr1:1498-1499".to_string(), Strand::Plus);

    tosa::output::write_boundary_bulk(
        &prefix,
        &boundary_totals,
        &boundary_types,
        &boundary_strands,
    )
    .unwrap();

    let output_path = format!("{}_boundary.tsv.gz", prefix);
    assert!(
        std::path::Path::new(&output_path).exists(),
        "Boundary output file should exist"
    );

    use flate2::read::GzDecoder;
    use std::io::Read;
    let file = std::fs::File::open(&output_path).unwrap();
    let mut decoder = GzDecoder::new(file);
    let mut content = String::new();
    decoder.read_to_string(&mut content).unwrap();

    assert!(
        content.contains("Boundary\tType\tStrand\tCount"),
        "Should have header"
    );
    assert!(
        content.contains("chr1:1200-1201\t5p\t+\t5"),
        "Should have 5' boundary entry"
    );
    assert!(
        content.contains("chr1:1498-1499\t3p\t+\t3"),
        "Should have 3' boundary entry"
    );
}

// ===========================================================================
// 12. Single-cell boundary output round-trip
// ===========================================================================
#[test]
fn test_single_boundary_output_write() {
    use std::collections::HashMap;
    use tosa::types::{BoundaryType, Strand};

    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir
        .path()
        .join("sc_boundary")
        .to_str()
        .unwrap()
        .to_string();

    let mut boundary_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    boundary_counts
        .entry("chr1:1200-1201".to_string())
        .or_default()
        .insert("AAAA-1".to_string(), 4);
    boundary_counts
        .entry("chr1:1498-1499".to_string())
        .or_default()
        .insert("BBBB-1".to_string(), 2);

    let mut cell_barcodes = HashSet::new();
    cell_barcodes.insert("AAAA-1".to_string());
    cell_barcodes.insert("BBBB-1".to_string());

    let mut boundary_types = HashMap::new();
    boundary_types.insert("chr1:1200-1201".to_string(), BoundaryType::FivePrime);
    boundary_types.insert("chr1:1498-1499".to_string(), BoundaryType::ThreePrime);

    let mut boundary_strands = HashMap::new();
    boundary_strands.insert("chr1:1200-1201".to_string(), Strand::Plus);
    boundary_strands.insert("chr1:1498-1499".to_string(), Strand::Plus);

    tosa::output::write_boundary_single(
        &prefix,
        &boundary_counts,
        &cell_barcodes,
        &boundary_types,
        &boundary_strands,
    )
    .unwrap();

    use flate2::read::GzDecoder;
    use std::io::Read;

    // Verify boundary_matrix.mtx.gz
    let mtx_path = format!("{}_boundary_matrix.mtx.gz", prefix);
    assert!(std::path::Path::new(&mtx_path).exists());
    let mut content = String::new();
    GzDecoder::new(std::fs::File::open(&mtx_path).unwrap())
        .read_to_string(&mut content)
        .unwrap();
    assert!(content.contains("%%MatrixMarket"));
    assert!(
        content.contains("2 2 2"),
        "Matrix dimensions should be 2×2 with 2 entries"
    );

    // Verify boundary_barcodes.tsv.gz
    let bc_path = format!("{}_boundary_barcodes.tsv.gz", prefix);
    let mut bc_content = String::new();
    GzDecoder::new(std::fs::File::open(&bc_path).unwrap())
        .read_to_string(&mut bc_content)
        .unwrap();
    assert!(bc_content.contains("AAAA-1"));
    assert!(bc_content.contains("BBBB-1"));

    // Verify boundary_features.tsv.gz
    let feat_path = format!("{}_boundary_features.tsv.gz", prefix);
    let mut feat_content = String::new();
    GzDecoder::new(std::fs::File::open(&feat_path).unwrap())
        .read_to_string(&mut feat_content)
        .unwrap();
    assert!(feat_content.contains("chr1:1200-1201\t5p\t+"));
    assert!(feat_content.contains("chr1:1498-1499\t3p\t+"));

    // Verify boundary_barcodes_detail.tsv.gz
    let detail_path = format!("{}_boundary_barcodes_detail.tsv.gz", prefix);
    let mut detail_content = String::new();
    GzDecoder::new(std::fs::File::open(&detail_path).unwrap())
        .read_to_string(&mut detail_content)
        .unwrap();
    assert!(detail_content.contains("Boundary\tType\tStrand\tBarcode\tCount"));
    assert!(detail_content.contains("AAAA-1\t4"));
    assert!(detail_content.contains("BBBB-1\t2"));
}

// ===========================================================================
// 13. Run pipeline – bulk mode (end-to-end via lib::run)
// ===========================================================================
#[test]
fn test_run_bulk_mode() {
    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir.path().join("run_bulk").to_str().unwrap().to_string();

    let config = tosa::types::RunConfig {
        mode: tosa::types::Mode::Bulk,
        bam_file: test_bam_path(),
        output_prefix: prefix.clone(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: None,
        strand_mode: tosa::types::StrandMode::Unstranded,
        gtf_file: Some(test_gtf_path()),
        verbose: false,
        threads: 1,
    };

    tosa::run(&config).unwrap();

    // Verify junction output
    let junction_path = format!("{}_junction.tsv.gz", prefix);
    assert!(std::path::Path::new(&junction_path).exists());

    // Verify boundary output
    let boundary_path = format!("{}_boundary.tsv.gz", prefix);
    assert!(std::path::Path::new(&boundary_path).exists());
}

// ===========================================================================
// 14. Run pipeline – bulk mode without GTF
// ===========================================================================
#[test]
fn test_run_bulk_no_gtf() {
    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir
        .path()
        .join("run_bulk_no_gtf")
        .to_str()
        .unwrap()
        .to_string();

    let config = tosa::types::RunConfig {
        mode: tosa::types::Mode::Bulk,
        bam_file: test_bam_path(),
        output_prefix: prefix.clone(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: None,
        strand_mode: tosa::types::StrandMode::Unstranded,
        gtf_file: None,
        verbose: false,
        threads: 1,
    };

    tosa::run(&config).unwrap();

    let junction_path = format!("{}_junction.tsv.gz", prefix);
    assert!(std::path::Path::new(&junction_path).exists());

    // No boundary output when GTF is not provided
    let boundary_path = format!("{}_boundary.tsv.gz", prefix);
    assert!(!std::path::Path::new(&boundary_path).exists());
}

// ===========================================================================
// 15. Run pipeline – single-cell mode
// ===========================================================================
#[test]
fn test_run_single_mode() {
    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir
        .path()
        .join("run_single")
        .to_str()
        .unwrap()
        .to_string();

    let config = tosa::types::RunConfig {
        mode: tosa::types::Mode::Single,
        bam_file: test_bam_path(),
        output_prefix: prefix.clone(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: Some(test_barcodes_path()),
        strand_mode: tosa::types::StrandMode::Unstranded,
        gtf_file: Some(test_gtf_path()),
        verbose: false,
        threads: 1,
    };

    tosa::run(&config).unwrap();

    // Verify junction output files
    assert!(std::path::Path::new(&format!("{}_matrix.mtx.gz", prefix)).exists());
    assert!(std::path::Path::new(&format!("{}_barcodes.tsv.gz", prefix)).exists());
    assert!(std::path::Path::new(&format!("{}_features.tsv.gz", prefix)).exists());

    // Verify boundary output files
    assert!(std::path::Path::new(&format!("{}_boundary_matrix.mtx.gz", prefix)).exists());
    assert!(std::path::Path::new(&format!("{}_boundary_barcodes.tsv.gz", prefix)).exists());
    assert!(std::path::Path::new(&format!("{}_boundary_features.tsv.gz", prefix)).exists());
}

// ===========================================================================
// 16. Run pipeline – single-cell mode without barcodes file
// ===========================================================================
#[test]
fn test_run_single_no_barcode_file() {
    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir
        .path()
        .join("run_single_nobc")
        .to_str()
        .unwrap()
        .to_string();

    let config = tosa::types::RunConfig {
        mode: tosa::types::Mode::Single,
        bam_file: test_bam_path(),
        output_prefix: prefix.clone(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: None,
        strand_mode: tosa::types::StrandMode::Unstranded,
        gtf_file: None,
        verbose: false,
        threads: 1,
    };

    tosa::run(&config).unwrap();

    assert!(std::path::Path::new(&format!("{}_matrix.mtx.gz", prefix)).exists());
}

// ===========================================================================
// CRAM tests: verify CRAM files produce identical results to BAM
// ===========================================================================

/// Check that CRAM example files exist (generated by generate_example_bam).
fn cram_files_exist() -> bool {
    std::path::Path::new(&test_cram_path()).exists()
}

// ---------------------------------------------------------------------------
// Helper: build a RunConfig for CRAM bulk mode
// ---------------------------------------------------------------------------
fn bulk_config_cram(
    strand: tosa::types::StrandMode,
    gtf: Option<String>,
) -> tosa::types::RunConfig {
    tosa::types::RunConfig {
        mode: tosa::types::Mode::Bulk,
        bam_file: test_cram_path(),
        output_prefix: "/dev/null".to_string(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: None,
        strand_mode: strand,
        gtf_file: gtf,
        verbose: false,
        threads: 1,
    }
}

// ===========================================================================
// 17. CRAM: Total read count
// ===========================================================================
#[test]
fn test_cram_count_total_reads() {
    if !cram_files_exist() {
        eprintln!("Skipping CRAM test: example.cram not found. Run `cargo run --example generate_example_bam` first.");
        return;
    }
    let total = tosa::bam_reader::count_total_reads(&test_cram_path(), 1).unwrap();
    assert_eq!(total, 26, "CRAM should have same 26 records as BAM");
}

// ===========================================================================
// 18. CRAM: Unstranded junction counting
// ===========================================================================
#[test]
fn test_cram_bulk_junction_unstranded() {
    if !cram_files_exist() {
        return;
    }
    let config = bulk_config_cram(tosa::types::StrandMode::Unstranded, None);
    let result = tosa::bam_reader::process_bam_records(&config, &HashSet::new(), None).unwrap();

    assert_eq!(
        result.junction_totals.len(),
        2,
        "CRAM: Should find exactly 2 junctions"
    );

    let j1 = result
        .junction_totals
        .iter()
        .find(|(k, _)| k.starts_with("chr1:1201-1499"))
        .map(|(_, &v)| v);
    let j2 = result
        .junction_totals
        .iter()
        .find(|(k, _)| k.starts_with("chr1:1701-1999"))
        .map(|(_, &v)| v);

    assert_eq!(j1, Some(13), "CRAM: chr1:1201-1499 should have 13 reads");
    assert_eq!(j2, Some(5), "CRAM: chr1:1701-1999 should have 5 reads");
}

// ===========================================================================
// 19. CRAM: XS strand junction counting
// ===========================================================================
#[test]
fn test_cram_bulk_junction_xs_strand() {
    if !cram_files_exist() {
        return;
    }
    let config = bulk_config_cram(tosa::types::StrandMode::XS, None);
    let result = tosa::bam_reader::process_bam_records(&config, &HashSet::new(), None).unwrap();

    assert_eq!(
        result.junction_totals.len(),
        3,
        "CRAM XS mode should yield 3 junction rows"
    );

    let get = |key: &str| result.junction_totals.get(key).copied().unwrap_or(0);
    assert_eq!(get("chr1:1201-1499:+"), 10);
    assert_eq!(get("chr1:1201-1499:-"), 3);
    assert_eq!(get("chr1:1701-1999:+"), 5);
}

// ===========================================================================
// 20. CRAM: Boundary counting with GTF
// ===========================================================================
#[test]
fn test_cram_bulk_boundary_with_gtf() {
    if !cram_files_exist() {
        return;
    }
    let config = bulk_config_cram(tosa::types::StrandMode::Unstranded, Some(test_gtf_path()));
    let boundary_index = tosa::gtf::parse_gtf(&test_gtf_path(), 1).unwrap();
    let result =
        tosa::bam_reader::process_bam_records(&config, &HashSet::new(), Some(&boundary_index))
            .unwrap();

    assert_eq!(result.junction_totals.len(), 2);
    assert!(!result.boundary_totals.is_empty());

    let get_b = |key: &str| result.boundary_totals.get(key).copied().unwrap_or(0);
    assert_eq!(get_b("chr1:1199-1201"), 2, "CRAM: 5' boundary at intron 1");
    assert_eq!(get_b("chr1:1498-1500"), 2, "CRAM: 3' boundary at intron 1");
    assert_eq!(get_b("chr1:1699-1701"), 1, "CRAM: 5' boundary at intron 2");
    assert_eq!(get_b("chr1:1998-2000"), 1, "CRAM: 3' boundary at intron 2");
}

// ===========================================================================
// 21. CRAM: Run pipeline – bulk mode (end-to-end)
// ===========================================================================
#[test]
fn test_cram_run_bulk_mode() {
    if !cram_files_exist() {
        return;
    }
    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir
        .path()
        .join("cram_bulk")
        .to_str()
        .unwrap()
        .to_string();

    let config = tosa::types::RunConfig {
        mode: tosa::types::Mode::Bulk,
        bam_file: test_cram_path(),
        output_prefix: prefix.clone(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: None,
        strand_mode: tosa::types::StrandMode::Unstranded,
        gtf_file: Some(test_gtf_path()),
        verbose: false,
        threads: 1,
    };

    tosa::run(&config).unwrap();

    assert!(std::path::Path::new(&format!("{}_junction.tsv.gz", prefix)).exists());
    assert!(std::path::Path::new(&format!("{}_boundary.tsv.gz", prefix)).exists());
}

// ===========================================================================
// 22. CRAM: Single-cell mode
// ===========================================================================
#[test]
fn test_cram_run_single_mode() {
    if !cram_files_exist() {
        return;
    }
    let tmpdir = tempfile::tempdir().unwrap();
    let prefix = tmpdir
        .path()
        .join("cram_single")
        .to_str()
        .unwrap()
        .to_string();

    let config = tosa::types::RunConfig {
        mode: tosa::types::Mode::Single,
        bam_file: test_cram_path(),
        output_prefix: prefix.clone(),
        min_anchor_length: 8,
        min_boundary_anchor_length: 8,
        min_intron_length: 20,
        max_intron_length: 500000,
        max_loci: 1,
        cell_barcode_file: Some(test_barcodes_path()),
        strand_mode: tosa::types::StrandMode::Unstranded,
        gtf_file: Some(test_gtf_path()),
        verbose: false,
        threads: 1,
    };

    tosa::run(&config).unwrap();

    assert!(std::path::Path::new(&format!("{}_matrix.mtx.gz", prefix)).exists());
    assert!(std::path::Path::new(&format!("{}_barcodes.tsv.gz", prefix)).exists());
    assert!(std::path::Path::new(&format!("{}_features.tsv.gz", prefix)).exists());
    assert!(std::path::Path::new(&format!("{}_boundary_matrix.mtx.gz", prefix)).exists());
}
