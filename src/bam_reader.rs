//! BAM/CRAM file processing: read iteration, CIGAR parsing, junction extraction, and boundary counting.

use log::{debug, info};
use rayon::prelude::*;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::{self, Read};
use rust_htslib::htslib;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::boundary::{self, BoundaryIndex};
use crate::junction;
use crate::types::{hash_read_name, JunctionKey, Mode, RunConfig, Strand, StrandMode};

/// Results from processing a BAM file.
pub struct ProcessingResult {
    /// Per-cell junction counts: junction_key -> (barcode -> count). Used in single mode.
    pub junction_counts: HashMap<String, HashMap<String, u32>>,
    /// Total junction counts: junction_key -> count. Used in bulk mode.
    pub junction_totals: HashMap<String, u32>,
    /// Strand assigned to each junction key.
    pub junction_strands: HashMap<String, Strand>,
    /// All observed cell barcodes.
    pub cell_barcodes: HashSet<String>,
    /// Per-cell boundary counts: boundary_key -> (barcode -> count). Used in single mode.
    pub boundary_counts: HashMap<String, HashMap<String, u32>>,
    /// Total boundary counts: boundary_key -> count. Used in bulk mode.
    pub boundary_totals: HashMap<String, u32>,
    /// Boundary type for each boundary key.
    pub boundary_types: HashMap<String, crate::types::BoundaryType>,
    /// Strand assigned to each boundary key.
    pub boundary_strands: HashMap<String, Strand>,
}

/// SAM fields Tosa needs from each record (everything except SEQ/QUAL).
/// By declaring these, CRAM can be decoded without a reference FASTA.
const TOSA_REQUIRED_FIELDS: u32 = htslib::sam_fields_SAM_QNAME
    | htslib::sam_fields_SAM_FLAG
    | htslib::sam_fields_SAM_RNAME
    | htslib::sam_fields_SAM_POS
    | htslib::sam_fields_SAM_CIGAR
    | htslib::sam_fields_SAM_AUX;

/// Count total mapped reads using the BAM/CRAM index.
pub fn count_total_reads(
    bam_file: &str,
    threads: usize,
) -> Result<u64, Box<dyn std::error::Error + Send + Sync>> {
    let mut bam_index_reader = IndexedReader::from_path(bam_file)?;
    if threads > 1 {
        bam_index_reader.set_threads(threads - 1)?;
    }
    let stats = bam_index_reader.index_stats()?;
    debug!("stats: {:?}", stats);
    let total_mapped_reads: u64 = stats.iter().map(|(_, _, mapped, _)| mapped).sum();
    Ok(total_mapped_reads)
}

/// Determine strand of a read based on the strand mode.
pub fn determine_strand(record: &bam::Record, strand_mode: &StrandMode) -> Strand {
    match strand_mode {
        StrandMode::Unstranded => Strand::Unknown,
        StrandMode::XS => match record.aux(b"XS") {
            Ok(Aux::Char(c)) => match c {
                b'+' => Strand::Plus,
                b'-' => Strand::Minus,
                _ => Strand::Unknown,
            },
            _ => Strand::Unknown,
        },
        StrandMode::RF => {
            // First-strand: read1+reverse => +, read1+forward => -
            //               read2+reverse => -, read2+forward => +
            let is_reverse = record.is_reverse();
            if record.is_first_in_template() {
                if is_reverse {
                    Strand::Plus
                } else {
                    Strand::Minus
                }
            } else if record.is_last_in_template() {
                if is_reverse {
                    Strand::Minus
                } else {
                    Strand::Plus
                }
            } else {
                // Single-end read in RF mode: reverse => +, forward => -
                if is_reverse {
                    Strand::Plus
                } else {
                    Strand::Minus
                }
            }
        }
        StrandMode::FR => {
            // Second-strand: read1+forward => +, read1+reverse => -
            //                read2+forward => -, read2+reverse => +
            let is_reverse = record.is_reverse();
            if record.is_first_in_template() {
                if is_reverse {
                    Strand::Minus
                } else {
                    Strand::Plus
                }
            } else if record.is_last_in_template() {
                if is_reverse {
                    Strand::Plus
                } else {
                    Strand::Minus
                }
            } else {
                // Single-end read in FR mode: forward => +, reverse => -
                if is_reverse {
                    Strand::Minus
                } else {
                    Strand::Plus
                }
            }
        }
    }
}

/// Extract aligned segments from a CIGAR string.
/// Returns a list of (start, end) reference coordinate intervals for M/=/X operations.
pub fn extract_aligned_segments(record: &bam::Record) -> Vec<(i64, i64)> {
    let mut segments = Vec::new();
    let mut current_pos = record.pos();

    let cigar_view = record.cigar();
    for op in cigar_view.iter() {
        match op {
            Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
                let len = *l as i64;
                segments.push((current_pos, current_pos + len));
                current_pos += len;
            }
            Cigar::Del(l) | Cigar::RefSkip(l) => {
                current_pos += *l as i64;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {
                // These do not consume reference bases
            }
        }
    }
    segments
}

/// Per-chromosome processing result (used internally for merging).
struct ChromResult {
    junction_counts: HashMap<JunctionKey, HashMap<String, u32>>,
    junction_totals: HashMap<JunctionKey, u32>,
    junction_strands: HashMap<JunctionKey, Strand>,
    junction_has_left_anchor: HashMap<JunctionKey, bool>,
    junction_has_right_anchor: HashMap<JunctionKey, bool>,
    cell_barcodes: HashSet<String>,
    boundary_counts: HashMap<String, HashMap<String, u32>>,
    boundary_totals: HashMap<String, u32>,
    boundary_types: HashMap<String, crate::types::BoundaryType>,
    boundary_strands: HashMap<String, Strand>,
}

/// Process one chromosome's records from an IndexedReader.
#[allow(clippy::too_many_arguments)]
fn process_chromosome(
    bam_file: &str,
    tid: u32,
    chrom: &str,
    config: &RunConfig,
    cell_barcodes_of_interest: &HashSet<String>,
    boundary_index: Option<&BoundaryIndex>,
    reference_names: &[String],
    progress_counter: &AtomicU64,
    total_mapped_reads: u64,
) -> Result<ChromResult, Box<dyn std::error::Error + Send + Sync>> {
    let mut reader = IndexedReader::from_path(bam_file)?;
    reader.set_cram_options(
        htslib::hts_fmt_option_CRAM_OPT_REQUIRED_FIELDS,
        TOSA_REQUIRED_FIELDS,
    )?;
    // Each thread has its own reader; no need for per-reader htslib IO threads
    reader.fetch(rust_htslib::bam::FetchDefinition::RegionString(
        chrom.as_bytes(),
        0,
        i64::MAX,
    ))?;

    let mut junction_state = junction::JunctionState::new();
    let mut junction_strands: HashMap<JunctionKey, Strand> = HashMap::new();
    let mut junction_has_left_anchor: HashMap<JunctionKey, bool> = HashMap::new();
    let mut junction_has_right_anchor: HashMap<JunctionKey, bool> = HashMap::new();
    let mut cell_barcodes: HashSet<String> = HashSet::new();

    let mut boundary_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut boundary_totals: HashMap<String, u32> = HashMap::new();
    let mut boundary_types: HashMap<String, crate::types::BoundaryType> = HashMap::new();
    let mut boundary_strands: HashMap<String, Strand> = HashMap::new();
    let mut processed_boundary_reads: HashMap<String, HashSet<u64>> = HashMap::new();
    let mut processed_boundary_umis: HashMap<String, HashSet<u64>> = HashMap::new();

    let mut local_read_count: u64 = 0;
    let is_single = config.mode == Mode::Single;

    for result in reader.records() {
        let record = result?;
        local_read_count += 1;

        // Progress logging (atomic counter shared across threads)
        if total_mapped_reads > 0 && local_read_count.is_multiple_of(10000) {
            let global_count = progress_counter.fetch_add(10000, Ordering::Relaxed) + 10000;
            let progress_percentage = (global_count * 100) / total_mapped_reads;
            if progress_percentage <= 100 {
                info!(
                    "Progress: {}% ({} / {})",
                    progress_percentage, global_count, total_mapped_reads
                );
            }
        }

        // Skip read if NH tag exceeds max_loci
        if let Ok(Aux::U8(nh)) = record.aux(b"NH") {
            if nh > config.max_loci as u8 {
                continue;
            }
        } else if let Ok(Aux::I32(nh)) = record.aux(b"NH") {
            if nh > config.max_loci as i32 {
                continue;
            }
        }

        // Extract reference name (chromosome) and start position
        let rec_tid = record.tid();
        if rec_tid < 0 {
            continue; // Unmapped read
        }

        // Borrow reference name — no per-read heap allocation
        let ref_name = &reference_names[rec_tid as usize];
        let mut current_pos = record.pos();

        // Extract Cell Barcode (CB) from tags if in single mode
        let cell_barcode = if is_single {
            match record.aux(b"CB") {
                Ok(Aux::String(cb_str)) => Some(cb_str.to_string()),
                _ => None,
            }
        } else {
            None
        };

        // Extract UMI (UB tag) if in single mode
        let umi = if is_single {
            match record.aux(b"UB") {
                Ok(Aux::String(ub_str)) => Some(ub_str.to_string()),
                _ => None,
            }
        } else {
            None
        };

        // Skip read if its barcode is not in the list of interest
        if let Some(cb) = &cell_barcode {
            if config.cell_barcode_file.is_some()
                && !cell_barcodes_of_interest.is_empty()
                && !cell_barcodes_of_interest.contains(cb)
            {
                continue;
            }
        }

        // If a cell barcode is present (for single mode), or always process for bulk mode
        if !is_single || cell_barcode.is_some() {
            if let Some(cb_str) = &cell_barcode {
                cell_barcodes.insert(cb_str.clone());
            }

            // Determine strand for this read
            let strand = determine_strand(&record, &config.strand_mode);
            let read_name_hash = hash_read_name(record.qname());

            // --- Junction extraction from CIGAR ---
            let cigar_view = record.cigar();
            let cigars: Vec<_> = cigar_view.iter().collect();
            for i in 0..cigars.len() {
                if let Cigar::RefSkip(len) = cigars[i] {
                    let intron_length = *len as i64;
                    if intron_length < config.min_intron_length
                        || intron_length > config.max_intron_length
                    {
                        current_pos += intron_length;
                        continue;
                    }

                    // Calculate left anchor length
                    let mut left_anchor_length: i64 = 0;
                    let mut j = i;
                    while j > 0 {
                        j -= 1;
                        match cigars[j] {
                            Cigar::Match(l) | Cigar::Equal(l) => {
                                left_anchor_length += *l as i64;
                                if left_anchor_length >= config.min_anchor_length {
                                    break;
                                }
                            }
                            Cigar::RefSkip(_) => continue,
                            _ => break,
                        }
                    }
                    let has_left_anchor = left_anchor_length >= config.min_anchor_length;

                    // Calculate right anchor length
                    let mut right_anchor_length: i64 = 0;
                    let mut k = i + 1;
                    while k < cigars.len() {
                        match cigars[k] {
                            Cigar::Match(r) | Cigar::Equal(r) => {
                                right_anchor_length += *r as i64;
                                if right_anchor_length >= config.min_anchor_length {
                                    break;
                                }
                            }
                            Cigar::RefSkip(_) => {
                                k += 1;
                                continue;
                            }
                            _ => break,
                        }
                        k += 1;
                    }
                    let has_right_anchor = right_anchor_length >= config.min_anchor_length;

                    // Use 1-based inclusive coordinates for junction IDs
                    let start = current_pos + 1; // 1-based intron start
                    let end = current_pos + intron_length; // 1-based intron end
                    let jkey = JunctionKey {
                        tid,
                        start,
                        end,
                        strand,
                    };

                    // Per-junction anchor tracking (regtools-style):
                    // OR-accumulate left/right anchor flags across all reads.
                    {
                        let left_flag = junction_has_left_anchor.entry(jkey).or_insert(false);
                        *left_flag = *left_flag || has_left_anchor;
                        let right_flag = junction_has_right_anchor.entry(jkey).or_insert(false);
                        *right_flag = *right_flag || has_right_anchor;
                    }

                    // Always count the read; filtering happens at output time.
                    junction_strands.entry(jkey).or_insert(strand);
                    junction::process_junction(
                        jkey,
                        cell_barcode.as_ref(),
                        umi.as_ref(),
                        &mut junction_state,
                        read_name_hash,
                        config.mode,
                    );
                    current_pos += intron_length;
                } else if let Cigar::SoftClip(_) = cigars[i] {
                    // SoftClip does not consume reference bases
                    continue;
                } else {
                    // FIX: Ins does NOT consume reference bases
                    current_pos += match cigars[i] {
                        Cigar::Match(l) | Cigar::Del(l) => *l as i64,
                        _ => 0,
                    };
                }
            }

            // --- Boundary counting ---
            if let Some(bi) = boundary_index {
                let segments = extract_aligned_segments(&record);
                boundary::count_boundaries(
                    ref_name,
                    &segments,
                    bi,
                    cell_barcode.as_ref(),
                    umi.as_ref(),
                    strand,
                    &mut boundary_counts,
                    &mut boundary_totals,
                    &mut boundary_types,
                    &mut boundary_strands,
                    &mut processed_boundary_reads,
                    &mut processed_boundary_umis,
                    read_name_hash,
                    config.mode,
                );
            }
        }
    }

    // Account for remaining reads not yet reported to progress counter
    let remainder = local_read_count % 10000;
    if remainder > 0 {
        progress_counter.fetch_add(remainder, Ordering::Relaxed);
    }

    Ok(ChromResult {
        junction_counts: junction_state.junction_counts,
        junction_totals: junction_state.junction_totals,
        junction_strands,
        junction_has_left_anchor,
        junction_has_right_anchor,
        cell_barcodes,
        boundary_counts,
        boundary_totals,
        boundary_types,
        boundary_strands,
    })
}

/// Process all records in a BAM file, extracting junction and boundary counts.
///
/// When `config.threads > 1`, processing is parallelized per-chromosome using rayon.
pub fn process_bam_records(
    config: &RunConfig,
    cell_barcodes_of_interest: &HashSet<String>,
    boundary_index: Option<&BoundaryIndex>,
) -> Result<ProcessingResult, Box<dyn std::error::Error + Send + Sync>> {
    let total_mapped_reads = count_total_reads(&config.bam_file, config.threads)?;
    info!("Total number of reads: {}", total_mapped_reads);

    // Get reference names (chromosome names) from BAM/CRAM header
    let bam_reader = bam::Reader::from_path(&config.bam_file)?;
    let header = bam_reader.header().to_owned();
    let reference_names: Vec<String> = header
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect();
    drop(bam_reader);

    // Build list of (tid, chrom_name) for all chromosomes
    let chroms: Vec<(u32, String)> = reference_names
        .iter()
        .enumerate()
        .map(|(i, name)| (i as u32, name.clone()))
        .collect();

    // Shared atomic progress counter
    let progress_counter = AtomicU64::new(0);

    // Process chromosomes in parallel using rayon
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build()?;

    let chrom_results: Vec<Result<ChromResult, Box<dyn std::error::Error + Send + Sync>>> = pool
        .install(|| {
            chroms
                .par_iter()
                .map(|(tid, chrom)| {
                    process_chromosome(
                        &config.bam_file,
                        *tid,
                        chrom,
                        config,
                        cell_barcodes_of_interest,
                        boundary_index,
                        &reference_names,
                        &progress_counter,
                        total_mapped_reads,
                    )
                })
                .collect()
        });

    // Merge per-chromosome results
    let mut junction_counts: HashMap<JunctionKey, HashMap<String, u32>> = HashMap::new();
    let mut junction_totals: HashMap<JunctionKey, u32> = HashMap::new();
    let mut junction_strands: HashMap<JunctionKey, Strand> = HashMap::new();
    let mut junction_has_left_anchor: HashMap<JunctionKey, bool> = HashMap::new();
    let mut junction_has_right_anchor: HashMap<JunctionKey, bool> = HashMap::new();
    let mut cell_barcodes: HashSet<String> = HashSet::new();
    let mut boundary_counts: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut boundary_totals: HashMap<String, u32> = HashMap::new();
    let mut boundary_types: HashMap<String, crate::types::BoundaryType> = HashMap::new();
    let mut boundary_strands: HashMap<String, Strand> = HashMap::new();

    for chrom_result in chrom_results {
        let cr = chrom_result?;

        // Merge junction data — keys are unique per chromosome so extend() is safe
        junction_counts.extend(cr.junction_counts);
        junction_totals.extend(cr.junction_totals);
        junction_strands.extend(cr.junction_strands);

        // Merge anchor flags with OR logic
        for (k, v) in cr.junction_has_left_anchor {
            let flag = junction_has_left_anchor.entry(k).or_insert(false);
            *flag = *flag || v;
        }
        for (k, v) in cr.junction_has_right_anchor {
            let flag = junction_has_right_anchor.entry(k).or_insert(false);
            *flag = *flag || v;
        }

        cell_barcodes.extend(cr.cell_barcodes);
        boundary_counts.extend(cr.boundary_counts);
        boundary_totals.extend(cr.boundary_totals);
        boundary_types.extend(cr.boundary_types);
        boundary_strands.extend(cr.boundary_strands);
    }

    info!(
        "Progress: 100% ({} / {})",
        total_mapped_reads, total_mapped_reads
    );

    // Filter junctions: only emit those where at least one read provided a
    // sufficient left anchor AND at least one read provided a sufficient right
    // anchor (regtools-style per-junction filtering).
    let anchor_pass = |k: &JunctionKey| -> bool {
        *junction_has_left_anchor.get(k).unwrap_or(&false)
            && *junction_has_right_anchor.get(k).unwrap_or(&false)
    };

    // Convert JunctionKey-keyed maps to String-keyed maps for output compatibility
    let junction_totals_out: HashMap<String, u32> = junction_totals
        .into_iter()
        .filter(|(k, _)| anchor_pass(k))
        .map(|(k, v)| (k.to_string_key(&reference_names), v))
        .collect();

    let junction_strands_out: HashMap<String, Strand> = junction_strands
        .into_iter()
        .filter(|(k, _)| anchor_pass(k))
        .map(|(k, v)| (k.to_string_key(&reference_names), v))
        .collect();

    let junction_counts_out: HashMap<String, HashMap<String, u32>> = junction_counts
        .into_iter()
        .filter(|(k, _)| anchor_pass(k))
        .map(|(k, v)| (k.to_string_key(&reference_names), v))
        .collect();

    Ok(ProcessingResult {
        junction_counts: junction_counts_out,
        junction_totals: junction_totals_out,
        junction_strands: junction_strands_out,
        cell_barcodes,
        boundary_counts,
        boundary_totals,
        boundary_types,
        boundary_strands,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Aux;

    /// Helper: build a minimal BAM record with the requested SAM flags.
    fn record_with_flags(flags: u16) -> bam::Record {
        let mut rec = bam::Record::new();
        // set() needs: name, cigar, seq, qualcargo
        // Provide a tiny dummy alignment so the record is valid.
        rec.set(
            b"dummy",
            Some(&bam::record::CigarString(vec![Cigar::Match(4)])),
            b"ACGT",
            &[30, 30, 30, 30],
        );
        rec.set_flags(flags);
        rec
    }

    // ---------------------------------------------------------------
    // Unstranded
    // ---------------------------------------------------------------
    #[test]
    fn test_strand_unstranded() {
        let rec = record_with_flags(0);
        assert_eq!(
            determine_strand(&rec, &StrandMode::Unstranded),
            Strand::Unknown
        );
    }

    // ---------------------------------------------------------------
    // XS tag
    // ---------------------------------------------------------------
    #[test]
    fn test_strand_xs_plus() {
        let mut rec = record_with_flags(0);
        rec.push_aux(b"XS", Aux::Char(b'+')).unwrap();
        assert_eq!(determine_strand(&rec, &StrandMode::XS), Strand::Plus);
    }

    #[test]
    fn test_strand_xs_minus() {
        let mut rec = record_with_flags(0);
        rec.push_aux(b"XS", Aux::Char(b'-')).unwrap();
        assert_eq!(determine_strand(&rec, &StrandMode::XS), Strand::Minus);
    }

    #[test]
    fn test_strand_xs_unknown_char() {
        let mut rec = record_with_flags(0);
        rec.push_aux(b"XS", Aux::Char(b'?')).unwrap();
        assert_eq!(determine_strand(&rec, &StrandMode::XS), Strand::Unknown);
    }

    #[test]
    fn test_strand_xs_no_tag() {
        let rec = record_with_flags(0);
        assert_eq!(determine_strand(&rec, &StrandMode::XS), Strand::Unknown);
    }

    // ---------------------------------------------------------------
    // RF mode – paired-end
    // ---------------------------------------------------------------
    // Flags: 0x1 = paired, 0x10 = reverse, 0x40 = read1, 0x80 = read2
    #[test]
    fn test_strand_rf_read1_reverse() {
        // read1 + reverse => Plus
        let rec = record_with_flags(0x1 | 0x40 | 0x10);
        assert_eq!(determine_strand(&rec, &StrandMode::RF), Strand::Plus);
    }

    #[test]
    fn test_strand_rf_read1_forward() {
        // read1 + forward => Minus
        let rec = record_with_flags(0x1 | 0x40);
        assert_eq!(determine_strand(&rec, &StrandMode::RF), Strand::Minus);
    }

    #[test]
    fn test_strand_rf_read2_reverse() {
        // read2 + reverse => Minus
        let rec = record_with_flags(0x1 | 0x80 | 0x10);
        assert_eq!(determine_strand(&rec, &StrandMode::RF), Strand::Minus);
    }

    #[test]
    fn test_strand_rf_read2_forward() {
        // read2 + forward => Plus
        let rec = record_with_flags(0x1 | 0x80);
        assert_eq!(determine_strand(&rec, &StrandMode::RF), Strand::Plus);
    }

    // RF – single-end (neither read1 nor read2 flag)
    #[test]
    fn test_strand_rf_single_end_reverse() {
        let rec = record_with_flags(0x10);
        assert_eq!(determine_strand(&rec, &StrandMode::RF), Strand::Plus);
    }

    #[test]
    fn test_strand_rf_single_end_forward() {
        let rec = record_with_flags(0);
        assert_eq!(determine_strand(&rec, &StrandMode::RF), Strand::Minus);
    }

    // ---------------------------------------------------------------
    // FR mode – paired-end
    // ---------------------------------------------------------------
    #[test]
    fn test_strand_fr_read1_forward() {
        // read1 + forward => Plus
        let rec = record_with_flags(0x1 | 0x40);
        assert_eq!(determine_strand(&rec, &StrandMode::FR), Strand::Plus);
    }

    #[test]
    fn test_strand_fr_read1_reverse() {
        // read1 + reverse => Minus
        let rec = record_with_flags(0x1 | 0x40 | 0x10);
        assert_eq!(determine_strand(&rec, &StrandMode::FR), Strand::Minus);
    }

    #[test]
    fn test_strand_fr_read2_forward() {
        // read2 + forward => Minus
        let rec = record_with_flags(0x1 | 0x80);
        assert_eq!(determine_strand(&rec, &StrandMode::FR), Strand::Minus);
    }

    #[test]
    fn test_strand_fr_read2_reverse() {
        // read2 + reverse => Plus
        let rec = record_with_flags(0x1 | 0x80 | 0x10);
        assert_eq!(determine_strand(&rec, &StrandMode::FR), Strand::Plus);
    }

    // FR – single-end
    #[test]
    fn test_strand_fr_single_end_forward() {
        let rec = record_with_flags(0);
        assert_eq!(determine_strand(&rec, &StrandMode::FR), Strand::Plus);
    }

    #[test]
    fn test_strand_fr_single_end_reverse() {
        let rec = record_with_flags(0x10);
        assert_eq!(determine_strand(&rec, &StrandMode::FR), Strand::Minus);
    }

    // ---------------------------------------------------------------
    // extract_aligned_segments
    // ---------------------------------------------------------------
    #[test]
    fn test_extract_segments_simple_match() {
        let mut rec = bam::Record::new();
        let seq = vec![b'A'; 50];
        let qual = vec![30u8; 50];
        rec.set(
            b"r1",
            Some(&bam::record::CigarString(vec![Cigar::Match(50)])),
            &seq,
            &qual,
        );
        rec.set_pos(100);
        let segs = extract_aligned_segments(&rec);
        assert_eq!(segs, vec![(100, 150)]);
    }

    #[test]
    fn test_extract_segments_with_intron() {
        // 10M 100N 10M  => two segments with a 100bp gap
        let seq = b"ACGTACGTACGTACGTACGT"; // 20bp
        let mut rec = bam::Record::new();
        rec.set(
            b"r2",
            Some(&bam::record::CigarString(vec![
                Cigar::Match(10),
                Cigar::RefSkip(100),
                Cigar::Match(10),
            ])),
            seq,
            &vec![30u8; 20],
        );
        rec.set_pos(1000);
        let segs = extract_aligned_segments(&rec);
        assert_eq!(segs, vec![(1000, 1010), (1110, 1120)]);
    }

    #[test]
    fn test_extract_segments_with_insertion_deletion() {
        // 5M 3I 5M 2D 5M
        let seq_len: usize = 5 + 3 + 5 + 5;
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        let mut rec = bam::Record::new();
        rec.set(
            b"r3",
            Some(&bam::record::CigarString(vec![
                Cigar::Match(5),
                Cigar::Ins(3),
                Cigar::Match(5),
                Cigar::Del(2),
                Cigar::Match(5),
            ])),
            &seq,
            &qual,
        );
        rec.set_pos(200);
        let segs = extract_aligned_segments(&rec);
        // 5M: [200, 205), 5M: [205, 210), 2D skips to 212, 5M: [212, 217)
        assert_eq!(segs, vec![(200, 205), (205, 210), (212, 217)]);
    }

    #[test]
    fn test_extract_segments_softclip_hardclip_pad() {
        // 3S 10M 2H — soft clip, hard clip, pad don't consume reference
        let seq_len: usize = 3 + 10; // softclip consumes query
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        let mut rec = bam::Record::new();
        rec.set(
            b"r4",
            Some(&bam::record::CigarString(vec![
                Cigar::SoftClip(3),
                Cigar::Match(10),
                Cigar::HardClip(2),
            ])),
            &seq,
            &qual,
        );
        rec.set_pos(500);
        let segs = extract_aligned_segments(&rec);
        assert_eq!(segs, vec![(500, 510)]);
    }

    #[test]
    fn test_extract_segments_equal_and_diff() {
        // 5= 3X — both consume reference like Match
        let seq_len: usize = 8;
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        let mut rec = bam::Record::new();
        rec.set(
            b"r5",
            Some(&bam::record::CigarString(vec![
                Cigar::Equal(5),
                Cigar::Diff(3),
            ])),
            &seq,
            &qual,
        );
        rec.set_pos(0);
        let segs = extract_aligned_segments(&rec);
        assert_eq!(segs, vec![(0, 5), (5, 8)]);
    }

    // ---------------------------------------------------------------
    // Helper: create a BAM record with specific fields
    // ---------------------------------------------------------------
    fn make_record(
        name: &[u8],
        flags: u16,
        tid: i32,
        pos: i64,
        cigar: Vec<Cigar>,
        seq_len: usize,
    ) -> bam::Record {
        let mut rec = bam::Record::new();
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        rec.set(name, Some(&bam::record::CigarString(cigar)), &seq, &qual);
        rec.set_flags(flags);
        rec.set_tid(tid);
        rec.set_pos(pos);
        rec.set_mapq(60);
        rec
    }

    /// Write a coordinate-sorted BAM file with given records and build its index.
    fn write_indexed_bam(path: &str, records: &[bam::Record]) {
        use rust_htslib::bam::{header::HeaderRecord, Format, Header, Writer};

        let mut header = Header::new();
        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "coordinate");
        header.push_record(&hd);

        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", 10_000_000);
        header.push_record(&sq);

        let mut writer = Writer::from_path(path, &header, Format::Bam).unwrap();
        for rec in records {
            writer.write(rec).unwrap();
        }
        drop(writer);

        bam::index::build(path, None, bam::index::Type::Bai, 0).unwrap();
    }

    // ---------------------------------------------------------------
    // process_bam_records: bulk mode – edge-case records
    // ---------------------------------------------------------------
    #[test]
    fn test_process_bam_bulk_edge_cases() {
        let _ = env_logger::Builder::from_default_env()
            .filter_level(log::LevelFilter::Info)
            .try_init();

        let tmpdir = tempfile::tempdir().unwrap();
        let bam_path = tmpdir.path().join("edge.bam");
        let bam_str = bam_path.to_str().unwrap();

        // Build records (sorted by tid, then pos; unmapped at end)
        let mut records: Vec<bam::Record> = Vec::new();

        // 1) Record with NH:i:5 as I32 → skipped (max_loci = 1)
        //    Covers L169-171
        let mut r_nh = make_record(b"nh_skip", 0, 0, 100, vec![Cigar::Match(10)], 10);
        r_nh.push_aux(b"NH", Aux::I32(5)).unwrap();
        records.push(r_nh);

        // 2) Record with NH:i:1 as I32 → kept (exercises L169-170 without skip)
        let mut r_nh_keep = make_record(b"nh_keep", 0, 0, 110, vec![Cigar::Match(10)], 10);
        r_nh_keep.push_aux(b"NH", Aux::I32(1)).unwrap();
        records.push(r_nh_keep);

        // 3) Short intron (50 < 70 = min_intron_length) → covers L220-221
        records.push(make_record(
            b"short_intron",
            0,
            0,
            400,
            vec![Cigar::Match(10), Cigar::RefSkip(50), Cigar::Match(10)],
            20,
        ));

        // 4) Multi-junction with RefSkip in LEFT anchor calculation
        //    CIGAR: 10M 100N 3M 200N 10M
        //    Processing the 200N: left anchor goes backwards through 3M (3<8),
        //    then encounters RefSkip(100N) → continue (L236), then 10M → 13>=8
        records.push(make_record(
            b"refskip_left",
            0,
            0,
            500,
            vec![
                Cigar::Match(10),
                Cigar::RefSkip(100),
                Cigar::Match(3),
                Cigar::RefSkip(200),
                Cigar::Match(10),
            ],
            23,
        ));

        // 5) Ins in LEFT anchor → _ => break (L237)
        //    CIGAR: 10M 2I 3M 200N 10M
        //    Processing 200N: left anchor → 3M (3<8), Ins(2) → break
        records.push(make_record(
            b"ins_left",
            0,
            0,
            600,
            vec![
                Cigar::Match(10),
                Cigar::Ins(2),
                Cigar::Match(3),
                Cigar::RefSkip(200),
                Cigar::Match(10),
            ],
            25,
        ));

        // 6) Multi-junction with RefSkip in RIGHT anchor calculation
        //    CIGAR: 10M 200N 3M 100N 10M
        //    Processing 200N: right anchor → 3M (3<8), RefSkip(100N) → continue (L253), 10M→13>=8
        records.push(make_record(
            b"refskip_right",
            0,
            0,
            800,
            vec![
                Cigar::Match(10),
                Cigar::RefSkip(200),
                Cigar::Match(3),
                Cigar::RefSkip(100),
                Cigar::Match(10),
            ],
            23,
        ));

        // 7) Ins in RIGHT anchor → _ => break (L254)
        //    CIGAR: 10M 200N 3M 2I 10M
        //    Processing 200N: right anchor → 3M (3<8), Ins(2) → break
        records.push(make_record(
            b"ins_right",
            0,
            0,
            1000,
            vec![
                Cigar::Match(10),
                Cigar::RefSkip(200),
                Cigar::Match(3),
                Cigar::Ins(2),
                Cigar::Match(10),
            ],
            25,
        ));

        // 8) SoftClip + Ins in CIGAR loop → covers L309 (SoftClip continue), L314 (_ => 0 for Ins)
        records.push(make_record(
            b"softclip_ins",
            0,
            0,
            2000,
            vec![
                Cigar::SoftClip(3),
                Cigar::Match(10),
                Cigar::Ins(2),
                Cigar::Match(10),
            ],
            25,
        ));

        // 9-10) Per-junction anchor tracking (regtools-style)
        //   bad_anchor (pos=5000): 10M 100N 3M → left=10>=8, right=3<8
        //   good_anchor (pos=5000): 10M 100N 10M → left=10>=8, right=10>=8
        //   Both share junction chr1:5011-5110
        //   Per-junction: left OK (both reads), right OK (good_anchor) → junction reported
        //   Both reads are counted.
        records.push(make_record(
            b"bad_anchor",
            0,
            0,
            5000,
            vec![Cigar::Match(10), Cigar::RefSkip(100), Cigar::Match(3)],
            13,
        ));
        records.push(make_record(
            b"good_anchor",
            0,
            0,
            5000,
            vec![Cigar::Match(10), Cigar::RefSkip(100), Cigar::Match(10)],
            20,
        ));

        // 11) Unmapped read (tid=-1, flag=0x4) → L178
        //     Goes at the end of sorted BAM
        {
            let mut rec = bam::Record::new();
            rec.set(b"unmapped", None, b"", &[]);
            rec.set_flags(0x4);
            rec.set_tid(-1);
            rec.set_pos(-1);
            records.push(rec);
        }

        write_indexed_bam(bam_str, &records);

        // Configure: bulk mode, min_anchor=8, min_intron=20, max_intron=500000
        let config = crate::types::RunConfig {
            mode: crate::types::Mode::Bulk,
            bam_file: bam_str.to_string(),
            output_prefix: tmpdir.path().join("out").to_str().unwrap().to_string(),
            min_anchor_length: 8,
            min_boundary_anchor_length: 8,
            min_intron_length: 20,
            max_intron_length: 500000,
            max_loci: 1,
            cell_barcode_file: None,
            strand_mode: StrandMode::Unstranded,
            gtf_file: None,
            verbose: false,
            threads: 1,
        };

        let result = process_bam_records(&config, &HashSet::new(), None).unwrap();

        // The good_anchor + bad_anchor pair → junction chr1:5011-5110 should be found
        // Junction keys include strand suffix ":." for Unstranded mode
        assert!(
            result.junction_totals.contains_key("chr1:5011-5110:."),
            "Expected junction chr1:5011-5110:., got: {:?}",
            result.junction_totals.keys().collect::<Vec<_>>()
        );
        // Per-junction anchor: bad_anchor provides left anchor, good_anchor provides both.
        // Both reads are always counted (filtering is at junction level, not read level).
        assert_eq!(*result.junction_totals.get("chr1:5011-5110:.").unwrap(), 2);

        // refskip_left produces junction at chr1:511-610 and chr1:614-813
        // refskip_right produces: pos=800, Match(10)→810, RefSkip(200): junction at chr1:811-1010
        // right anchor: 3+skip+10=13>=8 (RefSkip spans into next exon)
        assert!(
            result.junction_totals.contains_key("chr1:811-1010:."),
            "Expected junction chr1:811-1010:. from refskip_right"
        );
    }

    // ---------------------------------------------------------------
    // Per-junction anchor: different reads provide left / right
    // ---------------------------------------------------------------
    #[test]
    fn test_per_junction_anchor_from_different_reads() {
        let tmpdir = tempfile::tempdir().unwrap();
        let bam_path = tmpdir.path().join("perjunc.bam");
        let bam_str = bam_path.to_str().unwrap();

        let mut records: Vec<bam::Record> = Vec::new();

        // Read A: left=10>=8, right=3<8  → provides left anchor only
        records.push(make_record(
            b"read_a",
            0,
            0,
            6000,
            vec![Cigar::Match(10), Cigar::RefSkip(100), Cigar::Match(3)],
            13,
        ));

        // Read B: left=3<8, right=10>=8  → provides right anchor only
        // Same junction: chr1:6011-6110
        // left anchor: 3M only (no preceding RefSkip to span)
        records.push(make_record(
            b"read_b",
            0,
            0,
            6007,
            vec![Cigar::Match(3), Cigar::RefSkip(100), Cigar::Match(10)],
            13,
        ));

        write_indexed_bam(bam_str, &records);

        let config = crate::types::RunConfig {
            mode: crate::types::Mode::Bulk,
            bam_file: bam_str.to_string(),
            output_prefix: tmpdir.path().join("out").to_str().unwrap().to_string(),
            min_anchor_length: 8,
            min_boundary_anchor_length: 8,
            min_intron_length: 20,
            max_intron_length: 500000,
            max_loci: 1,
            cell_barcode_file: None,
            strand_mode: StrandMode::Unstranded,
            gtf_file: None,
            verbose: false,
            threads: 1,
        };

        let result = process_bam_records(&config, &HashSet::new(), None).unwrap();

        // Per-junction: Read A provides left>=8, Read B provides right>=8 → reported
        assert!(
            result.junction_totals.contains_key("chr1:6011-6110:."),
            "Expected per-junction anchor to report chr1:6011-6110:., got: {:?}",
            result.junction_totals.keys().collect::<Vec<_>>()
        );
        // Both reads counted
        assert_eq!(*result.junction_totals.get("chr1:6011-6110:.").unwrap(), 2);
    }

    // ---------------------------------------------------------------
    // Cigar::Diff (X = mismatch) breaks anchor accumulation
    // ---------------------------------------------------------------
    #[test]
    fn test_diff_breaks_anchor() {
        let tmpdir = tempfile::tempdir().unwrap();
        let bam_path = tmpdir.path().join("diff.bam");
        let bam_str = bam_path.to_str().unwrap();

        let mut records: Vec<bam::Record> = Vec::new();

        // CIGAR: 10M 1X 3M 200N 10M
        // Left anchor for 200N: backwards → 3M (3<8), Diff(1X) → break
        // Left anchor = 3 < 8 → left anchor fails
        // Right anchor = 10 >= 8 → right anchor passes
        // With only one read, the junction lacks a left anchor → NOT reported
        records.push(make_record(
            b"diff_left",
            0,
            0,
            7000,
            vec![
                Cigar::Match(10),
                Cigar::Diff(1),
                Cigar::Match(3),
                Cigar::RefSkip(200),
                Cigar::Match(10),
            ],
            24,
        ));

        // CIGAR: 10M 200N 3M 1X 10M
        // Right anchor for 200N: forwards → 3M (3<8), Diff(1X) → break
        // Right anchor = 3 < 8 → right anchor fails
        // Left anchor = 10 >= 8 → left anchor passes
        // With only one read, the junction lacks a right anchor → NOT reported
        records.push(make_record(
            b"diff_right",
            0,
            0,
            7500,
            vec![
                Cigar::Match(10),
                Cigar::RefSkip(200),
                Cigar::Match(3),
                Cigar::Diff(1),
                Cigar::Match(10),
            ],
            24,
        ));

        write_indexed_bam(bam_str, &records);

        let config = crate::types::RunConfig {
            mode: crate::types::Mode::Bulk,
            bam_file: bam_str.to_string(),
            output_prefix: tmpdir.path().join("out").to_str().unwrap().to_string(),
            min_anchor_length: 8,
            min_boundary_anchor_length: 8,
            min_intron_length: 20,
            max_intron_length: 500000,
            max_loci: 1,
            cell_barcode_file: None,
            strand_mode: StrandMode::Unstranded,
            gtf_file: None,
            verbose: false,
            threads: 1,
        };

        let result = process_bam_records(&config, &HashSet::new(), None).unwrap();

        // diff_left: junction chr1:7014-7213 — only right anchor → not reported (no left)
        assert!(
            !result.junction_totals.contains_key("chr1:7014-7213:."),
            "Diff should break left anchor; junction chr1:7014-7213 should NOT be reported"
        );

        // diff_right: junction chr1:7511-7710 — only left anchor → not reported (no right)
        assert!(
            !result.junction_totals.contains_key("chr1:7511-7710:."),
            "Diff should break right anchor; junction chr1:7511-7710 should NOT be reported"
        );
    }

    // ---------------------------------------------------------------
    // process_bam_records: single mode – CB tag and barcode filtering
    // ---------------------------------------------------------------
    #[test]
    fn test_process_bam_single_mode_barcode_filtering() {
        let tmpdir = tempfile::tempdir().unwrap();
        let bam_path = tmpdir.path().join("single.bam");
        let bam_str = bam_path.to_str().unwrap();

        let mut records: Vec<bam::Record> = Vec::new();

        // 1) Read WITHOUT CB tag in single mode → L187 (_ => None)
        //    Will be skipped (mode="single", cell_barcode=None → skip)
        records.push(make_record(
            b"no_cb",
            0,
            0,
            100,
            vec![Cigar::Match(10), Cigar::RefSkip(100), Cigar::Match(10)],
            20,
        ));

        // 2) Read with CB:Z:UNKNOWN-1 (not in interest) → L199 (continue)
        let mut r_unknown = make_record(
            b"unknown_bc",
            0,
            0,
            200,
            vec![Cigar::Match(10), Cigar::RefSkip(100), Cigar::Match(10)],
            20,
        );
        r_unknown.push_aux(b"CB", Aux::String("UNKNOWN-1")).unwrap();
        records.push(r_unknown);

        // 3) Read with CB:Z:KNOWN-1 (in interest) → processed normally
        let mut r_known = make_record(
            b"known_bc",
            0,
            0,
            300,
            vec![Cigar::Match(10), Cigar::RefSkip(100), Cigar::Match(10)],
            20,
        );
        r_known.push_aux(b"CB", Aux::String("KNOWN-1")).unwrap();
        records.push(r_known);

        write_indexed_bam(bam_str, &records);

        let mut barcodes_of_interest = HashSet::new();
        barcodes_of_interest.insert("KNOWN-1".to_string());

        let config = crate::types::RunConfig {
            mode: crate::types::Mode::Single,
            bam_file: bam_str.to_string(),
            output_prefix: tmpdir.path().join("out").to_str().unwrap().to_string(),
            min_anchor_length: 8,
            min_boundary_anchor_length: 8,
            min_intron_length: 20,
            max_intron_length: 500000,
            max_loci: 1,
            cell_barcode_file: Some("dummy_path".to_string()),
            strand_mode: StrandMode::Unstranded,
            gtf_file: None,
            verbose: false,
            threads: 1,
        };

        let result = process_bam_records(&config, &barcodes_of_interest, None).unwrap();

        // Only KNOWN-1 should appear in cell_barcodes
        assert_eq!(result.cell_barcodes.len(), 1);
        assert!(result.cell_barcodes.contains("KNOWN-1"));

        // Only one junction from the known_bc read
        assert_eq!(result.junction_totals.len(), 0); // junction goes through junction_counts not junction_totals in single mode
                                                     // In single mode, junctions are in junction_counts
        assert!(
            !result.junction_counts.is_empty(),
            "Expected junction from known_bc read"
        );
    }
}
