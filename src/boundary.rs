//! Exon-intron boundary read counting.
//!
//! Boundary coordinates straddle the exon-intron splice site so that a read must
//! cover bases on **both** the exonic and intronic sides to be counted. The width
//! of each boundary interval is `2 × anchor_length` (anchor bases on each side).
//!
//! With the default anchor length of 1, intron `chr2:6545675-6547042` (0-based
//! half-open) produces:
//! - 5' boundary: `chr2:6545674-6545676`  (1 exon base + 1 intron base)
//! - 3' boundary: `chr2:6547041-6547043`  (1 intron base + 1 exon base)

use crate::types::{BoundaryType, Mode, Strand};
use std::collections::{BTreeMap, HashMap, HashSet};

/// A single boundary entry derived from GTF annotation.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BoundaryEntry {
    /// Boundary coordinate string, e.g. "chr2:6545674-6545676"
    pub boundary_id: String,
    /// Start position of the boundary interval.
    pub start: i64,
    /// End position of the boundary interval (half-open).
    pub end: i64,
    /// Type of boundary (5' or 3').
    pub boundary_type: BoundaryType,
    /// Strand from GTF annotation.
    pub strand: Strand,
}

/// Index of all exon-intron boundaries, organized by chromosome.
///
/// For each chromosome, boundaries are stored in a BTreeMap keyed by start position
/// for efficient range lookups.
#[derive(Debug, Clone, Default)]
pub struct BoundaryIndex {
    /// chrom -> (start_pos -> Vec<BoundaryEntry>)
    pub boundaries: HashMap<String, BTreeMap<i64, Vec<BoundaryEntry>>>,
}

impl BoundaryIndex {
    /// Create an empty BoundaryIndex.
    pub fn new() -> Self {
        BoundaryIndex {
            boundaries: HashMap::new(),
        }
    }

    /// Add a boundary entry to the index.
    pub fn add(&mut self, chrom: &str, entry: BoundaryEntry) {
        self.boundaries
            .entry(chrom.to_string())
            .or_default()
            .entry(entry.start)
            .or_default()
            .push(entry);
    }

    /// Find all boundaries on a given chromosome whose interval is completely
    /// contained within [seg_start, seg_end).
    pub fn find_overlapping(
        &self,
        chrom: &str,
        seg_start: i64,
        seg_end: i64,
    ) -> Vec<&BoundaryEntry> {
        let mut results = Vec::new();
        if let Some(chrom_boundaries) = self.boundaries.get(chrom) {
            // Use BTreeMap range to efficiently find candidates
            for (_pos, entries) in chrom_boundaries.range(seg_start..seg_end) {
                for entry in entries {
                    // Boundary [entry.start, entry.end) must be fully within [seg_start, seg_end)
                    if entry.start >= seg_start && entry.end <= seg_end {
                        results.push(entry);
                    }
                }
            }
        }
        results
    }
}

/// Derive 5' and 3' boundary entries from an intron coordinate.
///
/// The boundary interval straddles the splice site with `anchor_length` bases
/// on each side:
/// - 5' boundary: `[intron_start - anchor_length, intron_start + anchor_length)`
/// - 3' boundary: `[intron_end - anchor_length, intron_end + anchor_length)`
pub fn intron_to_boundaries(
    chrom: &str,
    intron_start: i64,
    intron_end: i64,
    strand: Strand,
    anchor_length: i64,
) -> (BoundaryEntry, BoundaryEntry) {
    let five_prime = BoundaryEntry {
        boundary_id: format!(
            "{}:{}-{}",
            chrom,
            intron_start - anchor_length,
            intron_start + anchor_length
        ),
        start: intron_start - anchor_length,
        end: intron_start + anchor_length,
        boundary_type: BoundaryType::FivePrime,
        strand,
    };
    let three_prime = BoundaryEntry {
        boundary_id: format!(
            "{}:{}-{}",
            chrom,
            intron_end - anchor_length,
            intron_end + anchor_length
        ),
        start: intron_end - anchor_length,
        end: intron_end + anchor_length,
        boundary_type: BoundaryType::ThreePrime,
        strand,
    };
    (five_prime, three_prime)
}

/// Count boundary overlaps for a read's aligned segments.
///
/// For each aligned segment, finds boundaries in the index that are fully contained,
/// and increments counts accordingly.
#[allow(clippy::too_many_arguments)]
pub fn count_boundaries(
    chrom: &str,
    aligned_segments: &[(i64, i64)],
    boundary_index: &BoundaryIndex,
    cell_barcode: Option<&String>,
    umi: Option<&String>,
    strand: Strand,
    boundary_counts: &mut HashMap<String, HashMap<String, u32>>,
    boundary_totals: &mut HashMap<String, u32>,
    boundary_types: &mut HashMap<String, BoundaryType>,
    boundary_strands: &mut HashMap<String, Strand>,
    processed_boundary_reads: &mut HashMap<String, HashSet<u64>>,
    processed_boundary_umis: &mut HashMap<String, HashSet<u64>>,
    read_name_hash: u64,
    mode: Mode,
) {
    for (seg_start, seg_end) in aligned_segments {
        let overlapping = boundary_index.find_overlapping(chrom, *seg_start, *seg_end);
        for entry in overlapping {
            let key = &entry.boundary_id;

            // Check dedup: same read should not count same boundary twice (via u64 hash)
            let reads = processed_boundary_reads.entry(key.clone()).or_default();
            if !reads.insert(read_name_hash) {
                continue;
            }

            // Record type and strand
            boundary_types
                .entry(key.clone())
                .or_insert(entry.boundary_type);
            boundary_strands.entry(key.clone()).or_insert(strand);

            match mode {
                Mode::Single => {
                    if let Some(cb_str) = cell_barcode {
                        // UMI deduplication: if a UMI is present, ensure that
                        // (barcode, UMI) is unique per boundary before counting.
                        if let Some(umi_str) = umi {
                            let umi_hash = crate::types::hash_barcode_umi(cb_str, umi_str);
                            let umis = processed_boundary_umis.entry(key.clone()).or_default();
                            if !umis.insert(umi_hash) {
                                continue; // Same barcode+UMI already counted for this boundary
                            }
                        }
                        *boundary_counts
                            .entry(key.clone())
                            .or_default()
                            .entry(cb_str.clone())
                            .or_insert(0) += 1;
                    }
                }
                Mode::Bulk => {
                    *boundary_totals.entry(key.clone()).or_insert(0) += 1;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_intron_to_boundaries() {
        let (five_p, three_p) = intron_to_boundaries("chr2", 6545675, 6547042, Strand::Plus, 1);

        assert_eq!(five_p.boundary_id, "chr2:6545674-6545676");
        assert_eq!(five_p.start, 6545674);
        assert_eq!(five_p.end, 6545676);
        assert_eq!(five_p.boundary_type, BoundaryType::FivePrime);

        assert_eq!(three_p.boundary_id, "chr2:6547041-6547043");
        assert_eq!(three_p.start, 6547041);
        assert_eq!(three_p.end, 6547043);
        assert_eq!(three_p.boundary_type, BoundaryType::ThreePrime);
    }

    #[test]
    fn test_intron_to_boundaries_example2() {
        let (five_p, three_p) = intron_to_boundaries("chr1", 1000, 2000, Strand::Minus, 1);

        assert_eq!(five_p.boundary_id, "chr1:999-1001");
        assert_eq!(three_p.boundary_id, "chr1:1999-2001");
        assert_eq!(five_p.strand, Strand::Minus);
    }

    #[test]
    fn test_boundary_index_find_overlapping() {
        let mut index = BoundaryIndex::new();

        let (five_p, three_p) = intron_to_boundaries("chr2", 6545675, 6547042, Strand::Plus, 1);
        index.add("chr2", five_p);
        index.add("chr2", three_p);

        // Aligned segment that spans the 5' boundary
        let overlapping = index.find_overlapping("chr2", 6545600, 6545700);
        assert_eq!(overlapping.len(), 1);
        assert_eq!(overlapping[0].boundary_id, "chr2:6545674-6545676");
        assert_eq!(overlapping[0].boundary_type, BoundaryType::FivePrime);

        // Aligned segment that does NOT span the 5' boundary
        let overlapping = index.find_overlapping("chr2", 6545677, 6545800);
        assert_eq!(overlapping.len(), 0);

        // Aligned segment that spans the 3' boundary
        let overlapping = index.find_overlapping("chr2", 6547000, 6547100);
        assert_eq!(overlapping.len(), 1);
        assert_eq!(overlapping[0].boundary_id, "chr2:6547041-6547043");

        // Different chromosome
        let overlapping = index.find_overlapping("chr1", 6545600, 6545700);
        assert_eq!(overlapping.len(), 0);
    }

    #[test]
    fn test_count_boundaries_bulk() {
        let mut index = BoundaryIndex::new();
        let (five_p, three_p) = intron_to_boundaries("chr1", 1000, 2000, Strand::Plus, 1);
        index.add("chr1", five_p);
        index.add("chr1", three_p);

        let mut boundary_counts = HashMap::new();
        let mut boundary_totals = HashMap::new();
        let mut boundary_types = HashMap::new();
        let mut boundary_strands = HashMap::new();
        let mut processed = HashMap::new();
        let mut processed_umis = HashMap::new();

        // Segment that spans the 5' boundary at 999-1001
        let segments = vec![(900, 1100)];
        count_boundaries(
            "chr1",
            &segments,
            &index,
            None,
            None,
            Strand::Plus,
            &mut boundary_counts,
            &mut boundary_totals,
            &mut boundary_types,
            &mut boundary_strands,
            &mut processed,
            &mut processed_umis,
            12345, // read_name_hash
            Mode::Bulk,
        );

        assert_eq!(boundary_totals.get("chr1:999-1001"), Some(&1));
        assert_eq!(boundary_totals.get("chr1:1999-2001"), None); // Not overlapping
    }

    #[test]
    fn test_count_boundaries_dedup() {
        let mut index = BoundaryIndex::new();
        let (five_p, _three_p) = intron_to_boundaries("chr1", 1000, 2000, Strand::Plus, 1);
        index.add("chr1", five_p);

        let mut boundary_counts = HashMap::new();
        let mut boundary_totals = HashMap::new();
        let mut boundary_types = HashMap::new();
        let mut boundary_strands = HashMap::new();
        let mut processed = HashMap::new();
        let mut processed_umis = HashMap::new();

        let segments = vec![(900, 1100)];

        // Same read hash twice
        for _ in 0..2 {
            count_boundaries(
                "chr1",
                &segments,
                &index,
                None,
                None,
                Strand::Plus,
                &mut boundary_counts,
                &mut boundary_totals,
                &mut boundary_types,
                &mut boundary_strands,
                &mut processed,
                &mut processed_umis,
                12345,
                Mode::Bulk,
            );
        }

        assert_eq!(boundary_totals.get("chr1:999-1001"), Some(&1));
    }

    #[test]
    fn test_count_boundaries_umi_dedup() {
        // Same barcode + same UMI + same boundary → counted once
        let mut index = BoundaryIndex::new();
        let (five_p, _) = intron_to_boundaries("chr1", 1000, 2000, Strand::Plus, 1);
        index.add("chr1", five_p);

        let mut boundary_counts = HashMap::new();
        let mut boundary_totals = HashMap::new();
        let mut boundary_types = HashMap::new();
        let mut boundary_strands = HashMap::new();
        let mut processed = HashMap::new();
        let mut processed_umis = HashMap::new();

        let segments = vec![(900, 1100)];
        let barcode = "AAAA-1".to_string();
        let umi = "ACGT".to_string();

        // Two different reads with same barcode+UMI
        for hash in [11111u64, 22222u64] {
            count_boundaries(
                "chr1",
                &segments,
                &index,
                Some(&barcode),
                Some(&umi),
                Strand::Plus,
                &mut boundary_counts,
                &mut boundary_totals,
                &mut boundary_types,
                &mut boundary_strands,
                &mut processed,
                &mut processed_umis,
                hash,
                Mode::Single,
            );
        }

        assert_eq!(
            *boundary_counts
                .get("chr1:999-1001")
                .unwrap()
                .get("AAAA-1")
                .unwrap(),
            1,
            "Same barcode+UMI should be counted only once"
        );
    }

    #[test]
    fn test_count_boundaries_different_umis() {
        // Same barcode, different UMIs → each counted
        let mut index = BoundaryIndex::new();
        let (five_p, _) = intron_to_boundaries("chr1", 1000, 2000, Strand::Plus, 1);
        index.add("chr1", five_p);

        let mut boundary_counts = HashMap::new();
        let mut boundary_totals = HashMap::new();
        let mut boundary_types = HashMap::new();
        let mut boundary_strands = HashMap::new();
        let mut processed = HashMap::new();
        let mut processed_umis = HashMap::new();

        let segments = vec![(900, 1100)];
        let barcode = "AAAA-1".to_string();
        let umi1 = "ACGT".to_string();
        let umi2 = "TGCA".to_string();

        count_boundaries(
            "chr1",
            &segments,
            &index,
            Some(&barcode),
            Some(&umi1),
            Strand::Plus,
            &mut boundary_counts,
            &mut boundary_totals,
            &mut boundary_types,
            &mut boundary_strands,
            &mut processed,
            &mut processed_umis,
            11111,
            Mode::Single,
        );
        count_boundaries(
            "chr1",
            &segments,
            &index,
            Some(&barcode),
            Some(&umi2),
            Strand::Plus,
            &mut boundary_counts,
            &mut boundary_totals,
            &mut boundary_types,
            &mut boundary_strands,
            &mut processed,
            &mut processed_umis,
            22222,
            Mode::Single,
        );

        assert_eq!(
            *boundary_counts
                .get("chr1:999-1001")
                .unwrap()
                .get("AAAA-1")
                .unwrap(),
            2,
            "Different UMIs should be counted separately"
        );
    }
}
