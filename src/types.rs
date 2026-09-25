//! Common type definitions for Tosa.

use std::collections::hash_map::DefaultHasher;
use std::fmt;
use std::hash::{Hash, Hasher};

/// Mode of operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Bulk,
    Single,
}

impl fmt::Display for Mode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Mode::Bulk => write!(f, "bulk"),
            Mode::Single => write!(f, "single"),
        }
    }
}

/// Strand specificity mode for RNA-seq library preparation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StrandMode {
    /// No strand specificity. Uses XS tag if available, otherwise Unknown.
    Unstranded,
    /// Use XS tags provided by aligner.
    XS,
    /// First-strand (RF): read1 reverse = +, read1 forward = -.
    RF,
    /// Second-strand (FR): read1 forward = +, read1 reverse = -.
    FR,
}

impl StrandMode {
    /// Parse strand mode from a string argument.
    pub fn from_str_opt(s: Option<&String>) -> Self {
        match s.map(|v| v.as_str()) {
            Some("RF") => StrandMode::RF,
            Some("FR") => StrandMode::FR,
            Some("XS") => StrandMode::XS,
            _ => StrandMode::Unstranded,
        }
    }
}

impl fmt::Display for StrandMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StrandMode::Unstranded => write!(f, "unstranded"),
            StrandMode::XS => write!(f, "XS"),
            StrandMode::RF => write!(f, "RF"),
            StrandMode::FR => write!(f, "FR"),
        }
    }
}

/// Strand of a read or feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

/// Compact junction key for internal use. Zero heap allocation (24 bytes, stack-only).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct JunctionKey {
    pub tid: u32,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
}

impl JunctionKey {
    /// Convert to the string format "chrom:start-end:strand" used in output files.
    pub fn to_string_key(&self, reference_names: &[String]) -> String {
        format!(
            "{}:{}-{}:{}",
            reference_names[self.tid as usize], self.start, self.end, self.strand
        )
    }
}

/// Hash a read name (as bytes) to u64 for memory-efficient deduplication.
///
/// With 64-bit hashes the collision probability is negligible
/// (~1e-3 for 100M reads, i.e. n²/2⁶⁴).
pub fn hash_read_name(name: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    name.hash(&mut hasher);
    hasher.finish()
}

/// Hash a (barcode, UMI) pair to u64 for memory-efficient UMI deduplication.
///
/// Used in single-cell mode to ensure that reads with the same cell barcode
/// and UMI mapping to the same junction/boundary are counted only once.
pub fn hash_barcode_umi(barcode: &str, umi: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    barcode.hash(&mut hasher);
    umi.hash(&mut hasher);
    hasher.finish()
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

/// Type of exon-intron boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum BoundaryType {
    /// 5' splice site (exon → intron boundary).
    FivePrime,
    /// 3' splice site (intron → exon boundary).
    ThreePrime,
}

impl fmt::Display for BoundaryType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BoundaryType::FivePrime => write!(f, "5p"),
            BoundaryType::ThreePrime => write!(f, "3p"),
        }
    }
}

/// Configuration for a Tosa run, parsed from CLI arguments.
#[derive(Debug, Clone)]
pub struct RunConfig {
    /// Mode of operation: Bulk or Single.
    pub mode: Mode,
    /// Path to the BAM/CRAM file.
    pub bam_file: String,
    /// Output prefix for output files.
    pub output_prefix: String,
    /// Minimum anchor length for both sides of junctions.
    pub min_anchor_length: i64,
    /// Minimum anchor length on each side of exon-intron boundaries.
    pub min_boundary_anchor_length: i64,
    /// Minimum intron length for junctions.
    pub min_intron_length: i64,
    /// Maximum intron length for junctions.
    pub max_intron_length: i64,
    /// Maximum number of loci the read maps to (NH tag).
    pub max_loci: u32,
    /// Optional path to cell barcode file.
    pub cell_barcode_file: Option<String>,
    /// Strand specificity mode.
    pub strand_mode: StrandMode,
    /// Optional path to GTF annotation file.
    pub gtf_file: Option<String>,
    /// Enable verbose (debug) logging.
    pub verbose: bool,
    /// Number of threads for parallel processing.
    pub threads: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_mode_from_str_opt() {
        let rf = "RF".to_string();
        let fr = "FR".to_string();
        let xs = "XS".to_string();

        assert_eq!(StrandMode::from_str_opt(Some(&rf)), StrandMode::RF);
        assert_eq!(StrandMode::from_str_opt(Some(&fr)), StrandMode::FR);
        assert_eq!(StrandMode::from_str_opt(Some(&xs)), StrandMode::XS);
        assert_eq!(StrandMode::from_str_opt(None), StrandMode::Unstranded);
    }

    #[test]
    fn test_boundary_type_display() {
        assert_eq!(format!("{}", BoundaryType::FivePrime), "5p");
        assert_eq!(format!("{}", BoundaryType::ThreePrime), "3p");
    }

    #[test]
    fn test_strand_display() {
        assert_eq!(format!("{}", Strand::Plus), "+");
        assert_eq!(format!("{}", Strand::Minus), "-");
        assert_eq!(format!("{}", Strand::Unknown), ".");
    }

    #[test]
    fn test_strand_mode_display() {
        assert_eq!(format!("{}", StrandMode::Unstranded), "unstranded");
        assert_eq!(format!("{}", StrandMode::XS), "XS");
        assert_eq!(format!("{}", StrandMode::RF), "RF");
        assert_eq!(format!("{}", StrandMode::FR), "FR");
    }

    #[test]
    fn test_mode_display() {
        assert_eq!(format!("{}", Mode::Bulk), "bulk");
        assert_eq!(format!("{}", Mode::Single), "single");
    }

    #[test]
    fn test_junction_key() {
        let key = JunctionKey {
            tid: 0,
            start: 100,
            end: 200,
            strand: Strand::Plus,
        };
        let ref_names = vec!["chr1".to_string(), "chr2".to_string()];
        assert_eq!(key.to_string_key(&ref_names), "chr1:100-200:+");
    }

    #[test]
    fn test_hash_read_name() {
        let h1 = hash_read_name(b"read1");
        let h2 = hash_read_name(b"read2");
        let h1_again = hash_read_name(b"read1");
        assert_eq!(h1, h1_again);
        assert_ne!(h1, h2);
    }
}
