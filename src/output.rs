//! Output file writing for junction and boundary results.

use flate2::write::GzEncoder;
use flate2::Compression;
use itertools::Itertools;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use crate::junction;
use crate::types::BoundaryType;

/// Write junction results in bulk mode (TSV format).
pub fn write_junction_bulk(
    output_prefix: &str,
    junction_totals: &HashMap<String, u32>,
    junction_strands: &HashMap<String, crate::types::Strand>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut output_file = GzEncoder::new(
        File::create(format!("{}_junction.tsv.gz", output_prefix))?,
        Compression::default(),
    );
    writeln!(output_file, "Junction\tStrand\tCount")?;
    for (key, count) in junction_totals.iter().sorted_by(|(a, _), (b, _)| a.cmp(b)) {
        let (coords, strand_str) = junction::parse_junction_key(key);
        let _ = junction_strands; // strand is encoded in the key
        writeln!(output_file, "{}\t{}\t{}", coords, strand_str, count)?;
    }
    Ok(())
}

/// Write junction results in single-cell mode (MatrixMarket + TSV format).
pub fn write_junction_single(
    output_prefix: &str,
    junction_counts: &HashMap<String, HashMap<String, u32>>,
    cell_barcodes: &std::collections::HashSet<String>,
    junction_strands: &HashMap<String, crate::types::Strand>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut matrix_file = GzEncoder::new(
        File::create(format!("{}_matrix.mtx.gz", output_prefix))?,
        Compression::default(),
    );
    let mut barcodes_file = GzEncoder::new(
        File::create(format!("{}_barcodes.tsv.gz", output_prefix))?,
        Compression::default(),
    );
    let mut features_file = GzEncoder::new(
        File::create(format!("{}_features.tsv.gz", output_prefix))?,
        Compression::default(),
    );
    let mut output_tsv = GzEncoder::new(
        File::create(format!("{}_junction_barcodes.tsv.gz", output_prefix))?,
        Compression::default(),
    );

    // Write barcodes.tsv.gz
    let barcode_list: Vec<_> = cell_barcodes.iter().sorted().collect();
    for barcode in &barcode_list {
        writeln!(barcodes_file, "{}", barcode)?;
    }

    // Write features.tsv.gz (with strand info)
    let feature_list: Vec<_> = junction_counts.keys().sorted().collect();
    for feature_key in &feature_list {
        let (coords, strand_str) = junction::parse_junction_key(feature_key);
        writeln!(features_file, "{}\t{}", coords, strand_str)?;
    }

    // Build barcode index map (1-based for MatrixMarket)
    let barcode_map: HashMap<_, _> = barcode_list
        .iter()
        .enumerate()
        .map(|(i, b)| (b.as_str(), i + 1))
        .collect();

    // Count total non-zero entries
    let total_entries: usize = junction_counts.values().map(|c| c.len()).sum();

    // Matrix Market header
    let mut matrix_buffer: Vec<String> = Vec::new();
    matrix_buffer.push("%%MatrixMarket matrix coordinate integer general".to_string());
    matrix_buffer.push("%".to_string());
    matrix_buffer.push(format!(
        "{} {} {}",
        feature_list.len(),
        barcode_list.len(),
        total_entries
    ));

    // TSV header
    let mut tsv_buffer: Vec<String> = Vec::new();
    tsv_buffer.push("Feature\tStrand\tBarcode\tCount".to_string());

    let _ = junction_strands; // strand info encoded in keys

    for (i, feature_key) in feature_list.iter().enumerate() {
        let (coords, strand_str) = junction::parse_junction_key(feature_key);
        if let Some(cell_counts) = junction_counts.get(*feature_key) {
            for (barcode, count) in cell_counts
                .iter()
                .sorted_by_key(|(b, _)| barcode_map.get(b.as_str()).copied().unwrap_or(0))
            {
                if let Some(&j) = barcode_map.get(barcode.as_str()) {
                    // FIX: j is already 1-based, use it directly for MatrixMarket
                    matrix_buffer.push(format!("{} {} {}", i + 1, j, count));
                    tsv_buffer.push(format!(
                        "{}\t{}\t{}\t{}",
                        coords, strand_str, barcode, count
                    ));
                }
            }
        }
    }

    for line in matrix_buffer {
        writeln!(matrix_file, "{}", line)?;
    }
    for line in tsv_buffer {
        writeln!(output_tsv, "{}", line)?;
    }

    Ok(())
}

/// Write boundary results in bulk mode (TSV format).
pub fn write_boundary_bulk(
    output_prefix: &str,
    boundary_totals: &HashMap<String, u32>,
    boundary_types: &HashMap<String, BoundaryType>,
    boundary_strands: &HashMap<String, crate::types::Strand>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut output_file = GzEncoder::new(
        File::create(format!("{}_boundary.tsv.gz", output_prefix))?,
        Compression::default(),
    );
    writeln!(output_file, "Boundary\tType\tStrand\tCount")?;
    for (key, count) in boundary_totals.iter().sorted_by(|(a, _), (b, _)| a.cmp(b)) {
        let btype = boundary_types
            .get(key)
            .map(|t| format!("{}", t))
            .unwrap_or_else(|| ".".to_string());
        let strand = boundary_strands
            .get(key)
            .map(|s| format!("{}", s))
            .unwrap_or_else(|| ".".to_string());
        writeln!(output_file, "{}\t{}\t{}\t{}", key, btype, strand, count)?;
    }
    Ok(())
}

/// Write boundary results in single-cell mode (MatrixMarket + TSV format).
pub fn write_boundary_single(
    output_prefix: &str,
    boundary_counts: &HashMap<String, HashMap<String, u32>>,
    cell_barcodes: &std::collections::HashSet<String>,
    boundary_types: &HashMap<String, BoundaryType>,
    boundary_strands: &HashMap<String, crate::types::Strand>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut matrix_file = GzEncoder::new(
        File::create(format!("{}_boundary_matrix.mtx.gz", output_prefix))?,
        Compression::default(),
    );
    let mut barcodes_file = GzEncoder::new(
        File::create(format!("{}_boundary_barcodes.tsv.gz", output_prefix))?,
        Compression::default(),
    );
    let mut features_file = GzEncoder::new(
        File::create(format!("{}_boundary_features.tsv.gz", output_prefix))?,
        Compression::default(),
    );
    let mut output_tsv = GzEncoder::new(
        File::create(format!("{}_boundary_barcodes_detail.tsv.gz", output_prefix))?,
        Compression::default(),
    );

    // Write barcodes.tsv.gz
    let barcode_list: Vec<_> = cell_barcodes.iter().sorted().collect();
    for barcode in &barcode_list {
        writeln!(barcodes_file, "{}", barcode)?;
    }

    // Write features.tsv.gz (boundary_id + type + strand)
    let feature_list: Vec<_> = boundary_counts.keys().sorted().collect();
    for feature_key in &feature_list {
        let btype = boundary_types
            .get(*feature_key)
            .map(|t| format!("{}", t))
            .unwrap_or_else(|| ".".to_string());
        let strand = boundary_strands
            .get(*feature_key)
            .map(|s| format!("{}", s))
            .unwrap_or_else(|| ".".to_string());
        writeln!(features_file, "{}\t{}\t{}", feature_key, btype, strand)?;
    }

    // Build barcode index map (1-based)
    let barcode_map: HashMap<_, _> = barcode_list
        .iter()
        .enumerate()
        .map(|(i, b)| (b.as_str(), i + 1))
        .collect();

    let total_entries: usize = boundary_counts.values().map(|c| c.len()).sum();

    let mut matrix_buffer: Vec<String> = Vec::new();
    matrix_buffer.push("%%MatrixMarket matrix coordinate integer general".to_string());
    matrix_buffer.push("%".to_string());
    matrix_buffer.push(format!(
        "{} {} {}",
        feature_list.len(),
        barcode_list.len(),
        total_entries
    ));

    let mut tsv_buffer: Vec<String> = Vec::new();
    tsv_buffer.push("Boundary\tType\tStrand\tBarcode\tCount".to_string());

    for (i, feature_key) in feature_list.iter().enumerate() {
        let btype = boundary_types
            .get(*feature_key)
            .map(|t| format!("{}", t))
            .unwrap_or_else(|| ".".to_string());
        let strand = boundary_strands
            .get(*feature_key)
            .map(|s| format!("{}", s))
            .unwrap_or_else(|| ".".to_string());
        if let Some(cell_counts) = boundary_counts.get(*feature_key) {
            for (barcode, count) in cell_counts
                .iter()
                .sorted_by_key(|(b, _)| barcode_map.get(b.as_str()).copied().unwrap_or(0))
            {
                if let Some(&j) = barcode_map.get(barcode.as_str()) {
                    matrix_buffer.push(format!("{} {} {}", i + 1, j, count));
                    tsv_buffer.push(format!(
                        "{}\t{}\t{}\t{}\t{}",
                        feature_key, btype, strand, barcode, count
                    ));
                }
            }
        }
    }

    for line in matrix_buffer {
        writeln!(matrix_file, "{}", line)?;
    }
    for line in tsv_buffer {
        writeln!(output_tsv, "{}", line)?;
    }

    Ok(())
}
