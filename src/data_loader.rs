//! Data loading utilities for Tosa.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Load cell barcodes from a file into a HashSet.
///
/// Returns an empty set if no file path is provided.
pub fn load_cell_barcodes(
    file_path: Option<&String>,
) -> Result<HashSet<String>, Box<dyn std::error::Error>> {
    let mut barcodes = HashSet::new();
    if let Some(path) = file_path {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        for line in reader.lines() {
            let barcode = line?.trim().to_string();
            if !barcode.is_empty() {
                barcodes.insert(barcode);
            }
        }
    }
    Ok(barcodes)
}
