//! # Tosa
//!
//! Fast junction and exon-intron boundary read counting from RNA-seq/scRNA-seq BAM/CRAM files.
//!
//! Tosa processes mapped BAM/CRAM files to extract:
//! - **Junction read counts**: Reads spanning splice junctions (identified by `N` CIGAR operations)
//! - **Boundary read counts**: Reads overlapping exon-intron boundaries (when GTF annotation is provided)
//!
//! Supports both bulk RNA-seq and single-cell RNA-seq (10x Genomics-style barcodes),
//! with configurable strand specificity (unstranded, XS, RF, FR).

pub mod bam_reader;
pub mod boundary;
pub mod cli;
pub mod data_loader;
pub mod gtf;
pub mod junction;
pub mod output;
pub mod types;

use log::info;
use std::collections::HashSet;

use types::{Mode, RunConfig};

/// Run the Tosa pipeline with the given configuration.
///
/// This is the main entry point for the library, executing the full pipeline:
/// loading barcodes, parsing GTF, processing BAM, and writing output.
pub fn run(config: &RunConfig) -> Result<(), Box<dyn std::error::Error>> {
    // Log configuration
    info!("Running tosa v{}", env!("CARGO_PKG_VERSION"));
    info!("Mode: {}", config.mode);
    info!("Input file: {}", config.bam_file);
    info!("Output prefix: {}", config.output_prefix);
    info!("Minimum anchor length: {}", config.min_anchor_length);
    info!(
        "Boundary anchor length: {}",
        config.min_boundary_anchor_length
    );
    info!("Minimum intron length: {}", config.min_intron_length);
    info!("Maximum intron length: {}", config.max_intron_length);
    info!("Maximum loci (NH): {}", config.max_loci);
    info!("Strand mode: {}", config.strand_mode);
    info!("Threads: {}", config.threads);

    // Load cell barcodes of interest (single mode only)
    let cell_barcodes_of_interest = if config.mode == Mode::Single {
        let barcodes = data_loader::load_cell_barcodes(config.cell_barcode_file.as_ref())?;
        info!(
            "Cell barcodes of interest: {}",
            if barcodes.is_empty() {
                "None (processing all reads)".to_string()
            } else {
                format!("{} barcodes", barcodes.len())
            }
        );
        barcodes
    } else {
        HashSet::new()
    };

    // Parse GTF for boundary counting (if provided)
    let boundary_index = if let Some(ref gtf_path) = config.gtf_file {
        info!("GTF file: {}", gtf_path);
        Some(gtf::parse_gtf(gtf_path, config.min_boundary_anchor_length)?)
    } else {
        None
    };

    // Process BAM file
    let result = bam_reader::process_bam_records(
        config,
        &cell_barcodes_of_interest,
        boundary_index.as_ref(),
    )
    .map_err(|e| -> Box<dyn std::error::Error> { e })?;

    // Write output files
    info!("Writing output files");
    if config.mode == Mode::Single {
        output::write_junction_single(
            &config.output_prefix,
            &result.junction_counts,
            &result.cell_barcodes,
            &result.junction_strands,
        )?;

        if boundary_index.is_some() {
            output::write_boundary_single(
                &config.output_prefix,
                &result.boundary_counts,
                &result.cell_barcodes,
                &result.boundary_types,
                &result.boundary_strands,
            )?;
        }
    } else {
        output::write_junction_bulk(
            &config.output_prefix,
            &result.junction_totals,
            &result.junction_strands,
        )?;

        if boundary_index.is_some() {
            output::write_boundary_bulk(
                &config.output_prefix,
                &result.boundary_totals,
                &result.boundary_types,
                &result.boundary_strands,
            )?;
        }
    }

    info!("Finished processing");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_bam_path() -> String {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples/example.bam")
            .to_str()
            .unwrap()
            .to_string()
    }

    fn test_gtf_path() -> String {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples/annotation.gtf")
            .to_str()
            .unwrap()
            .to_string()
    }

    fn test_barcodes_path() -> String {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples/barcodes.tsv")
            .to_str()
            .unwrap()
            .to_string()
    }

    /// Ensure info! arguments are evaluated (covers L47-48: empty barcodes path)
    #[test]
    fn test_run_bulk_with_logger() {
        let _ = env_logger::Builder::from_default_env()
            .filter_level(log::LevelFilter::Info)
            .try_init();

        let tmpdir = tempfile::tempdir().unwrap();
        let prefix = tmpdir.path().join("lib_bulk").to_str().unwrap().to_string();

        let config = types::RunConfig {
            mode: types::Mode::Bulk,
            bam_file: test_bam_path(),
            output_prefix: prefix,
            min_anchor_length: 8,
            min_boundary_anchor_length: 8,
            min_intron_length: 20,
            max_intron_length: 500000,
            max_loci: 1,
            cell_barcode_file: None,
            strand_mode: types::StrandMode::Unstranded,
            gtf_file: Some(test_gtf_path()),
            verbose: false,
            threads: 1,
        };

        run(&config).unwrap();
    }

    /// Covers L50: the else branch of barcodes.is_empty() (non-empty barcodes)
    #[test]
    fn test_run_single_with_logger() {
        let _ = env_logger::Builder::from_default_env()
            .filter_level(log::LevelFilter::Info)
            .try_init();

        let tmpdir = tempfile::tempdir().unwrap();
        let prefix = tmpdir
            .path()
            .join("lib_single")
            .to_str()
            .unwrap()
            .to_string();

        let config = types::RunConfig {
            mode: types::Mode::Single,
            bam_file: test_bam_path(),
            output_prefix: prefix,
            min_anchor_length: 8,
            min_boundary_anchor_length: 8,
            min_intron_length: 20,
            max_intron_length: 500000,
            max_loci: 1,
            cell_barcode_file: Some(test_barcodes_path()),
            strand_mode: types::StrandMode::Unstranded,
            gtf_file: Some(test_gtf_path()),
            verbose: false,
            threads: 1,
        };

        run(&config).unwrap();
    }
}
