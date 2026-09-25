//! CLI argument definition for Tosa.

use crate::types::{Mode, RunConfig, StrandMode};
use clap::{Arg, Command};

/// Build the CLI command definition.
pub fn build_cli() -> Command {
    Command::new("tosa")
        .version(env!("CARGO_PKG_VERSION"))
        .author("NaotoKubota")
        .about("Extract junction and boundary reads from RNA-seq/scRNA-seq BAM/CRAM files")
        .arg(Arg::new("mode")
            .required(true)
            .value_parser(["bulk", "single"])
            .help("Mode of operation: 'bulk' or 'single'"))
        .arg(Arg::new("bam_file")
            .required(true)
            .help("Path to the BAM/CRAM file"))
        .arg(Arg::new("output_prefix")
            .required(true)
            .help("Output prefix for the output files"))
        .arg(Arg::new("anchor_length")
            .short('a')
            .long("anchor-length")
            .default_value("8")
            .value_parser(clap::value_parser!(i64))
            .help("Minimum anchor length for both sides of junctions"))
        .arg(Arg::new("boundary_anchor_length")
            .short('b')
            .long("boundary-anchor-length")
            .default_value("8")
            .value_parser(clap::value_parser!(i64))
            .help("Minimum anchor length on each side of exon-intron boundaries"))
        .arg(Arg::new("min_intron_length")
            .short('m')
            .long("min-intron-length")
            .default_value("20")
            .value_parser(clap::value_parser!(i64))
            .help("Minimum intron length for junctions"))
        .arg(Arg::new("max_intron_length")
            .short('M')
            .long("max-intron-length")
            .default_value("500000")
            .value_parser(clap::value_parser!(i64))
            .help("Maximum intron length for junctions"))
        .arg(Arg::new("max_loci")
            .short('l')
            .long("max-loci")
            .default_value("1")
            .value_parser(clap::value_parser!(u32))
            .help("Maximum number of loci the read maps to"))
        .arg(Arg::new("cell_barcode_file")
            .short('c')
            .long("cell-barcodes")
            .value_parser(clap::value_parser!(String))
            .help("Optional file specifying cell barcodes of interest"))
        .arg(Arg::new("strand")
            .short('s')
            .long("strand")
            .value_parser(["RF", "FR", "XS"])
            .help("Strand specificity of RNA library: RF (first-strand), FR (second-strand), XS (use XS tags). Omit for unstranded"))
        .arg(Arg::new("gtf_file")
            .short('g')
            .long("gtf")
            .value_parser(clap::value_parser!(String))
            .help("GTF annotation file for exon-intron boundary read counting"))
        .arg(Arg::new("verbose")
            .short('v')
            .long("verbose")
            .action(clap::ArgAction::SetTrue)
            .help("Enable verbose output to print all arguments"))
        .arg(Arg::new("threads")
            .short('p')
            .long("threads")
            .default_value("1")
            .value_parser(clap::value_parser!(usize))
            .help("Number of threads for parallel processing"))
}

/// Parse CLI matches into a RunConfig.
pub fn parse_config(matches: &clap::ArgMatches) -> RunConfig {
    RunConfig {
        mode: match matches.get_one::<String>("mode").unwrap().as_str() {
            "single" => Mode::Single,
            _ => Mode::Bulk,
        },
        bam_file: matches.get_one::<String>("bam_file").unwrap().clone(),
        output_prefix: matches.get_one::<String>("output_prefix").unwrap().clone(),
        min_anchor_length: *matches.get_one::<i64>("anchor_length").unwrap(),
        min_boundary_anchor_length: *matches.get_one::<i64>("boundary_anchor_length").unwrap(),
        min_intron_length: *matches.get_one::<i64>("min_intron_length").unwrap(),
        max_intron_length: *matches.get_one::<i64>("max_intron_length").unwrap(),
        max_loci: *matches.get_one::<u32>("max_loci").unwrap(),
        cell_barcode_file: matches.get_one::<String>("cell_barcode_file").cloned(),
        strand_mode: StrandMode::from_str_opt(matches.get_one::<String>("strand")),
        gtf_file: matches.get_one::<String>("gtf_file").cloned(),
        verbose: matches.get_flag("verbose"),
        threads: *matches.get_one::<usize>("threads").unwrap(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_cli_defaults() {
        let matches = build_cli().get_matches_from(["tosa", "bulk", "test.bam", "out_prefix"]);
        let config = parse_config(&matches);

        assert_eq!(config.mode, Mode::Bulk);
        assert_eq!(config.bam_file, "test.bam");
        assert_eq!(config.output_prefix, "out_prefix");
        assert_eq!(config.min_anchor_length, 8);
        assert_eq!(config.min_boundary_anchor_length, 8);
        assert_eq!(config.min_intron_length, 20);
        assert_eq!(config.max_intron_length, 500000);
        assert_eq!(config.max_loci, 1);
        assert_eq!(config.cell_barcode_file, None);
        assert_eq!(config.strand_mode, StrandMode::Unstranded);
        assert_eq!(config.gtf_file, None);
        assert!(!config.verbose);
        assert_eq!(config.threads, 1);
    }

    #[test]
    fn test_build_cli_all_options() {
        let matches = build_cli().get_matches_from([
            "tosa",
            "-a",
            "10",
            "-b",
            "3",
            "-m",
            "50",
            "-M",
            "1000000",
            "-l",
            "3",
            "-c",
            "barcodes.tsv",
            "-s",
            "RF",
            "-g",
            "annotation.gtf",
            "-v",
            "-p",
            "8",
            "single",
            "input.bam",
            "output",
        ]);
        let config = parse_config(&matches);

        assert_eq!(config.mode, Mode::Single);
        assert_eq!(config.bam_file, "input.bam");
        assert_eq!(config.output_prefix, "output");
        assert_eq!(config.min_anchor_length, 10);
        assert_eq!(config.min_boundary_anchor_length, 3);
        assert_eq!(config.min_intron_length, 50);
        assert_eq!(config.max_intron_length, 1000000);
        assert_eq!(config.max_loci, 3);
        assert_eq!(config.cell_barcode_file, Some("barcodes.tsv".to_string()));
        assert_eq!(config.strand_mode, StrandMode::RF);
        assert_eq!(config.gtf_file, Some("annotation.gtf".to_string()));
        assert!(config.verbose);
        assert_eq!(config.threads, 8);
    }

    #[test]
    fn test_parse_config_strand_modes() {
        // FR
        let matches = build_cli().get_matches_from(["tosa", "-s", "FR", "bulk", "t.bam", "o"]);
        assert_eq!(parse_config(&matches).strand_mode, StrandMode::FR);

        // XS
        let matches = build_cli().get_matches_from(["tosa", "-s", "XS", "bulk", "t.bam", "o"]);
        assert_eq!(parse_config(&matches).strand_mode, StrandMode::XS);

        // Unstranded (no -s flag)
        let matches = build_cli().get_matches_from(["tosa", "bulk", "t.bam", "o"]);
        assert_eq!(parse_config(&matches).strand_mode, StrandMode::Unstranded);
    }

    #[test]
    fn test_build_cli_invalid_mode() {
        let result = build_cli().try_get_matches_from(["tosa", "invalid", "test.bam", "out"]);
        assert!(result.is_err());
    }
}
