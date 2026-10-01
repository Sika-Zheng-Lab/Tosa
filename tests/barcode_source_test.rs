use rust_htslib::bam::{
    self,
    record::{Aux, Cigar, CigarString},
};
use std::{collections::HashSet, process::Command};
use tosa::{
    bam_reader::{process_bam_records, BarcodeStats},
    types::{BarcodeSource, RunConfig},
};

fn fixture(dir: &std::path::Path, cb: bool, ub: bool) -> RunConfig {
    let path = dir.join(format!("reads_{cb}_{ub}.bam"));
    let mut header = bam::Header::new();
    for chrom in ["chr1", "chr2"] {
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", chrom)
                .push_tag(b"LN", 10000),
        );
    }
    let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    for tid in 0..2 {
        for cell in ["SRR1", "SRR2"] {
            for junction in [true, false] {
                // Two mates with the same QNAME, plus a distinct read with the same UMI.
                for read in [1, 1, 2] {
                    let name = format!("{cell}.{tid}_{junction}_{read}");
                    let cigar = if junction {
                        vec![Cigar::Match(20), Cigar::RefSkip(300), Cigar::Match(20)]
                    } else {
                        vec![Cigar::Match(40)]
                    };
                    let mut record = bam::Record::new();
                    record.set(
                        name.as_bytes(),
                        Some(&CigarString(cigar)),
                        &[b'A'; 40],
                        &[30; 40],
                    );
                    record.set_tid(tid);
                    record.set_pos(1180);
                    record.set_flags(0);
                    record.push_aux(b"NH", Aux::U8(1)).unwrap();
                    if cb {
                        record.push_aux(b"CB", Aux::String(cell)).unwrap();
                    }
                    if ub {
                        record.push_aux(b"UB", Aux::String("ACGT")).unwrap();
                    }
                    writer.write(&record).unwrap();
                }
            }
        }
    }
    drop(writer);
    bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
    let gtf = dir.join("annotation.gtf");
    let mut annotation = String::new();
    for chrom in ["chr1", "chr2"] {
        for (start, end) in [(1001, 1200), (1501, 1700)] {
            annotation.push_str(&format!("{chrom}\ttest\texon\t{start}\t{end}\t.\t+\t.\tgene_id \"g{chrom}\"; transcript_id \"t{chrom}\";\n"));
        }
    }
    std::fs::write(&gtf, annotation).unwrap();
    let matches = tosa::cli::build_cli().get_matches_from([
        "tosa",
        "single",
        path.to_str().unwrap(),
        dir.join("out").to_str().unwrap(),
        "-g",
        gtf.to_str().unwrap(),
    ]);
    tosa::cli::parse_config(&matches)
}

#[test]
fn cb_and_qname_counts_match_with_and_without_umi() {
    let dir = tempfile::tempdir().unwrap();
    for ub in [false, true] {
        let cb = fixture(dir.path(), true, ub);
        let index = tosa::gtf::parse_gtf(cb.gtf_file.as_ref().unwrap(), 8).unwrap();
        let expected = process_bam_records(&cb, &HashSet::new(), Some(&index)).unwrap();
        assert_eq!(expected.junction_counts.len(), 2);
        assert_eq!(expected.boundary_counts.len(), 2);
        for counts in expected
            .junction_counts
            .values()
            .chain(expected.boundary_counts.values())
        {
            assert_eq!(counts.len(), 2);
            assert!(counts.values().all(|&n| n == if ub { 1 } else { 2 }));
        }
        let mut qname = fixture(dir.path(), false, ub);
        qname.barcode_source = BarcodeSource::Qname;
        qname.barcode_regex = Some(r"^(SRR[0-9]+)\.".into());
        for threads in [1, 2] {
            qname.threads = threads;
            let actual = process_bam_records(&qname, &HashSet::new(), Some(&index)).unwrap();
            assert_eq!(actual.junction_counts, expected.junction_counts);
            assert_eq!(actual.boundary_counts, expected.boundary_counts);
            assert_eq!(
                actual.barcode_stats,
                BarcodeStats {
                    processed_records: 24,
                    eligible_records: 24,
                    extracted: 24,
                    missing: 0,
                    whitelist_excluded: 0
                }
            );
        }
        qname.cell_barcode_file = Some("whitelist".into());
        let filtered =
            process_bam_records(&qname, &HashSet::from(["SRR1".into()]), Some(&index)).unwrap();
        assert_eq!(filtered.barcode_stats.whitelist_excluded, 12);
        assert_eq!(filtered.cell_barcodes, HashSet::from(["SRR1".into()]));
        assert!(filtered
            .junction_counts
            .values()
            .chain(filtered.boundary_counts.values())
            .all(|counts| counts.len() == 1 && counts.contains_key("SRR1")));
        let excluded =
            process_bam_records(&qname, &HashSet::from(["absent".into()]), Some(&index)).unwrap();
        assert_eq!(excluded.barcode_stats.whitelist_excluded, 24);
        assert!(excluded.cell_barcodes.is_empty());
        qname.cell_barcode_file = None;
        qname.barcode_regex = Some(r"^(SRR1)\.".into());
        let partial = process_bam_records(&qname, &HashSet::new(), None).unwrap();
        assert_eq!(partial.barcode_stats.extracted, 12);
        assert_eq!(partial.barcode_stats.missing, 12);
    }
}

#[test]
fn missing_ids_fail_before_output_and_report_diagnostics() {
    let dir = tempfile::tempdir().unwrap();
    let config = fixture(dir.path(), false, false);
    for extra in [
        vec![],
        vec!["--barcode-source", "qname", "--barcode-regex", "^(absent)"],
        vec!["--barcode-source", "qname", "--barcode-regex", "^()"],
        vec!["--barcode-source", "qname", "--barcode-regex", "^(absent)?"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_tosa"))
            .args(["single", &config.bam_file, &config.output_prefix])
            .args(extra)
            .output()
            .unwrap();
        assert!(!output.status.success());
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains("extracted=0, missing=24"), "{stderr}");
        assert!(
            stderr.contains("No cell IDs could be extracted"),
            "{stderr}"
        );
        assert!(!std::path::Path::new(&format!("{}_matrix.mtx.gz", config.output_prefix)).exists());
    }
    // No extraction candidates after NH filtering must not be misdiagnosed as missing IDs.
    let mut filtered = config.clone();
    filtered.max_loci = 0;
    let result = process_bam_records(&filtered, &HashSet::new(), None).unwrap();
    assert_eq!(result.barcode_stats.processed_records, 24);
    assert_eq!(result.barcode_stats.eligible_records, 0);
}

#[test]
fn invalid_barcode_options_are_rejected_before_opening_bam() {
    for extra in [
        vec!["--barcode-source", "qname"],
        vec!["--barcode-regex", "(x)"],
        vec!["--barcode-source", "qname", "--barcode-regex", "["],
        vec!["--barcode-source", "qname", "--barcode-regex", "^SRR"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_tosa"))
            .args(["single", "nonexistent.bam", "unused"])
            .args(extra)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8(output.stderr)
            .unwrap()
            .contains("--barcode-"));
    }
    let matches = tosa::cli::build_cli().get_matches_from([
        "tosa",
        "bulk",
        "missing.bam",
        "unused",
        "--barcode-source",
        "qname",
        "--barcode-regex",
        "(x)",
    ]);
    assert!(tosa::cli::parse_config(&matches)
        .compile_barcode_regex()
        .unwrap_err()
        .to_string()
        .contains("single mode"));
}

#[test]
fn empty_input_succeeds_and_qname_never_falls_back_to_cb() {
    let dir = tempfile::tempdir().unwrap();
    let mut config = fixture(dir.path(), true, false);
    config.barcode_source = BarcodeSource::Qname;
    config.barcode_regex = Some("^(absent)".into());
    assert!(process_bam_records(&config, &HashSet::new(), None).is_err());
    let path = dir.path().join("empty.bam");
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 10000),
    );
    drop(bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap());
    bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
    config.bam_file = path.to_str().unwrap().into();
    let result = process_bam_records(&config, &HashSet::new(), None).unwrap();
    assert_eq!(result.barcode_stats, BarcodeStats::default());
}
