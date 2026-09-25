//! Generates a small synthetic BAM file (and optionally a CRAM file) for trying out Tosa.
//!
//! Usage:
//!   cargo run --example generate_example_bam [output_dir]
//!
//! Creates the following files in the specified directory (defaults to `examples/`):
//!   - `example.bam` and `example.bam.bai`
//!   - `example.cram` and `example.cram.crai`
//!
//! The BAM/CRAM contains reads mapped to chr1 with:
//! - 13 unique junction reads spanning intron 1 (chr1:1201-1499)
//!   - 10 with XS:A:+ and 3 with XS:A:-
//!   - 1 duplicate pair (same QNAME) for dedup testing
//! - 5 junction reads spanning intron 2 (chr1:1701-1999), all XS:A:+
//! - 6 non-spliced reads overlapping exon-intron boundaries
//! - 1 multi-mapped read (NH:i:2, filtered by default)
//! - All reads carry a CB (cell barcode) tag for single-cell mode testing:
//!   AAAA-1 (junction 1 + strand), BBBB-1 (junction 1 - strand),
//!   CCCC-1 (junction 2 + boundary reads)
//! - All reads carry a UB (UMI) tag for UMI deduplication testing:
//!   Some reads within the same cell/junction share UMIs to verify dedup.
//!   After UMI dedup: AAAA-1 j1=6, BBBB-1 j1=2, CCCC-1 j2=3
//!
//! Matching GTF annotation is in `examples/annotation.gtf`.

use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::{
    self, record::Aux, record::Cigar, record::CigarString, Read, Record, Writer,
};
use std::io::Write;

fn main() {
    let out_dir = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "examples".to_string());
    let bam_path = format!("{}/example.bam", out_dir);

    // --- BAM header: one reference sequence chr1 (10 kb) ---
    let mut header = Header::new();
    let mut sq = HeaderRecord::new(b"SQ");
    sq.push_tag(b"SN", "chr1");
    sq.push_tag(b"LN", "10000");
    header.push_record(&sq);

    let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam)
        .expect("Failed to create BAM writer");

    // --- Gene model (GTF 1-based, + strand) ---
    //   Exon 1: 1000–1200   Exon 2: 1500–1700   Exon 3: 2000–2300
    //   Intron 1: 1201–1499 (0-based 1200..1499, len 299)
    //   Intron 2: 1701–1999 (0-based 1700..1999, len 299)
    //
    // Junction keys (1-based): chr1:1201-1499, chr1:1701-1999
    // Boundaries (0-based pairs, anchor_length=1):
    //   Intron 1 → 5': chr1:1199-1201  3': chr1:1498-1500
    //   Intron 2 → 5': chr1:1699-1701  3': chr1:1998-2000

    struct ReadSpec {
        qname: &'static str,
        pos: i64,
        cigar: Vec<Cigar>,
        xs: Option<u8>,
        nh: u8,
        cb: &'static str,
        ub: &'static str,
    }

    // All reads are 100 bp, sorted by position for coordinate-sorted BAM.
    // + and - strand reads are interleaved by position to maintain sort order.
    // Cell barcodes: AAAA-1 = junction1 + strand, BBBB-1 = junction1 - strand, CCCC-1 = junction2 + boundaries
    //
    // UMI (UB tag) assignments – designed so that UMI deduplication in single-cell
    // mode reduces counts:
    //   AAAA-1 junction1: 10 reads, 6 unique UMIs → deduped count = 6
    //   BBBB-1 junction1:  3 reads, 2 unique UMIs → deduped count = 2
    //   CCCC-1 junction2:  5 reads, 3 unique UMIs → deduped count = 3
    //   CCCC-1 5' bdry intron1: 2 reads, 1 UMI → deduped count = 1
    //   CCCC-1 3' bdry intron1: 2 reads, 1 UMI → deduped count = 1
    //   CCCC-1 5' bdry intron2: 1 read  → deduped count = 1
    //   CCCC-1 3' bdry intron2: 1 read  → deduped count = 1
    let reads: Vec<ReadSpec> = vec![
        // ── Junction 1 reads (intron at 0-based 1200..1499) ─────────
        ReadSpec {
            qname: "read01",
            pos: 1140,
            cigar: vec![Cigar::Match(60), Cigar::RefSkip(299), Cigar::Match(40)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "AAAA",
        },
        ReadSpec {
            qname: "read02",
            pos: 1145,
            cigar: vec![Cigar::Match(55), Cigar::RefSkip(299), Cigar::Match(45)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "AAAA",
        },
        ReadSpec {
            qname: "read03",
            pos: 1148,
            cigar: vec![Cigar::Match(52), Cigar::RefSkip(299), Cigar::Match(48)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "AAAA",
        },
        ReadSpec {
            qname: "read10",
            pos: 1148,
            cigar: vec![Cigar::Match(52), Cigar::RefSkip(299), Cigar::Match(48)],
            xs: Some(b'-'),
            nh: 1,
            cb: "BBBB-1",
            ub: "GGGG",
        },
        ReadSpec {
            qname: "read04",
            pos: 1150,
            cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "BBBB",
        },
        ReadSpec {
            qname: "read05",
            pos: 1150,
            cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "BBBB",
        },
        // Duplicate pair: same QNAME → counted only once (read-name dedup)
        ReadSpec {
            qname: "read_dup",
            pos: 1150,
            cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "CCCC",
        },
        ReadSpec {
            qname: "read_dup",
            pos: 1150,
            cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "CCCC",
        },
        ReadSpec {
            qname: "read06",
            pos: 1150,
            cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "DDDD",
        },
        ReadSpec {
            qname: "read11",
            pos: 1150,
            cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)],
            xs: Some(b'-'),
            nh: 1,
            cb: "BBBB-1",
            ub: "GGGG",
        },
        ReadSpec {
            qname: "read07",
            pos: 1152,
            cigar: vec![Cigar::Match(48), Cigar::RefSkip(299), Cigar::Match(52)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "DDDD",
        },
        ReadSpec {
            qname: "read12",
            pos: 1152,
            cigar: vec![Cigar::Match(48), Cigar::RefSkip(299), Cigar::Match(52)],
            xs: Some(b'-'),
            nh: 1,
            cb: "BBBB-1",
            ub: "HHHH",
        },
        ReadSpec {
            qname: "read08",
            pos: 1155,
            cigar: vec![Cigar::Match(45), Cigar::RefSkip(299), Cigar::Match(55)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "EEEE",
        },
        ReadSpec {
            qname: "read09",
            pos: 1160,
            cigar: vec![Cigar::Match(40), Cigar::RefSkip(299), Cigar::Match(60)],
            xs: Some(b'+'),
            nh: 1,
            cb: "AAAA-1",
            ub: "FFFF",
        },
        // ── Boundary reads (non-spliced, overlap exon-intron boundaries) ─
        ReadSpec {
            qname: "read13",
            pos: 1190,
            cigar: vec![Cigar::Match(100)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "LLLL",
        }, // 5' boundary intron 1
        ReadSpec {
            qname: "read14",
            pos: 1195,
            cigar: vec![Cigar::Match(100)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "LLLL",
        }, // 5' boundary intron 1
        ReadSpec {
            qname: "read15",
            pos: 1450,
            cigar: vec![Cigar::Match(100)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "MMMM",
        }, // 3' boundary intron 1
        ReadSpec {
            qname: "read16",
            pos: 1455,
            cigar: vec![Cigar::Match(100)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "MMMM",
        }, // 3' boundary intron 1
        // ── Junction 2 reads (intron at 0-based 1700..1999) ─────────
        ReadSpec {
            qname: "read17",
            pos: 1645,
            cigar: vec![Cigar::Match(55), Cigar::RefSkip(299), Cigar::Match(45)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "IIII",
        },
        ReadSpec {
            qname: "read18",
            pos: 1648,
            cigar: vec![Cigar::Match(52), Cigar::RefSkip(299), Cigar::Match(48)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "IIII",
        },
        ReadSpec {
            qname: "read19",
            pos: 1650,
            cigar: vec![Cigar::Match(50), Cigar::RefSkip(299), Cigar::Match(50)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "JJJJ",
        },
        ReadSpec {
            qname: "read20",
            pos: 1652,
            cigar: vec![Cigar::Match(48), Cigar::RefSkip(299), Cigar::Match(52)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "JJJJ",
        },
        ReadSpec {
            qname: "read21",
            pos: 1655,
            cigar: vec![Cigar::Match(45), Cigar::RefSkip(299), Cigar::Match(55)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "KKKK",
        },
        // ── More boundary reads ─────────────────────────────────────
        ReadSpec {
            qname: "read22",
            pos: 1690,
            cigar: vec![Cigar::Match(100)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "NNNN",
        }, // 5' boundary intron 2
        ReadSpec {
            qname: "read23",
            pos: 1950,
            cigar: vec![Cigar::Match(100)],
            xs: Some(b'+'),
            nh: 1,
            cb: "CCCC-1",
            ub: "OOOO",
        }, // 3' boundary intron 2
        // ── Multi-mapped read (filtered by default NH ≤ 1) ──────────
        ReadSpec {
            qname: "read_multi",
            pos: 5000,
            cigar: vec![Cigar::Match(100)],
            xs: Some(b'+'),
            nh: 2,
            cb: "AAAA-1",
            ub: "PPPP",
        },
    ];

    for rd in &reads {
        let mut rec = Record::new();
        let cigar_str = CigarString(rd.cigar.clone());

        // Query length = sum of query-consuming CIGAR ops (M, I, S, =, X)
        let qlen: u32 = rd
            .cigar
            .iter()
            .map(|c| match c {
                Cigar::Match(l)
                | Cigar::Ins(l)
                | Cigar::SoftClip(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l) => *l,
                _ => 0,
            })
            .sum();

        let seq = vec![b'A'; qlen as usize];
        let qual = vec![40u8; qlen as usize];

        rec.set(rd.qname.as_bytes(), Some(&cigar_str), &seq, &qual);
        rec.set_flags(0); // mapped, forward strand, unpaired
        rec.set_tid(0); // chr1
        rec.set_pos(rd.pos);
        rec.set_mapq(255);
        rec.push_aux(b"NH", Aux::U8(rd.nh)).unwrap();
        if let Some(xs) = rd.xs {
            rec.push_aux(b"XS", Aux::Char(xs)).unwrap();
        }
        rec.push_aux(b"CB", Aux::String(rd.cb)).unwrap();
        rec.push_aux(b"UB", Aux::String(rd.ub)).unwrap();

        writer.write(&rec).unwrap();
    }

    drop(writer);

    // Build BAI index
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 1).unwrap();

    println!("Created {bam_path} and {bam_path}.bai");
    println!("  26 records on chr1 (25 unique names + 1 duplicate pair)");
    println!("  Junction 1: chr1:1201-1499 (13 unique reads: 10 +strand, 3 -strand)");
    println!("  Junction 2: chr1:1701-1999 (5 reads, all +strand)");
    println!("  Boundary reads: 6 (2 per intron-end × 2 introns + 2 extra)");
    println!("  Multi-mapped (NH=2, filtered): 1");
    println!("  UMI dedup in single-cell mode:");
    println!("    AAAA-1 junction1: 10 reads → 6 unique UMIs");
    println!("    BBBB-1 junction1:  3 reads → 2 unique UMIs");
    println!("    CCCC-1 junction2:  5 reads → 3 unique UMIs");

    // --- Generate a temporary reference FASTA for CRAM encoding ---
    let tmp_dir = tempfile::tempdir().expect("Failed to create temp dir");
    let ref_path = tmp_dir.path().join("reference.fa");
    let ref_index_path = tmp_dir.path().join("reference.fa.fai");
    {
        let mut fa_file =
            std::fs::File::create(&ref_path).expect("Failed to create reference FASTA");
        // chr1: 10000 bp of 'N' (same length as the SQ header)
        writeln!(fa_file, ">chr1").unwrap();
        let seq = "N".repeat(10000);
        // Write in 80-char lines (standard FASTA)
        for chunk in seq.as_bytes().chunks(80) {
            fa_file.write_all(chunk).unwrap();
            fa_file.write_all(b"\n").unwrap();
        }
    }
    // Write a simple .fai index: name\tlength\toffset\tlinebases\tlinewidth
    {
        let mut fai_file = std::fs::File::create(&ref_index_path).expect("Failed to create .fai");
        // offset=6 (">chr1\n" = 6 bytes), linebases=80, linewidth=81 (80 + newline)
        writeln!(fai_file, "chr1\t10000\t6\t80\t81").unwrap();
    }

    // --- Generate CRAM from BAM ---
    let cram_path = format!("{}/example.cram", out_dir);
    {
        let mut bam_reader =
            bam::Reader::from_path(&bam_path).expect("Failed to open BAM for CRAM conversion");
        let header = bam::Header::from_template(bam_reader.header());
        let mut cram_writer = Writer::from_path(&cram_path, &header, bam::Format::Cram)
            .expect("Failed to create CRAM writer");
        cram_writer
            .set_reference(ref_path.to_str().unwrap())
            .expect("Failed to set CRAM reference");

        let mut record = Record::new();
        while let Some(result) = bam_reader.read(&mut record) {
            result.unwrap();
            cram_writer
                .write(&record)
                .expect("Failed to write CRAM record");
        }
    }

    // Build CRAI index
    bam::index::build(&cram_path, None, bam::index::Type::Csi(14), 1).unwrap();

    println!("Created {cram_path} and {cram_path}.crai");
}
