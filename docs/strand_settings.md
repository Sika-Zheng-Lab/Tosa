# Strand Settings — Detailed Reference

This document provides a comprehensive reference for strand-related settings across RNA-seq analysis tools.

The strand setting tables and library kit information below are adapted from the excellent
[Strand Settings](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/) page
of the [RNA-seq Bioinformatics](https://rnabio.org/) course, developed and maintained by the
[Griffith Lab](http://www.griffithlab.org/) at Washington University.
We are grateful for their effort in compiling and sharing this valuable resource with the community.

---

## Tosa strand options

The `-s`/`--strand` option specifies the strand specificity of your RNA-seq library.

| Tosa `-s` option | Library type | Description |
|---|---|---|
| *(omit)* | Unstranded | No strand information used |
| `XS` | Any | Use the `XS` auxiliary tag set by the aligner (e.g., HISAT2, STAR) |
| `RF` | First-strand (dUTP) | Read 1 maps to the reverse complement of the transcript |
| `FR` | Second-strand (ligation) | Read 1 maps to the transcript strand |

## How Tosa determines strand from RF / FR

The table below shows how Tosa infers the **gene strand** (+ or −) from the mapping orientation of each read. For example, if your library is RF and Read 1 maps to the reverse strand, Tosa assigns the junction to a **+** strand gene.

| Mode | Read 1 reverse | Read 1 forward | Read 2 reverse | Read 2 forward |
|---|---|---|---|---|
| **RF** (first-strand) | **+** | **−** | **−** | **+** |
| **FR** (second-strand) | **−** | **+** | **+** | **−** |

For single-end reads, Tosa uses the same rule as Read 1 in the corresponding mode.

## Correspondence with other tools

The table below shows equivalent strand settings across commonly used RNA-seq tools.

| Tool | First-strand (RF) | Second-strand (FR) | Unstranded |
|---|---|---|---|
| **Tosa** | `-s RF` | `-s FR` | *(omit `-s`)* or `-s XS` |
| **RegTools** | `-s RF` | `-s FR` | `-s XS` |
| TopHat | `--library-type fr-firststrand` | `--library-type fr-secondstrand` | `--library-type fr-unstranded` |
| HISAT2 | `--rna-strandness R` / `RF` | `--rna-strandness F` / `FR` | *(none)* |
| HTSeq | `--stranded reverse` | `--stranded yes` | `--stranded no` |
| featureCounts | `-s 2` | `-s 1` | `-s 0` |
| Salmon | `--libType ISR` | `--libType ISF` | `--libType IU` |
| Kallisto | `--rf-stranded` | `--fr-stranded` | *(none)* |
| StringTie | `--rf` | `--fr` | *(none)* |
| Picard | `SECOND_READ_TRANSCRIPTION_STRAND` | `FIRST_READ_TRANSCRIPTION_STRAND` | `NONE` |
| Trinity | `--SS_lib_type RF` | `--SS_lib_type FR` | *(none)* |
| RSEM | `--forward-prob 0` | `--forward-prob 1` | `--forward-prob 0.5` |
| IGV (5'→3' read orientation) | F2R1 | F1R2 | F2R1 or F1R2 |

## Common library kits

| Strand type | Example kits / methods |
|---|---|
| **First-strand (RF)** | dUTP, NSR, NNSR, Illumina TruSeq Stranded Total RNA, NEBNext Ultra II Directional, Watchmaker RNA Library Prep Kit with Polaris Depletion |
| **Second-strand (FR)** | Ligation, Standard SOLiD, NuGEN Encore, 10X Genomics 5' scRNA |
| **Unstranded** | Standard Illumina (non-stranded), NuGEN OvationV2, SMARTer Universal Low Input RNA Kit (TaKara), GDC normalized TCGA data |

> [!WARNING]
> The list above assumes kits are used as specified by their manufacturer. Sequencing providers may substitute adapters or modify protocols, which can flip the effective strandedness (e.g., substituting IDT xGen UDI-UMI adapters for NEB hairpin adapters changes RF to FR). Always confirm your data's strandedness empirically.

## Determining strandedness of your data

If you are unsure about your library's strand specificity, two practical approaches are:

1. **[check_strandedness](https://github.com/betsig/how_are_we_stranded_here)** — An automated tool that infers strandedness from your raw FASTQ data and a reference annotation.

2. **Visual inspection in [IGV](https://igv.org/)** — Color alignments by *First-of-pair strand* and check whether Read 1 / Read 2 orientations match the expected pattern for your library type.

## References

- Griffith Lab, *RNA-seq Bioinformatics — Strand Settings*: <https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/>
- Signal et al., *how_are_we_stranded_here*: <https://github.com/betsig/how_are_we_stranded_here>
