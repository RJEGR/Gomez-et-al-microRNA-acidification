# microRNA-abalone-acidification

Reproducibility repository for the manuscript:

> **Stage-specific microRNA turnover and the respiratory cost of early development in red abalone (*Haliotis rufescens*) larvae under ocean acidification**
> R. Gómez-Reyes et al. *Marine Biotechnology* (under revision).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

This repository contains the code and the analysis workflow needed to reproduce every table and figure in the manuscript, from the archived intermediate data objects through to the final figures. It is a curated, publication-focused release derived from the working repository [`RJEGR/Small-RNASeq-data-analysis`](https://github.com/RJEGR/Small-RNASeq-data-analysis). Public repository: https://github.com/RJEGR/microRNA-abalone-acidification

---

## 1. Data availability

| Data | Location |
|---|---|
| Raw small RNA sequencing reads (12 libraries) | NCBI SRA, BioProject **PRJNA1354232** |
| Reference transcriptome libraries (target-expression check) | NCBI SRA, BioProject **PRJNA488641** (Masonbrink et al. 2019) |
| Nuclear genome | **GCA_023055435.1** |
| Mitochondrial genome | **JALGQA010000616.1** |
| Intermediate data objects, count matrices, target tables, water-chemistry and respirometry measurements | Archived with this release on Zenodo (DOI above) |

The Zenodo archive bundles the ShortStack count matrix, the annotated miRNA sequences with their genomic coordinates (FASTA + GFF3), the pre- and post-intersection target-prediction tables, the DESeq2 result tables, and the raw water-chemistry and respirometry measurements.

---

## 2. Repository layout

```
microRNA-abalone-acidification/
├── README.md
├── LICENSE                     # MIT
├── config.R                    # central paths & parameters (single source of truth)
├── .gitignore
├── R/                          # shared helper functions
│   ├── FUNCTIONS.R
│   └── FUNCTIONS_pub.R
├── 00_upstream/                # QC + ShortStack + reference-transcriptome assembly (documented; run on HPC)
├── 01_loci/                    # genomic location DB of small-RNA loci
├── 02_targets/                 # RNAhybrid ∩ TargetScan target prediction & filtering
├── 03_diffexp/                 # DESeq2 across the four contrasts
├── 04_annotation_go/           # relational annotation DB + topGO enrichment
├── 05_respirometry/            # two-way ANOVA of oxygen consumption
├── 06_water_chemistry/         # CO2SYS carbonate-chemistry statistics (see README inside)
├── figures/                    # scripts that render Figures 1–3 (and panels)
├── tables/                     # scripts that build Tables 1–3
└── data/                       # inputs (populated from the Zenodo archive — see README inside)
```

---

## 3. What each manuscript output maps to

| Output | Script(s) | Key inputs |
|---|---|---|
| **Table 1** — respiratory metabolism | `05_respirometry/respirometry_two_way_anova.R` | `Abulon_{24,48,60,110}_{8,7.6}_Oxygen.xlsx`, `MTD_RESPIROMETRY.csv` |
| **Table 2** — sequencing information per sample | `figures/ALIGNMENT_DETAILS_EXPLORATORY.R` + QC (`00_upstream/`) | `alignment_details.tsv`, sample metadata |
| **Table 3** — DE miRNAs × 4 contrasts + targets + themes | `tables/build_table3.py` | `SEQUENCES_MERGED_DESEQ_RES.tsv`, `LONGER_RELATIONAL_DB.tsv`, `Biological_themes.tsv` |
| **Tables S1–S2** — seawater chemistry | `06_water_chemistry/chemistries_TA_DIC.R` | `quimicas.csv` (CO2SYS v2.0 values), `CO2_sys_inputs.tsv` |
| **Figure 1** — miRNA identification (length, precision, MFE, miR–miR\*, chromosomal map) | `figures/FIGURE_1.R`, `figures/LOCI_PLOT.R`, `figures/ALIGNMENT_DETAILS_EXPLORATORY.R`, `01_loci/*` | ShortStack `Results.txt`, `mir.fasta`, `Results.gff3`; `MIRGENEDB_2.1.tsv`; `RNA_LOCATION_MIR_DB.rds` |
| **Figure 2B** — PCA of VST counts | `figures/PCA.R`, `figures/PCA_variance_stabilized.R` | `IDENTICAL_SEQUENCES_MERGED_COUNT.rds` |
| **Figure 2C** — intergenic/intragenic proportion | `figures/INTRAGENIC_ANALYSIS.R` | `RNA_LOCATION_DB.tsv` |
| **Figure 2D** — GO enrichment | `04_annotation_go/TopGO_by_CONTRAST.R` (+ `_viz`) | `Gene_ontologies_DB.tsv`, `LONGER_RELATIONAL_DB.tsv` |
| **Figure 3** — miRNA degree & density | `figures/Figure_3.R`, `figures/OA_MIR_degree_analysis.R` | `LONGER_RELATIONAL_DB.tsv`, `SEQUENCES_MERGED_DESEQ_RES.tsv`, `Biological_themes.tsv` |

Shared prerequisites: `03_diffexp/1_DESEQ2_RUN_AND_PREP.R` produces the DESeq2 result table used by Table 3 and Figure 3; `04_annotation_go/JOIN_RELATIONAL_DB.R` builds the annotation table (`LONGER_RELATIONAL_DB.tsv`) and `Gene_ontologies_DB.tsv`.

---

## 4. Suggested execution order

1. `00_upstream/` — read QC (miRTrace) and ShortStack annotation/quantification *(run on HPC; produces the ShortStack outputs)*
2. `01_loci/` — build the genomic-location database of loci
3. `02_targets/` — RNAhybrid and TargetScan predictions, then their intersection
4. `03_diffexp/` — DESeq2 for the four contrasts
5. `04_annotation_go/` — annotation join and GO enrichment
6. `figures/` and `tables/` — render figures and tables
7. `05_respirometry/` and `06_water_chemistry/` — physiology and environment (independent of the sequencing branch)

---

## 5. Software

Analyses were run in **R** (respirometry under R 4.5.3) and **Python 3** (`tables/build_table3.py`). Principal tools and versions as reported in the manuscript:

- Upstream: miRTrace v1.0.1, RepeatMasker v4.1.4, **ShortStack v4.0**
- Target prediction: **RNAhybrid v2.1.2**, **TargetScan 7.0**, seqkit
- Annotation/networks: eggNOG-Mapper v2, DIAMOND v2.1.8, STRING v12.0
- R/Bioconductor: **DESeq2 v1.38.3**, **topGO v2.62**, GenomicRanges v1.62.1, Biostrings v2.78, rvest v1.0.5
- Reference transcriptome: HISAT2 + StringTie (Pertea et al. 2016)
- Respirometry: car, emmeans, rstatix, effectsize
- Carbonate chemistry: CO2SYS v2.0

Pinned package versions and how to regenerate a lockfile live in [`environment/`](environment/) (`requirements.md` + `sessionInfo.R`); run `renv::init()` locally to produce a `renv.lock` for the tagged release.

---

## 6. Externally-sourced (FAIR) inputs — not redeposited

The following inputs are public/community resources or are regenerable, so they are **documented here rather than archived** with this release (findable via their own persistent sources):

- **Reference genome & gene models.** *H. rufescens* nuclear (GCA_023055435.1) and mitochondrial (JALGQA010000616.1) assemblies; gene models / GTF and the sequences used for 3′UTR extraction obtained from **Ensembl Metazoa** (https://metazoa.ensembl.org/info/data/ftp/). The 3′UTR set (`utr_rmdup.fa`) and the concatenated `multi_genome.newid.fa` are produced from these with **seqkit**.
- **STRING interactome.** `protein.links.full.v12.0.txt` and `protein.info.v12.0.txt` from **STRING v12.0** (https://string-db.org/cgi/download), processed with the R code in `04_annotation_go/` (`3_STRING-db.R`); interactions are transferred by best-hit homology to *H. sapiens* (taxon 9606) and *C. elegans* (taxon 6239).
- **Clade-specific miRNA families.** `clade-specific_miRNA_families_of_animal_clades.txt` is a secondary file generated when running **miRTrace v1.0.1** (QC / taxonomic screen); regenerate it by running miRTrace.
- **Regenerable intermediates.** Objects such as `SEQUENCES_MERGED_DDS_DESEQ2.rds`, `SRNA_FUNCTION_PREDICTED_LONG_EXPRESSED.rds` and `LOCI2TARGETDB.tsv` are produced by the scripts in this repository and are recomputed on demand rather than archived.

## 7. Provenance of the deposited inputs

**Target prediction** (see the working repo's `RAW_TUTORIAL_BKP/MIRS_FUNCTIONAL_ANNOT.md`):

- 3′UTRs (54,432 sequences) were extracted from the Ensembl gene models with seqkit (`utr_rmdup.fa`).
- **RNAhybrid v2.1.2** was run with the mature and star sequences as queries against the de-duplicated 3′UTR set, seed constrained to positions 2–8 (`-f 2,8`), `-m 20000 -n 50`, EVD parameters `-s 3utr_human`; hits retained at Benjamini–Hochberg *q* < 0.05 → `mir_vs_utr_rmdup_RNAhybrid.out.tsv` (raw output ≈ **2.6 GB** — archived gzip-compressed on Zenodo, **not** tracked in git).
- **TargetScan 7.0** (7-nt seed, positions 2–8) was applied to the same UTR set; its overlap with the RNAhybrid significant hits is `mature_star_mir_vs_mir_vs_utr_rmdup_RNAhybrid.out.psig_targetscan.out` (in `data/target_prediction/`).
- The **intersection** of both tools is the high-confidence target set → `SRNA_REGULATORY_FUNCTION_DB.tsv` (post-intersection).

**Water chemistry** (Tables S1–S2):

- Raw discrete measurements (Total Alkalinity, DIC, salinity, temperature, pressure) are in `data/water_chemistry/` (`quimicas*.tsv/csv`, `CO2_sys_inputs.tsv`).
- The carbonate-system variables in `quimicas.csv` (DIC, TA, CO₃²⁻, Ω_ara, *p*CO₂, pH) were computed with **CO2SYS v2.0** — the method reported in the manuscript — from the measured Total Alkalinity and Dissolved Inorganic Carbon.
- The curated script `06_water_chemistry/chemistries_TA_DIC.R` runs only the **downstream statistics** on those values: outlier / normality (Shapiro–Wilk) / homoscedasticity (Levene) tests, pairwise *t*-tests and ANOVA (`rstatix`), salinity normalisation following Friis et al. (2003), a pH↔carbonate regression, and correlation matrices. The original exploratory `seacarb` recomputation has been removed, so the code and the manuscript now agree on CO2SYS v2.0.

**Reference transcriptome (target-expression check).** The public *H. rufescens* RNA-seq libraries (BioProject **PRJNA488641**; Masonbrink et al. 2019) were aligned to the genome with HISAT2 and assembled with StringTie in reference-guided mode (Pertea et al. 2016), producing `gene_count_matrix.csv` / `transcript_count_matrix.csv` — used to confirm that predicted targets are expressed during larval development. The full command-line protocol is in [`00_upstream/REFBASED_TRANSCRIPTOME_ASSEMBLY.md`](00_upstream/REFBASED_TRANSCRIPTOME_ASSEMBLY.md).

---

## 8. Status

- ✅ Scripts for every manuscript table and figure are included and organized by module — including Figure S1 (`figures/FIG_S1_redesign.R`), the WGCNA-free reconstructed GO visualisation (`04_annotation_go/TopGO_by_CONTRAST_viz_RECONSTRUCTED.R`) and the STRING-network script (`04_annotation_go/3_STRING-db.R`).
- ✅ `05_respirometry/` runs out of the box (inputs bundled + repository-relative paths).
- ✅ Deposited inputs now included in `data/`: read-length profiling (`prof_by_read_length_summary.tsv`, converted from `.rds` for interoperability), water chemistry (`data/water_chemistry/`) and the TargetScan-overlap table (`data/target_prediction/`).
- ✅ Path centralisation into `config.R` → `data/` is **complete**: every module (`00_upstream`–`06` + `figures/`) resolves paths through `DATA_DIR` / `ANNOT_DIR` / `RESULTS_DIR`. No absolute per-user paths remain in any script.
- ↗ Externally-sourced inputs (§6) are not redeposited; the ≈2.6 GB raw RNAhybrid output is archived compressed on Zenodo.

---

## 9. Citation & contact

If you use this code, please cite the manuscript and the Zenodo archive (DOI above).

**Ricardo Gómez-Reyes** — rgomez41@uabc.edu.mx — Universidad Autónoma de Baja California (UABC)

## 10. License

Released under the MIT License — see [`LICENSE`](LICENSE). © 2026 Ricardo Gómez-Reyes.
