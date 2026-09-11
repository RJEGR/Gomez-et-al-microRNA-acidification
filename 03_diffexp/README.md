# 03_diffexp — differential expression (DESeq2)

Produces the miRNA differential-expression results that feed **Table 3** and **Figure 3**,
plus the processed count object used for the **Figure 2B** PCA.

## Run
```
Rscript 03_diffexp/1_DESEQ2_RUN_AND_PREP.R      # from the repository root
```
Paths are read from `config.R` (`DATA_DIR`, `RESULTS_DIR`); the original absolute
`~/Documents/MIRNA_HALIOTIS/…` paths have been removed.

## Inputs (`data/`)
| File | Description | In bundle? |
|---|---|---|
| `Counts.txt` | ShortStack per-sample count matrix | large (~5 MB) — copy from your `MicroRNA_target_validation/Shortstacks_outputs/` or the Zenodo archive |
| `RNA_LOCATION_DB.tsv` | locus database; used to keep only `SRNAtype == "miR"` and map `MajorRNA` | large (~15 MB) — copy from `MicroRNA_target_validation/` or Zenodo |
| `METADATA.tsv` | library design / `colData` with the `CONTRAST_*` factor columns | ✅ bundled |
| `SRNA2MIRGENEDB.tsv` | `MajorRNA` → MirGeneDB id map (arm/family) | ✅ bundled |

## Outputs (`outputs/`)
- `IDENTICAL_SEQUENCES_MERGED_COUNT.rds` — 117 × 12 miRNA matrix (identical MajorRNA sequences merged).
- `SEQUENCES_MERGED_DESEQ_RES.tsv` — DE results across the four contrasts.
- `SEQUENCES_MERGED_DDS_DESEQ2.rds` — the DESeqDataSet object.

## Method
Four **independent pairwise** DESeq2 models (`design = ~ Design`, Wald test), reference level = control /
early stage; *p*-values adjusted with Benjamini–Hochberg. Differentially expressed = `padj < 0.05` **and**
`|log2FC| ≥ 1`. Contrasts: **A** low vs control pH @ 24 hpf · **B** low vs control pH @ 110 hpf ·
**C** 110 vs 24 hpf @ pH 8.0 · **D** 110 vs 24 hpf @ pH 7.6.

> Reproducibility note (kept from the source script): identical-`MajorRNA` miRNAs are summed before
> testing, and the `Name` recode is verified to give an exact **468/468** match against the published
> `SEQUENCES_MERGED_DESEQ_RES.tsv` (vs 196/468 without it).
