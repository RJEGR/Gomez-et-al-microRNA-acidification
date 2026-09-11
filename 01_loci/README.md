# 01_loci — genomic-location database of small-RNA loci

Builds the locus databases used by **Fig 1D**, **Fig 2C** and the target module:
`RNA_LOCATION_DB.tsv` and `RNA_LOCATION_MIR_DB.rds` (biotype and intergenic-vs-intragenic
classification, scaffold coordinates and inter-loci distances).

## Scripts (run in order)
```
1_SRNA_LOCATION_DB_PREP.R  →  2_SRNA_LOCATION_DB_PREP.R  →  3_SRNA_LOCATION_DB_VIZ.R  →  4_SRNA_LOCATION_MIRS.R
```
`consolidate_mirna_tibble.R` is a helper that assembles the consolidated miRNA tibble
(`MIRNA_CONSOLIDATED`).

## Paths
Resolved through `config.R` (`DATA_DIR` / `ANNOT_DIR` / `RESULTS_DIR`). No absolute per-user paths remain.

## Inputs
ShortStack outputs (`Results.txt`, `Results.gff3`, `mir.fasta`) from `data/`; ENSEMBL genome annotation
and RepeatMasker output from `data/annotation/` (FAIR-external — see top-level README §6).

> These are **provenance/build** scripts. The ready-made `RNA_LOCATION_DB.tsv` and
> `RNA_LOCATION_MIR_DB.rds` are provided as deposited inputs and consumed directly by the downstream
> figure modules, so the figures reproduce without re-running this stage.
