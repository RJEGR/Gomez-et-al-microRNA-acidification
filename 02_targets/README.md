# 02_targets — miRNA→mRNA target prediction & filtering

Builds the target database (RNAhybrid ∩ TargetScan → `SRNA_REGULATORY_FUNCTION_DB.tsv` /
`SEQUENCES_MERGED_MIRNA2TARGET_EXP_DB.tsv`) used by **Table 3** and **Figure 3**.

## Scripts
- `1_…`–`4_SRNA_REGULATORY_FUNCTION_DB_PREP.R` — build and filter the target DB.
- `5_SRNA_REGULATORY_FUNCTION_HOST_MIRS.R` — host-gene (intragenic) miRNAs.
- `CREATE_MIR_REGULATORY_DB.R` — assemble the regulatory DB.
- `RNAHYBRID_EXPLORATORY.R`, `TARGETSCAN_EXPLORATORY.R` — the two predictors (see top-level README §7 for how the raw outputs are produced).

## Paths
Resolved through `config.R` (`DATA_DIR` / `ANNOT_DIR` / `RESULTS_DIR`). No absolute per-user paths remain.

## Inputs
**Provenance/build** scripts. They consume: the 3′UTR set and predictor outputs (FAIR-external, see
README §6–§7 — `utr_rmdup.fa` from ENSEMBL + seqkit; the ≈2.6 GB raw RNAhybrid output; the TargetScan
overlap in `data/target_prediction/`), the StringTie quantification (`data/`), and eggNOG annotation
(`data/annotation/`).

> The high-confidence post-intersection database (`SRNA_REGULATORY_FUNCTION_DB.tsv`) is provided as a
> deposited product and consumed directly downstream, so Table 3 / Figure 3 reproduce without re-running
> the predictors.
