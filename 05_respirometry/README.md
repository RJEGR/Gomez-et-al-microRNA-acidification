# 05_respirometry — Table 1 (respiratory metabolism)

**Self-contained, runnable module.** Inputs live in `data_input/`; the script
uses a repository-relative path (no external `~/Documents/...` path).

## Run
From the repository root:
```
Rscript 05_respirometry/respirometry_two_way_anova.R
```
or from this folder:
```
cd 05_respirometry && Rscript respirometry_two_way_anova.R
```

## Inputs (`data_input/`)
- `Abulon_<hpf>_<pH>_Oxygen.xlsx` — SDR SensorDish oxygen traces, 24/48/60/108→110 hpf × pH 8.0 / 7.6 (the files carry `108`; the script recodes to 110 hpf, matching the manuscript).
- `MTD_RESPIROMETRY.csv` — chamber metadata.
- `conteos_abs_respirometrias.csv` — larval counts used for normalisation.

## Outputs (`out_respirometry/`)
`02_summary_by_cell.tsv` (Table 1 cells), `04_anova_table.tsv` (two-way ANOVA — reviewer point 3), `07_cumulative_oxygen.tsv`, `10_relative_expenditure.tsv`, `Fig_RM_stage_pH.png`, `08_results_paragraph.md`, plus assumption/robustness/post-hoc tables.

## R packages
`readxl, readr, dplyr, tidyr, purrr, stringr, ggplot2, lubridate, car, emmeans, rstatix, effectsize` (R ≥ 4.1; manuscript used R 4.5.3).

> Note: `data_input/` is git-ignored (raw data are distributed via the Zenodo archive) but is included in the delivered bundle so the module runs out of the box.
