# 06_water_chemistry — carbonate chemistry (Tables S1–S2)

**Inputs** in `data/water_chemistry/`:
- `quimicas.csv` — per-sample carbonate-system table (DIC, TA, CO₃²⁻, Ara = Ω_ara, salinity, …). **These carbonate variables were computed with CO2SYS v2.0** from measured Total Alkalinity and Dissolved Inorganic Carbon.
- `quimicas_filtered.tsv`, `quimicas_v2.tsv` — filtered / alternative tabulations.
- `CO2_sys_inputs.tsv` — the TA + DIC inputs handed to CO2SYS v2.0.

**Script:** `chemistries_TA_DIC.R` (curated copy; adapted from the `OA-research` repo,
https://github.com/RJEGR/OA-research/blob/main/chemistries_TA_DIC.R). It reads
`quimicas.csv` and runs the downstream statistics only: outlier detection, Shapiro–Wilk
normality, Levene homoscedasticity, pairwise *t*-tests / Tukey / ANOVA (`rstatix`),
salinity normalisation (nDIC/nTA) after Friis et al. (2003), a pH↔carbonate regression,
and correlation matrices. **The exploratory `seacarb` block was removed** — CO2SYS v2.0
is the method of record, so code and manuscript agree.

Optional helper dependencies (used only by parts of the script):
- `stats.R` — provides an `IC()` confidence-interval helper; the script now falls back to
  an inline `IC()` if `stats.R` is absent.
- `pH_aLLdatasets_stats.rds` — pH time-series summary used by the (optional) pH-regression
  section; supply it to run that part.

Run from the repository root (`Rscript 06_water_chemistry/chemistries_TA_DIC.R`) or from
this folder. Reported configuration: pH 7.6 ± 0.122, Ω_ara 1.14 ± 0.489, *p*CO₂ 987 ± 510
µatm (acidification) vs. pH 8.0 (control).

**Outputs:** Supplementary Table S1 (acidification, pH 7.6) and Table S2 (control, pH 8.0).
