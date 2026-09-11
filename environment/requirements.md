# Environment & versions

Analyses ran under **R 4.5.3** (respirometry; other R modules on R ≥ 4.1) and **Python 3**
(`tables/build_table3.py`). Versions as reported in the manuscript Methods.

## Command-line tools
- miRTrace v1.0.1 · RepeatMasker v4.1.4 · **ShortStack v4.0**
- **RNAhybrid v2.1.2** · **TargetScan 7.0** · seqkit
- HISAT2 + StringTie (Pertea et al. 2016 protocol — see `00_upstream/REFBASED_TRANSCRIPTOME_ASSEMBLY.md`)
- DIAMOND v2.1.8 · eggNOG-mapper v2 · STRING v12.0
- CO2SYS v2.0 (carbonate chemistry)

## R / Bioconductor packages
- **DESeq2 v1.38.3**, **topGO v2.62**, GenomicRanges v1.62.1, Biostrings v2.78, rvest v1.0.5, edgeR
- tidyverse (readr, dplyr, tidyr, purrr, stringr, ggplot2), readxl, lubridate, scales
- car, emmeans, rstatix, effectsize  (respirometry two-way ANOVA)
- GOSemSim + org.Hs.eg.db  (topGO semantic similarity)
- ggraph, tidygraph, igraph, visNetwork  (STRING network)

## Reproducing a lockfile
This repository does not ship a binary `renv.lock` (it must be resolved against your own library).
To pin versions on your machine, from the repo root:

```r
install.packages("renv")
renv::init()      # snapshots the packages the scripts load -> renv.lock
```

Commit the resulting `renv.lock`. For a plain, human-readable record instead, run:

```bash
Rscript environment/sessionInfo.R   # writes environment/sessionInfo.txt
```
