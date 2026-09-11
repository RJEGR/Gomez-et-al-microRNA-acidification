#!/usr/bin/env Rscript
# Reproducibility helper: capture the R session used to run this repository.
# Usage (from the repository root):  Rscript environment/sessionInfo.R
# Writes environment/sessionInfo.txt with sessionInfo() + attached package versions.

pkgs <- c("DESeq2", "topGO", "GenomicRanges", "Biostrings", "rvest", "edgeR",
          "tidyverse", "readxl", "car", "emmeans", "rstatix", "effectsize",
          "GOSemSim", "ggraph", "tidygraph", "igraph")

for (p in pkgs) {
  suppressWarnings(suppressMessages(
    if (requireNamespace(p, quietly = TRUE)) require(p, character.only = TRUE)
  ))
}

outdir <- if (dir.exists("environment")) "environment" else "."
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Wrote ", file.path(outdir, "sessionInfo.txt"))
