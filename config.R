# =============================================================================
# config.R — central paths and parameters for microRNA-red-abalone-pub
# Single source of truth. Scripts should `source("config.R")` and read paths
# from here instead of hard-coding absolute paths.
#
# NOTE: every module script resolves its paths through the variables below
# (DATA_DIR / ANNOT_DIR / RESULTS_DIR); no absolute per-user paths remain.
# =============================================================================

## Repository root (works when the .Rproj / repo is the working directory)
if (requireNamespace("here", quietly = TRUE)) {
  ROOT <- here::here()
} else {
  ROOT <- normalizePath(getwd())
}

## Input / output directories (relative to the repo)
DATA_DIR    <- file.path(ROOT, "data")          # inputs from the Zenodo archive
RESULTS_DIR <- file.path(ROOT, "outputs")       # intermediate reproducible outputs
FIGS_DIR    <- file.path(ROOT, "figures", "out")
TABLES_DIR  <- file.path(ROOT, "tables", "out")

for (d in c(RESULTS_DIR, FIGS_DIR, TABLES_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## Key input files (fill in as the data/ folder is populated)
SHORTSTACK_COUNTS <- file.path(DATA_DIR, "Results.txt")                 # ShortStack count matrix
MIRNA_FASTA       <- file.path(DATA_DIR, "mir.fasta")                   # miRNA sequences
MIRNA_GFF3        <- file.path(DATA_DIR, "Results.gff3")                # genomic coordinates
MIR_COUNT_RDS     <- file.path(DATA_DIR, "IDENTICAL_SEQUENCES_MERGED_COUNT.rds")
LOCATION_DB       <- file.path(DATA_DIR, "RNA_LOCATION_DB.tsv")
TARGET_DB         <- file.path(DATA_DIR, "SRNA_REGULATORY_FUNCTION_DB.tsv")
DESEQ_RES         <- file.path(DATA_DIR, "SEQUENCES_MERGED_DESEQ_RES.tsv")
RELATIONAL_DB     <- file.path(DATA_DIR, "LONGER_RELATIONAL_DB.tsv")
GENE_COUNTS       <- file.path(DATA_DIR, "gene_count_matrix.csv")      # StringTie mRNA matrix
PROF_READLEN      <- file.path(DATA_DIR, "prof_by_read_length_summary.tsv")  # Fig 1A / Fig 2 / Table 2
WATER_CHEM_DIR    <- file.path(DATA_DIR, "water_chemistry")            # Tables S1-S2 (seacarb inputs)
TARGET_PRED_DIR   <- file.path(DATA_DIR, "target_prediction")         # RNAhybrid∩TargetScan overlap
METADATA          <- file.path(DATA_DIR, "METADATA.tsv")              # library design / colData (CONTRAST_* columns)
SRNA2MIRGENEDB    <- file.path(DATA_DIR, "SRNA2MIRGENEDB.tsv")        # MajorRNA -> MirGeneDB id map
COUNTS_RAW        <- file.path(DATA_DIR, "Counts.txt")                # ShortStack per-sample counts (DESeq2 input)
ANNOT_DIR         <- file.path(DATA_DIR, "annotation")               # eggNOG + STRING + ENSEMBL annotation (some FAIR-external, see README §6)

## Analysis parameters (as reported in the manuscript)
DE_PADJ   <- 0.05     # DESeq2 adjusted-p threshold
DE_LFC    <- 1        # |log2FC| threshold
GO_PADJ   <- 0.05     # topGO elim-KS threshold
SEED      <- 2026
