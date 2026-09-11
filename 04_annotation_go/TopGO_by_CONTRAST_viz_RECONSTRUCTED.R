# =============================================================================
# Adapted visualization for the RECONSTRUCTED Boostrap table
# Input: Boostrap_topGO_intra_inter_reconstructed.tsv
# (rebuilt in Python because R/topGO could not be run; see summary .md)
# The reconstructed table already carries f1 (stage), f2 (pH), biotype_best_rank,
# term, parentTerm, size (REVIGO), Significant, n_miRNAs, classicFisher, p.adj.fdr.
# So the original recode_factor()/separate() preprocessing is NOT needed.
# =============================================================================
library(tidyverse)

if (file.exists("config.R")) source("config.R") else if (file.exists("../config.R")) source("../config.R")
if (!exists("DATA_DIR"))    DATA_DIR    <- if (dir.exists("data")) "data" else "../data"
if (!exists("RESULTS_DIR")) RESULTS_DIR <- file.path(dirname(DATA_DIR), "outputs")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir <- DATA_DIR   # reconstructed TSV lives in data/
DATA <- read_tsv(file.path(dir, "Boostrap_topGO_intra_inter_reconstructed.tsv"))

col_recode <- c(Intergenic = "#FFC107", Intragenic = "#2196F3")

# ---- Original-style lollipop: BP theme (parentTerm) x enrichment ratio -------
# 'Enrichment ratio' reproduces the original figure = REVIGO size / max(size)
# within each stage x biotype facet (the bars are driven by term generality, as
# in the published figure; they are NOT a p-value).
plot_df <- DATA %>%
  mutate(f1 = factor(f1, levels = c("24 hpf", "110 hpf"))) %>%
  group_by(f1, biotype_best_rank, parentTerm) %>%
  summarise(size = max(size), .groups = "drop") %>%
  group_by(f1, biotype_best_rank) %>%
  mutate(ratio = size / max(size)) %>%
  arrange(desc(ratio), .by_group = TRUE) %>%
  slice_head(n = 15) %>%          # top themes per facet; adjust as needed
  ungroup()

p <- plot_df %>%
  ggplot(aes(y = reorder(parentTerm, ratio), x = ratio, color = biotype_best_rank)) +
  geom_segment(aes(xend = 0, yend = parentTerm), linewidth = 1.3) +
  geom_point(size = 2) +
  facet_grid(f1 ~ biotype_best_rank, scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual("", values = col_recode) +
  labs(y = "microRNA features (Biological process)", x = "Enrichment ratio (REVIGO size / max)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black"),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7))

ggsave(p, filename = "Fig_reconstructed_target_themes.png",
       path = RESULTS_DIR, width = 8, height = 7, dpi = 300)

# ---- Alternative: order/filter by the honest hypergeometric test -------------
# NOTE: with the target-gene background the over-representation is NOT significant
# (foreground ~ whole target space -> tautological background). Use counts/size
# for the descriptive figure; use classicFisher/p.adj.fdr only if you switch to a
# genome-wide background later.
# DATA %>% filter(p.adj.fdr < 0.05)   # -> currently empty
