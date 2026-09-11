#!/usr/bin/env Rscript
# =============================================================================
# FIG_S1_redesign.R
# Rediseño de la Figura S1 (Marine Biotechnology, Haliotis rufescens miRNAs)
#
# Reemplaza el heatmap z-score de una sola condicion (pH 8.0) por 3 paneles
# que muestran AMBOS regimenes (pH 8.0 y pH 7.6) y la particion
# intergenico / intragenico de los loci de miRNA DE del desarrollo.
#
#   Panel A  Composicion (proporcion) intergenico/intragenico por direccion
#            del desarrollo y por pH  -> reproduce los conteos del texto.
#   Panel B  MA plot: log2FC vs log10(CPM medio de las 6 librerias del
#            regimen), coloreado por clase, por pH; nube gris = no significativos.
#   Panel C  Concordancia log2FC_C (pH 8.0) vs log2FC_D (pH 7.6) con linea
#            identidad -> robustez del programa de desarrollo frente al pH.
#
# Contrastes (formato ancho de ALL_117_miRNAs_4contrasts.tsv):
#   CONTRAST_C = 24 vs 110 hpf a pH 8.0   (lfc_C, padj_C) -> librerias pH=="Control"
#   CONTRAST_D = 24 vs 110 hpf a pH 7.6   (lfc_D, padj_D) -> librerias pH=="Low"
# Convencion de signo (S-dev): lfc > 0 = over-expressed a 24 hpf;
#                              lfc < 0 = over-expressed a 110 hpf.
# Significancia por contraste: padj < 0.05 & |log2FC| > 1.
# Clasificacion: el unico locus ambiguo ("Intergenic/Intragenic") cuenta
#                como Intergenico (para cuadrar 20/18 con el texto).
# CPM = counts / colSums(matriz de 117 miRNAs) * 1e6 (edgeR::cpm sobre esta matriz);
#       eje x del MA = log10 de la media de esas CPM en las 6 librerias del regimen.
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
})
options(stringsAsFactors = FALSE, readr.show_col_types = FALSE, bitmapType = "cairo")

# ---- rutas (sobreescribibles por variables de entorno) ----------------------
if (file.exists("config.R")) source("config.R") else if (file.exists("../config.R")) source("../config.R")
if (!exists("DATA_DIR"))    DATA_DIR    <- if (dir.exists("data")) "data" else "../data"
if (!exists("RESULTS_DIR")) RESULTS_DIR <- file.path(dirname(DATA_DIR), "outputs")
infile     <- Sys.getenv("ALL117", file.path(DATA_DIR, "ALL_117_miRNAs_4contrasts.tsv"))
counts_rds <- Sys.getenv("COUNTS", file.path(DATA_DIR, "IDENTICAL_SEQUENCES_MERGED_COUNT.rds"))
meta_file  <- Sys.getenv("META",   file.path(DATA_DIR, "METADATA_MICRORNAS.tsv"))
outdir     <- Sys.getenv("OUTDIR", RESULTS_DIR)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- datos DE ---------------------------------------------------------------
dat <- read_tsv(infile) %>%
  mutate(class = factor(ifelse(biotype_best_rank == "Intragenic",
                               "Intragenic", "Intergenic"),
                        levels = c("Intergenic", "Intragenic")))

sig <- function(lfc, padj) !is.na(padj) & padj < 0.05 & abs(lfc) > 1

dir_levels <- c("Over-expressed\nat 24 hpf", "Over-expressed\nat 110 hpf")
ph_levels  <- c("pH 8.0", "pH 7.6")
class_cols <- structure(c("#FFC107", "#2196F3"),
                        names = c("Intergenic", "Intragenic"))

# ---- CPM media por regimen (6 librerias c/u) --------------------------------
counts <- readRDS(counts_rds)                 # 117 miRNAs (MajorRNA) x 12 librerias
meta   <- read.delim(meta_file, check.names = FALSE)
stopifnot(all(meta$LIBRARY_ID %in% colnames(counts)))
cpm <- sweep(counts, 2, colSums(counts), "/") * 1e6
libs_pH8  <- meta$LIBRARY_ID[meta$pH == "Control"]   # pH 8.0  (CONTRAST_C)
libs_pH76 <- meta$LIBRARY_ID[meta$pH == "Low"]       # pH 7.6  (CONTRAST_D)
cpm_df <- tibble(MajorRNA     = rownames(counts),
                 meanCPM_pH8  = rowMeans(cpm[, libs_pH8,  drop = FALSE]),
                 meanCPM_pH76 = rowMeans(cpm[, libs_pH76, drop = FALSE]))
dat <- left_join(dat, cpm_df, by = "MajorRNA")
cat(sprintf("CPM unida a %d/%d loci (%d 8.0-libs, %d 7.6-libs)\n",
            sum(!is.na(dat$meanCPM_pH8)), nrow(dat), length(libs_pH8), length(libs_pH76)))

# =============================================================================
# Tema comun
# =============================================================================
theme_pub <- theme_bw(base_size = 11, base_family = "GillSans") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text       = element_text(face = "bold", size = 11),
        plot.title       = element_text(face = "bold", size = 11.5),
        plot.subtitle    = element_text(size = 9, colour = "grey30"),
        plot.tag         = element_text(face = "bold", size = 15),
        legend.position  = "right",
        axis.title       = element_text(size = 10))

# =============================================================================
# Panel A -- composicion (proporcion) por direccion y pH, con conteos
# =============================================================================
devA <- bind_rows(
  dat %>% filter(sig(lfc_C, padj_C)) %>% transmute(class, pH = "pH 8.0", up = lfc_C > 0),
  dat %>% filter(sig(lfc_D, padj_D)) %>% transmute(class, pH = "pH 7.6", up = lfc_D > 0)
) %>%
  mutate(pH        = factor(pH, levels = ph_levels),
         direction = factor(ifelse(up, dir_levels[1], dir_levels[2]), levels = dir_levels))

cntA <- devA %>% count(pH, direction, class, name = "n") %>%
  group_by(pH, direction) %>% mutate(tot = sum(n)) %>% ungroup()
totA <- cntA %>% distinct(pH, direction, tot)

# verificacion
cat("\n== Conteos (padj<0.05 & |log2FC|>1) ==\n")
print(as.data.frame(cntA %>% select(pH, direction, class, n) %>%
        pivot_wider(names_from = class, values_from = n, values_fill = 0)))

pA <- ggplot(cntA, aes(direction, n, fill = class)) +
  geom_col(position = "fill", width = 0.72, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_fill(vjust = 0.5),
            colour = "grey10", fontface = "bold", size = 3.6) +
  geom_text(data = totA, aes(x = direction, y = 1.04, label = paste0("n = ", tot)),
            inherit.aes = FALSE, size = 3.1, vjust = 0, colour = "grey20") +
  facet_wrap(~pH) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.10))) +
  scale_fill_manual(values = class_cols) +
  labs(tag = "A", x = NULL, y = "Proportion of DE loci", fill = "Locus class",
       title = "Intergenic vs intragenic composition of developmentally DE miRNA loci",
       subtitle = "Segment labels = number of loci; the intergenic (amber) share falls from the 24-hpf to the 110-hpf loci, more so at pH 7.6") +
  theme_pub

# =============================================================================
# Panel B -- MA plot: x = log10(CPM media del regimen); gris = no significativos
# =============================================================================
ma_all <- bind_rows(
  dat %>% transmute(mirna, class, pH = "pH 8.0", lfc = lfc_C, padj = padj_C, cpm = meanCPM_pH8),
  dat %>% transmute(mirna, class, pH = "pH 7.6", lfc = lfc_D, padj = padj_D, cpm = meanCPM_pH76)
) %>%
  filter(!is.na(lfc), !is.na(cpm), cpm > 0) %>%
  mutate(pH = factor(pH, levels = ph_levels),
         signif = !is.na(padj) & padj < 0.05 & abs(lfc) > 1)

pB <- ggplot(ma_all, aes(log10(cpm), lfc)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey45", linewidth = 0.35) +
  geom_point(data = subset(ma_all, !signif), colour = "grey80", size = 1.5, alpha = 0.6) +
  geom_point(data = subset(ma_all,  signif), aes(colour = class), size = 1.9, alpha = 0.9) +
  facet_wrap(~pH) +
  scale_colour_manual(values = class_cols) +
  labs(tag = "B",
       x = expression(log[10]~"mean CPM (6 libraries per pH regime)"),
       y = expression(log[2]~"fold change"),
       colour = "Locus class",
       title = "Expression-effect (MA) structure of DE loci, by pH regime",
       subtitle = "Grey = not significant. Positive = over-expressed at 24 hpf; negative = over-expressed at 110 hpf. Dashed: |log2FC| = 1") +
  theme_pub

# =============================================================================
# Panel C -- concordancia lfc_C (pH 8.0) vs lfc_D (pH 7.6)
# =============================================================================
conc <- dat %>%
  filter(sig(lfc_C, padj_C) | sig(lfc_D, padj_D)) %>%
  transmute(mirna, class, lfc_C, lfc_D)
r_sub <- cor(conc$lfc_C, conc$lfc_D, use = "complete.obs")
r_all <- with(dat, cor(lfc_C, lfc_D, use = "complete.obs"))
cat(sprintf("Pearson r (union DE, n=%d) = %.3f | (todos 117) = %.3f\n",
            nrow(conc), r_sub, r_all))

pC <- ggplot(conc, aes(lfc_C, lfc_D, colour = class)) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey35", linewidth = 0.45) +
  geom_point(size = 1.9, alpha = 0.85) +
  scale_colour_manual(values = class_cols) +
  coord_equal() +
  labs(tag = "C",
       x = expression(log[2]*"FC at pH 8.0  (CONTRAST_C)"),
       y = expression(log[2]*"FC at pH 7.6  (CONTRAST_D)"),
       colour = "Locus class",
       title = "Concordance of the developmental program across pH",
       subtitle = sprintf("Loci DE in either regime (n = %d); identity line; Pearson r = %.3f (all 117 loci: %.3f)",
                          nrow(conc), r_sub, r_all)) +
  theme_pub

# =============================================================================
# Guardado (PDF vectorial via cairo_pdf + PNG 300 dpi) y figura combinada
# =============================================================================
save_one <- function(p, base, w, h, res = 300) {
  ggsave(file.path(outdir, paste0(base, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  png(file.path(outdir, paste0(base, ".png")), width = w, height = h,
      units = "in", res = res, type = "cairo"); print(p); invisible(dev.off())
}
save_one(pA, "FIG_S1_panelA_composition", 7.6, 3.9)
save_one(pB, "FIG_S1_panelB_MA",          7.6, 3.9)
save_one(pC, "FIG_S1_panelC_concordance", 5.4, 5.2)

# combinada: apilado vertical A / B / C con grid base (sin patchwork)
save_stack <- function(plots, heights, file, w, h, res = 300) {
  if (grepl("\\.png$", file)) png(file, width = w, height = h, units = "in", res = res, type = "cairo")
  else cairo_pdf(file, width = w, height = h)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(length(plots), 1, heights = grid::unit(heights, "null"))))
  for (i in seq_along(plots)) {
    grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
    print(plots[[i]], newpage = FALSE)
    grid::popViewport()
  }
  grid::popViewport(); invisible(dev.off())
}
save_stack(list(pA, pB, pC), c(3.9, 3.9, 5.0),
           file.path(outdir, "FIG_S1_combined.pdf"), 8.6, 12.2)
save_stack(list(pA, pB, pC), c(3.9, 3.9, 5.0),
           file.path(outdir, "FIG_S1_combined.png"), 8.6, 12.2)

cat("\nListo. Archivos en:", normalizePath(outdir), "\n")
