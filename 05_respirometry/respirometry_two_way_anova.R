# =============================================================================
#  respirometry_two_way_anova.R
#
#  Haliotis rufescens larvae — respiratory metabolism (RM) under ocean
#  acidification.  Reproducible pipeline:
#     raw SDR oxygen traces -> blank- and individual-corrected rates
#     -> two-way ANOVA (developmental stage x pH, both fixed)
#     -> post-hoc contrasts -> cumulative oxygen consumption.
#
#  This script replaces the exploratory `respirometries.R`.  It reproduces the
#  published cell means (24 hpf: 51.19 +/- 10.21 vs 73.30 +/- 15.61; 110 hpf:
#  19.27 +/- 2.38 vs 30.11 +/- 4.48 pmol O2 ind-1 h-1 at pH 8.0 and 7.6) but
#     * replaces the stratified one-way Welch / Kruskal tests by the two-way
#       fixed-effects model requested by the Reviewer, so that the stage x pH
#       interaction can actually be tested;
#     * uses Type II / Type III sums of squares, because the design is
#       unbalanced and Type I SS would be order-dependent;
#     * repairs a key-construction bug in the metadata (see section 5);
#     * computes cumulative oxygen demand under BOTH definitions (unweighted
#       sum of stage means vs duration-weighted trapezoidal integration), which
#       is the point the Reviewer's "53 %" comment turns on.
#
#  Outputs (written to CFG$out_dir):
#    00_metadata_key_conflicts.tsv  metadata rows whose stored ID is wrong
#    01_rates_by_chamber.tsv        per-chamber rates and all intermediates
#    02_summary_by_cell.tsv         mean/SD/SEM/n per stage x pH cell (Table 2)
#    03_assumptions.tsv             Shapiro-Wilk, Levene
#    04_anova_table.tsv             F, df, p, adjusted p, eta2   <-- Reviewer pt 3
#    05_robustness.tsv              log, HC3, rank, permutation, subset models
#    06_posthoc_emmeans.tsv         Tukey-adjusted contrasts
#    07_cumulative_oxygen.tsv       cumulative demand under both definitions
#    10_relative_expenditure.tsv    share of the cumulative total per stage
#    08_results_paragraph.md        manuscript-ready wording
#    Fig_RM_stage_pH.png, Fig_diagnostics.png
#
#  Requires R >= 4.1 and: readxl readr dplyr tidyr purrr stringr ggplot2
#                         car emmeans rstatix effectsize
# =============================================================================

## ---- 0. Session ---------------------------------------------------------- ##

rm(list = ls())
options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)
set.seed(20260817)                      # permutation / bootstrap reproducibility

suppressPackageStartupMessages({
  library(readxl);  library(readr);   library(dplyr);   library(tidyr)
  library(purrr);   library(stringr); library(ggplot2); library(lubridate)
  library(car)          # Type II / III SS, heteroscedasticity-consistent tests
  library(emmeans)      # estimated marginal means and contrasts
  library(rstatix)      # Levene, Games-Howell, Tukey outlier rule
  library(effectsize)   # eta squared
})

## ---- 1. Configuration ---------------------------------------------------- ##
# Everything the user may need to change lives here.

CFG <- list(
  data_dir    = { d <- "data_input"; if (!dir.exists(d)) d <- "05_respirometry/data_input"; d },  # repo-relative (was an absolute ~/Documents path)
  file_regex  = "^Abulon_[0-9]+_[0-9.]+_Oxygen\\.xlsx$",   # excludes *_bkp.xlsx
  mtd_file    = "MTD_RESPIROMETRY.csv",
  out_dir     = "./out_respirometry",

  header_skip = 12,               # SDR export: 12 metadata lines before header
  spot_regex  = "^[A-D][0-9]$",   # optode chambers, laid out A1..D6
  chamber_L   = 0.0017,           # 1700 uL chamber volume, in litres
  umol_to_pmol= 1e6,              # umol -> pmol

  # Stage labels: the files and the metadata carry 108, the manuscript reports
  # 110 hpf.  Set to NULL to keep the raw values.
  stage_recode = c("108" = "110"),

  # Rate estimator:
  #   "scope" = (max - min) / elapsed time between those two points
  #             (the estimator used in the original analysis; kept as default
  #              so that the published means are reproduced exactly)
  #   "slope" = |OLS slope| of O2 on time over the whole trace
  #             (sensitivity analysis: uses every point, immune to single spikes)
  rate_method   = "scope",

  drop_outliers = TRUE,   # Tukey 1.5 x IQR rule, applied WITHIN each cell
  padjust       = "BH",   # adjustment of the ANOVA term p-values
  n_perm        = 5000,   # permutations for the distribution-free F test
  n_boot        = 2000,   # bootstrap resamples for the cumulative-demand CI
  alpha         = 0.05,
  control_pH    = "8",    # reference level of pH
  contrast_pH   = "7.6"   # "acidification" treatment reported in the abstract
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat("\n== ", sprintf(...), " ==\n", sep = "")
out <- function(f) file.path(CFG$out_dir, f)

#' Format a number exactly the way the file names and the metadata do:
#' 8 -> "8", 7.6 -> "7.6", 108 -> "108".
#' `format()` must NOT be used for this: it pads to a common width across the
#' vector (c(7.6, 8) becomes "7.6" and "8.0") and silently breaks every key
#' built from it.  This bug is why the pH factor must be built with fmt_num().
fmt_num <- function(x) {
  sub("\\.?0+$", "", formatC(as.numeric(x), format = "f", digits = 3))
}

## ---- 2. Readers ---------------------------------------------------------- ##

#' Read one SDR SensorDish oxygen export.
#'
#' The stage (hpf) and the nominal pH are encoded in the file name
#' (`Abulon_<hpf>_<pH>_Oxygen.xlsx`).  Parsing them with an explicit regex,
#' instead of `strsplit(...)[[2]]`, makes a misnamed file fail loudly.
#'
#' @param file path to one .xlsx export
#' @return tibble, one row per time point, one column per chamber
read_sdr <- function(file) {
  meta <- str_match(basename(file), "^Abulon_([0-9]+)_([0-9.]+)_Oxygen")
  if (anyNA(meta)) stop("Unparsable file name: ", basename(file))

  # `.name_repair = "unique_quiet"` names the trailing empty columns that the
  # SDR export leaves behind; the select() then drops them.
  readxl::read_excel(file, col_names = TRUE, skip = CFG$header_skip,
                     .name_repair = "unique_quiet") %>%
    select(`Date/Time`, `Time/Min.`, matches(CFG$spot_regex)) %>%
    rename(date_raw = `Date/Time`, min = `Time/Min.`) %>%
    mutate(date = suppressWarnings(lubridate::dmy_hms(date_raw)),
           hpf  = as.numeric(meta[, 2]),
           pH   = as.numeric(meta[, 3]),
           run  = basename(file)) %>%      # one file = one plate = one run
    filter(!is.na(date)) %>%
    select(-date_raw)
}

#' All traces, in long format, one row per chamber x time point.
read_all_traces <- function(dir) {
  files <- list.files(dir, pattern = CFG$file_regex, full.names = TRUE)
  if (!length(files)) stop("No SDR files matched in ", dir)
  message("Reading ", length(files), " SDR files")

  map_dfr(files, read_sdr) %>%
    pivot_longer(matches(CFG$spot_regex), names_to = "Spot", values_to = "Ox") %>%
    filter(!is.na(Ox)) %>%
    mutate(g  = str_sub(Spot, 1, 1),                       # chamber row A..D
           ID = paste(fmt_num(hpf), fmt_num(pH), Spot, sep = "-"))
}

## ---- 3. Rate estimation -------------------------------------------------- ##
# Oxygen is recorded in umol L-1.  Each chamber yields a depletion rate in
# umol L-1 h-1, later corrected for the blank and normalised per individual.

#' Oxygen depletion rate of a single chamber trace.
#'
#' @param ox numeric vector, oxygen concentration (umol L-1)
#' @param t  POSIXct vector, acquisition times
#' @param method "scope" or "slope" (see CFG$rate_method)
#' @return one-row tibble: scope, hour, rate, r2
#' The linear fit is computed for EVERY chamber whatever the estimator in use,
#' because two quality-control facts depend on it and both must be reportable:
#'   r2         how linear the trace is (the Methods claim "checked for
#'              linearity" is only defensible if this is actually reported);
#'   net_rise   TRUE when the fitted slope is positive, i.e. oxygen INCREASED
#'              over the run.  The (max - min)/dt estimator takes an absolute
#'              difference and therefore returns a positive "consumption" rate
#'              for such a chamber.  These must be counted, not hidden.
chamber_rate <- function(ox, t, method = CFG$rate_method) {
  h  <- as.numeric(difftime(t, min(t), units = "hours"))
  fm <- stats::lm(ox ~ h)
  sl <- unname(coef(fm)[2])
  qc <- tibble(r2 = summary(fm)$r.squared, slope = sl, net_rise = sl > 0,
               run_h = max(h))

  if (method == "scope") {
    i_max <- which.max(ox); i_min <- which.min(ox)
    dt <- abs(as.numeric(difftime(t[i_min], t[i_max], units = "hours")))
    bind_cols(tibble(scope = max(ox) - min(ox), hour = dt,
                     rate = if (dt > 0) (max(ox) - min(ox)) / dt else NA_real_,
                     # fraction of the record actually spanned by the estimator
                     dt_frac = if (max(h) > 0) dt / max(h) else NA_real_), qc)
  } else {
    bind_cols(tibble(scope = max(ox) - min(ox), hour = max(h),
                     rate = abs(sl), dt_frac = 1), qc)
  }
}

## ---- 4. Metadata --------------------------------------------------------- ##

#' Larval counts per chamber.
#'
#' IMPORTANT — the stored `ID` column of MTD_RESPIROMETRY.csv contains
#' transcription errors: at 60 hpf the blank chambers A5 and A6 of pH 7.8 and
#' of pH 8 are labelled "60-7.8-A1/A2" and "60-8-A1/A2".  Joining on that
#' column duplicates the A1/A2 blanks, drops A5/A6, and therefore shifts the
#' blank correction of the whole 60 hpf plate.  The key is rebuilt here from
#' the atomic columns hpf / pH / Lane; every disagreement is written out.
#'
#' @param legacy TRUE reproduces the original (buggy) join, for provenance
read_metadata <- function(legacy = FALSE) {
  raw <- read_csv(file.path(CFG$data_dir, CFG$mtd_file)) %>%
    drop_na(Lane, N) %>%
    filter(Lane != "NANA") %>%
    mutate(ID_key = paste(fmt_num(hpf), fmt_num(pH), Lane, sep = "-"))

  conflicts <- raw %>% filter(ID != ID_key) %>% select(hpf, pH, Lane, ID, ID_key)
  write_tsv(conflicts, out("00_metadata_key_conflicts.tsv"))
  if (nrow(conflicts) && !legacy)
    message("Metadata: ", nrow(conflicts),
            " rows whose stored ID disagrees with hpf/pH/Lane were repaired")

  raw %>%
    transmute(ID     = if (legacy) ID else ID_key,
              N      = as.numeric(N),
              Design = if_else(Design == "Control", "Blank", Design)) %>%
    # a "blank" in which larvae were later found is not a blank
    mutate(N = if_else(Design == "Blank" & N > 0, NA_real_, N)) %>%
    drop_na(N) %>%
    { if (legacy) . else distinct(., ID, .keep_all = TRUE) }
}

## ---- 5. Build the analysis table ----------------------------------------- ##

#' Chamber-level rates, blank-corrected, per individual, in pmol O2 ind-1 h-1.
#'
#' @param traces output of read_all_traces()
#' @param legacy passed to read_metadata()
#' @return list(rates = all chambers with outlier flags, blanks = blank summary)
build_rates <- function(traces, legacy = FALSE, method = CFG$rate_method) {
  mtd <- read_metadata(legacy)

  raw <- traces %>%
    group_by(run, hpf, pH, g, Spot, ID) %>%
    reframe(chamber_rate(Ox, date, method)) %>%
    inner_join(mtd, by = "ID", relationship = "many-to-many")

  # Background correction: mean blank rate of the SAME plate (stage x pH).
  blanks <- raw %>%
    filter(Design == "Blank") %>%
    group_by(hpf, pH) %>%
    summarise(rate_blank = mean(rate, na.rm = TRUE),
              sd_blank   = sd(rate,   na.rm = TRUE),
              n_blank    = n(), .groups = "drop")

  rates <- raw %>%
    filter(Design == "Experimental") %>%
    left_join(blanks, by = c("hpf", "pH")) %>%
    mutate(rate_adj  = rate - rate_blank,                     # umol L-1 h-1
           # per chamber volume and per individual -> pmol O2 ind-1 h-1
           r_ind_adj = rate_adj / N * CFG$chamber_L * CFG$umol_to_pmol) %>%
    group_by(hpf, pH) %>%
    mutate(is_outlier = rstatix::is_outlier(r_ind_adj),
           is_extreme = rstatix::is_extreme(r_ind_adj)) %>%
    ungroup()

  list(rates = rates, blanks = blanks)
}

#' Add the modelling factors.  `contr.sum` coding is set globally further down,
#' which is what makes Type III SS interpretable.
add_factors <- function(d) {
  relabel <- function(x, map = CFG$stage_recode) {
    x <- fmt_num(x)
    if (!is.null(map)) { hit <- x %in% names(map); x[hit] <- unname(map[x[hit]]) }
    x
  }
  d %>% mutate(
    stage = factor(relabel(hpf), levels = relabel(sort(unique(hpf)))),
    # levels are intersected with what is actually present, so that dropping a
    # treatment from data_dir (e.g. re-running on pH 8.0 vs 7.6 only) does not
    # leave an empty factor level that would break aov()/emmeans
    pH    = factor(fmt_num(pH),
                   levels = intersect(c("8", "7.8", "7.6"), fmt_num(sort(unique(pH), decreasing = TRUE)))),
    treat = factor(if_else(fmt_num(pH) == CFG$control_pH, "control", "acidified"),
                   levels = c("control", "acidified")))
}

traces <- read_all_traces(CFG$data_dir)
built  <- build_rates(traces, legacy = FALSE)
rates  <- add_factors(built$rates)
dat    <- if (CFG$drop_outliers) filter(rates, !is_outlier) else rates

write_tsv(rates, out("01_rates_by_chamber.tsv"))

say("Blank (row A) chambers retained per plate")
print(as.data.frame(built$blanks), digits = 4)

## ---- 5b. Trace quality control ------------------------------------------- ##
# Three facts about the raw traces that the Methods section must be able to
# state, and that decide whether a plate is usable at all.

trace_qc <- built$rates %>%          # experimental chambers only
  group_by(hpf, pH) %>%
  summarise(n = n(),
            run_h        = round(max(run_h), 2),          # length of the record
            dt_median_h  = round(median(hour), 2),        # window used by "scope"
            dt_frac_med  = round(median(dt_frac), 2),
            r2_median    = round(median(r2), 3),
            pct_r2_lt_09 = round(100 * mean(r2 < 0.9)),
            pct_net_rise = round(100 * mean(net_rise)),    # oxygen went UP
            .groups = "drop")

say("Trace quality control (experimental chambers)")
print(as.data.frame(trace_qc))
write_tsv(trace_qc, out("12_trace_qc.tsv"))

if (any(trace_qc$pct_net_rise > 25))
  warning("Plates where >25 % of chambers show a NET RISE in oxygen: ",
          paste(with(filter(trace_qc, pct_net_rise > 25), paste0(hpf, " hpf / pH ", pH)),
                collapse = "; "),
          ". The (max-min)/dt estimator reports these as positive consumption.")

## ---- 6. Descriptive table (Table 2) -------------------------------------- ##

cell_summary <- dat %>%
  group_by(stage, pH) %>%
  summarise(n = n(), mean = mean(r_ind_adj), sd = sd(r_ind_adj),
            sem = sd / sqrt(n), median = median(r_ind_adj), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

say("Table 2 — RM (pmol O2 ind-1 h-1) per stage x pH cell")
print(as.data.frame(cell_summary))
write_tsv(cell_summary, out("02_summary_by_cell.tsv"))

say("Design balance (n chambers per cell)")
print(table(stage = dat$stage, pH = dat$pH))
cat("\nChambers analysed N =", nrow(dat),
    "| before outlier removal =", nrow(rates),
    "| plate runs =", n_distinct(dat$run), "\n")

## ---- 7. Assumptions ------------------------------------------------------ ##
# Normality is a property of the residuals, not of the raw cell values, so it
# is tested on the residuals of the fitted model.

fit0 <- aov(r_ind_adj ~ stage * pH, data = dat)
sw   <- shapiro.test(residuals(fit0))
lv   <- rstatix::levene_test(dat, r_ind_adj ~ stage * pH, center = "median")

assump <- tibble(
  test      = c("Shapiro-Wilk (model residuals)", "Levene / Brown-Forsythe (cells)"),
  statistic = c(unname(sw$statistic), lv$statistic),
  df1       = c(NA, lv$df1), df2 = c(NA, lv$df2),
  p         = c(sw$p.value, lv$p)) %>%
  mutate(assumption_met = p > CFG$alpha)

say("Assumption checks")
print(as.data.frame(assump), digits = 4)
write_tsv(assump, out("03_assumptions.tsv"))

## ---- 8. TWO-WAY ANOVA — the Reviewer's request --------------------------- ##
# Model:  RM ~ stage + pH + stage:pH,  both factors fixed.
# The design is unbalanced (5-11 chambers per cell), so:
#   * Type I  SS are order-dependent  -> not used;
#   * Type II SS test each main effect after the other main effect but ignoring
#     the interaction -> most powerful when the interaction is weak;
#   * Type III SS test each term after all others, respecting marginality, and
#     require sum-to-zero contrasts -> reported alongside.
# Both are given because they answer slightly different questions and, in an
# unbalanced design, they do not coincide.

options(contrasts = c("contr.sum", "contr.poly"))
fit  <- aov(r_ind_adj ~ stage * pH, data = dat)

#' Normalise the several shapes car::Anova() returns (classical SS table, or a
#' Wald table when white.adjust is used) into one tidy data frame.
tidy_anova <- function(a, label) {
  d  <- as.data.frame(a) %>% tibble::rownames_to_column("term")
  nm <- names(d)
  d %>%
    mutate(SS = if ("Sum Sq" %in% nm) .data[["Sum Sq"]] else NA_real_,
           df = .data[["Df"]],
           F  = if ("F value" %in% nm) .data[["F value"]] else .data[["F"]],
           p  = .data[["Pr(>F)"]]) %>%
    filter(term != "(Intercept)") %>%
    select(term, SS, df, F, p) %>%
    mutate(model = label, .before = 1)
}

anova_tbl <- bind_rows(
  tidy_anova(car::Anova(fit, type = 2), "Type II SS"),
  tidy_anova(car::Anova(fit, type = 3), "Type III SS")) %>%
  group_by(model) %>%
  # p-values of the three model terms are adjusted together: they are three
  # hypotheses tested on the same response in the same experiment.
  mutate(p.adj = p.adjust(p, method = CFG$padjust),
         signif = case_when(is.na(p.adj)  ~ NA_character_,
                            p.adj < 0.001 ~ "***",
                            p.adj < 0.01  ~ "**",
                            p.adj < 0.05  ~ "*", TRUE ~ "ns")) %>%
  ungroup() %>%
  mutate(df_resid = df.residual(fit))

# Effect sizes on the Type II fit: partial eta2 (share of variance left after
# the other terms) and eta2 (share of total variance).
es_p <- effectsize::eta_squared(car::Anova(fit, type = 2), partial = TRUE, ci = .95) %>%
  as.data.frame() %>%
  select(term = Parameter, eta2_partial = Eta2_partial,
         eta2p_CI_low = CI_low, eta2p_CI_high = CI_high)
es_t <- effectsize::eta_squared(car::Anova(fit, type = 2), partial = FALSE) %>%
  as.data.frame() %>% select(term = Parameter, eta2 = Eta2)

anova_tbl <- anova_tbl %>% left_join(es_p, by = "term") %>% left_join(es_t, by = "term")

say("TWO-WAY ANOVA — r_ind_adj ~ stage * pH")
print(as.data.frame(anova_tbl), digits = 4)
write_tsv(anova_tbl, out("04_anova_table.tsv"))

## ---- 9. Robustness ------------------------------------------------------- ##
# Levene's test fails, so the OLS F-ratios above are reported together with
# four analyses that relax the assumption they lean on.  Agreement across them
# is what licenses the claim; disagreement must be declared in the Methods.

rob <- list()

# (a) log scale — the classical variance-stabilising transform.  CAVEAT: blank
#     correction pushes several chambers below zero, so a large additive shift
#     is required; with that shift the transform is nearly linear over most of
#     the range and mostly compresses the two largest cells.  Its null result
#     for pH should therefore NOT be read as evidence of no effect; HC3 and the
#     permutation test are the robust analyses to quote.
shift     <- if (min(dat$r_ind_adj) <= 0) abs(min(dat$r_ind_adj)) + 1 else 0
dat$log_r <- log(dat$r_ind_adj + shift)
rob$log   <- tidy_anova(car::Anova(aov(log_r ~ stage * pH, data = dat), type = 2),
                        sprintf("log(RM + %.2f), Type II", shift))

# (b) heteroscedasticity-consistent (HC3) Wald tests — same model, standard
#     errors that do not assume equal cell variances.
rob$hc3 <- tidy_anova(
  car::Anova(lm(r_ind_adj ~ stage * pH, data = dat), type = 2, white.adjust = "hc3"),
  "HC3 robust Wald, Type II")

# (c) Scheirer-Ray-Hare: the nonparametric two-way analogue (rank transform).
rob$srh <- local({
  d <- dat; d$rk <- rank(d$r_ind_adj)
  a  <- car::Anova(aov(rk ~ stage * pH, data = d), type = 2)
  MS <- var(d$rk)
  as.data.frame(a) %>% tibble::rownames_to_column("term") %>%
    filter(term != "Residuals") %>%
    transmute(model = "Scheirer-Ray-Hare", term,
              SS = `Sum Sq`, df = Df, F = `Sum Sq` / MS,      # F column holds H
              p = pchisq(`Sum Sq` / MS, Df, lower.tail = FALSE))
})

# (d) permutation test of the observed F statistics: no distributional
#     assumption at all.  Labels are permuted under the full null.
rob$perm <- local({
  obs <- tidy_anova(car::Anova(fit, type = 2), "obs") %>% filter(term != "Residuals")
  null <- replicate(CFG$n_perm, {
    d <- dat; d$r_ind_adj <- sample(d$r_ind_adj)
    tidy_anova(car::Anova(aov(r_ind_adj ~ stage * pH, data = d), type = 2), "n") %>%
      filter(term != "Residuals") %>% pull(F)
  })
  tibble(model = sprintf("Permutation (%d)", CFG$n_perm), term = obs$term,
         SS = NA_real_, df = obs$df, F = obs$F,
         p = rowMeans(null >= obs$F))
})

# (e) subset models: (i) without the 60 hpf plate, whose blank-corrected rates
#     collapse to ~0 with a variance an order of magnitude above the rest;
#     (ii) the two-level contrast actually reported in the manuscript.
rob$no60 <- tidy_anova(
  car::Anova(aov(r_ind_adj ~ stage * pH, data = filter(dat, stage != "60")), type = 2),
  "Type II, 60 hpf excluded")

# (only informative when a third pH level is present; otherwise it is the model
#  already fitted above)
if (nlevels(dat$pH) > 2) {
  dat2 <- dat %>% filter(pH %in% c(CFG$control_pH, CFG$contrast_pH)) %>% droplevels()
  rob$two_level <- tidy_anova(
    car::Anova(aov(r_ind_adj ~ stage * pH, data = dat2), type = 2),
    sprintf("Type II, pH %s vs %s only", CFG$control_pH, CFG$contrast_pH))
}

robust_tbl <- bind_rows(rob) %>% filter(term != "Residuals") %>%
  mutate(p.adj = p.adjust(p, CFG$padjust))

say("Robustness — same terms under five alternative analyses")
print(as.data.frame(robust_tbl), digits = 4)
write_tsv(robust_tbl, out("05_robustness.tsv"))

## ---- 9a. Cumulative-demand helpers (used by 9b and 11) ------------------- ##

#' Trapezoidal integral of y over x.
trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

#' Cumulative oxygen demand of the control and the acidified treatment under
#' the two competing definitions.
#' @return named vector: sums, integrals and the two percentage differences
cum_from <- function(d) {
  w <- d %>%
    filter(pH %in% c(CFG$control_pH, CFG$contrast_pH)) %>%
    group_by(stage, pH) %>% summarise(m = mean(r_ind_adj), .groups = "drop") %>%
    mutate(t = as.numeric(as.character(stage))) %>% arrange(t) %>%
    pivot_wider(names_from = pH, values_from = m, names_prefix = "p")
  ctl <- w[[paste0("p", CFG$control_pH)]]; acd <- w[[paste0("p", CFG$contrast_pH)]]
  c(sum_ctl = sum(ctl), sum_acd = sum(acd),
    sum_pct = 100 * (sum(acd) / sum(ctl) - 1),
    int_ctl = trapz(w$t, ctl), int_acd = trapz(w$t, acd),
    int_pct = 100 * (trapz(w$t, acd) / trapz(w$t, ctl) - 1))
}

## ---- 9b. Sensitivity grid: rate estimator x 60 hpf plate ----------------- ##
# Two analytical choices are not innocuous and must be reported:
#   * the rate estimator.  "scope" = (max - min)/dt takes the two extremes of a
#     noisy trace, so it is biased upward by an amount that depends on the
#     noise of each plate; "slope" uses every point of the trace.
#   * the 60 hpf plate, whose blank-corrected rates collapse towards zero (and
#     become negative under the slope estimator, i.e. the larval chambers
#     depleted oxygen more slowly than their own blanks — physically
#     impossible, and a signal that this run is not usable).
# The grid below shows how the pH effect and the cumulative percentage move
# across the four combinations.

sens_grid <- expand_grid(method = c("scope", "slope"),
                         stages = c("all stages", "60 hpf excluded")) %>%
  pmap_dfr(function(method, stages) {
    d <- add_factors(build_rates(traces, legacy = FALSE, method = method)$rates)
    if (CFG$drop_outliers) d <- filter(d, !is_outlier)
    if (stages == "60 hpf excluded") d <- filter(d, stage != "60") %>% droplevels()
    a <- tidy_anova(car::Anova(aov(r_ind_adj ~ stage * pH, data = d), type = 2), "II")
    cm <- cum_from(d)
    tibble(rate_estimator = method, subset = stages, n = nrow(d),
           F_pH   = a$F[a$term == "pH"],       p_pH   = a$p[a$term == "pH"],
           F_int  = a$F[a$term == "stage:pH"], p_int  = a$p[a$term == "stage:pH"],
           pct_sum = cm["sum_pct"], pct_integrated = cm["int_pct"])
  })

say("Sensitivity grid — rate estimator x inclusion of the 60 hpf plate")
print(as.data.frame(sens_grid), digits = 4)
write_tsv(sens_grid, out("09_sensitivity_grid.tsv"))

## ---- 10. Post-hoc contrasts ---------------------------------------------- ##
# The interaction is the biologically interesting term (is the acidification
# cost constant across development?), so pH is contrasted WITHIN each stage,
# Tukey-adjusted over the three pH levels; marginal contrasts are also given.

emm_by_stage <- emmeans(fit, ~ pH | stage)
ph_within  <- contrast(emm_by_stage, "pairwise", adjust = "tukey") %>%
  summary(infer = TRUE) %>% as.data.frame()
ph_marg    <- contrast(emmeans(fit, ~ pH),    "pairwise", adjust = "tukey") %>%
  summary(infer = TRUE) %>% as.data.frame()
stage_marg <- contrast(emmeans(fit, ~ stage), "pairwise", adjust = "tukey") %>%
  summary(infer = TRUE) %>% as.data.frame()

posthoc <- bind_rows(
  mutate(ph_within,  family = "pH within stage"),
  mutate(ph_marg,    family = "pH, averaged over stages",    stage = NA),
  mutate(stage_marg, family = "stage, averaged over pH", stage = NA))

say("Post-hoc — pH contrasts within each stage (Tukey-adjusted)")
print(as.data.frame(ph_within), digits = 4)
write_tsv(posthoc, out("06_posthoc_emmeans.tsv"))

# Games-Howell does not assume equal variances; it is the contrast to quote
# given that Levene's test failed.
gh <- dat %>% group_by(stage) %>% games_howell_test(r_ind_adj ~ pH)
say("Post-hoc — Games-Howell (unequal variances)")
print(as.data.frame(gh), digits = 4)

# Percentage elevation at pH 7.6 relative to pH 8.0, per stage.
elev <- cell_summary %>%
  filter(pH %in% c(CFG$control_pH, CFG$contrast_pH)) %>%
  select(stage, pH, mean) %>%
  pivot_wider(names_from = pH, values_from = mean, names_prefix = "pH") %>%
  mutate(pct_elevation = 100 * (.data[[paste0("pH", CFG$contrast_pH)]] /
                                .data[[paste0("pH", CFG$control_pH)]] - 1))
say("Elevation at pH %s relative to pH %s, per stage", CFG$contrast_pH, CFG$control_pH)
print(as.data.frame(elev), digits = 4)

## ---- 11. Cumulative oxygen consumption ----------------------------------- ##
# Two definitions, computed explicitly because they do NOT agree:
#   (a) unweighted sum of the four stage means — this is what "sum" means in
#       the current Table 2 legend, and it is what reproduces the previously
#       reported 53 %.  Its units are pmol O2 ind-1 h-1: it is a mean elevation
#       in rate, NOT a cumulative consumption.
#   (b) trapezoidal integration over 24-110 hpf, each rate weighted by the
#       duration of the interval it represents (24-48 h, 48-60 h, 60-110 h).
#       Units pmol O2 ind-1: this is the quantity the word "cumulative" names,
#       and it is the definition adopted in the manuscript (+57.3 %).

obs_cum <- cum_from(dat)

# Bootstrap: chambers are resampled within each cell, preserving the design.
boot <- replicate(CFG$n_boot, {
  d <- dat %>% group_by(stage, pH) %>% slice_sample(prop = 1, replace = TRUE) %>% ungroup()
  cum_from(d)[c("sum_pct", "int_pct")]
})

cumO2 <- tibble(
  definition = c("(a) unweighted sum of the 4 stage means",
                 "(b) trapezoidal integration, 24-110 hpf"),
  units      = c("pmol O2 ind-1 h-1", "pmol O2 ind-1"),
  control    = c(obs_cum["sum_ctl"], obs_cum["int_ctl"]),
  acidified  = c(obs_cum["sum_acd"], obs_cum["int_acd"]),
  pct_increase = c(obs_cum["sum_pct"], obs_cum["int_pct"]),
  ci_low     = c(quantile(boot["sum_pct", ], .025), quantile(boot["int_pct", ], .025)),
  ci_high    = c(quantile(boot["sum_pct", ], .975), quantile(boot["int_pct", ], .975)))

# Provenance: the same two quantities computed with the original, unrepaired
# metadata join, to document where the published 53 % came from.
legacy_dat <- add_factors(build_rates(traces, legacy = TRUE)$rates) %>%
  { if (CFG$drop_outliers) filter(., !is_outlier) else . }
legacy_cum <- cum_from(legacy_dat)
cumO2 <- bind_rows(cumO2, tibble(
  definition = c("(a) as published: unweighted sum, unrepaired metadata key",
                 "(b) as published: integration, unrepaired metadata key"),
  units      = c("pmol O2 ind-1 h-1", "pmol O2 ind-1"),
  control    = c(legacy_cum["sum_ctl"], legacy_cum["int_ctl"]),
  acidified  = c(legacy_cum["sum_acd"], legacy_cum["int_acd"]),
  pct_increase = c(legacy_cum["sum_pct"], legacy_cum["int_pct"]),
  ci_low = NA_real_, ci_high = NA_real_))

say("Cumulative oxygen demand at pH %s relative to pH %s", CFG$contrast_pH, CFG$control_pH)
print(as.data.frame(cumO2), digits = 4)
write_tsv(cumO2, out("07_cumulative_oxygen.tsv"))

## ---- 11b. Relative energy expenditure per stage --------------------------- ##
# The trapezoidal integral can be rewritten exactly as a weighted sum of the
# stage means, the weight of each stage being half the duration of the
# intervals adjacent to it:
#     24 hpf -> 12 h, 48 hpf -> 18 h, 60 hpf -> 31 h, 110 hpf -> 25 h  (86 h)
# Written this way, "the fraction of the cumulative total attributable to each
# stage" is well defined, the fractions sum to 100 %, and the total is
# identical to the trapezoidal integral reported above.  This is the column
# that the Table 2 legend calls relative energy expenditure.

#' Trapezoidal node weights (hours) for a vector of sampling times.
trapz_weights <- function(t) {
  n <- length(t); w <- numeric(n)
  w[1] <- (t[2] - t[1]) / 2
  w[n] <- (t[n] - t[n - 1]) / 2
  if (n > 2) w[2:(n - 1)] <- (t[3:n] - t[1:(n - 2)]) / 2
  w
}

rel_exp <- dat %>%
  group_by(pH, stage) %>%
  summarise(mean_rate = mean(r_ind_adj), .groups = "drop") %>%
  mutate(t = as.numeric(as.character(stage))) %>%
  arrange(pH, t) %>%
  group_by(pH) %>%
  mutate(weight_h         = trapz_weights(t),
         contribution     = mean_rate * weight_h,       # pmol O2 ind-1
         cumulative_total = sum(contribution),
         rel_expenditure  = 100 * contribution / cumulative_total) %>%
  ungroup() %>%
  select(pH, stage, mean_rate, weight_h, contribution, cumulative_total,
         rel_expenditure)

stopifnot(all.equal(
  unique(rel_exp$cumulative_total[rel_exp$pH == CFG$control_pH]),
  unname(obs_cum["int_ctl"])))          # the two routes must agree exactly

say("Relative energy expenditure (share of the cumulative total, per stage)")
print(as.data.frame(rel_exp), digits = 4)
write_tsv(rel_exp, out("10_relative_expenditure.tsv"))

## ---- 12. Manuscript-ready numbers ---------------------------------------- ##

pick   <- function(tm, md = "Type II SS") filter(anova_tbl, model == md, term == tm)[1, ]
fmt_p  <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)
S <- pick("stage"); H <- pick("pH"); I <- pick("stage:pH")

results_md <- sprintf(
"### Two-way ANOVA — respiratory metabolism

Response: blank-corrected, individual-normalised oxygen consumption
(pmol O2 ind-1 h-1). Factors: developmental stage (%s hpf) and pH (%s), both
fixed. Unbalanced design, Type II sums of squares, N = %d chambers from %d
plate runs (%d-%d per stage x pH cell). p-values adjusted across the three
model terms (%s).

| Term | df | F | p | p adj. | partial eta2 |
|---|---|---|---|---|---|
| Developmental stage | %d, %d | %.2f | %s | %s | %.3f |
| pH | %d, %d | %.2f | %s | %s | %.3f |
| Stage x pH | %d, %d | %.2f | %s | %s | %.3f |

Cumulative oxygen demand at pH %s relative to pH %s:
  (a) unweighted sum of the four stage means: %+.1f %% (bootstrap 95%% CI %.1f to %.1f)
  (b) duration-weighted trapezoidal integration: %+.1f %% (bootstrap 95%% CI %.1f to %.1f)

Assumptions: Shapiro-Wilk on residuals W = %.3f, p = %.3f; Levene F = %.2f,
p = %.2g -- variances are heterogeneous, see 05_robustness.tsv.
",
paste(levels(dat$stage), collapse = ", "), paste(levels(dat$pH), collapse = ", "),
nrow(dat), n_distinct(dat$run), min(table(dat$stage, dat$pH)), max(table(dat$stage, dat$pH)),
CFG$padjust,
S$df, S$df_resid, S$F, fmt_p(S$p), fmt_p(S$p.adj), S$eta2_partial,
H$df, H$df_resid, H$F, fmt_p(H$p), fmt_p(H$p.adj), H$eta2_partial,
I$df, I$df_resid, I$F, fmt_p(I$p), fmt_p(I$p.adj), I$eta2_partial,
CFG$contrast_pH, CFG$control_pH,
cumO2$pct_increase[1], cumO2$ci_low[1], cumO2$ci_high[1],
cumO2$pct_increase[2], cumO2$ci_low[2], cumO2$ci_high[2],
assump$statistic[1], assump$p[1], assump$statistic[2], assump$p[2])

cat(results_md)
writeLines(results_md, out("08_results_paragraph.md"))

## ---- 13. Figures ---------------------------------------------------------- ##

pHpalette <- c(`8` = "#4575b4", `7.8` = "#abdda4", `7.6` = "#d73027")

p_main <- ggplot(dat, aes(stage, r_ind_adj, colour = pH, fill = pH)) +
  stat_boxplot(geom = "errorbar", width = .25, position = position_dodge(.8)) +
  geom_boxplot(width = .6, alpha = .25, outlier.shape = NA,
               position = position_dodge(.8)) +
  geom_point(position = position_jitterdodge(jitter.width = .12, dodge.width = .8),
             size = 1.2, alpha = .6) +
  scale_colour_manual("pH", values = pHpalette) +
  scale_fill_manual("pH", values = pHpalette) +
  labs(x = "Developmental stage (hpf)",
       y = expression("RM ("*O[2]~"pmol"~ind^-1~h^-1*")")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.border = element_blank(),
        panel.grid.minor = element_blank())

ggsave(out("Fig_RM_stage_pH.png"), p_main, width = 7, height = 4, dpi = 300)

png(out("Fig_diagnostics.png"), width = 1600, height = 1600, res = 180)
par(mfrow = c(2, 2)); plot(fit); dev.off()

# Software versions, so that the Methods section can be filled in from the run
# that produced the numbers rather than from memory.
writeLines(capture.output(sessionInfo()), out("11_sessionInfo.txt"))

say("Done — outputs in %s", normalizePath(CFG$out_dir))

## ---- 14. Caveat that must travel with these F statistics ------------------ ##
# Each stage x pH cell corresponds to ONE SDR plate run on larvae drawn from
# ONE culture vessel.  The chambers within a cell are therefore technical
# replicates of a single experimental unit: the pH F-ratio is tested against
# among-chamber variation, and the pH effect is confounded with plate/run.
# The analysis above is the one the Reviewer asked for, but the Methods must
# state the level of replication.  If independent vessels exist per treatment,
# the correct model is
#     lme4::lmer(r_ind_adj ~ stage * pH + (1 | vessel), data = dat)
# and the pH effect should be tested against the vessel variance component.
