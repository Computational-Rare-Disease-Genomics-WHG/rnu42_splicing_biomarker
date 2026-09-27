#!/usr/bin/env Rscript
# Layout: patchwork only -- never gridExtra arrangeGrob/grid.arrange
# =============================================================================
# figure2.R
# Paper Figure 2 (RNU4-2 splicing biomarker manuscript).
#
# Restructured 2026-07-26: panels A-B (LASSO coef bar, insT P(case) by
# cohort) move in from figure1.R (old panels E-F). The heatmap (old panel D)
# moves out to figure3.R, alongside the promoted Stem III classifier figure
# (old figure_s5.R). insT variant-group P(case)/annotation table/STX16 TPM
# (old panels A-C) shift down to C-E; the old figure1.R unmatched-control
# P(case) violin (old panel G) is folded into panel E as an extra category
# rather than its own panel. All computation code is preserved below.
#
# Re-laid-out 2026-07-27 into 3 rows x 2 panels (row 3 is a single
# full-width panel):
#   Panel A (row 1, 33%) — Bar chart of the K insT-signature LASSO
#             coefficients (beta at lambda.1se), ordered by magnitude
#             descending.
#   Panel B (row 1, 66%) — insT-LASSO predicted P(case) by cohort and group
#             (matched controls / ReNU insT train / ReNU insT holdout),
#             faceted by cohort (Whole-Blood / Lymphocytes).
#   Panel C (row 2, 66%) — Annotation table for the 6 insT signature events:
#             per-track transcript tier (MANE / Alternate / Unannotated, one
#             label each for the orange insT-associated track and the grey
#             other track, to the left of the diagram) | Gene | Event (alt
#             donor +/-N nt) | mean PSI (WB controls, %) | mean PSI (WB
#             n.64_65insT carriers, %) | Protein consequence | NMD. PSI and
#             NMD values aligned to the orange (insT-associated) track's
#             row; PSI oriented to that track too, not always the raw
#             rMATS long-form PSI. WB only (no lymphocyte column). Data
#             from a5ss_candidate_full_table.csv.
#   Panel D (row 2, 33%) — log2(TPM+1) boxplot for STX16 (the only predicted
#             frameshift among the insT events), GEL cases vs controls.
#   Panel E (row 3, full width) — insT LASSO P(case) by RNU4-2 variant group
#             (insT T-loop, non-insT T-loop, Stem III, n.72_73del, biallelic
#             RNU4-2, matched controls; boxplot + jitter) with an added
#             "WB unmatched controls" category (violin, no jitter -- ~5984
#             QC-pass whole-blood controls scored by
#             score_gel_controls_insT.R, $ctrl_prob_hist in
#             gel_controls_insT_robustness.rds).
#
# Inputs:
#   outputs/loo_results_insT_lasso.rds
#   outputs/gel_controls_insT_robustness.rds
#   outputs/psi_matrix.rds
#   outputs/nava_psi_matrix.rds
#   outputs/gel_controls_psi_cand_pool.rds
#   outputs/a5ss_candidate_full_table.csv
#   data/salmon_quant_paths.csv
#   data/biallelic_rnu42_platekeys.txt
#
# Outputs:
#   outputs/figures/figure2.pdf
#   outputs/figure2_summary.rds  (airlock-safe aggregate stats)
#
# Run from workspace root: pixi run Rscript scripts/gel/figure2.R
#
# DATA PROTECTION: all panels derive from GEL protected patient data +
# Nava cases. Selected event_ids + gene symbols + coefs remain inside the
# RE. The summary RDS carries only counts / aggregate stats.
# =============================================================================

# ---- libpath setup (pixi-aware) --------------------------------------------
# Under pixi, use only the pixi env's packages -- mixing in the GEL system
# libpath or the ~/re_gecip user lib (built against a different R) causes ABI
# mismatches: undefined-symbol crashes loading rlang / Rcpp. Outside pixi,
# fall back to GEL system R + the user lib.
invisible(local({
  user_lib <- path.expand("~/re_gecip/machine_learning/R_packages")
  if (nzchar(Sys.getenv("PIXI_ENVIRONMENT_NAME"))) {
    keep <- !normalizePath(.libPaths(), mustWork = FALSE) %in%
            normalizePath(user_lib,    mustWork = FALSE)
    if (any(!keep)) {
      message("[pixi] removing ", sum(!keep),
              " GEL-system libpath entries to avoid ABI mismatch")
      .libPaths(.libPaths()[keep])
    }
    return(invisible(NULL))
  }
  v <- substr(as.character(getRversion()), 1, 3)
  candidate <- switch(v,
    "4.0" = "/tools/aws-workspace-apps/ce/R/4.0.2",
    "4.3" = "/tools/aws-workspace-ubuntu-apps/ce/R/4.3.3",
    "4.5" = "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2",
    NULL)
  if (!is.null(candidate) && dir.exists(candidate)) {
    .libPaths(c(.libPaths(), candidate))
  }
}))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

# ---- Configuration ----------------------------------------------------------
GEL_PSI_PATH  <- Sys.getenv("GEL_PSI",  "outputs/psi_matrix.rds")
NAVA_PSI_PATH <- Sys.getenv("NAVA_PSI", "outputs/nava_psi_matrix.rds")

OUT_DIR_FIG  <- "outputs/figures"
OUT_FIG_2    <- file.path(OUT_DIR_FIG, "figure2.pdf")
OUT_SUMMARY  <- "outputs/figure2_summary.rds"

WB_RMATS_DIR   <- Sys.getenv("GEL_RMATS_DIR",   "data/rmats")
FASTA_PATH     <- Sys.getenv("FASTA_PATH",
  "/public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa")
CACHE_DIR_F2   <- "outputs/figure2_cache"

INST_HGVS          <- "n.64_65insT"
HGVS_DEL72       <- "n.72_73del"
INST_LOO_PATH    <- Sys.getenv("INST_LOO_RDS",    "outputs/loo_results_insT_lasso.rds")
INST_ROB_PATH    <- Sys.getenv("INST_ROB",       "outputs/gel_controls_insT_robustness.rds")
CAND_POOL_PATH   <- Sys.getenv("CAND_POOL",      "outputs/gel_controls_psi_cand_pool.rds")
BIALLELIC_FILE   <- Sys.getenv("BIALLELIC_FILE", "data/biallelic_rnu42_platekeys.txt")
ANNOT_PATH       <- Sys.getenv("ANNOT_CSV",      "outputs/a5ss_candidate_full_table.csv")
QUANT_PATHS_CSV  <- Sys.getenv("QUANT_PATHS",    "data/salmon_quant_paths.csv")
GENO_FILE        <- Sys.getenv("GENO_CSV",        "data/sample_genotype.csv")
NAVA_SEV_FILE    <- Sys.getenv("NAVA_SEV",        "outputs/nava_severity_labels.csv")
MANIFEST_PATH    <- Sys.getenv("MANIFEST_PATH",       "data/control_manifest.csv")
QCFAIL_MANIFEST  <- Sys.getenv("QCFAIL_MANIFEST_PATH","data/qcfail_manifest.csv")
MANIFEST_KEY     <- Sys.getenv("MANIFEST_KEY_COL", "rnaseq_platekey")
MANIFEST_YOB     <- Sys.getenv("MANIFEST_YOB_COL", "age_at_rnaseq")
DEV_CASE_BAMS    <- Sys.getenv("DEV_CASE_BAMS",
                               "data/rmats/inputs/renu_11_cases_paths.txt")
DEV_CTRL_BAMS    <- Sys.getenv("DEV_CTRL_BAMS",
                               "data/rmats/inputs/renu_11_cases_control_paths.txt")

# Cohort plot labels (paper text uses Whole-Blood / Lymphocytes in place of
# the database tokens GEL / Nava). Internal logic still keys on GEL / Nava.
WB_LAB <- "whole blood"
LY_LAB <- "lymphocytes"
COHORT_LABELS <- c(GEL = WB_LAB, Nava = LY_LAB)

# Paper figure palette (mirrors figure1.R; set 2026-04-29).
PAL_WB     <- "#2C5AA0"  # Whole-Blood navy
PAL_LY     <- "#2AA198"  # Lymphocytes teal
PAL_CTRL   <- "#888888"  # Ctrl grey
# RNU4-2 structural region (panels C + E). Oranges are reserved for the
# region split and not reused on the panel E clinical-severity axis.
PAL_TLOOP  <- "#D55E00"  # T-loop   (panel C box; panel E T-loop>Stem III)
PAL_STEM3  <- "#F4A582"  # Stem III (panel C box; panel E Stem III>T-loop)
# Nava clinical severity (panel E). Purple pair -- the ColorBrewer PuOr
# complement of the region oranges, so the two severity axes stay distinct.
PAL_SEVERE <- "#5E3C99"  # Severe phenotype (panel E)
PAL_MILD   <- "#B2ABD2"  # Mild phenotype   (panel E)

dir.create(CACHE_DIR_F2, recursive = TRUE, showWarnings = FALSE)

stopifnot("GEL PSI matrix not found"           = file.exists(GEL_PSI_PATH),
          "Nava PSI matrix not found"          = file.exists(NAVA_PSI_PATH),
          "insT LOO RDS not found"             = file.exists(INST_LOO_PATH),
          "Cand pool matrix not found"         = file.exists(CAND_POOL_PATH),
          "Biallelic platekeys file not found" = file.exists(BIALLELIC_FILE),
          "Sample genotype CSV not found"      = file.exists(GENO_FILE),
          "Nava severity labels file not found" = file.exists(NAVA_SEV_FILE))
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Load inputs ------------------------------------------------------------
cat("Loading GEL PSI:    ", GEL_PSI_PATH,  "\n", sep = "")
gel <- readRDS(GEL_PSI_PATH)
cat("Loading Nava PSI:   ", NAVA_PSI_PATH, "\n", sep = "")
nava <- readRDS(NAVA_PSI_PATH)

# HGVS + Nava severity lookup -- the authoritative sources for per-case
# genotype identity (independent of any LASSO fit's training population;
# per_case below is built from these, not from any *_loo$per_case_meta).
geno     <- read.csv(GENO_FILE, stringsAsFactors = FALSE)
nava_lbl <- local({
  first <- readLines(NAVA_SEV_FILE, n = 1L)
  sep   <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
  read.table(NAVA_SEV_FILE, header = TRUE, sep = sep, stringsAsFactors = FALSE,
             quote = "\"", fill = TRUE)
})
nava_lbl$sample_id <- paste0("nava:case:", nava_lbl$sample_n)
hgvs_lookup <- setNames(c(geno$hgvs,      nava_lbl$hgvs),
                         c(geno$sample_id, nava_lbl$sample_id))
cat(sprintf("HGVS lookup loaded: %d GEL + %d Nava entries\n",
            nrow(geno), nrow(nava_lbl)))

# Shared HGVS genotype palette (Set3, frequency-ranked) — used by Panel E
# so n.72_73del gets a consistent colour with figure3.R's heatmap.
hgvs_all_vals   <- hgvs_lookup[!is.na(hgvs_lookup)]
hgvs_levels_all <- names(sort(table(hgvs_all_vals), decreasing = TRUE))
n_hgvs_all      <- length(hgvs_levels_all)
brewer_cols_all  <- RColorBrewer::brewer.pal(max(3L, min(n_hgvs_all, 12L)), "Set3")
PAL_HGVS_SHARED  <- if (n_hgvs_all <= 12L) {
  setNames(brewer_cols_all[seq_len(n_hgvs_all)], hgvs_levels_all)
} else {
  setNames(colorRampPalette(brewer_cols_all)(n_hgvs_all), hgvs_levels_all)
}

# ---- Load insT LASSO model --------------------------------------------------
cat("Loading insT LOO RDS: ", INST_LOO_PATH, "\n", sep = "")
insT_loo <- readRDS(INST_LOO_PATH)

insT_cand    <- insT_loo$a5ss_cand
insT_sel_idx <- insT_loo$selected_idx
K_insT       <- insT_loo$K_selected
insT_gel_ev  <- insT_cand$gel_ev_id[insT_sel_idx]
insT_nava_ev <- insT_cand$nava_ev_id[insT_sel_idx]
insT_dir_v   <- insT_cand$direction[insT_sel_idx]
insT_mu_gel  <- insT_loo$gel_ctrl_means[insT_gel_ev]
insT_mu_nava <- insT_loo$nava_ctrl_means[insT_gel_ev]
insT_beta    <- insT_loo$fit_at_1se$beta[insT_gel_ev]
insT_icept   <- insT_loo$fit_at_1se$intercept
cat(sprintf("insT model: K=%d selected sites\n", K_insT))

# Biallelic samples live in the 6k cand pool matrix (manifest doesn't exclude them)
cat("Loading cand pool matrix: ", CAND_POOL_PATH, "\n", sep = "")
ctrl_pool <- readRDS(CAND_POOL_PATH)
stopifnot("Cand pool missing ijc/sjc" = all(c("ijc", "sjc") %in% names(ctrl_pool)))

biallelic_keys <- trimws(readLines(BIALLELIC_FILE))
biallelic_keys <- biallelic_keys[nzchar(biallelic_keys)]
stopifnot("Biallelic platekeys file is empty" = length(biallelic_keys) > 0L)

# ---- bb-logit deviation matrix per cohort ----------------------------------
# orient_ij = direction_j * (bb_logit(IJC,SJC)_ij - mu_ctrl,cohort,j)
bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}

# Sample metadata
gel_case_ids  <- gel$samples$sample_id[gel$samples$group  == "Case"]
gel_ctrl_ids  <- gel$samples$sample_id[gel$samples$group  == "Control"]
nava_case_ids <- nava$samples$sample_id[nava$samples$group == "Case"]
nava_ctrl_ids <- nava$samples$sample_id[nava$samples$group == "Control"]

# Built from the raw HGVS/severity sources (independent of any LASSO fit's
# training population) so every case genotype is represented here, not just
# insT/stemIII carriers -- insT_loo$per_case_meta only covers insT's own
# training set (insT cases + controls) and silently drops every other
# genotype group (non-insT T-loop, Stem III, n.72_73del, ...).
per_case <- data.frame(
  sample_id = c(gel_case_ids, nava_case_ids),
  cohort    = c(rep("GEL",  length(gel_case_ids)),
                rep("Nava", length(nava_case_ids))),
  stringsAsFactors = FALSE)
per_case$hgvs <- unname(hgvs_lookup[per_case$sample_id])
stopifnot("Some cases missing HGVS from hgvs_lookup" = !anyNA(per_case$hgvs))
per_case$is_inst <- per_case$cohort == "GEL" & per_case$hgvs == INST_HGVS

pub_theme <- theme_classic(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 9, face = "bold"),
        axis.text        = element_text(size = 9),
        axis.title       = element_text(size = 9),
        plot.title       = element_text(size = 9),
        legend.title     = element_text(size = 9),
        legend.text      = element_text(size = 9),
        plot.margin      = margin(2, 2, 2, 2))

# ---- Helpers shared across panels ------------------------------------------
# Asterisk convention reported in the figure legend (paper text):
#   ns p>0.05, * p<=0.05, ** p<=0.01, *** p<=0.001, **** p<=0.0001
read_bam_list <- function(path) {
  if (!file.exists(path)) return(character(0))
  raw  <- paste(readLines(path), collapse = ",")
  bams <- trimws(unlist(strsplit(raw, ",")))
  bams <- bams[nzchar(bams)]
  sub("\\.bam$", "", basename(bams))
}

p_to_stars <- function(p) {
  if (is.na(p))            "ns"
  else if (p <= 0.0001)    "****"
  else if (p <= 0.001)     "***"
  else if (p <= 0.01)      "**"
  else if (p <= 0.05)      "*"
  else                     "ns"
}

# Plotmath-parseable p-value label: plain decimal for p >= 0.001, otherwise
# "3.3 x 10^-16" style (mantissa %*% 10^exponent) -- use with
# annotate/geom_text(..., parse = TRUE).
format_p_sci <- function(p) {
  if (p >= 0.001) return(sprintf("p == %.3f", p))
  exp10  <- floor(log10(p))
  mant_r <- round(p / 10^exp10, 1)
  if (mant_r >= 10) { mant_r <- mant_r / 10; exp10 <- exp10 + 1 }
  sprintf("p == %.1f %%*%% 10^%d", mant_r, exp10)
}

sigmoid <- function(z) 1 / (1 + exp(-z))

# Assign each case to an RNU4-2 structural region by the n. position parsed
# from its HGVS string. This replaces the SGE function_score split with a
# position-based one: it tracks SGE moderate/strong closely (T-loop ~ SGE
# strong, Stem III ~ SGE moderate) but is defined purely by variant
# location, so it does not depend on SGE assay coverage.
#   T-loop   : n. position 63-70   (recurrent n.64_65insT cluster)
#   Stem III : n. position 75-79
# Variants outside both ranges, or with an unparseable HGVS, return NA and
# are dropped from the region comparison.
RNU4_TLOOP <- c(63L, 70L)
RNU4_STEM3 <- c(75L, 79L)
classify_rnu4_region <- function(hgvs) {
  out <- rep(NA_character_, length(hgvs))
  for (i in seq_along(hgvs)) {
    h <- hgvs[i]
    if (is.na(h) || !nzchar(h)) next
    coord <- regmatches(h, regexpr("n\\.[0-9]+(_[0-9]+)?", h))
    if (length(coord) == 0L) next
    pos <- as.integer(regmatches(coord, gregexpr("[0-9]+", coord))[[1]])
    if (length(pos) == 0L) next
    # T-loop takes precedence for the rare variant spanning both regions.
    if (any(pos >= RNU4_TLOOP[1] & pos <= RNU4_TLOOP[2]))
      out[i] <- "T-loop"
    else if (any(pos >= RNU4_STEM3[1] & pos <= RNU4_STEM3[2]))
      out[i] <- "Stem III"
  }
  out
}

# Per-cohort attached labels for jittered case dots in panels A/B.
per_case$cohort_lab <- unname(COHORT_LABELS[per_case$cohort])
per_case$region_bucket <- classify_rnu4_region(per_case$hgvs)

cat("\nRNU4-2 region assignment over cases (cohort x region):\n")
print(addmargins(table(cohort = per_case$cohort,
                       region = factor(per_case$region_bucket,
                                       levels = c("Stem III", "T-loop")),
                       useNA  = "always")))

# ---- insT LASSO P(case) for all GEL + Nava samples -------------------------
stopifnot(
  "insT selected gel events missing from gel$ijc"   = all(insT_gel_ev  %in% colnames(gel$ijc)),
  "insT selected nava events missing from nava$ijc" = all(insT_nava_ev %in% colnames(nava$ijc))
)
insT_gel_x  <- bb_logit(gel$ijc[,  insT_gel_ev,  drop = FALSE],
                         gel$sjc[,  insT_gel_ev,  drop = FALSE])
insT_nava_x <- bb_logit(nava$ijc[, insT_nava_ev, drop = FALSE],
                         nava$sjc[, insT_nava_ev, drop = FALSE])
colnames(insT_nava_x) <- insT_gel_ev   # canonical gel_ev_id names

insT_orient_gel  <- sweep(sweep(insT_gel_x,  2, insT_mu_gel,  "-"), 2, insT_dir_v, "*")
insT_orient_nava <- sweep(sweep(insT_nava_x, 2, insT_mu_nava, "-"), 2, insT_dir_v, "*")

insT_lp_gel  <- as.numeric(insT_orient_gel  %*% insT_beta) + insT_icept
insT_lp_nava <- as.numeric(insT_orient_nava %*% insT_beta) + insT_icept
names(insT_lp_gel)  <- rownames(insT_orient_gel)
names(insT_lp_nava) <- rownames(insT_orient_nava)
insT_prob_gel  <- sigmoid(insT_lp_gel)
insT_prob_nava <- sigmoid(insT_lp_nava)

insT_score_gel  <- setNames(rowSums(insT_orient_gel,  na.rm = TRUE), rownames(insT_orient_gel))
insT_score_nava <- setNames(rowSums(insT_orient_nava, na.rm = TRUE), rownames(insT_orient_nava))

stopifnot(
  "insT selected events missing from cand pool ijc" =
    all(insT_gel_ev %in% colnames(ctrl_pool$ijc))
)
biallelic_in_pool <- intersect(biallelic_keys, rownames(ctrl_pool$ijc))
stopifnot("No biallelic samples found in cand pool matrix" =
            length(biallelic_in_pool) > 0L)
cat(sprintf("Biallelic: %d/%d platekeys found in cand pool\n",
            length(biallelic_in_pool), length(biallelic_keys)))

bial_gel_x  <- bb_logit(ctrl_pool$ijc[biallelic_in_pool, insT_gel_ev, drop = FALSE],
                         ctrl_pool$sjc[biallelic_in_pool, insT_gel_ev, drop = FALSE])
bial_orient <- sweep(sweep(bial_gel_x, 2, insT_mu_gel, "-"), 2, insT_dir_v, "*")
bial_lp        <- as.numeric(bial_orient %*% insT_beta) + insT_icept
bial_prob      <- setNames(sigmoid(bial_lp), biallelic_in_pool)
bial_score_insT <- setNames(rowSums(bial_orient, na.rm = TRUE), biallelic_in_pool)
cat(sprintf("Biallelic insT P(case) range [%.3f, %.3f]  score range [%.3f, %.3f]\n",
            min(bial_prob), max(bial_prob),
            min(bial_score_insT), max(bial_score_insT)))

# ============================================================================
# NEW PANEL A: insT-specific LASSO coefficient bar (moved from figure1.R
# panel E). K selected coefficients from loo_results_insT_lasso.rds, ordered
# by coef magnitude descending.
# ============================================================================
gene_insT_a    <- insT_cand$GeneSymbol[insT_sel_idx]
coef_insT_a    <- unname(insT_beta)
sel_ord_insT_a <- order(-coef_insT_a)
gene_insT_ord  <- gene_insT_a[sel_ord_insT_a]
coef_insT_ord  <- coef_insT_a[sel_ord_insT_a]
if (anyDuplicated(gene_insT_ord))
  gene_insT_ord <- make.unique(gene_insT_ord, sep = ".")

coef_insT_df <- data.frame(
  gene_lab = factor(gene_insT_ord, levels = gene_insT_ord),
  coef     = coef_insT_ord,
  stringsAsFactors = FALSE)

fig2_panelA <- ggplot(coef_insT_df, aes(x = gene_lab, y = coef)) +
  geom_col(width = 0.5, fill = PAL_CTRL) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "LASSO coef (lambda.1se)") +
  pub_theme +
  theme(axis.text.x  = element_text(size = 9, face = "italic",
                                    angle = 45, hjust = 1),
        plot.margin  = margin(12, 4, 12, 4))

# ============================================================================
# NEW PANEL B: insT-specific LASSO predicted P(case) separation (moved from
# figure1.R panel F). Full-train linpred at lambda.1se passed through
# sigmoid, faceted by cohort (Whole-Blood / Lymphocytes). Adds a fourth
# group: non-insT ReNU cases scored with the insT-specific fit to show
# whether the signature generalises across RNU4-2 variants.
# ============================================================================
X_train_insT    <- insT_loo$X_train
beta_insT_vec   <- insT_loo$fit_at_1se$beta
intercept_insT  <- insT_loo$fit_at_1se$intercept
lp_train_insT   <- as.numeric(X_train_insT %*% beta_insT_vec + intercept_insT)
names(lp_train_insT) <- insT_loo$sample_ids_train
prob_train_insT <- sigmoid(lp_train_insT)

y_insT           <- insT_loo$y_train
cohort_tr_insT   <- insT_loo$cohort_train
is_gel_insT      <- cohort_tr_insT == "GEL"
is_nava_insT     <- cohort_tr_insT == "Nava"
is_case_insT     <- y_insT == 1L
is_ctrl_insT     <- y_insT == 0L

pr_ctrl_gel_insT     <- prob_train_insT[is_gel_insT  & is_ctrl_insT]
pr_ctrl_nava_insT    <- prob_train_insT[is_nava_insT & is_ctrl_insT]
pr_case_tr_gel_insT  <- prob_train_insT[is_gel_insT  & is_case_insT]
pr_case_tr_nava_insT <- prob_train_insT[is_nava_insT & is_case_insT]

per_case_insT_b   <- insT_loo$per_case_meta
prob_ho_gel_insT  <- per_case_insT_b$prob_full[per_case_insT_b$cohort == "GEL"  &
                                                per_case_insT_b$set   == "holdout"]
prob_ho_nava_insT <- per_case_insT_b$prob_full[per_case_insT_b$cohort == "Nava" &
                                                per_case_insT_b$set   == "holdout"]

b_df <- rbind(
  data.frame(cohort = "GEL",  group = "matched controls",
             prob = pr_ctrl_gel_insT,    stringsAsFactors = FALSE),
  data.frame(cohort = "GEL",  group = "n.64_65insT (train)",
             prob = pr_case_tr_gel_insT, stringsAsFactors = FALSE),
  data.frame(cohort = "GEL",  group = "n.64_65insT (holdout)",
             prob = prob_ho_gel_insT,    stringsAsFactors = FALSE),
  data.frame(cohort = "Nava", group = "matched controls",
             prob = pr_ctrl_nava_insT,   stringsAsFactors = FALSE),
  data.frame(cohort = "Nava", group = "n.64_65insT (train)",
             prob = pr_case_tr_nava_insT,stringsAsFactors = FALSE),
  data.frame(cohort = "Nava", group = "n.64_65insT (holdout)",
             prob = prob_ho_nava_insT,   stringsAsFactors = FALSE)
)
b_df$group  <- factor(b_df$group,
                      levels = c("matched controls", "n.64_65insT (train)", "n.64_65insT (holdout)"))
b_df$cohort <- factor(b_df$cohort, levels = c("GEL", "Nava"),
                      labels = c(WB_LAB, LY_LAB))
b_df$colour_class <- ifelse(
  b_df$group == "matched controls", "matched controls",
  paste("ReNU", as.character(b_df$cohort)))
renu_wb_lab_b <- paste("ReNU", WB_LAB)
renu_ly_lab_b <- paste("ReNU", LY_LAB)

# Mann-Whitney U p-value + AUC (rank-based, threshold-free -- same formula
# as auc_two_sample() in loo_evaluate_insT_lasso.R) for train-vs-matched-
# controls and holdout-vs-matched-controls, per cohort. BH-adjusted jointly
# across all 4 tests (consistent with figure_s5.R/figure_s6.R convention).
auc_two_sample <- function(case_scores, ctrl_scores) {
  case_scores <- case_scores[!is.na(case_scores)]
  ctrl_scores <- ctrl_scores[!is.na(ctrl_scores)]
  if (length(case_scores) == 0L || length(ctrl_scores) == 0L) return(NA_real_)
  r <- rank(c(case_scores, ctrl_scores))
  U <- sum(r[seq_along(case_scores)]) - length(case_scores) * (length(case_scores) + 1L) / 2
  U / (length(case_scores) * length(ctrl_scores))
}

panelB_tests <- list(
  list(cohort = WB_LAB, set = "train",   x2 = 2, case = pr_case_tr_gel_insT,  ctrl = pr_ctrl_gel_insT),
  list(cohort = WB_LAB, set = "holdout", x2 = 3, case = prob_ho_gel_insT,     ctrl = pr_ctrl_gel_insT),
  list(cohort = LY_LAB, set = "train",   x2 = 2, case = pr_case_tr_nava_insT, ctrl = pr_ctrl_nava_insT),
  list(cohort = LY_LAB, set = "holdout", x2 = 3, case = prob_ho_nava_insT,    ctrl = pr_ctrl_nava_insT)
)
panelB_p   <- vapply(panelB_tests, function(t)
  wilcox.test(t$case, t$ctrl, exact = FALSE)$p.value, numeric(1L))
panelB_auc <- vapply(panelB_tests, function(t) auc_two_sample(t$case, t$ctrl), numeric(1L))
panelB_p_adj <- p.adjust(panelB_p, method = "BH")

for (i in seq_along(panelB_tests)) {
  t <- panelB_tests[[i]]
  cat(sprintf("Panel B %s %s vs matched controls: p=%.3e (BH=%.3e) AUC=%.3f n_case=%d n_ctrl=%d\n",
              t$cohort, t$set, panelB_p[i], panelB_p_adj[i], panelB_auc[i],
              length(t$case[!is.na(t$case)]), length(t$ctrl[!is.na(t$ctrl)])))
}

format_p_auc <- function(p, auc) sprintf('%s * "," ~ AUC == %.2f', format_p_sci(p), auc)

# Two brackets per facet (ctrl-vs-train, ctrl-vs-holdout), stacked at
# different heights so they don't overlap; y-axis extended above 1 to fit.
panelB_ann <- data.frame(
  cohort = factor(vapply(panelB_tests, `[[`, character(1L), "cohort"),
                  levels = c(WB_LAB, LY_LAB), labels = c(WB_LAB, LY_LAB)),
  x    = 1,
  xend = vapply(panelB_tests, `[[`, numeric(1L), "x2"),
  y    = ifelse(vapply(panelB_tests, `[[`, character(1L), "set") == "train", 1.06, 1.16),
  label = mapply(format_p_auc, panelB_p_adj, panelB_auc),
  stringsAsFactors = FALSE
)

fig2_panelB <- ggplot(b_df, aes(x = group, y = prob, colour = colour_class)) +
  # Holdout group is only n=2/3 per cohort -- show points with a median
  # crossbar instead of a boxplot (IQR undefined at n<=3).
  geom_boxplot(data = function(d) d[d$group != "n.64_65insT (holdout)", ],
               fill = NA, outlier.shape = NA,
               width = 0.45, linewidth = 0.3) +
  stat_summary(data = function(d) d[d$group == "n.64_65insT (holdout)", ],
               fun = median, geom = "crossbar",
               width = 0.45, fatten = 0, linewidth = 0.3) +
  geom_jitter(width = 0.13, size = 0.85, alpha = 0.75) +
  geom_segment(data = panelB_ann, aes(x = x, xend = xend, y = y, yend = y),
              colour = "grey30", linewidth = 0.35, inherit.aes = FALSE) +
  geom_text(data = panelB_ann, aes(x = (x + xend) / 2, y = y * 1.02, label = label),
           parse = TRUE, size = 2.2, colour = "grey20", hjust = 0.5, inherit.aes = FALSE) +
  scale_colour_manual(values = setNames(
                        c(PAL_CTRL, PAL_WB, PAL_LY),
                        c("matched controls", renu_wb_lab_b, renu_ly_lab_b)),
                      name = NULL, guide = "none") +
  scale_y_continuous(limits = c(0, 1.25), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  # Display-only relabel (line breaks) so the x-axis text can sit
  # horizontal instead of angled -- underlying group values (used above for
  # the Wilcoxon/AUC tests) are unchanged.
  scale_x_discrete(labels = c("matched controls"       = "matched\ncontrols",
                              "n.64_65insT (train)"    = "n.64_65insT\n(train)",
                              "n.64_65insT (holdout)"  = "n.64_65insT\n(holdout)")) +
  facet_wrap(~ cohort, nrow = 1L) +
  labs(x = NULL, y = "insT LASSO predicted P(case)") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
        strip.text   = element_text(size = 9, face = "plain"))

# ============================================================================
# NEW PANEL E: insT LASSO P(case) by RNU4-2 variant group, combined with the
# 6k whole-blood unmatched-control distribution (folded in as an extra
# category rather than its own panel; moved from figure1.R panel G).
# Group 1: non-insT T-loop carriers (T-loop HGVS, not n.64_65insT)
# Group 2: Stem III carriers (positions 75-79)
# Group 3: n.72_73del -- position 72-73, outside both T-loop [63-70] and Stem III
#           [75-79]; n=2 (confirmed via diagnostic print -- both carriers
#           share the exact same HGVS string, matched correctly)
# Group 4: Biallelic RNU4-2 (from gel_controls_psi_cand_pool.rds)
# Group 5: Unmatched controls -- ~5984 QC-pass whole-blood controls scored
#           with the insT-specific fit (score_gel_controls_insT.R); shown as
#           a violin (no jitter -- too many points) rather than boxplot.
# insT carriers (n.64_65insT) are the model training set -- excluded here.
# cat_other_het kept as "Stem III / other" for deprecated panel code that runs
# below but is not assembled into the output figure.
# ============================================================================
cat_unmatched  <- "WB unmatched controls"
cat_new_ctrl   <- "Matched controls"
cat_tloop_inst <- "n.64_65insT"
cat_tloop_ni   <- "other T-loop"
cat_stem3_only <- "Stem III"
cat_del72      <- HGVS_DEL72          # "n.72_73del"
cat_other_het  <- "Stem III / other"  # deprecated panels only
cat_biallelic  <- "Biallelic RNU4-2"
# Shared group order for this panel
shared_grp_lev <- c(cat_new_ctrl, cat_unmatched, cat_tloop_inst, cat_tloop_ni,
                    cat_stem3_only, cat_del72, cat_biallelic)

per_case$p_case_insT <- ifelse(
  per_case$cohort == "GEL",
  insT_prob_gel[per_case$sample_id],
  insT_prob_nava[per_case$sample_id])
per_case$score_insT <- ifelse(
  per_case$cohort == "GEL",
  insT_score_gel[per_case$sample_id],
  insT_score_nava[per_case$sample_id])

per_case$new_group <- with(per_case, {
  grp <- rep(NA_character_, nrow(per_case))
  grp[region_bucket == "T-loop"  & hgvs == INST_HGVS] <- cat_tloop_inst
  grp[region_bucket == "T-loop"  & hgvs != INST_HGVS] <- cat_tloop_ni
  grp[region_bucket == "Stem III" | is.na(region_bucket)] <- cat_other_het
  grp
})
ctrl_prob_insT  <- c(insT_prob_gel[gel_ctrl_ids],   insT_prob_nava[nava_ctrl_ids])
ctrl_score_insT <- c(insT_score_gel[gel_ctrl_ids], insT_score_nava[nava_ctrl_ids])
cat("\nPanel E group counts:\n")
print(table(group = per_case$new_group, useNA = "always"))

grp_sel <- !is.na(per_case$new_group)
# Map from the "Stem III / other" new_group label to the two display categories
new_a_cat_raw <- per_case$new_group[grp_sel]

# DIAGNOSTIC: is_del72_a below is an EXACT string match against HGVS_DEL72
# ("n.72_73del"). Print every raw HGVS string currently in the catch-all
# bucket (quoted, so leading/trailing whitespace or a differently-formatted
# variant -- e.g. "n.72_73delTG" specifying the deleted bases -- is visible)
# to check whether all n.72_73del carriers actually match this exact string.
cat_other_rows <- per_case[grp_sel, , drop = FALSE][new_a_cat_raw == cat_other_het, ]
cat(sprintf("\nDIAGNOSTIC: %d cases in the Stem III/other catch-all bucket (region_bucket == \"Stem III\" or NA):\n",
            nrow(cat_other_rows)))
print(data.frame(sample_id = cat_other_rows$sample_id,
                 cohort    = cat_other_rows$cohort,
                 hgvs_quoted = sprintf("\"%s\"", cat_other_rows$hgvs),
                 exact_match_HGVS_DEL72 = cat_other_rows$hgvs == HGVS_DEL72,
                 stringsAsFactors = FALSE),
      row.names = FALSE)
cat(sprintf("HGVS_DEL72 constant is: \"%s\"\n", HGVS_DEL72))

is_del72_a    <- new_a_cat_raw == cat_other_het &
                 !is.na(per_case$hgvs[grp_sel]) &
                 per_case$hgvs[grp_sel] == HGVS_DEL72
new_a_cat_raw[is_del72_a]                  <- cat_del72
new_a_cat_raw[new_a_cat_raw == cat_other_het] <- cat_stem3_only

new_a_df <- rbind(
  data.frame(category = cat_new_ctrl, value = unname(ctrl_prob_insT),
             hgvs = NA_character_, stringsAsFactors = FALSE),
  data.frame(category = new_a_cat_raw,
             value    = per_case$p_case_insT[grp_sel],
             hgvs     = per_case$hgvs[grp_sel],
             stringsAsFactors = FALSE),
  data.frame(category = cat_biallelic, value = unname(bial_prob),
             hgvs = NA_character_, stringsAsFactors = FALSE)
)

# Unmatched whole-blood controls (violin category) -- optional input.
if (file.exists(INST_ROB_PATH)) {
  insT_rob  <- readRDS(INST_ROB_PATH)
  stopifnot("gel_controls_insT_robustness.rds missing ctrl_prob field" =
              "ctrl_prob" %in% names(insT_rob))
  # Prefer ctrl_prob_hist (biallelic + HD excluded, matches figure_s6.R n);
  # fall back to ctrl_prob (dedup only) for older RDS without that field.
  prob_vec_unmatched <- if ("ctrl_prob_hist" %in% names(insT_rob))
    insT_rob$ctrl_prob_hist else insT_rob$ctrl_prob
  n_unmatched    <- length(prob_vec_unmatched)
  n_unmatched_na <- sum(!is.finite(prob_vec_unmatched))
  cat(sprintf("\nPanel E: %d whole-blood controls P(case) loaded from %s (%d NA dropped by ggplot)\n",
              n_unmatched, INST_ROB_PATH, n_unmatched_na))
  # Tukey-rule outlier count (matches geom_boxplot's default outlier draw):
  # any point outside [Q1 - 1.5 IQR, Q3 + 1.5 IQR].
  ub_vals <- prob_vec_unmatched[is.finite(prob_vec_unmatched)]
  ub_q    <- quantile(ub_vals, c(0.25, 0.75), na.rm = TRUE)
  ub_iqr  <- ub_q[2L] - ub_q[1L]
  ub_lo   <- ub_q[1L] - 1.5 * ub_iqr
  ub_hi   <- ub_q[2L] + 1.5 * ub_iqr
  n_out_lo <- sum(ub_vals < ub_lo)
  n_out_hi <- sum(ub_vals > ub_hi)
  cat(sprintf(
    "Panel E: unmatched-controls Tukey outliers: %d total (%d below %.3f, %d above %.3f; IQR=[%.3f, %.3f])\n",
    n_out_lo + n_out_hi, n_out_lo, ub_lo, n_out_hi, ub_hi, ub_q[1L], ub_q[2L]))
  new_a_df <- rbind(new_a_df,
    data.frame(category = cat_unmatched, value = unname(prob_vec_unmatched),
               hgvs = NA_character_, stringsAsFactors = FALSE))
} else {
  cat(sprintf("\n[panel E] %s not found; skipping unmatched-controls category.\n", INST_ROB_PATH))
  n_unmatched <- 0L
}

# Display groups: shared_grp_lev, minus "WB unmatched controls" if unavailable.
new_a_lev <- if (n_unmatched > 0L) shared_grp_lev else setdiff(shared_grp_lev, cat_unmatched)
new_a_df$category <- factor(new_a_df$category, levels = new_a_lev)

new_a_pal <- c(
  "WB unmatched controls" = PAL_WB,
  "Matched controls"    = PAL_CTRL,
  "n.64_65insT"         = PAL_TLOOP,
  "other T-loop"        = PAL_TLOOP,
  "Stem III"            = PAL_STEM3,
  "n.72_73del"          = if (!is.null(PAL_HGVS_SHARED) && "n.72_73del" %in% names(PAL_HGVS_SHARED))
                            PAL_HGVS_SHARED[["n.72_73del"]] else PAL_LY,
  "Biallelic RNU4-2"    = PAL_SEVERE
)

x_labels_e <- setNames(new_a_lev, new_a_lev)
x_labels_e[cat_new_ctrl]  <- "Matched\ncontrols"
x_labels_e[cat_biallelic] <- "Biallelic\nRNU4-2"
if (n_unmatched > 0L)
  x_labels_e[cat_unmatched] <- sprintf("WB controls\n(n=%s)",
                                       format(n_unmatched, big.mark = ","))

cat(sprintf("\nPanel E sample counts: %s\n",
            paste(levels(new_a_df$category),
                  as.integer(table(new_a_df$category)),
                  sep = "=n", collapse = ", ")))

# Mann-Whitney U p-value + AUC (auc_two_sample(), defined above for Panel B).
# Two comparisons, each rendered as its own fixed-height row pair (not
# stacked per-group above each box, which gets uneven/cluttered once a
# group needs more than one comparison) -- like a compact results table
# overlaid on the plot, with a row label on the left (outside the
# categorical x range) identifying each comparison. Each cell is two lines
# (p, then AUC) -- a single combined line (format_p_auc(), used for Panel
# B) was too wide and overlapped neighbouring groups' text here. All text
# plain grey/black (no colour-coding); matched-controls row sits above the
# n.64_65insT row, both well clear of the y=1.00 gridline/tick label.
#   Row 1: every group vs matched controls (existing comparison).
#   Row 2: other T-loop / Stem III / n.72_73del vs n.64_65insT -- does the
#     founder-variant-trained signature generalise to these other RNU4-2
#     variant classes (as opposed to "is it elevated over background",
#     which row 1 answers).
# n.72_73del is n=2 (see diagnostic print above) -- p/AUC are well-defined
# but reflect only 2 carriers, not a powered test. BH-adjusted separately
# within each row (they're answering different questions, not one family).
panelE_groups <- setdiff(new_a_lev, cat_new_ctrl)
panelE_p   <- vapply(panelE_groups, function(g)
  wilcox.test(new_a_df$value[new_a_df$category == g], ctrl_prob_insT, exact = FALSE)$p.value,
  numeric(1L))
panelE_auc <- vapply(panelE_groups, function(g)
  auc_two_sample(new_a_df$value[new_a_df$category == g], ctrl_prob_insT), numeric(1L))
panelE_p_adj <- setNames(p.adjust(panelE_p, method = "BH"), panelE_groups)
panelE_auc   <- setNames(panelE_auc, panelE_groups)

for (g in panelE_groups)
  cat(sprintf("Panel E %s vs matched controls: p=%.3e (BH=%.3e) AUC=%.3f n=%d\n",
              g, panelE_p[[g]], panelE_p_adj[[g]], panelE_auc[[g]],
              sum(new_a_df$category == g)))

insT_case_vals_e <- new_a_df$value[new_a_df$category == cat_tloop_inst]
panelE_vs_insT_groups <- c(cat_tloop_ni, cat_stem3_only, cat_del72)
panelE_vs_insT_p   <- vapply(panelE_vs_insT_groups, function(g)
  wilcox.test(new_a_df$value[new_a_df$category == g], insT_case_vals_e, exact = FALSE)$p.value,
  numeric(1L))
panelE_vs_insT_auc <- vapply(panelE_vs_insT_groups, function(g)
  auc_two_sample(new_a_df$value[new_a_df$category == g], insT_case_vals_e), numeric(1L))
panelE_vs_insT_p_adj <- setNames(p.adjust(panelE_vs_insT_p, method = "BH"), panelE_vs_insT_groups)
panelE_vs_insT_auc   <- setNames(panelE_vs_insT_auc, panelE_vs_insT_groups)

for (g in panelE_vs_insT_groups)
  cat(sprintf("Panel E %s vs n.64_65insT: p=%.3e (BH=%.3e) AUC=%.3f n=%d\n",
              g, panelE_vs_insT_p[[g]], panelE_vs_insT_p_adj[[g]], panelE_vs_insT_auc[[g]],
              sum(new_a_df$category == g)))

cat("\nPanel E median P(case) per group:\n")
for (g in new_a_lev)
  cat(sprintf("  %s: median=%.3f n=%d\n", g,
              median(new_a_df$value[new_a_df$category == g], na.rm = TRUE),
              sum(new_a_df$category == g)))

cat("\nPanel E extremes:\n")
if (n_unmatched > 0L) {
  max_unmatched_p <- max(new_a_df$value[new_a_df$category == cat_unmatched], na.rm = TRUE)
  cat(sprintf("  Highest-scoring WB unmatched control: P(case)=%.3f\n", max_unmatched_p))
} else {
  cat("  No unmatched WB controls loaded -- skipping.\n")
}
min_insT_p <- min(new_a_df$value[new_a_df$category == cat_tloop_inst], na.rm = TRUE)
cat(sprintf("  Lowest-scoring n.64_65insT individual: P(case)=%.3f\n", min_insT_p))

# Two fixed-height comparison "rows", matched-controls on top (higher y),
# vs-n.64_65insT below it, well clear of the y=1.00 gridline/tick label.
# p-value only (AUC dropped from the plot -- still computed/logged to
# console above). All text plain black/grey (no colour-coding). Row labels
# sit left of the first category; needs clip="off" (added to the plot
# below) since that's outside the categorical x range.
Y_CTRL <- 1.20
Y_INST <- 1.08

format_p_sci_vec <- function(p) vapply(p, format_p_sci, character(1L))

panelE_row1_ann <- data.frame(
  category = factor(panelE_groups, levels = new_a_lev), y = Y_CTRL,
  label = format_p_sci_vec(panelE_p_adj[panelE_groups]),
  stringsAsFactors = FALSE)
panelE_row2_ann <- data.frame(
  category = factor(panelE_vs_insT_groups, levels = new_a_lev), y = Y_INST,
  label = format_p_sci_vec(panelE_vs_insT_p_adj[panelE_vs_insT_groups]),
  stringsAsFactors = FALSE)

panelE_row_labels <- data.frame(
  x = 0.25, y = c(Y_CTRL, Y_INST),
  label = c("vs matched\ncontrols", "vs n.64_65insT"),
  stringsAsFactors = FALSE)

fig2_panelE <- ggplot(new_a_df, aes(x = category, y = value, colour = category)) +
  geom_violin(data = function(d) d[d$category == cat_unmatched, ],
              aes(fill = category), colour = NA, alpha = 0.7, width = 0.7,
              trim = TRUE) +
  # Narrow boxplot overlay: the violin alone tapers to near-invisible width
  # in the sparse tails, hiding lone points. geom_boxplot's default 1.5*IQR
  # outlier rule draws those as distinct dots (mirrors figure3.R Panel A3).
  geom_boxplot(data = function(d) d[d$category == cat_unmatched, ],
               width = 0.08, outlier.size = 0.9, outlier.alpha = 0.8,
               colour = "grey30", fill = NA, linewidth = 0.3) +
  # Small groups (n typically <=3) get dots + median crossbar; larger
  # categorical groups keep a boxplot. IQR is unreliable at n<=3, and a
  # narrow crossbar is easier to read than a degenerate box.
  geom_boxplot(data = function(d) d[!d$category %in% c(cat_unmatched, cat_tloop_ni,
                                                       cat_del72, cat_biallelic), ],
               fill = NA, outlier.shape = NA,
               width = 0.35, linewidth = 0.3) +
  stat_summary(data = function(d) d[d$category %in% c(cat_tloop_ni, cat_del72,
                                                      cat_biallelic), ],
               fun = median, geom = "crossbar",
               width = 0.35, fatten = 0, linewidth = 0.3) +
  geom_jitter(data = function(d) d[d$category != cat_unmatched, ],
              width = 0.10, size = 1.1, alpha = 0.75) +
  geom_text(data = panelE_row1_ann, aes(x = category, y = y, label = label),
           parse = TRUE, size = 2.0, colour = "grey20", hjust = 0.5,
           inherit.aes = FALSE) +
  geom_text(data = panelE_row2_ann, aes(x = category, y = y, label = label),
           parse = TRUE, size = 2.0, colour = "grey20", hjust = 0.5,
           inherit.aes = FALSE) +
  geom_text(data = panelE_row_labels, aes(x = x, y = y, label = label),
           size = 2.0, colour = "grey30", hjust = 1, lineheight = 0.85,
           inherit.aes = FALSE) +
  scale_x_discrete(limits = new_a_lev, labels = x_labels_e) +
  scale_colour_manual(values = new_a_pal, name = NULL, guide = "none") +
  scale_fill_manual(values = new_a_pal, name = NULL, guide = "none") +
  scale_y_continuous(limits = c(0, 1.28), breaks = seq(0, 1, 0.25)) +
  coord_cartesian(clip = "off") +
  labs(x = NULL,
       y = "insT LASSO predicted P(case)") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(4, 4, 4, 24))

# ============================================================================
# NEW PANEL C: insT event annotation table
# 6 insT signature events; one row per event.
# Each row shows two tracks (long and short A5SS isoform) coloured by which
# isoform is upregulated in ReNU insT carriers. Gene name (italic) is shown
# above the track. A small nt-distance bracket sits below the tracks
# (asterisk = naive frameshift: alt region fully in CDS, length not a
# multiple of 3). A per-track transcript annotation tier label (MANE /
# Alternate / Unannotated, grey track above orange track) sits to the left
# of the diagram. Columns to the right: PSI ctrl/insT (%), then a filled/
# empty circle for NMD -- all aligned to the orange (insT-associated)
# track's row.
# ============================================================================

make_event_schematic <- function(annot_path, sel_ev_ids, psi_ctrl, psi_case, tloop_col) {
  if (!file.exists(annot_path)) {
    cat("make_event_schematic: annotation file not found:", annot_path, "\n")
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, size = 3, colour = "grey40",
                      label = "Event annotations\nnot available") +
             theme_void())
  }
  annot <- read.csv(annot_path, stringsAsFactors = FALSE)
  tbl   <- annot[match(sel_ev_ids, annot$event_id), , drop = FALSE]
  tbl   <- tbl[!is.na(tbl$event_id), , drop = FALSE]
  if (nrow(tbl) == 0L) {
    cat("make_event_schematic: no matching event_ids found in annotation file\n")
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, size = 3, colour = "grey40",
                      label = "No events matched") +
             theme_void())
  }

  # Per-event binary annotations
  # NMD, one per track: is_nmd = towards transcript (protein_consequence is
  # already computed from towards_biotype only); is_away_nmd = away
  # transcript, computed directly from away_biotype with the same rule
  # (mirrors the towards-side check in export_a5ss_candidate_full_table.R).
  tbl$is_nmd      <- tbl$protein_consequence == "NMD predicted"
  tbl$is_away_nmd <- !is.na(tbl$away_biotype) &
    tbl$away_biotype %in% c("nonsense_mediated_decay", "protein_coding_CDS_not_defined")
  # Transcript annotation tier, per track: MANE Select / other annotated
  # transcript ("Alternate") / no matching transcript ("Unannotated"). Shown
  # as a full-word label next to each track rather than a single Ann. circle,
  # since towards_tier/away_tier describe two different transcripts (the
  # orange insT-associated isoform and the grey other isoform respectively).
  tier_bucket <- function(tier) {
    ifelse(is.na(tier), "Unannotated",
           ifelse(tier == "MANE_Select", "MANE", "Alternate"))
  }
  tbl$towards_tier_label <- tier_bucket(tbl$towards_tier)
  tbl$away_tier_label    <- tier_bucket(tbl$away_tier)
  # UTR-only: use cds_class if present, else fall back to protein_consequence string
  tbl$is_utr <- if ("cds_class" %in% names(tbl)) {
    !is.na(tbl$cds_class) & tbl$cds_class %in% c("5UTR", "3UTR")
  } else {
    grepl("UTR-only", tbl$protein_consequence, ignore.case = TRUE)
  }
  # Naive frameshift: alt region fully within CDS and its length isn't a
  # multiple of 3 -- flagged with an asterisk on the nt-distance bracket.
  tbl$is_frameshift <- if ("cds_class" %in% names(tbl)) {
    !is.na(tbl$cds_class) & tbl$cds_class == "CDS" & (tbl$alt_len_nt %% 3L != 0L)
  } else {
    grepl("frameshift", tbl$protein_delta, ignore.case = TRUE)
  }

  n_ev <- nrow(tbl)

  # Layout constants (half-height relative to original design)
  EXON_W     <- 0.30
  ALTW_MIN   <- 0.06
  ALTW_MAX   <- 0.18
  INTRON_W   <- 0.20
  TRACK_H    <- 0.07   # base half-height; UTR-only events use TRACK_H / 2
  TRACK_SEP  <- 0.16
  EVENT_H    <- 0.50
  BRACK_DROP <- 0.12
  BRACK_TICK <- 0.04
  NAME_GAP   <- 0.05
  # Column spacing. The transcript-tier column and the two NMD circles
  # (away/towards) both sit to the LEFT of the diagrams (x < 0): tier text
  # (right-aligned, needs full-word room -- "Unannotated" is the longest
  # label) first, then the NMD circle column, then the diagram itself at
  # x=0. PSI columns stay on the right, using the tighter COL_SEP spacing.
  COL_SEP     <- 0.056   # column spacing for the right-side PSI columns
  DIAGRAM_GAP <- 0.05    # gap between the NMD circle column and the diagram
  TIER_GAP    <- 0.05    # gap between the tier-label column and the NMD circles
  TIER_W      <- 0.24    # width reserved for the tier label column
  CIRC_SZ     <- 2.0

  alts      <- tbl$alt_len_nt
  alt_range <- max(alts) - min(alts)
  alt_w_ev  <- if (alt_range < 1) rep((ALTW_MIN + ALTW_MAX) / 2, n_ev) else
    ALTW_MIN + (ALTW_MAX - ALTW_MIN) * (alts - min(alts)) / alt_range

  NMD_X       <- -DIAGRAM_GAP                  # NMD circle column, x < 0
  TIER_X      <- NMD_X - TIER_GAP              # right-aligned anchor, x < 0
  X_LIM_L     <- TIER_X - TIER_W
  MAX_X       <- EXON_W + ALTW_MAX + INTRON_W
  PSI_CTRL_X  <- MAX_X + 0.04
  PSI_CASE_X  <- PSI_CTRL_X + COL_SEP
  X_LIM_R     <- PSI_CASE_X + 0.05

  tbl$cds_suffix <- ""   # populated per-event in loop below

  rects    <- vector("list", n_ev * 2L)
  segs     <- vector("list", n_ev * 5L)
  gene_tx  <- vector("list", n_ev)
  cds_tx   <- vector("list", n_ev)
  brack_tx <- vector("list", n_ev)
  dpsi_tx  <- vector("list", n_ev)
  tier_tx  <- vector("list", n_ev)
  circles  <- vector("list", n_ev)
  r_i <- 0L; s_i <- 0L

  for (i in seq_len(n_ev)) {
    yc     <- -(i - 1L) * EVENT_H
    y_long <- yc + TRACK_SEP / 2
    y_shor <- yc - TRACK_SEP / 2
    aw     <- alt_w_ev[i]
    dir_i  <- tbl$direction[i]
    x_long_don <- EXON_W + aw
    x_right    <- EXON_W + aw + INTRON_W

    # ReNU-associated isoform always on the bottom track (y_shor), grey on top.
    # dir_i == 1: long form (wider) is ReNU-up → bottom=long, top=short
    # dir_i == -1: short form is ReNU-up         → bottom=short, top=long
    track_h_i    <- if (tbl$is_utr[i]) TRACK_H / 2 else TRACK_H
    bot_xmax     <- if (dir_i == 1L) x_long_don else EXON_W
    top_xmax     <- if (dir_i == 1L) EXON_W     else x_long_don
    bot_intron_x <- bot_xmax
    top_intron_x <- top_xmax

    r_i <- r_i + 1L
    rects[[r_i]] <- data.frame(
      xmin = 0, xmax = top_xmax,
      ymin = y_long - track_h_i, ymax = y_long + track_h_i,
      fill = "grey75", stringsAsFactors = FALSE)

    r_i <- r_i + 1L
    rects[[r_i]] <- data.frame(
      xmin = 0, xmax = bot_xmax,
      ymin = y_shor - track_h_i, ymax = y_shor + track_h_i,
      fill = tloop_col, stringsAsFactors = FALSE)

    s_i <- s_i + 1L
    segs[[s_i]] <- data.frame(x = top_intron_x, xend = x_right,
                               y = y_long, yend = y_long, stringsAsFactors = FALSE)
    s_i <- s_i + 1L
    segs[[s_i]] <- data.frame(x = bot_intron_x, xend = x_right,
                               y = y_shor, yend = y_shor, stringsAsFactors = FALSE)

    bx1 <- EXON_W; bx2 <- x_long_don
    by  <- y_shor - BRACK_DROP
    s_i <- s_i + 1L
    segs[[s_i]] <- data.frame(x = bx1, xend = bx2, y = by, yend = by,
                               stringsAsFactors = FALSE)
    s_i <- s_i + 1L
    segs[[s_i]] <- data.frame(x = bx1, xend = bx1, y = by, yend = by + BRACK_TICK,
                               stringsAsFactors = FALSE)
    s_i <- s_i + 1L
    segs[[s_i]] <- data.frame(x = bx2, xend = bx2, y = by, yend = by + BRACK_TICK,
                               stringsAsFactors = FALSE)

    brack_tx[[i]] <- data.frame(
      x = (bx1 + bx2) / 2, y = by - 0.04,
      label = sprintf("%s\n%d nt", if (tbl$is_frameshift[i]) "*" else " ",
                      tbl$alt_len_nt[i]),
      stringsAsFactors = FALSE)

    # Gene name (italic) + CDS class suffix (non-italic) above the long track
    label_y   <- y_long + track_h_i + NAME_GAP
    gene_tx[[i]] <- data.frame(x = 0, y = label_y,
                                label = tbl$gene_symbol[i],
                                stringsAsFactors = FALSE)
    cds_suffix <- switch(
      if ("cds_class" %in% names(tbl)) tbl$cds_class[i] else "",
      "CDS"         = "CDS",
      "5UTR"        = "5'UTR",
      "3UTR"        = "3'UTR",
      "partial_CDS" = "CDS/UTR",
      "")
    tbl$cds_suffix[i] <- cds_suffix
    # CDS suffix placed at a fixed offset to clear the longest gene name (~7 chars)
    cds_tx[[i]] <- data.frame(x = 0.17, y = label_y,
                               label = cds_suffix,
                               stringsAsFactors = FALSE)

    # PSI text columns (centred on each column x, same y as circle row).
    # psi_ctrl/psi_case are the raw rMATS long-form (inclusion) PSI, which
    # isn't always the orange (insT-associated) track -- when dir_i == -1
    # the short form is the insT-associated isoform, so flip to 1-PSI here
    # so the printed numbers always match the orange track's PSI.
    ev_id_i    <- tbl$event_id[i]
    psi_ctrl_v <- psi_ctrl[ev_id_i]
    psi_case_v <- psi_case[ev_id_i]
    if (dir_i == -1L) {
      psi_ctrl_v <- 1 - psi_ctrl_v
      psi_case_v <- 1 - psi_case_v
    }
    dpsi_tx[[i]] <- data.frame(
      x     = c(PSI_CTRL_X, PSI_CASE_X),
      y     = c(y_shor, y_shor),
      label = c(
        if (!is.na(psi_ctrl_v)) sprintf("%.0f%%", 100 * psi_ctrl_v) else "NA",
        if (!is.na(psi_case_v)) sprintf("%.0f%%", 100 * psi_case_v) else "NA"),
      stringsAsFactors = FALSE)

    # Two circles per event, lined up with the two tracks/tier labels:
    # away (top, y_long) and towards (bottom, y_shor).
    circles[[i]] <- data.frame(
      x    = c(NMD_X, NMD_X),
      y    = c(y_long, y_shor),
      fill = c(if (tbl$is_away_nmd[i]) "grey20" else "white",
              if (tbl$is_nmd[i])      "grey20" else "white"),
      stringsAsFactors = FALSE)

    # Transcript tier labels, one per track: grey (top, "away") and the
    # insT-associated (bottom, "towards") track. Sit to the left of the
    # diagrams, right-aligned; same neutral colour for both (not tied to
    # track colour).
    tier_tx[[i]] <- data.frame(
      x     = c(TIER_X, TIER_X),
      y     = c(y_long, y_shor),
      label = c(tbl$away_tier_label[i], tbl$towards_tier_label[i]),
      stringsAsFactors = FALSE)
  }

  rects_df  <- do.call(rbind, rects[!sapply(rects,   is.null)])
  segs_df   <- do.call(rbind, segs[!sapply(segs,    is.null)])
  gene_df   <- do.call(rbind, gene_tx[!sapply(gene_tx,  is.null)])
  cds_df    <- do.call(rbind, cds_tx[!sapply(cds_tx,   is.null)])
  brack_df  <- do.call(rbind, brack_tx[!sapply(brack_tx, is.null)])
  dpsi_df   <- do.call(rbind, dpsi_tx[!sapply(dpsi_tx,  is.null)])
  circ_df   <- do.call(rbind, circles[!sapply(circles,  is.null)])
  tier_df   <- do.call(rbind, tier_tx[!sapply(tier_tx,  is.null)])

  # Column headers, all rotated 90 degrees
  first_label_y <- 0 + TRACK_SEP / 2 + TRACK_H + NAME_GAP
  hdr_y  <- first_label_y + 0.08
  hdr_df <- data.frame(
    x     = c(NMD_X, PSI_CTRL_X, PSI_CASE_X),
    y     = hdr_y,
    label = c("NMD", "PSI ctrl", "PSI insT"),
    stringsAsFactors = FALSE)

  y_range_top <- hdr_y + 0.55
  y_range_bot <- -(n_ev - 1L) * EVENT_H - TRACK_SEP / 2 - BRACK_DROP - 0.12

  ggplot() +
    geom_rect(data = rects_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = "grey40", linewidth = 0.25, inherit.aes = FALSE) +
    scale_fill_identity() +
    geom_segment(data = segs_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 colour = "grey40", linewidth = 0.35, inherit.aes = FALSE) +
    geom_text(data = gene_df,
              aes(x = x, y = y, label = label),
              colour = "black", size = 2.4, hjust = 0, fontface = "italic",
              inherit.aes = FALSE) +
    geom_text(data = cds_df[nchar(cds_df$label) > 0, , drop = FALSE],
              aes(x = x, y = y, label = label),
              colour = "grey30", size = 2.4, hjust = 0, fontface = "plain",
              inherit.aes = FALSE) +
    geom_text(data = brack_df,
              aes(x = x, y = y, label = label),
              colour = "grey30", size = 2.0, hjust = 0.5, lineheight = 0.85,
              inherit.aes = FALSE) +
    geom_text(data = dpsi_df,
              aes(x = x, y = y, label = label),
              colour = "grey20", size = 2.0, hjust = 0.5, inherit.aes = FALSE) +
    geom_text(data = tier_df,
              aes(x = x, y = y, label = label),
              colour = "grey20", size = 2.0, hjust = 1, inherit.aes = FALSE) +
    geom_point(data = circ_df,
               aes(x = x, y = y, fill = fill),
               shape = 21, size = CIRC_SZ, colour = "grey40", stroke = 0.5,
               inherit.aes = FALSE) +
    geom_text(data = hdr_df,
              aes(x = x, y = y, label = label),
              colour = "grey20", size = 2.4, hjust = 0, angle = 90, fontface = "bold",
              inherit.aes = FALSE) +
    scale_x_continuous(limits = c(X_LIM_L, X_LIM_R), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_range_bot, y_range_top), expand = c(0, 0)) +
    theme_void(base_size = 9) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# Order events by descending |LASSO beta| (highest-weight site first)
beta_order_idx   <- order(abs(insT_beta), decreasing = TRUE)
insT_ev_ordered  <- insT_gel_ev[beta_order_idx]
cat(sprintf("Panel C event order (by |beta|): %s\n",
            paste(seq_along(insT_ev_ordered), insT_ev_ordered, sep = "=", collapse = ", ")))

# Per-site mean PSI (raw, not delta) in WB controls vs WB n.64_65insT
# carriers. Uses raw PSI from psi_matrix.rds rather than rMATS IncLevelDiff
# (all-case average). WB only -- lymphocyte column dropped, this panel now
# rests on whole blood.
compute_carrier_psi_f2 <- function(psi_mat, sample_ids, ev_ids) {
  ev_ids     <- ev_ids[ev_ids %in% colnames(psi_mat)]
  sample_ids <- sample_ids[sample_ids %in% rownames(psi_mat)]
  if (length(sample_ids) == 0L)
    return(setNames(rep(NA_real_, length(ev_ids)), ev_ids))
  apply(psi_mat[sample_ids, ev_ids, drop = FALSE], 2L, median, na.rm = TRUE)
}

insT_gel_carrier_ids <- per_case$sample_id[per_case$is_inst]

insT_psi_ctrl_gel <- compute_carrier_psi_f2(gel$psi, gel_ctrl_ids,         insT_gel_ev)
insT_psi_case_gel <- compute_carrier_psi_f2(gel$psi, insT_gel_carrier_ids, insT_gel_ev)

fig2_panelC <- tryCatch(
  make_event_schematic(ANNOT_PATH, insT_ev_ordered,
                       psi_ctrl  = insT_psi_ctrl_gel,
                       psi_case  = insT_psi_case_gel,
                       tloop_col = PAL_TLOOP),
  error = function(e) {
    cat("Panel C schematic error:", conditionMessage(e), "\n")
    patchwork::plot_spacer()
  })

# ============================================================================
# NEW PANEL D: STX16 + RPS15A log2(TPM+1) boxplots -- cases vs controls,
# faceted (STX16 top, RPS15A bottom), sharing x-axis text (facet_grid rows
# only labels the bottom-most facet's x-axis by default). STX16 is the only
# predicted frameshift event among the 6 insT signature events; RPS15A is
# also significant and shown alongside it (moved out of figure_s5.R Panel A
# to sit here instead). Two Mann-Whitney U p-values, BH-adjusted jointly across
# both genes. TPM loaded from per-sample Salmon quant.genes.sf files.
# ============================================================================
compute_gene_tpm_f2 <- function(quant_paths_csv, ensg_target, sample_ids, label) {
  if (!file.exists(quant_paths_csv)) {
    cat("compute_gene_tpm_f2: quant paths CSV not found:", quant_paths_csv, "\n")
    return(NULL)
  }
  qp <- read.csv(quant_paths_csv, stringsAsFactors = FALSE)
  stopifnot("quant_paths_csv must have columns platekey and file_path" =
              all(c("platekey", "file_path") %in% colnames(qp)))
  qp <- qp[qp$platekey %in% sample_ids, , drop = FALSE]
  if (nrow(qp) == 0L) {
    cat("compute_gene_tpm_f2: no matching samples found for label:", label, "\n")
    return(NULL)
  }
  res <- do.call(rbind, lapply(seq_len(nrow(qp)), function(i) {
    sid <- qp$platekey[i]
    sf  <- qp$file_path[i]
    if (!file.exists(sf)) {
      cat("  missing:", sf, "\n")
      return(NULL)
    }
    d <- read.delim(sf, stringsAsFactors = FALSE)
    ensg_base <- sub("\\.\\d+$", "", d$Name)
    hit <- d[ensg_base == ensg_target, , drop = FALSE]
    if (nrow(hit) == 0L) return(NULL)
    tpm <- sum(hit$TPM, na.rm = TRUE)
    data.frame(sample_id = sid, tpm = tpm, stringsAsFactors = FALSE)
  }))
  if (is.null(res) || nrow(res) == 0L) return(NULL)
  res$group <- label
  res
}

# Find the STX16 + RPS15A events in the insT annotation to confirm ENSGs
fig2_panelD <- patchwork::plot_spacer()
tryCatch({
  if (file.exists(ANNOT_PATH)) {
    annot_c <- read.csv(ANNOT_PATH, stringsAsFactors = FALSE)
    insT_annot <- annot_c[match(insT_gel_ev, annot_c$event_id), , drop = FALSE]
    panelD_genes <- c("STX16", "RPS15A")
    gene_rows <- lapply(panelD_genes, function(g)
      insT_annot[!is.na(insT_annot$gene_symbol) &
                   toupper(insT_annot$gene_symbol) == g, , drop = FALSE])
    names(gene_rows) <- panelD_genes
    missing_genes <- panelD_genes[vapply(gene_rows, nrow, integer(1L)) == 0L]
    if (length(missing_genes) > 0L)
      cat(sprintf("Panel D: %s not found among insT events -- omitted\n",
                  paste(missing_genes, collapse = ", ")))
    panelD_genes <- setdiff(panelD_genes, missing_genes)

    # insT cases only (exclude Stem III and other variant carriers)
    inst_gel_ids_c <- per_case$sample_id[per_case$is_inst]
    cat(sprintf("Panel D: n_insT_cases = %d\n", length(inst_gel_ids_c)))

    gene_tpm_df <- setNames(vector("list", length(panelD_genes)), panelD_genes)
    gene_p_raw  <- setNames(rep(NA_real_, length(panelD_genes)), panelD_genes)
    for (g in panelD_genes) {
      ensg <- gene_rows[[g]]$ensg[1L]
      cat(sprintf("Panel D: %s ENSG = %s\n", g, ensg))
      tpm_case <- compute_gene_tpm_f2(QUANT_PATHS_CSV, ensg, inst_gel_ids_c, "WB n.64_65insT")
      tpm_ctrl <- compute_gene_tpm_f2(QUANT_PATHS_CSV, ensg, gel_ctrl_ids,   "WB control")
      if (is.null(tpm_case) || is.null(tpm_ctrl)) {
        cat(sprintf("Panel D: %s TPM data unavailable -- omitted\n", g))
        next
      }
      # Naive sanity check (not a model): does the raw % TPM change roughly
      # match the raw % PSI increase of the towards (insT-associated)
      # isoform, in the direction implied by WHICH SIDE is the NMD
      # transcript? No efficiency term, no per-sample pairing -- just two
      # point estimates compared directly. Runs for every gene in
      # panelD_genes, not just STX16.
      #
      # Two distinct scenarios, not interchangeable:
      #   towards == NMD (e.g. STX16): cases shift INTO the NMD isoform ->
      #     more transcript degraded -> expect TPM to DECREASE in cases.
      #   away == NMD (e.g. RPS15A): cases shift OUT of the NMD isoform
      #     (towards the normal one) -> less transcript degraded -> expect
      #     TPM to INCREASE in cases. Getting this backwards (as an earlier
      #     version of this check did) makes an increase look like an
      #     unexplained "decrease".
      g_ev_id       <- gene_rows[[g]]$event_id[1L]
      g_dir         <- gene_rows[[g]]$direction[1L]
      towards_biotype_g <- gene_rows[[g]]$towards_biotype[1L]
      away_biotype_g    <- gene_rows[[g]]$away_biotype[1L]
      nmd_biotypes  <- c("nonsense_mediated_decay", "protein_coding_CDS_not_defined")
      towards_is_nmd <- !is.na(towards_biotype_g) && towards_biotype_g %in% nmd_biotypes
      away_is_nmd    <- !is.na(away_biotype_g)    && away_biotype_g    %in% nmd_biotypes

      psi_ctrl_v <- insT_psi_ctrl_gel[g_ev_id]
      psi_case_v <- insT_psi_case_gel[g_ev_id]
      if (g_dir == -1) { psi_ctrl_v <- 1 - psi_ctrl_v; psi_case_v <- 1 - psi_case_v }
      delta_psi_pp <- 100 * (psi_case_v - psi_ctrl_v)   # +ve = towards-isoform PSI increase
      tpm_pc       <- 100 * (mean(tpm_case$tpm, na.rm = TRUE) /
                             mean(tpm_ctrl$tpm, na.rm = TRUE) - 1)  # +ve = TPM increase

      if (towards_is_nmd && !away_is_nmd) {
        cat(sprintf(
          "Panel D naive check (%s): towards-isoform (NMD) PSI %.1f%% -> %.1f%% (%+.1fpp); observed TPM change %+.1f%% (naive 1:1 expectation: ~%.1f%% decrease)\n",
          g, 100 * psi_ctrl_v, 100 * psi_case_v, delta_psi_pp, tpm_pc, -delta_psi_pp))
      } else if (away_is_nmd && !towards_is_nmd) {
        cat(sprintf(
          "Panel D naive check (%s): towards-isoform (non-NMD) PSI %.1f%% -> %.1f%% (%+.1fpp, shifting AWAY from the NMD isoform); observed TPM change %+.1f%% (naive 1:1 expectation: ~%.1f%% increase)\n",
          g, 100 * psi_ctrl_v, 100 * psi_case_v, delta_psi_pp, tpm_pc, delta_psi_pp))
      } else {
        cat(sprintf(
          "Panel D naive check (%s): skipped -- towards/away NMD status not a clean single side (towards_is_nmd=%s, away_is_nmd=%s)\n",
          g, towards_is_nmd, away_is_nmd))
      }

      tpm_df <- rbind(tpm_case, tpm_ctrl)
      tpm_df$log2tpm <- log2(tpm_df$tpm + 1)
      tpm_df$group   <- factor(tpm_df$group, levels = c("WB control", "WB n.64_65insT"))
      tpm_df$gene    <- g
      wt <- wilcox.test(tpm_df$log2tpm[tpm_df$group == "WB n.64_65insT"],
                        tpm_df$log2tpm[tpm_df$group == "WB control"], exact = FALSE)
      gene_p_raw[g]  <- wt$p.value
      gene_tpm_df[[g]] <- tpm_df
      cat(sprintf("Panel D: %s n_insT_case=%d n_ctrl=%d Mann-Whitney U p=%.3e\n",
                  g, sum(tpm_df$group == "WB n.64_65insT"),
                  sum(tpm_df$group == "WB control"), wt$p.value))

      # Raw-TPM effect-size reporting (paper-facing numbers). Median-based:
      # TPM is right-skewed and n is small, so median is more robust than mean.
      # 95% CI on the median ratio via bootstrap over paired resampling of case
      # and control TPMs (2000 reps), reported as % change = 100 * (ratio - 1).
      case_tpm <- tpm_case$tpm
      ctrl_tpm <- tpm_ctrl$tpm
      med_ctrl <- median(ctrl_tpm, na.rm = TRUE)
      med_case <- median(case_tpm,  na.rm = TRUE)
      tpm_pc_median <- 100 * (med_case / med_ctrl - 1)

      set.seed(20260915L)
      B <- 2000L
      boot_pc <- replicate(B, {
        cb <- sample(case_tpm, length(case_tpm), replace = TRUE)
        rb <- sample(ctrl_tpm, length(ctrl_tpm), replace = TRUE)
        100 * (median(cb, na.rm = TRUE) / median(rb, na.rm = TRUE) - 1)
      })
      ci_lo <- unname(quantile(boot_pc, 0.025, na.rm = TRUE))
      ci_hi <- unname(quantile(boot_pc, 0.975, na.rm = TRUE))
      direction_word <- if (tpm_pc_median < 0) "reduction" else "increase"

      # Per-sample PSI ranges (towards-isoform, oriented as in the naive
      # check above). Reports controls and insT cases separately.
      psi_ctrl_vec <- gel$psi[gel_ctrl_ids,   g_ev_id]
      psi_case_vec <- gel$psi[inst_gel_ids_c, g_ev_id]
      if (g_dir == -1) {
        psi_ctrl_vec <- 1 - psi_ctrl_vec
        psi_case_vec <- 1 - psi_case_vec
      }
      psi_ctrl_vec <- psi_ctrl_vec[is.finite(psi_ctrl_vec)]
      psi_case_vec <- psi_case_vec[is.finite(psi_case_vec)]
      cat(sprintf(
        "Panel D: %s towards-isoform PSI range -- controls [%.1f%%, %.1f%%] (median %.1f%%, n=%d); insT cases [%.1f%%, %.1f%%] (median %.1f%%, n=%d)\n",
        g,
        100 * min(psi_ctrl_vec),  100 * max(psi_ctrl_vec),
        100 * median(psi_ctrl_vec), length(psi_ctrl_vec),
        100 * min(psi_case_vec),  100 * max(psi_case_vec),
        100 * median(psi_case_vec), length(psi_case_vec)))

      cat(sprintf(
        "Panel D: %s raw TPM range -- controls [%.2f, %.2f] (median %.2f); insT cases [%.2f, %.2f] (median %.2f)\n",
        g,
        min(ctrl_tpm, na.rm = TRUE), max(ctrl_tpm, na.rm = TRUE), med_ctrl,
        min(case_tpm, na.rm = TRUE), max(case_tpm, na.rm = TRUE), med_case))
      cat(sprintf(
        "Panel D: %s median raw TPM %s: %+.1f%% (95%% CI %+.1f%% to %+.1f%%; ctrl %.2f -> case %.2f; Mann-Whitney U P=%.3e)\n",
        g, direction_word, tpm_pc_median, ci_lo, ci_hi,
        med_ctrl, med_case, wt$p.value))
    }

    panelD_genes <- panelD_genes[!is.na(gene_p_raw[panelD_genes])]
    if (length(panelD_genes) == 0L) {
      cat("Panel D: no genes with usable TPM data -- panel shows spacer\n")
    } else {
      tpm_df   <- do.call(rbind, gene_tpm_df[panelD_genes])
      tpm_df$gene <- factor(tpm_df$gene, levels = c("STX16", "RPS15A"))  # top, bottom

      p_adj <- p.adjust(gene_p_raw[panelD_genes], method = "BH")
      cat(sprintf("Panel D BH-adjusted p: %s\n",
                  paste(sprintf("%s=%.3e", panelD_genes, p_adj), collapse = ", ")))
      ann_df <- do.call(rbind, lapply(panelD_genes, function(g) {
        d <- tpm_df[tpm_df$gene == g, , drop = FALSE]
        data.frame(gene = g, y = max(d$tpm, na.rm = TRUE) * 1.05,
                  label = format_p_sci(p_adj[[g]]), stringsAsFactors = FALSE)
      }))
      ann_df$gene <- factor(ann_df$gene, levels = c("STX16", "RPS15A"))

      tpm_pal_c <- c("WB control" = PAL_CTRL, "WB n.64_65insT" = PAL_TLOOP)
      fig2_panelD <- ggplot(tpm_df, aes(x = group, y = tpm, colour = group)) +
        geom_boxplot(fill = NA, outlier.shape = NA, width = 0.45, linewidth = 0.35) +
        geom_jitter(width = 0.12, size = 1.2, alpha = 0.75) +
        geom_segment(data = ann_df, aes(x = 1, xend = 2, y = y, yend = y),
                    colour = "grey30", linewidth = 0.35, inherit.aes = FALSE) +
        geom_text(data = ann_df, aes(x = 1.5, y = y * 1.02, label = label),
                 parse = TRUE, size = 2.4, colour = "grey20", hjust = 0.5,
                 inherit.aes = FALSE) +
        facet_grid(rows = vars(gene), scales = "free_y") +
        scale_colour_manual(values = tpm_pal_c, guide = "none") +
        # Display-only relabel (line break after "WB") so the x-axis text
        # can sit horizontal instead of angled -- underlying group values
        # (used above for the Wilcoxon test / palette matching) are unchanged.
        scale_x_discrete(labels = c("WB control"     = "WB\ncontrol",
                                    "WB n.64_65insT" = "WB\nn.64_65insT")) +
        labs(x = NULL, y = "TPM") +
        pub_theme +
        theme(strip.text.y = element_text(size = 9, face = "italic"),
              axis.text.x = element_text(angle = 0, hjust = 0.5))
    }
  } else {
    cat("Panel D: ANNOT_PATH not found -- panel omitted\n")
  }
}, error = function(e) {
  cat("Panel D error:", conditionMessage(e), "\n")
})

# ---- Combine into Figure 2 PDF ---------------------------------------------
mm_to_in <- function(mm) mm / 25.4

# 5-panel layout, 3 rows x 2 panels (row 3 is a single full-width panel):
#   Row 1: A (insT coef bar, 20%) | B (insT P(case) WB/Lymph, 80%)
#   Row 2: C (A5SS event schematics, 66%) | D (STX16 TPM, 33%)
#   Row 3: E (P(case) by RNU4-2 variant group + unmatched controls)  full width
row1_f2 <- (fig2_panelA + labs(tag = "A")) |
           (fig2_panelB + labs(tag = "B"))
row1_f2 <- row1_f2 + plot_layout(widths = c(1, 4))
# Panel C is theme_void() with no y-axis; wrap_elements() stops patchwork
# from reserving left-margin space to align it with axis-bearing neighbours
# (was causing left whitespace + vertical squish).
row2_f2 <- (wrap_elements(full = fig2_panelC) + labs(tag = "C")) |
           (fig2_panelD + labs(tag = "D"))
row2_f2 <- row2_f2 + plot_layout(widths = c(2, 1))
row3_f2 <- fig2_panelE + labs(tag = "E")

fig2_combined <- row1_f2 / row2_f2 / row3_f2 +
  plot_layout(heights = c(0.65, 0.85, 0.8))

pdf(OUT_FIG_2, width = mm_to_in(170), height = mm_to_in(200))
print(fig2_combined & theme(plot.tag = element_text(size = 14, face = "bold"),
                            plot.tag.position = "topleft"))
dev.off()
cat(sprintf("Figure 2 saved to: %s\n", OUT_FIG_2))

# ---- Summary RDS (small, airlock-safe) -------------------------------------
# Counts + aggregate stats only. No event_ids, gene symbols, or coef values.
summary_out <- list(
  K_insT       = K_insT,
  # Panel E: insT LASSO group counts
  panel_e_n_ctrl          = length(c(gel_ctrl_ids, nava_ctrl_ids)),
  panel_e_n_unmatched     = n_unmatched,
  panel_e_n_insT_tloop    = sum(per_case$new_group == cat_tloop_inst, na.rm = TRUE),
  panel_e_n_noninsT_tloop = sum(per_case$new_group == cat_tloop_ni,   na.rm = TRUE),
  panel_e_n_stem3         = sum(!is.na(per_case$region_bucket) &
                                  per_case$region_bucket == "Stem III", na.rm = TRUE),
  panel_e_n_biallelic     = length(biallelic_in_pool))
saveRDS(summary_out, OUT_SUMMARY)
cat(sprintf("Summary saved to: %s\n", OUT_SUMMARY))

cat("\n=== Done. ===\n")
cat(sprintf("Panel E: %d insT T-loop, %d non-insT T-loop, %d Stem III, %d biallelic, %d unmatched controls\n",
            summary_out$panel_e_n_insT_tloop, summary_out$panel_e_n_noninsT_tloop,
            summary_out$panel_e_n_stem3, summary_out$panel_e_n_biallelic,
            summary_out$panel_e_n_unmatched))
cat("\nNOTE: GEL-derived figure stays in workspace unless airlocked.\n")
