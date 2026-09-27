#!/usr/bin/env Rscript
# Layout: patchwork only -- never gridExtra arrangeGrob/grid.arrange
# =============================================================================
# figure3.R
# Paper Figure 3 (RNU4-2 splicing biomarker manuscript): development of the
# Stem III signature.
#
# Promoted from figure_s3.R to a main paper figure 2026-07-26 (was
# Supplementary Figure S3). Panel D (heatmap of direction-oriented
# bb-logit-PSI deviation across all insT + stemIII signature sites) moves
# in from figure2.R's old panel D. Panels A-C otherwise unchanged from the
# 2026-07-06 restructure: all pre-restructure per-site characterisation
# panels (formerly A-J, none of which fed the output PDF) were deleted
# 2026-07-21, along with the donor-seq/TPM/NDD-enrichment computations that
# only ever fed those archived panels or the summary RDS.
#
#   Panel A - Stem III LASSO coefficient bar (by coef magnitude). Moved from
#             figure2.R panel B/C/D restructure.
#   Panel B - Stem III event annotation table: one row per Stem III-only
#             signature event. WB-only mean PSI (ctrl vs Stem III carriers,
#             %, oriented to the orange Stem III-associated track) | per-track
#             transcript annotation tier (MANE/Alternate/Unannotated, left of
#             the diagram) | Protein consequence | NMD. Analogous to
#             figure2.R Panel C but for the Stem III-specific events (events
#             in stemIII set not shared with insT set).
#   Panel C - Stem III LASSO P(case) boxplot by variant group (matched
#             controls / n.64_65insT / other T-loop / Stem III / n.72_73del /
#             biallelic RNU4-2), with the ~6k unmatched whole-blood controls
#             folded in as an extra violin category (no jitter) rather than
#             a separate panel -- matches figure2.R Panel E's category
#             set/order/labels/palette.
#   Panel D - Heatmap of direction-oriented bb-logit-PSI deviation across
#             all insT + stemIII signature sites, for all carriers; rows =
#             sites (T-loop | Shared | Stem III), cols = cases (T-loop |
#             Stem III | unassigned). Moved from figure2.R's old panel D.
#
# Primary inputs:
#   outputs/loo_results_insT_lasso.rds
#   outputs/loo_results_stemIII_lasso.rds
#   outputs/psi_matrix.rds
#   outputs/nava_psi_matrix.rds
#   outputs/gel_controls_psi_cand_pool.rds
#   outputs/a5ss_candidate_full_table.csv
#   data/biallelic_rnu42_platekeys.txt
#
# Outputs:
#   outputs/figures/figure3.pdf
#   outputs/figure3_summary.rds   (counts + stats; airlock-safe)
#
# Run from workspace root: pixi run Rscript scripts/gel/figure3.R
#
# DATA PROTECTION: event_ids, gene symbols, and LASSO coefs are GEL protected.
# The figure PDF + summary.rds carry only aggregate stats and are airlock-safe.
# =============================================================================

# ---- 0. libpath setup (pixi-aware) -----------------------------------------
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
  if (dir.exists(user_lib)) {
    .libPaths(c(user_lib, .libPaths()))
  }
}))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(pheatmap)
  library(RColorBrewer)
  library(gridExtra)
})

# ---- 1. Configuration ------------------------------------------------------
LOO_PATH        <- Sys.getenv("LOO_RDS",         "outputs/loo_results_insT_lasso.rds")
STEMIII_LOO_PATH <- Sys.getenv("STEMIII_LOO_RDS","outputs/loo_results_stemIII_lasso.rds")
GEL_PSI_PATH    <- Sys.getenv("GEL_PSI",         "outputs/psi_matrix.rds")
NAVA_PSI_PATH   <- Sys.getenv("NAVA_PSI",        "outputs/nava_psi_matrix.rds")
WB_RMATS_DIR    <- Sys.getenv("GEL_RMATS_DIR",   "data/rmats")
ANNOT_PATH       <- Sys.getenv("ANNOT_CSV",
  "outputs/a5ss_candidate_full_table.csv")
CAND_POOL_PATH   <- Sys.getenv("CAND_POOL",      "outputs/gel_controls_psi_cand_pool.rds")
BIALLELIC_FILE   <- Sys.getenv("BIALLELIC_FILE", "data/biallelic_rnu42_platekeys.txt")
STEMIII_ROB_PATH <- Sys.getenv("STEMIII_ROB_RDS", "outputs/gel_controls_stemIII_robustness.rds")
GENO_FILE        <- Sys.getenv("GENO_CSV",         "data/sample_genotype.csv")
NAVA_SEV_FILE    <- Sys.getenv("NAVA_SEV",         "outputs/nava_severity_labels.csv")

OUT_DIR_FIG <- "outputs/figures"
OUT_FIG_3   <- file.path(OUT_DIR_FIG, "figure3.pdf")
OUT_SUMMARY <- "outputs/figure3_summary.rds"

EVENT_TYPES   <- c("SE", "A5SS", "A3SS", "MXE", "RI")
FDR_THRESHOLD <- 0.05

WB_LAB <- "whole blood"
LY_LAB <- "lymphocytes"
COHORT_LABELS <- c(GEL = WB_LAB, Nava = LY_LAB)

# Palette (mirrors figure1.R / figure2.R)
PAL_WB      <- "#2C5AA0"
PAL_LY      <- "#2AA198"
PAL_INC     <- "#D55E00"   # increased-usage (insT signature colour)
PAL_DEC     <- "#0072B2"   # decreased-usage
PAL_BG      <- "#BBBBBB"   # grey background
PAL_CTRL    <- "#888888"
PAL_CASE    <- "#D55E00"
PAL_STEMIII <- "#F4A582"   # Stem III colour (matches figure2.R PAL_STEM3)
PAL_TLOOP   <- "#D55E00"   # T-loop colour  (matches figure2.R PAL_TLOOP)
PAL_SEVERE  <- "#5E3C99"   # Severe phenotype (matches figure2.R)
PAL_MILD    <- "#B2ABD2"   # Mild phenotype   (matches figure2.R)
PAL_SHARED  <- "#7B2D8B"   # shared insT+stemIII events

dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

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

placeholder_panel <- function(label, msg) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 3, colour = "grey40") +
    theme_void() +
    theme(panel.border = element_rect(fill = NA, colour = "grey70",
                                      linewidth = 0.3),
          plot.margin  = margin(2, 2, 2, 2))
}

# ============================================================================
# ---- Stem III classifier setup (Panels A, B, C) ----------------------------
# Loads PSI matrices + insT/stemIII LOO objects at script top level.
# ============================================================================

INST_HGVS   <- "n.64_65insT"
HGVS_DEL72  <- "n.72_73del"
RNU4_TLOOP  <- c(63L, 70L)
RNU4_STEM3  <- c(75L, 79L)

classify_rnu4_region_s3 <- function(hgvs) {
  out <- rep(NA_character_, length(hgvs))
  for (i in seq_along(hgvs)) {
    h <- hgvs[i]
    if (is.na(h) || !nzchar(h)) next
    coord <- regmatches(h, regexpr("n\\.[0-9]+(_[0-9]+)?", h))
    if (length(coord) == 0L) next
    pos <- as.integer(regmatches(coord, gregexpr("[0-9]+", coord))[[1]])
    if (length(pos) == 0L) next
    if (any(pos >= RNU4_TLOOP[1] & pos <= RNU4_TLOOP[2]))
      out[i] <- "T-loop"
    else if (any(pos >= RNU4_STEM3[1] & pos <= RNU4_STEM3[2]))
      out[i] <- "Stem III"
  }
  out
}

bb_logit_s3 <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}
sigmoid_s3  <- function(z) 1 / (1 + exp(-z))

# Load PSI matrices
gel_s3  <- readRDS(GEL_PSI_PATH)
nava_s3 <- readRDS(NAVA_PSI_PATH)

gel_case_ids_s3  <- gel_s3$samples$sample_id[gel_s3$samples$group  == "Case"]
gel_ctrl_ids_s3  <- gel_s3$samples$sample_id[gel_s3$samples$group  == "Control"]
nava_case_ids_s3 <- nava_s3$samples$sample_id[nava_s3$samples$group == "Case"]
nava_ctrl_ids_s3 <- nava_s3$samples$sample_id[nava_s3$samples$group == "Control"]

# Load insT LOO object
insT_loo_s3 <- readRDS(LOO_PATH)
insT_cand_s3    <- insT_loo_s3$a5ss_cand
insT_sel_idx_s3 <- insT_loo_s3$selected_idx
K_insT_s3       <- insT_loo_s3$K_selected
insT_gel_ev_s3  <- insT_cand_s3$gel_ev_id[insT_sel_idx_s3]
insT_nava_ev_s3 <- insT_cand_s3$nava_ev_id[insT_sel_idx_s3]
insT_dir_s3     <- insT_cand_s3$direction[insT_sel_idx_s3]
insT_mu_gel_s3  <- insT_loo_s3$gel_ctrl_means[insT_gel_ev_s3]
insT_mu_nava_s3 <- insT_loo_s3$nava_ctrl_means[insT_gel_ev_s3]
insT_beta_s3    <- insT_loo_s3$fit_at_1se$beta[insT_gel_ev_s3]
insT_icept_s3   <- insT_loo_s3$fit_at_1se$intercept

# Load stemIII LOO object
stemIII_loo_s3 <- readRDS(STEMIII_LOO_PATH)
stemIII_cand_s3    <- stemIII_loo_s3$a5ss_cand
stemIII_sel_idx_s3 <- stemIII_loo_s3$selected_idx
K_stemIII_s3       <- stemIII_loo_s3$K_selected
stemIII_sel_ev_s3  <- stemIII_cand_s3$gel_ev_id[stemIII_sel_idx_s3]
stemIII_dir_s3     <- stemIII_cand_s3$direction[stemIII_sel_idx_s3]
stemIII_gene_s3    <- stemIII_cand_s3$GeneSymbol[stemIII_sel_idx_s3]
stemIII_mu_gel_s3  <- stemIII_loo_s3$gel_ctrl_means[stemIII_sel_ev_s3]
stemIII_mu_nava_s3 <- stemIII_loo_s3$nava_ctrl_means[stemIII_sel_ev_s3]
stemIII_beta_s3    <- unname(stemIII_loo_s3$fit_at_1se$beta[stemIII_sel_ev_s3])
stemIII_icept_s3   <- stemIII_loo_s3$fit_at_1se$intercept
cat(sprintf("Stem III classifier: K=%d selected sites\n", K_stemIII_s3))

# Cand pool + biallelic keys
ctrl_pool_s3 <- if (file.exists(CAND_POOL_PATH)) readRDS(CAND_POOL_PATH) else NULL
biallelic_keys_s3 <- if (file.exists(BIALLELIC_FILE)) {
  k <- trimws(readLines(BIALLELIC_FILE)); k[nzchar(k)]
} else character(0L)
biallelic_in_pool_s3 <- if (!is.null(ctrl_pool_s3) && length(biallelic_keys_s3) > 0L)
  intersect(biallelic_keys_s3, rownames(ctrl_pool_s3$ijc)) else character(0L)
cat(sprintf("Biallelic: %d/%d keys in cand pool\n",
            length(biallelic_in_pool_s3), length(biallelic_keys_s3)))

# Build HGVS lookup for all cases from authoritative genotype files
# (insT_loo per_case_meta is insT-only; this covers all variant groups)
all_hgvs_s3 <- setNames(character(0), character(0))
if (file.exists(GENO_FILE)) {
  geno_s3 <- read.csv(GENO_FILE, stringsAsFactors = FALSE)
  if (all(c("sample_id", "hgvs") %in% names(geno_s3)))
    all_hgvs_s3 <- c(all_hgvs_s3, setNames(geno_s3$hgvs, geno_s3$sample_id))
  cat(sprintf("Genotype file loaded: %d GEL entries\n", nrow(geno_s3)))
} else {
  cat("WARNING: GENO_FILE not found -- non-insT groups will be empty\n")
}
if (file.exists(NAVA_SEV_FILE)) {
  nava_lbl_s3 <- local({
    first <- readLines(NAVA_SEV_FILE, n = 1L)
    sep   <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
    read.table(NAVA_SEV_FILE, header = TRUE, sep = sep, stringsAsFactors = FALSE,
               quote = "\"", fill = TRUE)
  })
  nava_lbl_s3$sample_id <- paste0("nava:case:", nava_lbl_s3$sample_n)
  all_hgvs_s3 <- c(all_hgvs_s3, setNames(nava_lbl_s3$hgvs, nava_lbl_s3$sample_id))
  cat(sprintf("Nava severity file loaded: %d Nava entries\n", nrow(nava_lbl_s3)))
} else {
  cat("WARNING: NAVA_SEV_FILE not found -- Nava non-insT groups will be empty\n")
}

# Shared HGVS genotype palette (Set3, frequency-ranked) -- used by Panel D
# (heatmap) so each genotype gets a consistent colour.
hgvs_all_vals_s3   <- all_hgvs_s3[!is.na(all_hgvs_s3)]
hgvs_levels_all_s3 <- names(sort(table(hgvs_all_vals_s3), decreasing = TRUE))
n_hgvs_all_s3      <- length(hgvs_levels_all_s3)
brewer_cols_all_s3  <- RColorBrewer::brewer.pal(max(3L, min(n_hgvs_all_s3, 12L)), "Set3")
PAL_HGVS_SHARED <- if (n_hgvs_all_s3 <= 12L) {
  setNames(brewer_cols_all_s3[seq_len(n_hgvs_all_s3)], hgvs_levels_all_s3)
} else {
  setNames(colorRampPalette(brewer_cols_all_s3)(n_hgvs_all_s3), hgvs_levels_all_s3)
}

# All-cases data frame with HGVS and region
all_cases_s3 <- rbind(
  data.frame(sample_id = gel_case_ids_s3,  cohort = "GEL",  stringsAsFactors = FALSE),
  data.frame(sample_id = nava_case_ids_s3, cohort = "Nava", stringsAsFactors = FALSE)
)
all_cases_s3$hgvs          <- all_hgvs_s3[all_cases_s3$sample_id]
all_cases_s3$region_bucket <- classify_rnu4_region_s3(all_cases_s3$hgvs)
cat("\nAll-cases HGVS region distribution:\n")
print(table(cohort = all_cases_s3$cohort,
            region = factor(all_cases_s3$region_bucket,
                            levels = c("T-loop", "Stem III")), useNA = "always"))

# insT P(case) for all samples
insT_gel_x_s3  <- bb_logit_s3(
  gel_s3$ijc[,  insT_gel_ev_s3,  drop = FALSE],
  gel_s3$sjc[,  insT_gel_ev_s3,  drop = FALSE])
insT_nava_x_s3 <- bb_logit_s3(
  nava_s3$ijc[, insT_nava_ev_s3, drop = FALSE],
  nava_s3$sjc[, insT_nava_ev_s3, drop = FALSE])
colnames(insT_nava_x_s3) <- insT_gel_ev_s3

insT_orient_gel_s3  <- sweep(sweep(insT_gel_x_s3,  2, insT_mu_gel_s3,  "-"), 2, insT_dir_s3, "*")
insT_orient_nava_s3 <- sweep(sweep(insT_nava_x_s3, 2, insT_mu_nava_s3, "-"), 2, insT_dir_s3, "*")

# stemIII linear predictor for all samples
sIII_gel_x_s3  <- bb_logit_s3(
  gel_s3$ijc[,  stemIII_sel_ev_s3, drop = FALSE],
  gel_s3$sjc[,  stemIII_sel_ev_s3, drop = FALSE])
sIII_nava_x_s3 <- bb_logit_s3(
  nava_s3$ijc[, stemIII_sel_ev_s3, drop = FALSE],
  nava_s3$sjc[, stemIII_sel_ev_s3, drop = FALSE])

sIII_orient_gel_s3  <- sweep(sweep(sIII_gel_x_s3,  2, stemIII_mu_gel_s3,  "-"), 2, stemIII_dir_s3, "*")
sIII_orient_nava_s3 <- sweep(sweep(sIII_nava_x_s3, 2, stemIII_mu_nava_s3, "-"), 2, stemIII_dir_s3, "*")

sIII_prob_gel_s3  <- sigmoid_s3(as.numeric(sIII_orient_gel_s3  %*% stemIII_beta_s3) + stemIII_icept_s3)
sIII_prob_nava_s3 <- sigmoid_s3(as.numeric(sIII_orient_nava_s3 %*% stemIII_beta_s3) + stemIII_icept_s3)
names(sIII_prob_gel_s3)  <- rownames(sIII_orient_gel_s3)
names(sIII_prob_nava_s3) <- rownames(sIII_orient_nava_s3)

# Biallelic stemIII P(case)
bial_sIII_prob_s3 <- numeric(0)
if (!is.null(ctrl_pool_s3) && length(biallelic_in_pool_s3) > 0L &&
    all(stemIII_sel_ev_s3 %in% colnames(ctrl_pool_s3$ijc))) {
  bx <- bb_logit_s3(ctrl_pool_s3$ijc[biallelic_in_pool_s3, stemIII_sel_ev_s3, drop = FALSE],
                    ctrl_pool_s3$sjc[biallelic_in_pool_s3, stemIII_sel_ev_s3, drop = FALSE])
  bo <- sweep(sweep(bx, 2, stemIII_mu_gel_s3, "-"), 2, stemIII_dir_s3, "*")
  bial_sIII_prob_s3 <- sigmoid_s3(as.numeric(bo %*% stemIII_beta_s3) + stemIII_icept_s3)
  cat(sprintf("Biallelic stemIII P(case) range [%.3f, %.3f]\n",
              min(bial_sIII_prob_s3), max(bial_sIII_prob_s3)))
}

# ============================================================================
# ---- 4c. Panel A: 2x2 composite of old figure 2 panels A-D ---------------
# Row 1: insT LASSO P(case) by variant group | Stem III LASSO coef bar
# Row 2: Stem III LASSO P(case) by variant group | ~6k WB controls violin
# ============================================================================

# A1: Coefficient bar
coef_s3_val     <- unname(stemIII_loo_s3$fit_at_1se$beta[stemIII_sel_ev_s3])
sel_ord_s3      <- order(-coef_s3_val)
gene_s3_ord     <- stemIII_gene_s3[sel_ord_s3]
coef_s3_ord     <- coef_s3_val[sel_ord_s3]
if (anyDuplicated(gene_s3_ord)) gene_s3_ord <- make.unique(gene_s3_ord, sep = ".")

# Signature membership per site (same classification as Panel D's heatmap
# row annotation): "Shared" = also selected by the insT model, "Stem III" =
# Stem III-only. Coloured with Panel D's Signature palette (PAL_WB/PAL_STEMIII).
coef_s3_sig <- ifelse(stemIII_sel_ev_s3[sel_ord_s3] %in% insT_gel_ev_s3,
                     "Shared", "Stem III-only")

coef_s3_df <- data.frame(
  gene_lab  = factor(gene_s3_ord, levels = gene_s3_ord),
  coef      = coef_s3_ord,
  signature = factor(coef_s3_sig, levels = c("Shared", "Stem III-only")),
  stringsAsFactors = FALSE)

fig_s3_A1 <- ggplot(coef_s3_df, aes(x = gene_lab, y = coef, fill = signature)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = c(Shared = PAL_WB, `Stem III-only` = PAL_STEMIII),
                    name = "Signature") +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "LASSO coef (lambda.1se)") +
  pub_theme +
  theme(axis.text.x    = element_text(size = 7, face = "italic", angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text     = element_text(size = 9),
        plot.margin     = margin(12, 4, 12, 4))

# A2: P(case) by group (+ unmatched controls, merged in like figure2.R
# Panel E -- same category set/order/labels/palette preference).
cat_unmatched_s3  <- "WB unmatched controls"
cat_new_ctrl_s3   <- "Matched controls"
cat_tloop_inst_s3 <- "n.64_65insT"
cat_tloop_ni_s3   <- "other T-loop"
cat_stem3_only_s3 <- "Stem III"
cat_del72_s3      <- HGVS_DEL72
cat_biallelic_s3  <- "Biallelic RNU4-2"
shared_grp_lev_s3 <- c(cat_new_ctrl_s3, cat_unmatched_s3, cat_tloop_inst_s3, cat_tloop_ni_s3,
                        cat_stem3_only_s3, cat_del72_s3, cat_biallelic_s3)

get_sIII_prob_s3 <- function(df)
  ifelse(df$cohort == "GEL", sIII_prob_gel_s3[df$sample_id], sIII_prob_nava_s3[df$sample_id])

inst_tloop_s3    <- all_cases_s3[!is.na(all_cases_s3$hgvs) &
                                   all_cases_s3$hgvs == INST_HGVS, ]
noninst_tloop_s3 <- all_cases_s3[!is.na(all_cases_s3$region_bucket) &
                                    all_cases_s3$region_bucket == "T-loop" &
                                    (is.na(all_cases_s3$hgvs) |
                                     all_cases_s3$hgvs != INST_HGVS), ]
stem3_pc_s3      <- all_cases_s3[!is.na(all_cases_s3$region_bucket) &
                                    all_cases_s3$region_bucket == "Stem III", ]
del_pc_s3        <- all_cases_s3[!is.na(all_cases_s3$hgvs) &
                                    all_cases_s3$hgvs == HGVS_DEL72, ]

ctrl_prob_s3   <- c(sIII_prob_gel_s3[gel_ctrl_ids_s3], sIII_prob_nava_s3[nava_ctrl_ids_s3])

a2_df <- data.frame(
  prob  = c(ctrl_prob_s3,
            get_sIII_prob_s3(inst_tloop_s3),
            get_sIII_prob_s3(noninst_tloop_s3),
            get_sIII_prob_s3(stem3_pc_s3),
            get_sIII_prob_s3(del_pc_s3),
            bial_sIII_prob_s3),
  group = c(
    rep(cat_new_ctrl_s3,   length(ctrl_prob_s3)),
    rep(cat_tloop_inst_s3, nrow(inst_tloop_s3)),
    rep(cat_tloop_ni_s3,   nrow(noninst_tloop_s3)),
    rep(cat_stem3_only_s3, nrow(stem3_pc_s3)),
    rep(cat_del72_s3,      nrow(del_pc_s3)),
    rep(cat_biallelic_s3,  length(bial_sIII_prob_s3))),
  stringsAsFactors = FALSE)

cat("\nStem III P(case) group counts:\n")
print(table(group = a2_df$group, useNA = "always"))

# A3: ~6k WB controls violin (stemIII model)
# Requires pre-computed robustness RDS -- no fallback scoring path.
if (!file.exists(STEMIII_ROB_PATH)) {
  stop(sprintf("Panel A3: STEMIII_ROB_PATH not found (%s) -- run score_gel_controls_stemIII.R first",
               STEMIII_ROB_PATH))
}
sIII_rob_s3 <- readRDS(STEMIII_ROB_PATH)
prob_vec_s3 <- if ("ctrl_prob_hist" %in% names(sIII_rob_s3))
  sIII_rob_s3$ctrl_prob_hist else sIII_rob_s3$ctrl_prob
n_ctrl_s3a3 <- length(prob_vec_s3)
cat(sprintf("Panel A: %d unmatched WB controls loaded from %s\n", n_ctrl_s3a3, STEMIII_ROB_PATH))
n_exceed_s3a3 <- sum(prob_vec_s3 > 0.5, na.rm = TRUE)
cat(sprintf("Panel A: %d / %d unmatched controls exceed Stem III LASSO P(case) 0.5\n",
            n_exceed_s3a3, n_ctrl_s3a3))
# Tukey-rule outlier count (matches geom_boxplot's default outlier draw):
# any point outside [Q1 - 1.5 IQR, Q3 + 1.5 IQR].
ub_vals_s3 <- prob_vec_s3[is.finite(prob_vec_s3)]
ub_q_s3    <- quantile(ub_vals_s3, c(0.25, 0.75), na.rm = TRUE)
ub_iqr_s3  <- ub_q_s3[2L] - ub_q_s3[1L]
ub_lo_s3   <- ub_q_s3[1L] - 1.5 * ub_iqr_s3
ub_hi_s3   <- ub_q_s3[2L] + 1.5 * ub_iqr_s3
n_out_lo_s3 <- sum(ub_vals_s3 < ub_lo_s3)
n_out_hi_s3 <- sum(ub_vals_s3 > ub_hi_s3)
cat(sprintf(
  "Panel C: unmatched-controls Tukey outliers: %d total (%d below %.3f, %d above %.3f; IQR=[%.3f, %.3f])\n",
  n_out_lo_s3 + n_out_hi_s3, n_out_lo_s3, ub_lo_s3, n_out_hi_s3, ub_hi_s3,
  ub_q_s3[1L], ub_q_s3[2L]))

# Merge unmatched controls into a2_df as an extra category (matches
# figure2.R Panel E: violin, no jitter -- too many points -- rather than a
# separate panel).
a2_df <- rbind(a2_df, data.frame(
  prob = unname(prob_vec_s3), group = cat_unmatched_s3, stringsAsFactors = FALSE))
a2_df$group <- factor(a2_df$group, levels = shared_grp_lev_s3)

c_pal_s3 <- c(
  "WB unmatched controls" = PAL_WB,
  "Matched controls"      = PAL_CTRL,
  "n.64_65insT"           = PAL_TLOOP,
  "other T-loop"          = PAL_TLOOP,
  "Stem III"              = PAL_STEMIII,
  "n.72_73del"            = if (!is.null(PAL_HGVS_SHARED) && "n.72_73del" %in% names(PAL_HGVS_SHARED))
                              PAL_HGVS_SHARED[["n.72_73del"]] else PAL_LY,
  "Biallelic RNU4-2"      = PAL_SEVERE)

x_labels_s3a2 <- setNames(shared_grp_lev_s3, shared_grp_lev_s3)
x_labels_s3a2[cat_new_ctrl_s3]  <- "Matched\ncontrols"
x_labels_s3a2[cat_biallelic_s3] <- "Biallelic\nRNU4-2"
x_labels_s3a2[cat_unmatched_s3] <- sprintf("WB controls\n(n=%s)",
                                           format(n_ctrl_s3a3, big.mark = ","))

# Mann-Whitney U p-value (AUC computed/logged but not plotted -- matches
# figure2.R Panel E's final design) for two comparisons, each its own
# fixed-height row (not stacked per-group, which gets uneven once a group
# needs more than one comparison): row 1 = every group vs matched controls;
# row 2 = n.64_65insT / other T-loop / n.72_73del vs Stem III (does the
# Stem III-trained signature generalise to these other RNU4-2 variant
# classes -- Stem III itself is the reference here, not insT, since this
# is the Stem III LASSO panel). BH-adjusted separately within each row
# (different questions, not one test family). Plain text, no colour-coding.
# Row labels are a SEPARATE narrow ggplot patchworked to the left of the
# main plot (not geom_text at a negative x with clip="off") -- the two
# plots share an identical y scale so patchwork aligns their panel grobs
# exactly, which is what actually keeps the labels level with the p-value
# rows regardless of how many x categories the main plot has. A same-plot
# clip="off" annotation instead measures its offset in data units, which
# get squeezed as more categories are added, so the same nominal margin
# no longer reaches (or overshoots) the intended position.
auc_two_sample <- function(case_scores, ctrl_scores) {
  case_scores <- case_scores[!is.na(case_scores)]
  ctrl_scores <- ctrl_scores[!is.na(ctrl_scores)]
  if (length(case_scores) == 0L || length(ctrl_scores) == 0L) return(NA_real_)
  r <- rank(c(case_scores, ctrl_scores))
  U <- sum(r[seq_along(case_scores)]) - length(case_scores) * (length(case_scores) + 1L) / 2
  U / (length(case_scores) * length(ctrl_scores))
}
format_p_sci <- function(p) {
  if (p >= 0.001) return(sprintf("p == %.3f", p))
  exp10  <- floor(log10(p))
  mant_r <- round(p / 10^exp10, 1)
  if (mant_r >= 10) { mant_r <- mant_r / 10; exp10 <- exp10 + 1 }
  sprintf("p == %.1f %%*%% 10^%d", mant_r, exp10)
}
format_p_sci_vec_s3 <- function(p) vapply(p, format_p_sci, character(1L))

panelC_groups_s3 <- setdiff(shared_grp_lev_s3, cat_new_ctrl_s3)
panelC_p_s3   <- vapply(panelC_groups_s3, function(g)
  wilcox.test(a2_df$prob[a2_df$group == g], ctrl_prob_s3, exact = FALSE)$p.value, numeric(1L))
panelC_auc_s3 <- vapply(panelC_groups_s3, function(g)
  auc_two_sample(a2_df$prob[a2_df$group == g], ctrl_prob_s3), numeric(1L))
panelC_p_adj_s3 <- setNames(p.adjust(panelC_p_s3, method = "BH"), panelC_groups_s3)
panelC_auc_s3    <- setNames(panelC_auc_s3, panelC_groups_s3)

for (g in panelC_groups_s3)
  cat(sprintf("Panel C %s vs matched controls: p=%.3e (BH=%.3e) AUC=%.3f n=%d\n",
              g, panelC_p_s3[[g]], panelC_p_adj_s3[[g]], panelC_auc_s3[[g]],
              sum(a2_df$group == g)))

stem3_vals_s3 <- a2_df$prob[a2_df$group == cat_stem3_only_s3]
panelC_vs_stem3_groups_s3 <- c(cat_tloop_inst_s3, cat_tloop_ni_s3, cat_del72_s3)
panelC_vs_stem3_p_s3   <- vapply(panelC_vs_stem3_groups_s3, function(g)
  wilcox.test(a2_df$prob[a2_df$group == g], stem3_vals_s3, exact = FALSE)$p.value, numeric(1L))
panelC_vs_stem3_auc_s3 <- vapply(panelC_vs_stem3_groups_s3, function(g)
  auc_two_sample(a2_df$prob[a2_df$group == g], stem3_vals_s3), numeric(1L))
panelC_vs_stem3_p_adj_s3 <- setNames(p.adjust(panelC_vs_stem3_p_s3, method = "BH"),
                                     panelC_vs_stem3_groups_s3)
panelC_vs_stem3_auc_s3   <- setNames(panelC_vs_stem3_auc_s3, panelC_vs_stem3_groups_s3)

for (g in panelC_vs_stem3_groups_s3)
  cat(sprintf("Panel C %s vs Stem III: p=%.3e (BH=%.3e) AUC=%.3f n=%d\n",
              g, panelC_vs_stem3_p_s3[[g]], panelC_vs_stem3_p_adj_s3[[g]],
              panelC_vs_stem3_auc_s3[[g]], sum(a2_df$group == g)))

cat("\nPanel C median P(case) per group:\n")
for (g in shared_grp_lev_s3)
  cat(sprintf("  %s: median=%.3f n=%d\n", g,
              median(a2_df$prob[a2_df$group == g], na.rm = TRUE),
              sum(a2_df$group == g)))

Y_CTRL_S3  <- 1.24
Y_STEM3_S3 <- 1.08

panelC_row1_ann_s3 <- data.frame(
  group = factor(panelC_groups_s3, levels = shared_grp_lev_s3), y = Y_CTRL_S3,
  label = format_p_sci_vec_s3(panelC_p_adj_s3[panelC_groups_s3]),
  stringsAsFactors = FALSE)
panelC_row2_ann_s3 <- data.frame(
  group = factor(panelC_vs_stem3_groups_s3, levels = shared_grp_lev_s3), y = Y_STEM3_S3,
  label = format_p_sci_vec_s3(panelC_vs_stem3_p_adj_s3[panelC_vs_stem3_groups_s3]),
  stringsAsFactors = FALSE)

PANELC_Y_LIM <- c(0, 1.32)
PANELC_Y_BRK <- seq(0, 1, 0.25)

# Row labels ("vs matched controls" / "vs Stem III") deliberately omitted --
# to be added manually in Illustrator. A same-plot clip="off" annotation
# broke alignment as category count changed, and a separate patchworked
# label column broke alignment between A and C row heights; not worth
# fighting further for two lines of static text.
fig_s3_A2 <- ggplot(a2_df, aes(x = group, y = prob, colour = group)) +
  geom_violin(data = function(d) d[d$group == cat_unmatched_s3, ],
              aes(fill = group), colour = NA, alpha = 0.7, width = 0.7, trim = TRUE) +
  # Narrow boxplot overlay: the violin alone tapers to near-invisible width in
  # the sparse tails, hiding lone points. geom_boxplot's default 1.5*IQR
  # outlier rule draws those as distinct dots.
  geom_boxplot(data = function(d) d[d$group == cat_unmatched_s3, ],
               width = 0.08, outlier.size = 0.9, outlier.alpha = 0.8,
               colour = "grey30", fill = NA, linewidth = 0.3) +
  # Small groups (n typically <=3) get dots + median crossbar; larger
  # categorical groups keep a boxplot. IQR is unreliable at n<=3, and a
  # narrow crossbar is easier to read than a degenerate box.
  geom_boxplot(data = function(d) d[!d$group %in% c(cat_unmatched_s3, cat_tloop_ni_s3,
                                                     cat_del72_s3, cat_biallelic_s3), ],
               fill = NA, outlier.shape = NA, width = 0.35, linewidth = 0.3) +
  stat_summary(data = function(d) d[d$group %in% c(cat_tloop_ni_s3, cat_del72_s3,
                                                    cat_biallelic_s3), ],
               fun = median, geom = "crossbar",
               width = 0.35, fatten = 0, linewidth = 0.3) +
  geom_jitter(data = function(d) d[d$group != cat_unmatched_s3, ],
              width = 0.10, size = 0.85, alpha = 0.75) +
  geom_text(data = panelC_row1_ann_s3, aes(x = group, y = y, label = label),
           parse = TRUE, size = 2.2, colour = "grey20", hjust = 0.5,
           inherit.aes = FALSE) +
  scale_x_discrete(limits = shared_grp_lev_s3, labels = x_labels_s3a2) +
  scale_colour_manual(values = c_pal_s3, guide = "none") +
  scale_fill_manual(values = c_pal_s3, guide = "none") +
  scale_y_continuous(limits = PANELC_Y_LIM, breaks = PANELC_Y_BRK) +
  labs(x = NULL, y = "Stem III LASSO P(case)") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(4, 4, 4, 4))

# ============================================================================
# ---- 4d. New Panel B: Stem III-only transcript schematic ------------------
# Same make_event_schematic as figure2.R Panel B; PAL_STEMIII for cases.
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
  # orange Stem III-associated isoform and the grey other isoform respectively).
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
    # isn't always the orange (Stem III-associated) track -- when dir_i == -1
    # the short form is the Stem III-associated isoform, so flip to 1-PSI
    # here so the printed numbers always match the orange track's PSI.
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
    # Stem III-associated (bottom, "towards") track. Sit to the left of the
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
    label = c("NMD", "PSI ctrl", "PSI stemIII"),
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

# stemIII-only event IDs = stemIII set minus insT set
s3only_mask_s3   <- !stemIII_sel_ev_s3 %in% insT_gel_ev_s3
s3only_gel_ev_s3 <- stemIII_sel_ev_s3[s3only_mask_s3]
cat(sprintf("Stem III-only events: %d / %d (excluding %d shared with insT)\n",
            length(s3only_gel_ev_s3), K_stemIII_s3,
            sum(!s3only_mask_s3)))

# Order Stem III-only events by descending |beta| (highest-weight first)
stemIII_nava_ev_s3  <- stemIII_cand_s3$nava_ev_id[stemIII_sel_idx_s3]
s3only_beta_val_s3  <- unname(stemIII_loo_s3$fit_at_1se$beta[s3only_gel_ev_s3])
s3only_beta_ord_s3  <- order(abs(s3only_beta_val_s3), decreasing = TRUE)
s3only_ev_ordered_s3     <- s3only_gel_ev_s3[s3only_beta_ord_s3]
s3only_nava_ev_ordered_s3 <- stemIII_nava_ev_s3[s3only_mask_s3][s3only_beta_ord_s3]

# Per-site mean PSI (raw, not delta) in WB controls vs WB Stem III carriers.
# Uses raw PSI from psi_matrix.rds, not rMATS IncLevelDiff (all-case
# average). WB only, matching figure2.R Panel C's rework.
compute_carrier_psi_s3 <- function(psi_mat, sample_ids, ev_ids) {
  ev_ids     <- ev_ids[ev_ids %in% colnames(psi_mat)]
  sample_ids <- sample_ids[sample_ids %in% rownames(psi_mat)]
  if (length(sample_ids) == 0L)
    return(setNames(rep(NA_real_, length(ev_ids)), ev_ids))
  apply(psi_mat[sample_ids, ev_ids, drop = FALSE], 2L, median, na.rm = TRUE)
}

stem3_gel_ids <- stem3_pc_s3$sample_id[stem3_pc_s3$cohort == "GEL"]
cat(sprintf("Panel B PSI: GEL Stem III n=%d\n", length(stem3_gel_ids)))

s3_psi_ctrl_gel <- compute_carrier_psi_s3(gel_s3$psi, gel_ctrl_ids_s3,   s3only_ev_ordered_s3)
s3_psi_case_gel <- compute_carrier_psi_s3(gel_s3$psi, stem3_gel_ids,     s3only_ev_ordered_s3)

fig_s3_panelB <- tryCatch(
  make_event_schematic(ANNOT_PATH, s3only_ev_ordered_s3,
                       psi_ctrl  = s3_psi_ctrl_gel,
                       psi_case  = s3_psi_case_gel,
                       tloop_col = PAL_STEMIII),
  error = function(e) {
    cat("Panel B schematic error:", conditionMessage(e), "\n")
    patchwork::plot_spacer()
  })

# ============================================================================
# ---- 9b. Panel D: heatmap of direction-oriented bb-logit-PSI deviation ----
# across all insT + stemIII-only signature sites, for all carriers (GEL +
# Nava). Moved from figure2.R's old panel D.
#   Rows  = sites  (T-loop | Shared | Stem III, fixed order)
#   Cols  = cases  (T-loop | Stem III | unassigned, within-group ward.D2)
# ============================================================================
fig_s3_panelD <- patchwork::plot_spacer()
case_region_s3 <- factor(character(0), levels = c("T-loop", "Stem III", "unassigned"))
tryCatch({
  # stemIII-only = events in the stemIII set not in the insT set
  shared_ev_s3     <- intersect(insT_gel_ev_s3, stemIII_sel_ev_s3)
  s3_only_mu_gel   <- stemIII_mu_gel_s3[s3only_mask_s3]
  s3_only_mu_nava  <- stemIII_mu_nava_s3[s3only_mask_s3]
  s3_only_dir      <- stemIII_dir_s3[s3only_mask_s3]
  s3_only_gel_ev   <- s3only_gel_ev_s3
  s3_only_nava_ev  <- stemIII_nava_ev_s3[s3only_mask_s3]
  K_s3_only        <- sum(s3only_mask_s3)
  cat(sprintf("Panel D: K_insT=%d, K_s3_only=%d -> %d combined sites\n",
              K_insT_s3, K_s3_only, K_insT_s3 + K_s3_only))

  stopifnot(
    "s3-only gel events missing from gel ijc"   = all(s3_only_gel_ev  %in% colnames(gel_s3$ijc)),
    "s3-only nava events missing from nava ijc" = all(s3_only_nava_ev %in% colnames(nava_s3$ijc)))

  s3only_gel_x  <- bb_logit_s3(gel_s3$ijc[,  s3_only_gel_ev,  drop = FALSE],
                                gel_s3$sjc[,  s3_only_gel_ev,  drop = FALSE])
  s3only_nava_x <- bb_logit_s3(nava_s3$ijc[, s3_only_nava_ev, drop = FALSE],
                                nava_s3$sjc[, s3_only_nava_ev, drop = FALSE])
  colnames(s3only_nava_x) <- s3_only_gel_ev  # canonical gel IDs

  orient_s3only_gel  <- sweep(sweep(s3only_gel_x,  2, s3_only_mu_gel,  "-"), 2, s3_only_dir, "*")
  orient_s3only_nava <- sweep(sweep(s3only_nava_x, 2, s3_only_mu_nava, "-"), 2, s3_only_dir, "*")

  # Gene labels for sites: insT first, then stemIII-only
  insT_gene_hm       <- insT_cand_s3$GeneSymbol[insT_sel_idx_s3]
  s3only_gene_hm     <- stemIII_gene_s3[s3only_mask_s3]
  all_gene_raw       <- c(insT_gene_hm, s3only_gene_hm)
  all_gene_lab       <- make.unique(all_gene_raw, sep = ".")

  # Combined orient: insT sites | stemIII-only sites, cases only
  # orient_all is [n_cases x n_sites]  (rows=cases, cols=sites)
  orient_gel_hm  <- cbind(
    insT_orient_gel_s3[gel_case_ids_s3,    , drop = FALSE],
    orient_s3only_gel[gel_case_ids_s3,  , drop = FALSE])
  orient_nava_hm <- cbind(
    insT_orient_nava_s3[nava_case_ids_s3,  , drop = FALSE],
    orient_s3only_nava[nava_case_ids_s3, , drop = FALSE])
  colnames(orient_gel_hm)  <- all_gene_lab
  colnames(orient_nava_hm) <- all_gene_lab
  orient_all <- rbind(orient_gel_hm, orient_nava_hm)
  all_case_ids <- c(gel_case_ids_s3, nava_case_ids_s3)
  stopifnot(identical(rownames(orient_all), all_case_ids))

  # Drop cases with any NA orient (never impute)
  n_na <- rowSums(is.na(orient_all))
  if (any(n_na > 0L)) {
    dropped <- names(n_na)[n_na > 0L]
    warning(sprintf("Panel D: dropping %d case(s) with NA orient: %s",
                    length(dropped), paste(dropped, collapse = ", ")))
    orient_all <- orient_all[n_na == 0L, , drop = FALSE]
  }

  # Per-site max-abs normalisation (each column scaled to [-1, 1])
  col_maxabs <- apply(orient_all, 2L, function(x) max(abs(x), na.rm = TRUE))
  col_maxabs[col_maxabs == 0] <- 1
  orient_norm <- sweep(orient_all, 2L, col_maxabs, "/")

  # Site signature annotation (rows of orient_plot after transpose)
  site_sig_vals <- c(
    ifelse(insT_gel_ev_s3 %in% shared_ev_s3, "Shared", "T-loop"),  # K_insT insT sites
    rep("Stem III", K_s3_only))                                     # stemIII-only sites

  # orient_plot: rows = sites, cols = cases  (transpose of orient_norm)
  orient_plot <- t(orient_norm)
  rownames(orient_plot) <- all_gene_lab

  # Ann_site: per-row (site) annotation
  ann_site <- data.frame(
    Signature = factor(site_sig_vals, levels = c("T-loop", "Shared", "Stem III")),
    row.names = all_gene_lab,
    stringsAsFactors = FALSE)

  # Apply fixed row order: T-loop, Shared, Stem III
  site_group_order <- order(match(site_sig_vals, c("T-loop", "Shared", "Stem III")))
  # Reverse row order within the "Shared" group only (cosmetic display
  # choice; identified by group membership, not by gene name, since
  # LASSO-selected gene identifiers are GEL protected and never hardcoded).
  shared_pos_hm <- which(site_sig_vals[site_group_order] == "Shared")
  site_group_order[shared_pos_hm] <- rev(site_group_order[shared_pos_hm])
  orient_plot <- orient_plot[site_group_order, , drop = FALSE]
  ann_site    <- ann_site[site_group_order, , drop = FALSE]

  # HGVS lookup from stable CSV sources (loaded at top of script)
  case_ids_hm <- colnames(orient_plot)
  ann_hgvs <- if (!is.null(all_hgvs_s3)) {
    unname(all_hgvs_s3[case_ids_hm])
  } else {
    rep(NA_character_, length(case_ids_hm))
  }
  ann_tissue <- unname(COHORT_LABELS[ifelse(case_ids_hm %in% gel_case_ids_s3, "GEL", "Nava")])
  ann_region_raw <- classify_rnu4_region_s3(ann_hgvs)
  ann_region_raw[is.na(ann_region_raw)] <- "unassigned"

  tissue_levels <- intersect(c(WB_LAB, LY_LAB), ann_tissue)
  region_levels <- intersect(c("Stem III", "T-loop", "unassigned"), ann_region_raw)
  hgvs_levels   <- names(sort(table(ann_hgvs), decreasing = TRUE))

  ann_row_hm <- data.frame(
    Tissue   = factor(ann_tissue,     levels = tissue_levels),
    Region   = factor(ann_region_raw, levels = region_levels),
    Genotype = factor(ann_hgvs,       levels = hgvs_levels),
    row.names = case_ids_hm,
    stringsAsFactors = FALSE, check.names = FALSE)

  # Within-group column clustering: T-loop | Stem III | unassigned. Each
  # block is clustered independently, so there is no single tree spanning
  # all cases -- cluster_within_hm returns both the leaf order (used to fix
  # the column order, as before) and the hclust object itself, so the three
  # block-level trees can be drawn as an overlay above the heatmap later.
  cluster_within_hm <- function(ids, mat) {
    if (length(ids) <= 1L) return(list(order = ids, hc = NULL))
    hc <- hclust(dist(mat[ids, , drop = FALSE], method = "euclidean"), method = "ward.D2")
    list(order = ids[hc$order], hc = hc)
  }
  region_lookup_hm <- setNames(as.character(ann_row_hm$Region), rownames(ann_row_hm))
  tloop_ids  <- case_ids_hm[region_lookup_hm[case_ids_hm] == "T-loop"]
  stem3_ids  <- case_ids_hm[region_lookup_hm[case_ids_hm] == "Stem III"]
  other_ids  <- case_ids_hm[!region_lookup_hm[case_ids_hm] %in% c("T-loop", "Stem III")]
  # cluster_within_hm takes case IDs and the case-oriented matrix (t(orient_plot))
  tloop_res  <- cluster_within_hm(tloop_ids, t(orient_plot))
  stem3_res  <- cluster_within_hm(stem3_ids, t(orient_plot))
  other_res  <- cluster_within_hm(other_ids, t(orient_plot))
  tloop_ord  <- tloop_res$order
  stem3_ord  <- stem3_res$order
  other_ord  <- other_res$order
  col_order  <- c(tloop_ord, stem3_ord, other_ord)
  gaps_col_pos <- cumsum(c(length(tloop_ord), length(stem3_ord)))
  gaps_col_pos <- gaps_col_pos[gaps_col_pos < length(col_order)]
  orient_plot  <- orient_plot[, col_order, drop = FALSE]
  case_region_s3 <<- factor(region_lookup_hm[col_order],
                             levels = c("T-loop", "Stem III", "unassigned"))

  cat(sprintf("Column blocks: %d T-loop | %d Stem III | %d unassigned\n",
              length(tloop_ord), length(stem3_ord), length(other_ord)))

  # Colour palettes
  PAL_TISSUE_HM <- setNames(c(PAL_WB, PAL_LY)[seq_along(tissue_levels)], tissue_levels)
  PAL_REGION_HM <- c(`Stem III`   = PAL_STEMIII,
                     `T-loop`     = PAL_TLOOP,
                     unassigned   = "grey85")[region_levels]
  # Use the shared HGVS palette (computed early) so this heatmap and Panel A
  # (P(case) by variant group) match exactly.
  PAL_HGVS_HM <- if (!is.null(PAL_HGVS_SHARED)) {
    PAL_HGVS_SHARED[intersect(hgvs_levels, names(PAL_HGVS_SHARED))]
  } else {
    n_hgvs_hm      <- length(hgvs_levels)
    brewer_cols_hm  <- RColorBrewer::brewer.pal(max(3L, min(n_hgvs_hm, 12L)), "Set3")
    if (n_hgvs_hm <= 12L) setNames(brewer_cols_hm[seq_len(n_hgvs_hm)], hgvs_levels) else
      setNames(colorRampPalette(brewer_cols_hm)(n_hgvs_hm), hgvs_levels)
  }
  PAL_SIGNATURE_HM <- c(`T-loop` = PAL_TLOOP, Shared = PAL_WB, `Stem III` = PAL_STEMIII)
  ann_colors_hm <- list(Tissue    = PAL_TISSUE_HM,
                        Region    = PAL_REGION_HM,
                        Genotype  = PAL_HGVS_HM,
                        Signature = PAL_SIGNATURE_HM)

  # Italic gene-symbol row labels
  clean_label_hm <- function(x) sub("\\.[0-9]+$", "", x)
  gene_labels_hm <- as.expression(lapply(
    rownames(orient_plot),
    function(x) bquote(italic(.(clean_label_hm(x))))))

  # Colour scale: #1f78b4 -> white -> #E64B35
  breaks_v <- seq(-1, 1, length.out = 102)
  colors_v <- colorRampPalette(c("#1f78b4", "white", "#E64B35"))(101)

  ph <- pheatmap::pheatmap(
    orient_plot,
    cluster_rows         = FALSE,
    cluster_cols         = FALSE,
    treeheight_row       = 0,
    treeheight_col       = 0,
    gaps_col             = gaps_col_pos,
    show_rownames        = TRUE,
    show_colnames        = FALSE,
    labels_row           = gene_labels_hm,
    annotation_col       = ann_row_hm,   # sample (column) annotations
    annotation_row       = ann_site,     # site (row) annotations
    annotation_colors    = ann_colors_hm,
    annotation_names_col = TRUE,
    annotation_names_row = TRUE,
    color                = colors_v,
    breaks               = breaks_v,
    border_color         = NA,
    fontsize             = 8,
    fontsize_row         = 7,
    legend               = TRUE,
    annotation_legend    = TRUE,
    silent               = TRUE)

  # Stacked legend post-processing
  gt <- ph$gtable

  # ---- Dendrogram overlay ---------------------------------------------
  # cluster_cols = FALSE above (column order is fixed manually via
  # gaps_col), so pheatmap draws no tree. Each region block was clustered
  # independently (cluster_within_hm), so there's no single tree spanning
  # all cases -- draw the three block-level trees as their own mini-trees,
  # positioned using pheatmap's own gap-aware column-position formula
  # (reimplemented here from pheatmap:::find_coordinates, pheatmap 1.0.13,
  # rather than depending on an unexported internal) so every leaf lines up
  # exactly with its heatmap column, gaps included.
  n_total_hm    <- ncol(orient_plot)
  cell_size_hm  <- (1 / n_total_hm) * (unit(1, "npc") - length(gaps_col_pos) * unit(4, "bigpts"))
  gaps2_of_hm   <- function(idx) sum(gaps_col_pos < idx)
  col_center_hm <- function(idx) idx * cell_size_hm + gaps2_of_hm(idx) * unit(4, "bigpts") - cell_size_hm / 2

  # Segment table for one hclust: 3 segments per merge (left riser, right
  # riser, horizontal bar), x in leaf-position units (1..n, fractional for
  # internal nodes), y in raw merge height (0 at the leaves).
  hclust_segments <- function(hc) {
    n <- length(hc$order)
    # Leaf j (original row index, encoded as -j in hc$merge) is plotted at
    # its RANK within hc$order, not at j itself -- the heatmap columns are
    # arranged by ids[hc$order] (cluster_within_hm above), so leaf j's
    # x-position must be match(j, hc$order) or risers connect to the wrong
    # columns and branches cross ones that aren't their actual children.
    leaf_pos <- match(seq_len(n), hc$order)
    x <- setNames(as.numeric(leaf_pos), as.character(-seq_len(n)))
    y <- setNames(rep(0, n), as.character(-seq_len(n)))
    out <- vector("list", nrow(hc$merge))
    for (i in seq_len(nrow(hc$merge))) {
      a <- as.character(hc$merge[i, 1]); b <- as.character(hc$merge[i, 2])
      xa <- x[[a]]; xb <- x[[b]]; ya <- y[[a]]; yb <- y[[b]]
      yi <- hc$height[i]; xi <- (xa + xb) / 2
      out[[i]] <- data.frame(x0 = c(xa, xb, xa), y0 = c(ya, yb, yi),
                              x1 = c(xa, xb, xb), y1 = c(yi, yi, yi))
      x[[as.character(i)]] <- xi
      y[[as.character(i)]] <- yi
    }
    do.call(rbind, out)
  }

  blocks_for_tree <- list(
    list(hc = tloop_res$hc, offset = 0L),
    list(hc = stem3_res$hc, offset = length(tloop_ord)),
    list(hc = other_res$hc, offset = length(tloop_ord) + length(stem3_ord)))

  seg_x0 <- seg_x1 <- vector("list", 0L)
  seg_y0 <- seg_y1 <- numeric(0)
  for (blk in blocks_for_tree) {
    if (is.null(blk$hc) || max(blk$hc$height) <= 0) next
    segs  <- hclust_segments(blk$hc)
    max_h <- max(blk$hc$height)
    seg_x0 <- c(seg_x0, lapply(blk$offset + segs$x0, col_center_hm))
    seg_x1 <- c(seg_x1, lapply(blk$offset + segs$x1, col_center_hm))
    seg_y0 <- c(seg_y0, segs$y0 / max_h)
    seg_y1 <- c(seg_y1, segs$y1 / max_h)
  }

  if (length(seg_x0) > 0L) {
    tree_grob_hm <- grid::segmentsGrob(
      x0 = do.call(grid::unit.c, seg_x0), x1 = do.call(grid::unit.c, seg_x1),
      y0 = unit(seg_y0, "npc"), y1 = unit(seg_y1, "npc"),
      gp = grid::gpar(lwd = 0.6, col = "grey20"))
    ann_col_pos_hm <- gt$layout[gt$layout$name == "col_annotation", ]
    stopifnot(nrow(ann_col_pos_hm) == 1L)
    gt <- gtable::gtable_add_rows(gt, unit(24, "bigpts"), pos = ann_col_pos_hm$t - 1L)
    gt <- gtable::gtable_add_grob(gt, tree_grob_hm,
      t = ann_col_pos_hm$t, l = ann_col_pos_hm$l, r = ann_col_pos_hm$r,
      name = "col_tree_overlay")
  }

  # Make the "Signature" row-annotation label horizontal and left-aligned.
  # pheatmap nests text grobs inside gTrees; walk recursively to find and fix.
  fix_rot_label <- function(g, target) {
    if (inherits(g, "text") && isTRUE(any(as.character(g$label) == target))) {
      g$rot   <- 0
      g$hjust <- 0
      g$vjust <- 0.5
      g$x     <- unit(0, "npc")
    } else if (!is.null(g$children)) {
      g$children <- do.call(grid::gList, lapply(g$children, fix_rot_label, target))
    }
    g
  }
  gt$grobs <- lapply(gt$grobs, fix_rot_label, "Signature")

  ann_idx <- which(gt$layout$name == "annotation_legend")
  leg_idx <- which(gt$layout$name == "legend")
  stopifnot(length(ann_idx) == 1L, length(leg_idx) == 1L)
  legend_title_hm <- grid::textGrob(
    "dir x (logit-PSI\n- mu_ctrl) / max",
    x = 0, hjust = 0,
    gp = grid::gpar(fontsize = 9, fontface = "bold"))
  stacked_legend_hm <- gridExtra::arrangeGrob(
    gt$grobs[[ann_idx]], legend_title_hm, gt$grobs[[leg_idx]], grid::nullGrob(),
    ncol = 1L, heights = grid::unit(c(4, 0.8, 1.6, 3), "null"))
  gt$grobs[[ann_idx]] <- stacked_legend_hm
  gt$grobs[[leg_idx]] <- grid::nullGrob()

  # Wrap for patchwork embedding; negative bottom margin clips pheatmap's excess space
  fig_s3_panelD <- patchwork::wrap_elements(full = gt) +
    theme(plot.margin = margin(0, 0, -8, 0, "mm"))
  cat("Panel D: heatmap built successfully\n")
}, error = function(e) {
  cat("Panel D heatmap error:", conditionMessage(e), "\n")
  fig_s3_panelD <<- patchwork::plot_spacer()
})

# ---- 10. Combine into Figure 3 PDF -----------------------------------------
# Layout (4 panels):
#   Row 1: A (Stem III LASSO coef bar) | B (Stem III-only event annotation table)
#   Row 2: C (Stem III LASSO P(case) by variant group, incl. unmatched controls) -- full width
#   Row 3: D (insT + stemIII heatmap)                                            -- full width
mm_to_in <- function(mm) mm / 25.4

# fig_s3_panelB is theme_void() with no y-axis; wrap_elements() stops
# patchwork reserving left-margin space to align it with axis-bearing
# fig_s3_A1 in the same row (mirrors figure2.R Panel C's fix).
row1_s3 <- (fig_s3_A1 + labs(tag = "A")) |
           (patchwork::wrap_elements(full = fig_s3_panelB) + labs(tag = "B"))
row1_s3 <- row1_s3 + plot_layout(widths = c(0.25, 0.75))
row2_s3 <- fig_s3_A2 + labs(tag = "C")
panelD_w <- fig_s3_panelD + labs(tag = "D")

fig_s3_combined <- row1_s3 / row2_s3 / panelD_w +
  plot_layout(heights = c(0.35, 0.45, 1.1))

pdf(OUT_FIG_3, width = mm_to_in(170), height = mm_to_in(210))
print(fig_s3_combined & theme(plot.tag = element_text(size = 14, face = "bold"),
                              plot.tag.position = "topleft"))
dev.off()
cat(sprintf("Figure 3 saved to: %s\n", OUT_FIG_3))

# ---- 11. Summary RDS (airlock-safe) ----------------------------------------
# Counts and group-level aggregate stats; no event_ids, gene symbols, or coefs.
summary_out <- list(
  K_stemIII = K_stemIII_s3,
  K_s3_only = sum(s3only_mask_s3),
  K_insT    = K_insT_s3,
  # Panel A: group sample counts
  panel_a_n_ctrl          = length(c(gel_ctrl_ids_s3, nava_ctrl_ids_s3)),
  panel_a_n_insT_tloop    = nrow(inst_tloop_s3),
  panel_a_n_noninsT_tloop = nrow(noninst_tloop_s3),
  panel_a_n_stem3         = nrow(stem3_pc_s3),
  panel_a_n_del72         = nrow(del_pc_s3),
  panel_a_n_biallelic     = length(biallelic_in_pool_s3),
  # Panel D: case counts by region in heatmap
  panel_d_n_tloop = sum(case_region_s3 == "T-loop"),
  panel_d_n_stem3 = sum(case_region_s3 == "Stem III"),
  panel_d_n_other = sum(case_region_s3 == "unassigned"))
saveRDS(summary_out, OUT_SUMMARY)
cat(sprintf("Summary saved to: %s\n", OUT_SUMMARY))

cat("\n=== Done. Figure 3 ===\n")
cat(sprintf("Panel A: K=%d Stem III sites; %d ctrl, %d insT T-loop, %d non-insT T-loop, %d Stem III, %d biallelic\n",
            K_stemIII_s3,
            summary_out$panel_a_n_ctrl, summary_out$panel_a_n_insT_tloop,
            summary_out$panel_a_n_noninsT_tloop, summary_out$panel_a_n_stem3,
            summary_out$panel_a_n_biallelic))
cat(sprintf("Panel B: %d Stem III-only events in annotation table\n", length(s3only_gel_ev_s3)))
cat(sprintf("Panel D heatmap: %d T-loop, %d Stem III, %d other cases\n",
            summary_out$panel_d_n_tloop, summary_out$panel_d_n_stem3,
            summary_out$panel_d_n_other))
cat("\nNOTE: GEL-derived figure stays in workspace unless airlocked.\n")
