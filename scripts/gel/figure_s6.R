#!/usr/bin/env Rscript
# Layout: patchwork only -- never gridExtra arrangeGrob/grid.arrange
# =============================================================================
# figure_s6.R
# Supplementary Figure S4 (RNU4-2 splicing biomarker manuscript).
# Phase-2 robustness panels for the insT-specific LASSO uniform-sum signature.
#
#   Panel A — Per-site PSI distribution, one facet per selected site (1 row),
#             ordered by LASSO coef magnitude. Three groups per facet:
#             6k Whole-Blood ctrls (light grey, violin only), Dev controls
#             (mid grey, box+jitter), Dev insT cases (navy = PAL_WB,
#             box+jitter). Two per-site Mann-Whitney U tests (BH-adjusted
#             jointly across all K*2 tests): Dev ctrl vs 6k ctrl
#             (consistency) and Dev insT case vs 6k ctrl (signal). Stars
#             sit directly above the box being tested. The 6k group is
#             rendered as a violin to keep the PDF small with n>>1.
#   Panel B — PCA on raw PSI at K sites: unmatched controls, matched
#             controls, insT cases. Shape = group. Three subpanels coloured
#             by sex, age bin (10-yr), mapped reads bin.
#   Panel C — FPR as function of number of sites aggregated (1..K-1).
#             All C(K,k) site combinations per k. Score = LASSO linear
#             predictor restricted to that subset (pre-sigmoid). Threshold =
#             min case score for each combination. Each combination's FPR is
#             its own jittered point; boxplot summarises the distribution
#             across combinations at that k.
#
# Inputs:
#   outputs/loo_results_insT_lasso.rds  (a5ss_cand, selected_idx,
#                                        gel_ctrl_means, fit_at_1se,
#                                        sample_ids_train, cohort_train,
#                                        y_train, ho_scores)
#   outputs/gel_controls_psi_cand_pool.rds  (6k whole-blood with $psi/$ijc/$sjc
#                                            for all 483 concordant A5SS events)
#   outputs/psi_matrix.rds              (dev cohort PSI/IJC/SJC)
#   data/biallelic_rnu42_platekeys.txt  (optional)
#   data/renu_highdepth_case_path.txt   (optional)
#   data/rmats/inputs/renu_11_cases_paths.txt
#   data/rmats/inputs/renu_11_cases_control_paths.txt
#
# Outputs:
#   outputs/figures/figure_s6.pdf
#   outputs/figure_s6_summary.rds
#
# Run from workspace root: pixi run Rscript scripts/gel/figure_s6.R
#
# DATA PROTECTION: 6k cohort orient values are GEL-protected. Airlock the
# PDF + summary RDS for paper figure use.
# =============================================================================

# ---- libpath setup (pixi-aware) --------------------------------------------
# Under pixi, drop the GEL user_lib entirely — pkgs there are built against
# system R and produce "undefined symbol: R_mkClosure" style ABI errors when
# loaded under the pixi r-base. Outside pixi, prepend user_lib as usual.
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
  library(grid)
  library(gridExtra)
  library(patchwork)
  library(ggpubr)
  library(scales)
})

# ---- Configuration ----------------------------------------------------------
LOO_PATH      <- Sys.getenv("LOO_RDS",  "outputs/loo_results_insT_lasso.rds")
CTRL_MAT_PATH <- Sys.getenv("CTRL_MAT", "outputs/gel_controls_psi_cand_pool.rds")
DEV_PSI_PATH  <- Sys.getenv("DEV_PSI",  "outputs/psi_matrix.rds")
BIALLELIC_FILE<- Sys.getenv("BIALLELIC_FILE",
                            "data/biallelic_rnu42_platekeys.txt")
HD_CASE_FILE  <- Sys.getenv("HD_CASE_FILE",
                            "data/renu_highdepth_case_path.txt")
DEV_CASE_BAMS <- Sys.getenv("DEV_CASE_BAMS",
                            "data/rmats/inputs/renu_11_cases_paths.txt")
DEV_CTRL_BAMS <- Sys.getenv("DEV_CTRL_BAMS",
                            "data/rmats/inputs/renu_11_cases_control_paths.txt")

OUT_DIR_FIG <- "outputs/figures"
OUT_FIG_S4  <- file.path(OUT_DIR_FIG, "figure_s6.pdf")
OUT_SUMMARY <- "outputs/figure_s6_summary.rds"

# Paper figure palette (mirrors figure1.R; set 2026-04-29).
PAL_WB        <- "#2C5AA0"  # Whole-Blood navy
PAL_CTRL      <- "#888888"  # Ctrl grey

WB_LAB        <- "whole blood"
LY_LAB        <- "lymphocytes"
DEV_CTRL_LAB  <- "Dev control"
DEV_CASE_LAB  <- "Dev insT case"

stopifnot("insT LASSO LOO RDS not found"  = file.exists(LOO_PATH),
          "6k pool matrix not found"     = file.exists(CTRL_MAT_PATH),
          "Dev PSI matrix not found"     = file.exists(DEV_PSI_PATH))
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers ----------------------------------------------------------------
# Asterisk convention (paper-wide): ns p>0.05, * p<=0.05, ** p<=0.01,
# *** p<=0.001, **** p<=0.0001. Mirrors p_to_stars in figure1.R/figure2.R;
# named with the _s2 suffix here so the helper can move into a shared
# utils file later without colliding.
p_to_stars_s2 <- function(p) {
  if (is.na(p))         "ns"
  else if (p <= 0.0001) "****"
  else if (p <= 0.001)  "***"
  else if (p <= 0.01)   "**"
  else if (p <= 0.05)   "*"
  else                  "ns"
}

# Plotmath-parseable p-value label: plain decimal for p >= 0.001, otherwise
# "3.3 x 10^-16" style (mantissa %*% 10^exponent) -- use with
# geom_text(..., parse = TRUE). Mirrors figure2.R/figure3.R's format_p_sci().
format_p_sci_s2 <- function(p) {
  if (is.na(p)) return("\"ns\"")
  if (p >= 0.001) return(sprintf("p == %.3f", p))
  exp10  <- floor(log10(p))
  mant_r <- round(p / 10^exp10, 1)
  if (mant_r >= 10) { mant_r <- mant_r / 10; exp10 <- exp10 + 1 }
  sprintf("p == %.1f %%*%% 10^%d", mant_r, exp10)
}

bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}

bblogit_at_events <- function(ijc_mat, sjc_mat) {
  out <- bb_logit(ijc_mat, sjc_mat)
  out[(ijc_mat + sjc_mat) == 0] <- NA_real_
  out
}

score_v4 <- function(bblogit_mat, mu_ctrl, direction) {
  oriented <- sweep(bblogit_mat, 2L, mu_ctrl,   "-")
  oriented <- sweep(oriented,    2L, direction, "*")
  data.frame(
    sample_id = rownames(bblogit_mat),
    score     = rowSums(oriented, na.rm = TRUE),
    n_events  = rowSums(!is.na(bblogit_mat)),
    stringsAsFactors = FALSE
  )
}

read_bam_list <- function(path) {
  if (!file.exists(path)) return(character(0))
  raw  <- paste(readLines(path), collapse = ",")
  bams <- trimws(unlist(strsplit(raw, ",")))
  bams <- bams[nzchar(bams)]
  sub("\\.bam$", "", basename(bams))
}

# ---- Load inputs ------------------------------------------------------------
cat("Loading insT LASSO LOO RDS: ", LOO_PATH,      "\n", sep = "")
ins <- readRDS(LOO_PATH)
cat("Loading 6k pool matrix:     ", CTRL_MAT_PATH, "\n", sep = "")
ctrl <- readRDS(CTRL_MAT_PATH)
cat("Loading dev PSI:             ", DEV_PSI_PATH,  "\n", sep = "")
dev <- readRDS(DEV_PSI_PATH)

stopifnot(
  "insT RDS missing required fields" =
    all(c("a5ss_cand", "selected_idx", "gel_ctrl_means",
          "fit_at_1se", "sample_ids_train", "cohort_train", "y_train",
          "ho_scores") %in% names(ins)),
  "6k pool matrix missing required fields" =
    all(c("psi", "ijc", "sjc", "samples") %in% names(ctrl)),
  "dev matrix missing required fields" =
    all(c("psi", "ijc", "sjc", "samples") %in% names(dev))
)

sel_meta_v      <- ins$a5ss_cand[ins$selected_idx, , drop = FALSE]
selected_events <- sel_meta_v$gel_ev_id
direction_sel   <- sel_meta_v$direction
mu_ctrl_sel     <- ins$gel_ctrl_means[selected_events]
K_sel           <- ins$K_selected
beta_sel        <- ins$fit_at_1se$beta[ins$selected_idx]
intercept_1se   <- ins$fit_at_1se$intercept
cat(sprintf("insT-LASSO: K=%d selected events\n", K_sel))

events_missing <- setdiff(selected_events, colnames(ctrl$psi))
if (length(events_missing) > 0L) {
  stop("6k pool matrix missing ", length(events_missing), " insT-selected event(s). ",
       "Re-run combine_rmats_batches_cand_pool.R.")
}

# Identify insT GEL cases (positional IDs) from the insT RDS
gel_insT_train_ids <- ins$sample_ids_train[ins$cohort_train == "GEL" &
                                            ins$y_train == 1L]
gel_insT_ho_ids    <- ins$ho_scores$gel$sample_id
gel_insT_all_ids   <- c(gel_insT_train_ids, gel_insT_ho_ids)
cat(sprintf("insT GEL cases: %d (train=%d, holdout=%d)\n",
            length(gel_insT_all_ids),
            length(gel_insT_train_ids), length(gel_insT_ho_ids)))

# Selected-site labels + ordering (by LASSO coef magnitude)
sel_coef_v   <- unname(beta_sel)
sel_gene     <- sel_meta_v$GeneSymbol
sel_gene_lab <- if (anyDuplicated(sel_gene))
  make.unique(sel_gene, sep = ".") else sel_gene
sel_order    <- order(-sel_coef_v, -sel_meta_v$g_auc)
sel_gene_ord <- sel_gene_lab[sel_order]

# ---- Score 6k controls + dedup dev cohort ----------------------------------
ctrl_ijc_sel <- ctrl$ijc[, selected_events, drop = FALSE]
ctrl_sjc_sel <- ctrl$sjc[, selected_events, drop = FALSE]
ctrl_bblogit <- bblogit_at_events(ctrl_ijc_sel, ctrl_sjc_sel)
ctrl_scores  <- score_v4(ctrl_bblogit, mu_ctrl_sel, direction_sel)
# Per-site orient matrix (N x K) — used for per-site FP analysis (panel C)
ctrl_orient_mat <- sweep(sweep(ctrl_bblogit, 2L, mu_ctrl_sel, "-"),
                         2L, direction_sel, "*")

dev_case_keys <- read_bam_list(DEV_CASE_BAMS)
dev_ctrl_keys <- read_bam_list(DEV_CTRL_BAMS)
dev_all_keys  <- c(dev_case_keys, dev_ctrl_keys)
overlap_keys  <- intersect(ctrl_scores$sample_id, dev_all_keys)
cat(sprintf("Dev/6k overlap: %d sample(s) — deduplicated from 6k cohort\n",
            length(overlap_keys)))

# `_dedup` removes only the dev cases/controls from the 6k pool (used for
# biallelic + HD lookup so we don't accidentally re-find them as dev
# samples). The grey-histogram cohort `_hist` (computed lower down once
# biallelic_in_6k + hd_in_6k are known) further excludes those overlay
# groups so the histogram count is "controls that aren't in any
# overlay strip".
ctrl_scores_dedup <- ctrl_scores[!ctrl_scores$sample_id %in% overlap_keys, ,
                                 drop = FALSE]
cat(sprintf("After dev-overlap dedup: %d 6k samples scored\n",
            nrow(ctrl_scores_dedup)))

# Dev P(case) — recompute from dev$ijc/$sjc using full-train coefs.
stopifnot(
  "dev case BAM list count" = length(dev_case_keys) == 11L,
  "dev ctrl BAM list count" = length(dev_ctrl_keys) == 33L,
  "dev$samples not in case-then-ctrl order" =
    identical(tolower(dev$samples$group),
              c(rep("case", 11), rep("control", 33))))
dev_real_ids        <- c(dev_case_keys, dev_ctrl_keys)
positional_to_plate <- setNames(dev_real_ids, rownames(dev$psi))
plate_to_positional <- setNames(rownames(dev$psi), dev_real_ids)

dev_ijc_sel <- dev$ijc[, selected_events, drop = FALSE]
dev_sjc_sel <- dev$sjc[, selected_events, drop = FALSE]
dev_bblogit <- bblogit_at_events(dev_ijc_sel, dev_sjc_sel)

# ---- Biallelic RNU4-2 carrier overlay --------------------------------------
biallelic_keys <- if (file.exists(BIALLELIC_FILE)) {
  trimws(readLines(BIALLELIC_FILE))
} else {
  cat(sprintf("Biallelic file not found at %s — skipping biallelic overlay\n",
              BIALLELIC_FILE))
  character(0)
}
biallelic_keys <- biallelic_keys[nzchar(biallelic_keys)]
biallelic_in_6k <- intersect(biallelic_keys, ctrl_scores_dedup$sample_id)
cat(sprintf("Biallelic RNU4-2 carriers in 6k cohort: %d/%d\n",
            length(biallelic_in_6k), length(biallelic_keys)))

# ---- High-depth ReNU case (a known case in the 6k pool) -------------------
hd_keys <- if (file.exists(HD_CASE_FILE)) {
  read_bam_list(HD_CASE_FILE)
} else {
  cat(sprintf("HD case file not found at %s — skipping HD overlay\n",
              HD_CASE_FILE))
  character(0)
}
hd_in_6k <- intersect(hd_keys, ctrl_scores_dedup$sample_id)
cat(sprintf("HD ReNU case(s) in 6k (post-dedup): %d\n", length(hd_in_6k)))

# Held-out insT LASSO P(case) for the HD ReNU sample(s). This case sits in
# the 6k WB pool but was excluded from the training set, so this is a true
# held-out prediction (n=1, no AUC). Uses the same orient matrix / beta /
# intercept as panel C's LASSO linear predictor.
if (length(hd_in_6k) > 0L) {
  hd_orient <- ctrl_orient_mat[hd_in_6k, , drop = FALSE]
  hd_lin    <- as.numeric(intercept_1se +
                          rowSums(sweep(hd_orient, 2L, beta_sel, "*"),
                                  na.rm = TRUE))
  hd_prob   <- 1 / (1 + exp(-hd_lin))
  for (i in seq_along(hd_in_6k)) {
    cat(sprintf(
      "HD ReNU (held-out) insT LASSO P(case): %.4f  (linpred=%+.4f)  sample=%s\n",
      hd_prob[i], hd_lin[i], hd_in_6k[i]))
  }
}

# ---- Clean control cohort for FP analysis: 6k MINUS dev overlap + biallelic + HD --
# Excludes any sample that could plausibly be a true positive so the FP rate
# is estimated on genuinely unaffected controls.
hist_excl_keys   <- unique(c(overlap_keys, biallelic_in_6k, hd_in_6k))
ctrl_scores_hist <- ctrl_scores[!ctrl_scores$sample_id %in% hist_excl_keys, ,
                                drop = FALSE]
cat(sprintf("Clean control cohort (dev + biallelic + HD excluded): %d 6k samples\n",
            nrow(ctrl_scores_hist)))

# ---- Plot theme ------------------------------------------------------------
pub_theme <- theme_classic(base_size = 9) +
  # All text standardised to 9 pt across figure_s6.

  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 9, face = "bold"),
        axis.text        = element_text(size = 9),
        axis.title       = element_text(size = 9),
        axis.line        = element_line(lineend = "square"),
        legend.title     = element_text(size = 9),
        legend.text      = element_text(size = 9),
        plot.margin      = margin(2, 2, 2, 2))

# ============================================================================
# Panel A: PSI distribution per selected site, split into 3 groups:
#   1) 6k whole-blood controls (post-dedup)
#   2) Dev controls (33 from psi_matrix.rds)
#   3) Dev ReNU cases (11 from psi_matrix.rds)
# One facet per site, ordered by LASSO coef magnitude (matches Fig 1F).
# Sanity check: dev cases should be visibly shifted away from both control
# groups at the high-coef sites; the two control groups should overlap.
# (Was Panel B before the 2026-05-03 restructure — old Panel A "% IJC+SJC>0"
# lollipop was dropped because every site detected at >=99% of the 6k cohort,
# rendering the panel uninformative.)
# ============================================================================
GROUP_6K       <- "unmatched controls"
GROUP_DEV_CTRL <- "matched controls"
GROUP_DEV_CASE <- DEV_CASE_LAB
b_levels <- c(GROUP_6K, GROUP_DEV_CTRL, GROUP_DEV_CASE)

ctrl_6k_psi  <- ctrl$psi[ctrl_scores_hist$sample_id, selected_events, drop = FALSE]
dev_case_psi <- dev$psi[gel_insT_all_ids, selected_events, drop = FALSE]
dev_ctrl_psi <- dev$psi[plate_to_positional[dev_ctrl_keys], selected_events,
                        drop = FALSE]

build_psi_long <- function(psi_mat, group_lab) {
  data.frame(
    gene_lab = rep(sel_gene_lab[match(colnames(psi_mat), selected_events)],
                   each = nrow(psi_mat)),
    psi      = as.numeric(psi_mat),
    group    = group_lab,
    stringsAsFactors = FALSE)
}

psi_long <- rbind(
  build_psi_long(ctrl_6k_psi,  GROUP_6K),
  build_psi_long(dev_ctrl_psi, GROUP_DEV_CTRL),
  build_psi_long(dev_case_psi, GROUP_DEV_CASE))
psi_long <- psi_long[!is.na(psi_long$psi), , drop = FALSE]
psi_long$gene_lab <- factor(psi_long$gene_lab, levels = sel_gene_ord)
psi_long$group    <- factor(psi_long$group,    levels = b_levels)

cat(sprintf("Panel A: %s n=%s, %s n=%d, %s n=%d\n",
            GROUP_6K,       format(nrow(ctrl_6k_psi),  big.mark = ","),
            GROUP_DEV_CTRL, nrow(dev_ctrl_psi),
            GROUP_DEV_CASE, nrow(dev_case_psi)))

# ---- Per-event quantification % (console diagnostic) -----------------------
# For each selected signature event, report the % of samples with non-NA PSI
# in each of the three groups. Complements the panel A distribution plot.
quant_pct <- data.frame(
  gene          = sel_gene_lab,
  event_id      = selected_events,
  pct_6k        = 100 * colMeans(!is.na(ctrl_6k_psi)),
  pct_dev_ctrl  = 100 * colMeans(!is.na(dev_ctrl_psi)),
  pct_dev_case  = 100 * colMeans(!is.na(dev_case_psi)),
  n_na_6k       = colSums(is.na(ctrl_6k_psi)),
  n_na_dev_ctrl = colSums(is.na(dev_ctrl_psi)),
  n_na_dev_case = colSums(is.na(dev_case_psi)),
  stringsAsFactors = FALSE)
cat(sprintf(
  "\nPer-event quantification %% (samples with non-NA PSI) -- 6k n=%d, dev ctrl n=%d, dev case n=%d:\n",
  nrow(ctrl_6k_psi), nrow(dev_ctrl_psi), nrow(dev_case_psi)))
print(quant_pct[order(quant_pct$pct_6k), , drop = FALSE],
      row.names = FALSE, digits = 4)

# ---- TPM at NA-PSI (dev cohort only, per signature event) -------------------
# For each selected event, look up the gene's TPM in the dev samples where
# the event was NOT quantified (PSI = NA) vs those where it was. Answers
# "does 'unquantifiable event' just mean the gene isn't expressed in that
# sample?". Uses the same per-sample Salmon quant paths as figure2.R /
# figure_s5.R; 6k cohort has no per-sample TPM available, so skipped there.
ANNOT_PATH_S4      <- Sys.getenv("ANNOT_CSV",   "outputs/a5ss_candidate_full_table.csv")
QUANT_PATHS_CSV_S4 <- Sys.getenv("QUANT_PATHS", "data/salmon_quant_paths.csv")

if (!file.exists(ANNOT_PATH_S4)) {
  cat(sprintf("\nTPM-at-NA-PSI check skipped: annotation table not found at %s\n",
              ANNOT_PATH_S4))
} else if (!file.exists(QUANT_PATHS_CSV_S4)) {
  cat(sprintf("\nTPM-at-NA-PSI check skipped: quant paths CSV not found at %s\n",
              QUANT_PATHS_CSV_S4))
} else {
  annot_s4 <- read.csv(ANNOT_PATH_S4, stringsAsFactors = FALSE)
  qp_s4    <- read.csv(QUANT_PATHS_CSV_S4, stringsAsFactors = FALSE)
  stopifnot("quant paths CSV must have platekey and file_path columns" =
              all(c("platekey", "file_path") %in% names(qp_s4)))
  ensg_lookup_s4 <- function(sym) {
    hit <- annot_s4$ensg[!is.na(annot_s4$gene_symbol) &
                          toupper(annot_s4$gene_symbol) == toupper(sym)]
    if (length(hit) == 0L) NA_character_ else hit[1L]
  }
  gene_ensg_sel <- vapply(sel_gene, ensg_lookup_s4, character(1L))
  needed_ensg   <- unique(gene_ensg_sel[!is.na(gene_ensg_sel)])

  dev_all_pos <- rownames(dev$psi)
  tpm_cache <- matrix(NA_real_,
                      nrow = length(dev_all_pos),
                      ncol = length(needed_ensg),
                      dimnames = list(dev_all_pos, needed_ensg))
  for (pid in dev_all_pos) {
    pk  <- positional_to_plate[pid]
    row <- qp_s4[qp_s4$platekey == pk, , drop = FALSE]
    if (nrow(row) == 0L || !file.exists(row$file_path[1L])) next
    sf  <- read.delim(row$file_path[1L], stringsAsFactors = FALSE)
    sf_ensg <- sub("\\..*$", "", sf$Name)
    for (e in needed_ensg) {
      hit <- sf$TPM[sf_ensg == e]
      if (length(hit)) tpm_cache[pid, e] <- sum(hit, na.rm = TRUE)
    }
  }

  gel_case_ids_s4 <- paste0("gel:case:",    seq_len(11L))
  cat("\nTPM in NA-PSI vs quantified dev samples, per signature event:\n")
  for (j in seq_along(selected_events)) {
    e   <- selected_events[j]
    sym <- sel_gene[j]
    ensg <- gene_ensg_sel[j]
    if (is.na(ensg)) {
      cat(sprintf("  %s (%s): ENSG lookup failed -- skipped\n", sym, e))
      next
    }
    dev_psi_e <- dev$psi[dev_all_pos, e]
    na_ids <- dev_all_pos[is.na(dev_psi_e)]
    ok_ids <- dev_all_pos[!is.na(dev_psi_e)]
    if (length(na_ids) == 0L) {
      cat(sprintf("  %s (%s, %s): quantified in all %d dev samples\n",
                  sym, e, ensg, length(dev_all_pos)))
      next
    }
    tpm_na <- tpm_cache[na_ids, ensg]
    tpm_ok <- tpm_cache[ok_ids, ensg]
    cat(sprintf("  %s (%s, %s): NA in %d/%d dev samples\n",
                sym, e, ensg, length(na_ids), length(dev_all_pos)))
    cat(sprintf("    TPM in NA-PSI samples:     n=%d median=%.2f min=%.2f max=%.2f\n",
                sum(!is.na(tpm_na)),
                median(tpm_na, na.rm = TRUE),
                min(tpm_na,    na.rm = TRUE),
                max(tpm_na,    na.rm = TRUE)))
    cat(sprintf("    TPM in quantified samples: n=%d median=%.2f min=%.2f max=%.2f\n",
                sum(!is.na(tpm_ok)),
                median(tpm_ok, na.rm = TRUE),
                min(tpm_ok,    na.rm = TRUE),
                max(tpm_ok,    na.rm = TRUE)))
    grp_lab <- ifelse(na_ids %in% gel_case_ids_s4, "case", "ctrl")
    cat(sprintf("    per-sample TPM in NA-PSI:  %s\n",
                paste(sprintf("%s(%s)=%.2f", na_ids, grp_lab, tpm_na),
                      collapse = ", ")))
  }

  # ---- 6k cohort TOP1 TPM: NA-PSI vs quantified -----------------------------
  # Derives per-6k-sample TPM from the control_manifest.csv BAM paths (same
  # .bam -> .quant.genes.sf substitution used by build_salmon_quant_paths.sh).
  # Cached to RDS since it's ~5984 file reads; invalidated on MANIFEST mtime.
  MANIFEST_S4  <- Sys.getenv("MANIFEST_PATH", "data/control_manifest.csv")
  TOP1_CACHE   <- Sys.getenv("TOP1_CACHE",    "outputs/top1_tpm_6k_cache.rds")
  top1_j       <- which(sel_gene == "TOP1")
  if (length(top1_j) == 0L) {
    cat("\n6k TOP1 TPM check skipped: TOP1 not among selected signature events\n")
  } else if (!file.exists(MANIFEST_S4)) {
    cat(sprintf("\n6k TOP1 TPM check skipped: manifest not found at %s\n", MANIFEST_S4))
  } else {
    top1_ev   <- selected_events[top1_j]
    top1_ensg <- gene_ensg_sel[top1_j]
    if (is.na(top1_ensg)) {
      cat("\n6k TOP1 TPM check skipped: TOP1 ENSG lookup failed\n")
    } else {
      cohort_ids_6k <- ctrl_scores_hist$sample_id
      top1_psi_6k   <- ctrl$psi[cohort_ids_6k, top1_ev]
      na_6k <- cohort_ids_6k[is.na(top1_psi_6k)]
      ok_6k <- cohort_ids_6k[!is.na(top1_psi_6k)]

      # Load or (re)build cache: named numeric vector platekey -> TOP1 TPM
      need_rebuild <- !file.exists(TOP1_CACHE) ||
        file.mtime(MANIFEST_S4) > file.mtime(TOP1_CACHE)
      if (need_rebuild) {
        cat(sprintf("\nBuilding TOP1 TPM cache for 6k cohort (n=%d samples)...\n",
                    length(cohort_ids_6k)))
        manifest_6k <- read.csv(MANIFEST_S4, stringsAsFactors = FALSE)
        stopifnot("manifest missing file_path" = "file_path" %in% names(manifest_6k))
        manifest_6k$platekey <- sub("\\.bam$", "", basename(manifest_6k$file_path))
        manifest_6k$sf_path  <- sub("\\.bam$", ".quant.genes.sf", manifest_6k$file_path)
        pk_to_sf <- setNames(manifest_6k$sf_path, manifest_6k$platekey)

        tpm_6k <- setNames(rep(NA_real_, length(cohort_ids_6k)), cohort_ids_6k)
        n_missing_file <- 0L
        for (i in seq_along(cohort_ids_6k)) {
          pid <- cohort_ids_6k[i]
          sf  <- pk_to_sf[pid]
          if (is.na(sf) || !file.exists(sf)) { n_missing_file <- n_missing_file + 1L; next }
          sf_df <- read.delim(sf, stringsAsFactors = FALSE)
          hit   <- sf_df$TPM[sub("\\..*$", "", sf_df$Name) == top1_ensg]
          if (length(hit)) tpm_6k[pid] <- sum(hit, na.rm = TRUE)
          if (i %% 500L == 0L) cat(sprintf("  ... %d/%d\n", i, length(cohort_ids_6k)))
        }
        cat(sprintf("  built TOP1 TPM for %d/%d samples (%d missing quant.sf files)\n",
                    sum(!is.na(tpm_6k)), length(cohort_ids_6k), n_missing_file))
        saveRDS(list(tpm = tpm_6k, ensg = top1_ensg,
                     built = Sys.time(), n = length(cohort_ids_6k)),
                TOP1_CACHE)
      }
      cache_obj <- readRDS(TOP1_CACHE)
      tpm_6k    <- cache_obj$tpm
      missing_pids <- setdiff(cohort_ids_6k, names(tpm_6k))
      if (length(missing_pids) > 0L) {
        cat(sprintf("Note: %d cohort samples not in cache -- delete %s to rebuild\n",
                    length(missing_pids), TOP1_CACHE))
      }
      tpm_na <- tpm_6k[intersect(na_6k, names(tpm_6k))]
      tpm_ok <- tpm_6k[intersect(ok_6k, names(tpm_6k))]

      cat(sprintf("\n6k cohort TOP1 event (%s, ENSG %s):\n", top1_ev, top1_ensg))
      cat(sprintf("  quantified (PSI non-NA): n=%d samples\n", length(ok_6k)))
      cat(sprintf("  not quantified (PSI NA): n=%d samples\n", length(na_6k)))
      cat(sprintf("  TPM range in quantified samples:     n=%d min=%.2f median=%.2f mean=%.2f max=%.2f\n",
                  sum(!is.na(tpm_ok)),
                  min(tpm_ok, na.rm = TRUE),  median(tpm_ok, na.rm = TRUE),
                  mean(tpm_ok, na.rm = TRUE), max(tpm_ok, na.rm = TRUE)))
      cat(sprintf("  TPM range in not-quantified samples: n=%d min=%.2f median=%.2f mean=%.2f max=%.2f\n",
                  sum(!is.na(tpm_na)),
                  min(tpm_na, na.rm = TRUE),  median(tpm_na, na.rm = TRUE),
                  mean(tpm_na, na.rm = TRUE), max(tpm_na, na.rm = TRUE)))
    }
  }
}

# Per-site Mann-Whitney U tests, two comparisons:
#   (1) Dev control vs 6k control     -> consistency check (should be ns)
#   (2) Dev ReNU case vs 6k control   -> signal check (should be sig at most)
# BH adjustment is applied jointly across all 30 tests (15 sites x 2
# comparisons). Stars are annotated above each facet (one strip near the
# top per comparison), colour-coded by comparison so the consistency vs
# signal layers stay visually distinct.
do_mwu_psi <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) >= 2L && length(b) >= 2L)
    suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
  else NA_real_
}
panel_c_tests <- do.call(rbind, lapply(seq_along(selected_events), function(j) {
  e <- selected_events[j]
  data.frame(
    gene_lab     = sel_gene_lab[j],
    p_ctrl_vs_6k = do_mwu_psi(dev_ctrl_psi[, e], ctrl_6k_psi[, e]),
    p_case_vs_6k = do_mwu_psi(dev_case_psi[, e], ctrl_6k_psi[, e]),
    stringsAsFactors = FALSE)
}))
all_p <- c(panel_c_tests$p_ctrl_vs_6k, panel_c_tests$p_case_vs_6k)
all_q <- p.adjust(all_p, method = "BH")
n_sites <- nrow(panel_c_tests)
panel_c_tests$q_ctrl_vs_6k <- all_q[seq_len(n_sites)]
panel_c_tests$q_case_vs_6k <- all_q[n_sites + seq_len(n_sites)]
panel_c_tests$gene_lab     <- factor(panel_c_tests$gene_lab,
                                     levels = sel_gene_ord)

cat(sprintf("\nPanel A per-site MWU (BH across all %d tests):\n",
            nrow(panel_c_tests) * 2L))
panel_c_tests_disp <- panel_c_tests[order(panel_c_tests$gene_lab), ]
panel_c_tests_disp$stars_ctrl <- vapply(panel_c_tests_disp$q_ctrl_vs_6k,
                                         p_to_stars_s2, character(1L))
panel_c_tests_disp$stars_case <- vapply(panel_c_tests_disp$q_case_vs_6k,
                                         p_to_stars_s2, character(1L))
print(data.frame(
  gene = panel_c_tests_disp$gene_lab,
  q_dev_ctrl_vs_6k_ctrl = signif(panel_c_tests_disp$q_ctrl_vs_6k, 3),
  stars_ctrl = panel_c_tests_disp$stars_ctrl,
  q_dev_case_vs_6k_ctrl = signif(panel_c_tests_disp$q_case_vs_6k, 3),
  stars_case = panel_c_tests_disp$stars_case,
  stringsAsFactors = FALSE), row.names = FALSE)

# Bracket-and-p-value annotation per facet, lifted above y=1 so PSIs that
# saturate at 1 don't hide the labels. Two comparisons stacked
# vertically:
#   Lower bracket (y=1.10): 6k ctrl (x=1) <-> Dev ctrl (x=2), grey
#   Upper bracket (y=1.22): 6k ctrl (x=1) <-> Dev case (x=3), navy
# Each bracket spans the two boxes being tested; the BH-adjusted p-value
# (not asterisks -- real numbers requested so magnitude is legible) sits
# above the bracket midpoint. Colour matches the dev group (so the bracket
# colour reads "this is the dev-ctrl test" / "this is the dev-case
# test" without needing a legend).
panel_c_ann <- data.frame(
  gene_lab   = rep(panel_c_tests$gene_lab, 2L),
  comparison = rep(c("vs Dev ctrl", "vs Dev case"),
                   each = nrow(panel_c_tests)),
  label      = c(vapply(panel_c_tests$q_ctrl_vs_6k, format_p_sci_s2,
                        character(1L)),
                 vapply(panel_c_tests$q_case_vs_6k, format_p_sci_s2,
                        character(1L))),
  x_start    = 1,                           # always vs 6k ctrl
  x_end      = rep(c(2, 3), each = nrow(panel_c_tests)),
  y_bar      = rep(c(1.10, 1.22), each = nrow(panel_c_tests)),
  colour     = rep(c(PAL_CTRL, PAL_WB), each = nrow(panel_c_tests)),
  stringsAsFactors = FALSE)
panel_c_ann$x_text <- (panel_c_ann$x_start + panel_c_ann$x_end) / 2

# 6k group has n>>1; rendering 6018 jittered points per facet bloats the
# PDF and makes the panel look noisy. Use a violin layer for the 6k group
# only (vector density curve, small file size) and keep boxplot+jitter
# for the small dev groups where individual points are informative.
fig_s2a <- ggplot(psi_long, aes(x = group, y = psi, colour = group)) +
  geom_violin(data = subset(psi_long, group == GROUP_6K),
              fill = NA, linewidth = 0.25, scale = "width",
              trim = FALSE, na.rm = TRUE) +
  geom_boxplot(fill = NA, outlier.shape = NA,
               width = 0.35, linewidth = 0.25) +
  geom_jitter(data = subset(psi_long, group != GROUP_6K),
              width = 0.15, size = 0.4, alpha = 0.6) +
  geom_segment(data = panel_c_ann,
               aes(x = x_start, xend = x_end,
                   y = y_bar, yend = y_bar),
               inherit.aes = FALSE,
               colour = panel_c_ann$colour, linewidth = 0.3) +
  geom_text(data = panel_c_ann,
            aes(x = x_text, y = y_bar, label = label),
            inherit.aes = FALSE,
            colour = panel_c_ann$colour,
            parse  = TRUE,
            vjust  = -0.2,
            size   = 9 / .pt) +   # 9 pt via geom_text's mm scale
  facet_wrap(~ gene_lab, nrow = 1L, strip.position = "bottom") +
  scale_y_continuous(limits = c(0, 1.32),
                     breaks = c(0, 0.5, 1)) +
  scale_colour_manual(
    values = setNames(c("#cccccc", PAL_CTRL, PAL_WB), b_levels),
    name   = NULL) +
  labs(x = NULL, y = "PSI") +
  pub_theme +
  theme(axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        strip.text          = element_text(size = 9, face = "italic"),
        strip.background    = element_blank(),
        strip.placement     = "outside",
        panel.spacing = unit(0.8, "mm"),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"))

# ---- Orient matrices for panel C (also used by Panel A MWU) ----------------
# Per-site orient matrix for dev insT cases
dev_case_orient <- sweep(
  sweep(dev_bblogit[gel_insT_all_ids, , drop = FALSE], 2L, mu_ctrl_sel, "-"),
  2L, direction_sel, "*")

# Per-site orient matrix for clean control cohort
ctrl_orient_clean <- ctrl_orient_mat[ctrl_scores_hist$sample_id, , drop = FALSE]
n_ctrl_clean      <- nrow(ctrl_orient_clean)

# ============================================================================
# Panel C: FPR as a function of number of sites aggregated (1..K_sel).
# For each k, enumerate all C(K_sel, k) site combinations.
# Score = sum_j(beta_j * orient_ij) — LASSO linear predictor, pre-sigmoid.
# Threshold = min case score for that combination.
# FPR = fraction of clean controls exceeding that threshold.
# At k = K_sel this equals the full-model linear predictor (cf. Fig 1G).
# Dot = mean FPR across combinations; bar = min-to-max range.
# Dot shape: filled circle = at least one combination with 0 FPs;
#            open circle   = all combinations have FPs.
# Min-count label at bottom of each bar shows exact lowest-FP count.
# ============================================================================
combo_fpr_raw_list <- vector("list", K_sel)
combo_fpr_list <- lapply(seq_len(K_sel), function(k) {
  combos   <- combn(K_sel, k, simplify = FALSE)
  n_fp_v   <- vapply(combos, function(idx) {
    b      <- beta_sel[idx]
    case_s <- rowSums(sweep(dev_case_orient[, idx, drop = FALSE],
                            2L, b, "*"), na.rm = TRUE)
    ctrl_s <- rowSums(sweep(ctrl_orient_clean[, idx, drop = FALSE],
                            2L, b, "*"), na.rm = TRUE)
    thresh <- min(case_s, na.rm = TRUE)
    sum(ctrl_s > thresh, na.rm = TRUE)
  }, numeric(1L))
  fpr_v <- n_fp_v / n_ctrl_clean * 100
  combo_fpr_raw_list[[k]] <<- data.frame(k = k, fpr = fpr_v,
                                         stringsAsFactors = FALSE)
  data.frame(
    k        = k,
    mean_fpr = mean(fpr_v),
    min_fpr  = min(fpr_v),
    max_fpr  = max(fpr_v),
    min_n_fp = min(n_fp_v),
    max_n_fp = max(n_fp_v),
    n_combos = length(fpr_v),
    stringsAsFactors = FALSE
  )
})
combo_fpr_df <- do.call(rbind, combo_fpr_list)
combo_fpr_df$any_zero <- combo_fpr_df$min_n_fp == 0
combo_fpr_raw_df <- do.call(rbind, combo_fpr_raw_list)

cat(sprintf("\nPanel C - FPR by number of sites (n_unmatched_ctrl=%d):\n", n_ctrl_clean))
print(combo_fpr_df, row.names = FALSE)

cat("\nZero-FP combinations (k=3 to K-1):\n")
for (k in seq(3L, K_sel - 1L)) {
  combos  <- combn(K_sel, k, simplify = FALSE)
  zero_combos <- Filter(function(idx) {
    b      <- beta_sel[idx]
    case_s <- rowSums(sweep(dev_case_orient[, idx, drop = FALSE], 2L, b, "*"), na.rm = TRUE)
    ctrl_s <- rowSums(sweep(ctrl_orient_clean[, idx, drop = FALSE], 2L, b, "*"), na.rm = TRUE)
    sum(ctrl_s > min(case_s, na.rm = TRUE), na.rm = TRUE) == 0L
  }, combos)
  cat(sprintf("  k=%d: %d zero-FP combination(s) out of %d\n",
              k, length(zero_combos), length(combos)))
  for (idx in zero_combos) {
    cat(sprintf("    %s\n", paste(sel_gene_lab[idx], collapse = " + ")))
  }
}

# Every C(K_sel, k) combination's FPR shown as its own jittered point
# (boxplot summarises the distribution across combinations at that k) --
# replaces the earlier mean-dot + min/max-range design, which collapsed
# each k down to 3 numbers and hid how the combinations were distributed.
fig_s2c <- ggplot(combo_fpr_raw_df[combo_fpr_raw_df$k < K_sel, ],
                  aes(x = factor(k), y = fpr)) +
  geom_boxplot(colour = PAL_WB, fill = NA, outlier.shape = NA,
               width = 0.5, linewidth = 0.4) +
  geom_jitter(colour = PAL_WB, width = 0.15, size = 0.7, alpha = 0.5) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0.02, 0.1))) +
  labs(x = "Number of sites",
       y = sprintf("FPR in unmatched controls (%%, n=%s)",
                   format(n_ctrl_clean, big.mark = ","))) +
  pub_theme +
  # Extra bottom padding so Panel C's plot region sits at roughly the same
  # y as Panel A's plot region -- A has a legend below its x-axis, C does
  # not; without this pad C ends up sitting visibly lower on the page.
  theme(plot.margin = margin(2, 2, 25, 2))

# ============================================================================
# Combine into figure_s6 PDF  (patchwork throughout for axis alignment)
# ============================================================================
mm_to_in <- function(mm) mm / 25.4

# ============================================================================
# Panel B: PCA on raw PSI at K selected sites — confounding assessment.
# Samples: unmatched controls (6k clean), matched controls (33 dev), insT
# cases. Shape = group. Three subpanels coloured by sex, age bin, mapped
# reads bin. Metadata from control_manifest.csv.
#
# NOTE: dev cohort and 6k controls come from different rMATS runs; raw PSI
# is run-dependent, so PC axes may partly reflect run/batch. Interpret
# within-cohort covariate gradients more than between-cohort separation.
#
# Expected manifest columns (adjust via env vars if names differ):
#   MANIFEST_KEY_COL  — sample identifier (default: rnaseq_platekey)
#   MANIFEST_SEX_COL  — sex (default: participant_phenotypic_sex)
#   MANIFEST_YOB_COL  — year of birth or age (default: year_of_birth)
#   MANIFEST_MAP_COL  — mapped unique reads (default: mapped_unique_reads)
# ============================================================================
MANIFEST_PATH      <- Sys.getenv("MANIFEST_PATH",       "data/control_manifest.csv")
QCFAIL_MANIFEST    <- Sys.getenv("QCFAIL_MANIFEST_PATH","data/qcfail_manifest.csv")
MANIFEST_KEY       <- Sys.getenv("MANIFEST_KEY_COL",    "rnaseq_platekey")
MANIFEST_SEX       <- Sys.getenv("MANIFEST_SEX_COL",    "participant_phenotypic_sex")
MANIFEST_YOB       <- Sys.getenv("MANIFEST_YOB_COL",    "age_at_rnaseq")
MANIFEST_MAP       <- Sys.getenv("MANIFEST_MAP_COL",    "mapped_unique_reads")

stopifnot("Control manifest not found" = file.exists(MANIFEST_PATH))
manifest <- read.csv(MANIFEST_PATH, stringsAsFactors = FALSE)
if (file.exists(QCFAIL_MANIFEST)) {
  qcfail  <- read.csv(QCFAIL_MANIFEST, stringsAsFactors = FALSE)
  shared  <- intersect(colnames(manifest), colnames(qcfail))
  combined <- rbind(manifest[, shared, drop = FALSE],
                    qcfail[,    shared, drop = FALSE])
  manifest <- combined[!duplicated(combined[[MANIFEST_KEY]]), , drop = FALSE]
  cat(sprintf("Merged control + qcfail manifests: %d unique samples\n", nrow(manifest)))
} else {
  cat(sprintf("qcfail manifest not found at %s — using control manifest only\n",
              QCFAIL_MANIFEST))
}
stopifnot("Manifest key column missing" = MANIFEST_KEY %in% colnames(manifest))

# ---- Build bb-logit PSI matrices per group ---------------------------------
# Use bblogit (with pseudocounts) so NA-free; consistent with scoring metric.
psi_unctrl  <- ctrl_bblogit[ctrl_scores_hist$sample_id, selected_events, drop = FALSE]
unctrl_ids  <- ctrl_scores_hist$sample_id

matched_pos <- plate_to_positional[dev_ctrl_keys]
psi_matched <- dev_bblogit[matched_pos, selected_events, drop = FALSE]
matched_ids <- dev_ctrl_keys

psi_cases   <- dev_bblogit[gel_insT_all_ids, selected_events, drop = FALSE]
cases_ids   <- positional_to_plate[gel_insT_all_ids]

GROUP_UC <- "Unmatched controls"
GROUP_MC <- "Matched controls"
GROUP_IC <- "insT cases"

pca_psi <- rbind(psi_unctrl, psi_matched, psi_cases)
pca_group <- factor(
  c(rep(GROUP_UC, nrow(psi_unctrl)),
    rep(GROUP_MC, nrow(psi_matched)),
    rep(GROUP_IC, nrow(psi_cases))),
  levels = c(GROUP_UC, GROUP_MC, GROUP_IC))
pca_ids <- c(unctrl_ids, matched_ids, cases_ids)

# ---- PCA (complete cases only) ---------------------------------------------
complete_rows <- complete.cases(pca_psi)
cat(sprintf("Panel B PCA: %d/%d samples complete at all K=%d sites\n",
            sum(complete_rows), nrow(pca_psi), K_sel))
pca_fit  <- prcomp(pca_psi[complete_rows, , drop = FALSE], scale. = TRUE)
pct_var  <- round(summary(pca_fit)$importance[2, 1:2] * 100, 1L)

pca_df <- data.frame(
  PC1       = pca_fit$x[, 1L],
  PC2       = pca_fit$x[, 2L],
  group     = pca_group[complete_rows],
  sample_id = pca_ids[complete_rows],
  stringsAsFactors = FALSE)

# ---- Join metadata ---------------------------------------------------------
meta_want <- intersect(c(MANIFEST_KEY, MANIFEST_SEX, MANIFEST_YOB, MANIFEST_MAP),
                       colnames(manifest))
meta <- manifest[, meta_want, drop = FALSE]
colnames(meta)[colnames(meta) == MANIFEST_KEY] <- "sample_id"
if (MANIFEST_SEX %in% colnames(meta)) colnames(meta)[colnames(meta) == MANIFEST_SEX] <- "sex"
if (MANIFEST_YOB %in% colnames(meta)) colnames(meta)[colnames(meta) == MANIFEST_YOB] <- "age"
if (MANIFEST_MAP %in% colnames(meta)) colnames(meta)[colnames(meta) == MANIFEST_MAP] <- "mapped_reads"
pca_df <- merge(pca_df, meta, by = "sample_id", all.x = TRUE)

# ---- Bin covariates (quartiles computed from data) -------------------------
make_quartile_bins <- function(x, suffix = "") {
  vals   <- x[!is.na(x)]
  breaks <- unique(quantile(vals, probs = seq(0, 1, 0.25), na.rm = TRUE))
  if (length(breaks) < 2L) return(rep(NA_character_, length(x)))
  labs <- sprintf("%g-%g%s",
                  round(breaks[-length(breaks)]),
                  round(breaks[-1]),
                  suffix)
  labs[length(labs)] <- sprintf(">%g%s", round(breaks[length(breaks) - 1L]), suffix)
  cut(x, breaks = breaks, labels = labs, include.lowest = TRUE)
}

if ("age" %in% colnames(pca_df)) {
  pca_df$age_bin <- make_quartile_bins(pca_df$age)
}
if ("mapped_reads" %in% colnames(pca_df)) {
  pca_df$reads_bin <- make_quartile_bins(pca_df$mapped_reads / 1e6, suffix = "M")
}
if ("sex" %in% colnames(pca_df)) {
  pca_df$sex <- factor(pca_df$sex)
}

# ---- Plot helper -----------------------------------------------------------
# Unmatched controls: filled square (15); matched: filled circle (16);
# insT cases: filled triangle (17). Sizes kept equal for matched/cases so
# triangles don't appear disproportionately large.
PCA_SHAPES <- c(15L, 16L, 17L)
names(PCA_SHAPES) <- c(GROUP_UC, GROUP_MC, GROUP_IC)
PCA_SIZES  <- c(0.3, 1.2, 1.2)
names(PCA_SIZES)  <- c(GROUP_UC, GROUP_MC, GROUP_IC)
PCA_ALPHAS <- c(0.25, 0.9, 1.0)
names(PCA_ALPHAS) <- c(GROUP_UC, GROUP_MC, GROUP_IC)

# n-count labels for the shared group legend.
n_per_group  <- table(pca_group[complete_rows])
group_labs_n <- setNames(
  sprintf("%s (n=%s)", c(GROUP_UC, GROUP_MC, GROUP_IC),
          format(as.integer(n_per_group[c(GROUP_UC, GROUP_MC, GROUP_IC)]),
                 big.mark = ",")),
  c(GROUP_UC, GROUP_MC, GROUP_IC))

# add_group_legend: applied to the last PCA plot to surface the shared group
# shape legend at the bottom via ggplot2 >= 3.5 per-guide positions.
add_group_legend <- function(p) {
  p +
    scale_shape_manual(values = PCA_SHAPES, labels = group_labs_n, name = NULL) +
    guides(shape = guide_legend(direction = "horizontal", nrow = 1L,
                                override.aes = list(size = 2.2, colour = "grey30")))
}

make_pca_plot <- function(df, colour_col, colour_lab, colour_scale = NULL) {
  # Drop rows with NA for this covariate so they don't appear as grey
  df_ok  <- df[!is.na(df[[colour_col]]), , drop = FALSE]
  # Render unmatched controls first (bottom layer)
  df_ord <- df_ok[order(df_ok$group, decreasing = TRUE), ]
  p <- ggplot(df_ord,
              aes(x = PC1, y = PC2,
                  shape = group, colour = .data[[colour_col]],
                  size  = group, alpha = group)) +
    geom_point() +
    scale_size_manual( values = PCA_SIZES,  guide = "none") +
    scale_alpha_manual(values = PCA_ALPHAS, guide = "none") +
    labs(x      = sprintf("PC1 (%.1f%%)", pct_var[1]),
         y      = sprintf("PC2 (%.1f%%)", pct_var[2]),
         colour = colour_lab) +
    pub_theme +
    theme(legend.key.size = unit(3, "mm"),
          legend.text     = element_text(size = 9),
          legend.title    = element_text(size = 9))
  if (!is.null(colour_scale)) p <- p + colour_scale
  p
}

# Build one plot per available covariate
pca_plot_list <- list()

if ("sex" %in% colnames(pca_df)) {
  pca_plot_list[["sex"]] <- make_pca_plot(
    pca_df, "sex", "Sex",
    scale_colour_manual(
      values = c(Female = "#E64B35", Male = "#4DBBD5", Other = "#7E6148")))
}
if ("age_bin" %in% colnames(pca_df)) {
  age_pal <- colorRampPalette(c("#fdae6b", "#e6550d", "#7f2704"))(
    nlevels(pca_df$age_bin))
  pca_plot_list[["age"]] <- make_pca_plot(
    pca_df, "age_bin", "Age (years)",
    scale_colour_manual(values = setNames(age_pal, levels(pca_df$age_bin))))
}
if ("reads_bin" %in% colnames(pca_df)) {
  reads_pal <- colorRampPalette(c("#9ecae1", "#3182bd", "#08519c"))(
    nlevels(pca_df$reads_bin))
  pca_plot_list[["reads"]] <- make_pca_plot(
    pca_df, "reads_bin", "Mapped reads",
    scale_colour_manual(values = setNames(reads_pal, levels(pca_df$reads_bin))))
}

n_pca_plots <- length(pca_plot_list)
if (n_pca_plots == 0L) {
  warning("No covariate columns found in manifest — Panel B skipped")
  panel_b_pw <- NULL
} else {
  for (nm in names(pca_plot_list)) {
    pca_plot_list[[nm]] <- pca_plot_list[[nm]] +
      scale_shape_manual(values = PCA_SHAPES, labels = group_labs_n, name = NULL) +
      if (nm == "age") guides(shape = guide_legend(
                                position  = "bottom",
                                direction = "horizontal",
                                nrow      = 1L,
                                override.aes = list(size = 2.2, colour = "grey30")))
      else guides(shape = "none")
  }
  pca_plot_list[[1L]] <- pca_plot_list[[1L]] + labs(tag = "B") +
    theme(plot.tag = element_text(size = 14, face = "bold"),
          plot.tag.position = "topleft")
  panel_b_pw <- Reduce(`+`, pca_plot_list) + plot_layout(nrow = 1L, guides = "keep")
  cat(sprintf("Panel B: PCA on %d samples (PC1=%.1f%%, PC2=%.1f%%); %d covariate plot(s)\n",
              sum(complete_rows), pct_var[1], pct_var[2], n_pca_plots))
}

# Layout: 170 x 160 mm — ggpubr::ggarrange for the outer assembly (this
# figure is an exception to the "patchwork only" rule: patchwork was
# forcing Panel C to align its axis against Panel B's leftmost PCA below,
# breaking the intended layout).
#   Row 1 = A (75%) | C (25%) — per-site PSI facets | FPR vs n sites
#   Row 2 = B (full width)    — PCA confounding panels (sex | age | reads)
tag_theme <- theme(plot.tag = element_text(size = 14, face = "bold"),
                   plot.tag.position = "topleft")
top_row <- ggpubr::ggarrange(
  fig_s2a + labs(tag = "A") + tag_theme,
  fig_s2c + labs(tag = "C") + tag_theme,
  ncol = 2L, widths = c(3, 1))

if (!is.null(panel_b_pw)) {
  full_fig <- ggpubr::ggarrange(top_row, panel_b_pw,
                                 ncol = 1L, heights = c(1.2, 0.7))
  fig_height <- 160
} else {
  full_fig <- top_row
  fig_height <- 80
}

pdf(OUT_FIG_S4, width = mm_to_in(170), height = mm_to_in(fig_height))
print(full_fig)
dev.off()
cat(sprintf("Supplementary Figure S4 saved to: %s\n", OUT_FIG_S4))

# ============================================================================
# Summary RDS (small numbers, airlock-safe)
# ============================================================================
summary_out <- list(
  K_selected             = K_sel,
  n_6k_after_dedup       = nrow(ctrl_scores_dedup),
  n_dev_overlap          = length(overlap_keys),
  n_biallelic_in_6k      = length(biallelic_in_6k),
  n_hd_in_6k             = length(hd_in_6k),
  # Panel A — per-site PSI by group
  panel_a_n              = c(n_6k_ctrl  = nrow(ctrl_6k_psi),
                             n_dev_ctrl = nrow(dev_ctrl_psi),
                             n_dev_case = nrow(dev_case_psi)),
  panel_a_median_psi     = vapply(split(psi_long, psi_long$group),
                                  function(d) median(d$psi, na.rm = TRUE),
                                  numeric(1L)),
  panel_a_q_dev_ctrl_vs_6k = setNames(panel_c_tests$q_ctrl_vs_6k,
                                      as.character(panel_c_tests$gene_lab)),
  panel_a_q_dev_case_vs_6k = setNames(panel_c_tests$q_case_vs_6k,
                                      as.character(panel_c_tests$gene_lab)),
  # Panel C — FPR by number of sites aggregated
  panel_c_n_ctrl         = n_ctrl_clean,
  panel_c_combo_fpr      = combo_fpr_df,
  panel_c_combo_fpr_raw  = combo_fpr_raw_df,
  # Panel B — PCA confounding assessment
  panel_b_n_complete     = if (exists("pca_df")) sum(complete.cases(pca_df[, c("PC1","PC2")])) else NA_integer_,
  panel_b_pct_var        = if (exists("pct_var")) pct_var else NULL,
  panel_b_n_groups       = if (exists("pca_group")) table(pca_group[complete_rows]) else NULL
)
saveRDS(summary_out, OUT_SUMMARY)
cat(sprintf("Summary saved to: %s\n", OUT_SUMMARY))

cat("\n=== Done. Headline numbers ===\n")
cat("Panel A (per-group n + median PSI across all sites):\n")
panel_a_n_named <- summary_out$panel_a_n
names(panel_a_n_named) <- b_levels
cat(sprintf("  %s\n",
            paste(names(panel_a_n_named),
                  sprintf("n=%s median_psi=%.2f",
                          format(panel_a_n_named, big.mark = ","),
                          summary_out$panel_a_median_psi[b_levels]),
                  sep = " ", collapse = " | ")))
n_tests <- length(summary_out$panel_a_q_dev_ctrl_vs_6k) * 2L
cat(sprintf("  Per-site MWU q (BH/%d): dev-ctrl-vs-6k %d/%d sig at q<0.05; dev-insT-case-vs-6k %d/%d sig\n",
            n_tests,
            sum(summary_out$panel_a_q_dev_ctrl_vs_6k < 0.05, na.rm = TRUE),
            length(summary_out$panel_a_q_dev_ctrl_vs_6k),
            sum(summary_out$panel_a_q_dev_case_vs_6k < 0.05, na.rm = TRUE),
            length(summary_out$panel_a_q_dev_case_vs_6k)))
cat(sprintf("Panel C: FPR by n sites (n_ctrl=%s) — k=1 mean=%.1f%% (range %.1f-%.1f%%), k=%d mean=%.1f%% (min=%d/%d)\n",
            format(n_ctrl_clean, big.mark = ","),
            combo_fpr_df$mean_fpr[1],
            combo_fpr_df$min_fpr[1],
            combo_fpr_df$max_fpr[1],
            K_sel,
            combo_fpr_df$mean_fpr[K_sel],
            as.integer(combo_fpr_df$min_n_fp[K_sel]),
            n_ctrl_clean))
cat("\nNOTE: 6k cohort orient values are GEL-protected; airlock the PDF + summary RDS.\n")
