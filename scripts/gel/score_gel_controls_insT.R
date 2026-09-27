#!/usr/bin/env Rscript
#
# score_gel_controls_insT.R
# Apply the insT-specific LASSO uniform-sum signature to the ~6k GEL whole-
# blood control cohort.
#
# Score formula (matches loo_evaluate_insT_lasso.R exactly):
#   bb_logit(ijc, sjc) = log( ((ijc+0.5) / (ijc+sjc+1)) /
#                              (1 - (ijc+0.5) / (ijc+sjc+1)) )
#   feature_ij = direction_j * (bb_logit_ij - mu_ctrl_GEL_j)
#   score_i    = sum_j feature_ij                  (j in selected events)
# mu_ctrl_GEL is taken from the insT RDS (computed on the 33 GEL training
# controls). NA cells (IJC+SJC == 0 or missing) are excluded; never imputed.
#
# Run from workspace root: Rscript scripts/gel/score_gel_controls_insT.R
#
# Inputs:
#   outputs/loo_results_insT_lasso.rds       - selected sites + gel_ctrl_means
#   outputs/gel_controls_psi_cand_pool.rds   - 6k matrix for full candidate pool
#                                              (build via combine_rmats_batches_cand_pool.R)
#   outputs/psi_matrix.rds                   - dev cohort PSI/IJC/SJC
#   data/control_manifest.csv
#   data/rmats/inputs/renu_11_cases_paths.txt
#   data/rmats/inputs/renu_11_cases_control_paths.txt
#   data/ndd_carriers.csv                    - optional
#   data/biallelic_rnu42_platekeys.txt       - optional
#   data/renu_highdepth_case_path.txt        - optional
#
# Outputs (subject to GEL airlock):
#   outputs/gel_controls_insT_robustness.pdf
#   outputs/gel_controls_insT_robustness.rds
#   outputs/top_ctrl_insT.csv
#
# DATA PROTECTION: insT selected event identifiers + directions are GEL-
# protected. The PDF uses anonymous FP-N labels; top-N CSV ships platekeys
# only after airlock review.

local({
  v <- substr(as.character(getRversion()), 1, 3)
  candidate <- switch(
    v,
    "4.0" = "/tools/aws-workspace-apps/ce/R/4.0.2",
    "4.3" = "/tools/aws-workspace-ubuntu-apps/ce/R/4.3.3",
    "4.5" = "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2",
    NULL
  )
  if (!is.null(candidate) && dir.exists(candidate)) {
    .libPaths(c(.libPaths(), candidate))
  }
})
suppressPackageStartupMessages({
  library(ggplot2)
})

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
INST_RDS       <- "outputs/loo_results_insT_lasso.rds"
CTRL_MAT       <- "outputs/gel_controls_psi_cand_pool.rds"
DEV_PSI        <- "outputs/psi_matrix.rds"
MANIFEST       <- "data/control_manifest.csv"
DEV_CASE_BAMS  <- "data/rmats/inputs/renu_11_cases_paths.txt"
DEV_CTRL_BAMS  <- "data/rmats/inputs/renu_11_cases_control_paths.txt"
HD_CASE_FILE   <- "data/renu_highdepth_case_path.txt"
NDD_FILE       <- "data/ndd_carriers.csv"
BIALLELIC_FILE <- "data/biallelic_rnu42_platekeys.txt"
OUT_PDF        <- "outputs/gel_controls_insT_robustness.pdf"
OUT_RDS        <- "outputs/gel_controls_insT_robustness.rds"
TOP_CSV        <- "outputs/top_ctrl_insT.csv"

SCORE_PARITY_TOL <- 1e-6

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}

bblogit_at_events <- function(ijc_mat, sjc_mat) {
  out <- bb_logit(ijc_mat, sjc_mat)
  out[(ijc_mat + sjc_mat) == 0] <- NA_real_
  out
}

score_insT <- function(bblogit_mat, mu_ctrl, direction) {
  oriented <- sweep(bblogit_mat, 2L, mu_ctrl,   "-")
  oriented <- sweep(oriented,    2L, direction, "*")
  data.frame(
    sample_id = rownames(bblogit_mat),
    score     = rowSums(oriented, na.rm = TRUE),
    n_events  = rowSums(!is.na(bblogit_mat)),
    stringsAsFactors = FALSE
  )
}

lasso_prob <- function(bblogit_mat, mu_ctrl, direction, beta, intercept) {
  X <- sweep(bblogit_mat, 2L, mu_ctrl,   "-")
  X <- sweep(X,           2L, direction, "*")
  X[is.na(X)] <- 0  # missing event -> no deviation from ctrl mean (consistent with rowSums na.rm=TRUE)
  lp <- as.numeric(X %*% beta + intercept)
  setNames(1 / (1 + exp(-lp)), rownames(bblogit_mat))
}

read_bam_list <- function(path) {
  raw  <- paste(readLines(path), collapse = ",")
  bams <- trimws(unlist(strsplit(raw, ",")))
  bams <- bams[nzchar(bams)]
  sub("\\.bam$", "", basename(bams))
}

# --------------------------------------------------------------------------
# Load insT RDS, extract selected events + directions + ctrl means
# --------------------------------------------------------------------------
ins <- readRDS(INST_RDS)
stopifnot(
  "insT RDS missing a5ss_cand"      = "a5ss_cand"      %in% names(ins),
  "insT RDS missing selected_idx"   = "selected_idx"    %in% names(ins),
  "insT RDS missing gel_ctrl_means" = "gel_ctrl_means"  %in% names(ins),
  "insT RDS missing gel_score"      = "gel_score"       %in% names(ins),
  "insT RDS missing K_selected"     = "K_selected"      %in% names(ins)
)

sel_meta        <- ins$a5ss_cand[ins$selected_idx, , drop = FALSE]
selected_events <- sel_meta$gel_ev_id
direction_sel   <- sel_meta$direction
mu_ctrl_sel     <- ins$gel_ctrl_means[selected_events]
K_sel           <- ins$K_selected

stopifnot(
  "selected_events / sel_meta mismatch"  = length(selected_events) == nrow(sel_meta),
  "Direction length mismatch"            = length(direction_sel)   == K_sel,
  "ctrl mean / selected length mismatch" = length(mu_ctrl_sel)     == K_sel,
  "ctrl means contain NA"               = !anyNA(mu_ctrl_sel),
  "Direction has zeros"                  = !any(direction_sel == 0)
)
message("insT-LASSO: K=", K_sel, " selected events")

# Identify insT GEL case positional IDs (train + holdout)
gel_insT_train_ids <- ins$sample_ids_train[ins$cohort_train == "GEL" &
                                            ins$y_train == 1L]
gel_insT_ho_ids    <- ins$ho_scores$gel$sample_id
gel_insT_all_ids   <- c(gel_insT_train_ids, gel_insT_ho_ids)
message("insT GEL cases: ", length(gel_insT_all_ids),
        " (train=", length(gel_insT_train_ids),
        ", holdout=", length(gel_insT_ho_ids), ")")

# --------------------------------------------------------------------------
# Load 6k pool matrix + dev matrix + manifest
# --------------------------------------------------------------------------
ctrl     <- readRDS(CTRL_MAT)
dev      <- readRDS(DEV_PSI)
manifest <- read.csv(MANIFEST, stringsAsFactors = FALSE)

stopifnot(
  "6k matrix missing $psi/$ijc/$sjc - run combine_rmats_batches_cand_pool.R" =
    all(c("psi", "ijc", "sjc", "samples", "events") %in% names(ctrl)),
  "Dev matrix missing $ijc/$sjc - rebuild via build_psi_matrix.R" =
    all(c("ijc", "sjc") %in% names(dev))
)

events_missing <- setdiff(selected_events, colnames(ctrl$psi))
if (length(events_missing) > 0L) {
  stop("6k pool matrix is missing ", length(events_missing),
       " insT-selected event(s).\n  Missing: ",
       paste(events_missing, collapse = "; "),
       "\n  Re-run combine_rmats_batches_cand_pool.R (check that psi_matrix.rds",
       " and nava_psi_matrix.rds are current).")
}
message("All ", K_sel, " selected events present in 6k pool matrix")

# --------------------------------------------------------------------------
# SANITY: recompute all GEL dev scores from psi_matrix.rds and verify they
# match ins$gel_score. Catches silent drift in feature construction.
# --------------------------------------------------------------------------
dev_ijc_sel <- dev$ijc[, selected_events, drop = FALSE]
dev_sjc_sel <- dev$sjc[, selected_events, drop = FALSE]
dev_bblogit <- bblogit_at_events(dev_ijc_sel, dev_sjc_sel)
dev_recomp  <- score_insT(dev_bblogit, mu_ctrl_sel, direction_sel)

dev_recomp_named <- setNames(dev_recomp$score, dev_recomp$sample_id)
shared_dev       <- intersect(names(ins$gel_score), names(dev_recomp_named))
stopifnot("No shared sample IDs between ins$gel_score and dev recompute" =
            length(shared_dev) >= 1L)
parity_diff <- max(abs(ins$gel_score[shared_dev] - dev_recomp_named[shared_dev]),
                   na.rm = TRUE)
message(sprintf("Dev cohort score parity: max |delta| = %.3g (n=%d shared)",
                parity_diff, length(shared_dev)))
if (parity_diff > SCORE_PARITY_TOL) {
  warning("Dev scores recomputed here drift from ins$gel_score by >",
          SCORE_PARITY_TOL,
          ". Likely cause: psi_matrix.rds was rebuilt after the insT RDS. ",
          "Re-run loo_evaluate_insT_lasso.R before this script.")
}

# --------------------------------------------------------------------------
# Relabel dev$samples with real platekeys
# --------------------------------------------------------------------------
dev_case_keys <- read_bam_list(DEV_CASE_BAMS)
dev_ctrl_keys <- read_bam_list(DEV_CTRL_BAMS)
stopifnot(
  "Expected 11 dev case BAMs"    = length(dev_case_keys) == 11L,
  "Expected 33 dev control BAMs" = length(dev_ctrl_keys) == 33L,
  "dev$samples not in case-then-control order" =
    identical(tolower(dev$samples$group), c(rep("case", 11), rep("control", 33)))
)
dev_real_ids        <- c(dev_case_keys, dev_ctrl_keys)
positional_to_plate <- setNames(dev_real_ids, rownames(dev$psi))

# --------------------------------------------------------------------------
# Score 6k matrix BEFORE dedup
# --------------------------------------------------------------------------
ctrl_ijc_sel    <- ctrl$ijc[, selected_events, drop = FALSE]
ctrl_sjc_sel    <- ctrl$sjc[, selected_events, drop = FALSE]
ctrl_bblogit    <- bblogit_at_events(ctrl_ijc_sel, ctrl_sjc_sel)
ctrl_scores_all <- score_insT(ctrl_bblogit, mu_ctrl_sel, direction_sel)
# LASSO P(case) for every 6k sample: sigmoid(X %*% beta + intercept)
beta_sel        <- ins$fit_at_1se$beta[selected_events]
intercept_insT  <- ins$fit_at_1se$intercept
ctrl_prob_all   <- lasso_prob(ctrl_bblogit, mu_ctrl_sel, direction_sel,
                              beta_sel, intercept_insT)

# --------------------------------------------------------------------------
# Paired-score reproducibility for dev/6k overlap
# --------------------------------------------------------------------------
all_dev_keys    <- c(dev_case_keys, dev_ctrl_keys)
overlap_case_pk <- intersect(rownames(ctrl_bblogit), dev_case_keys)
overlap_ctrl_pk <- intersect(rownames(ctrl_bblogit), dev_ctrl_keys)
overlap_all_pk  <- c(overlap_case_pk, overlap_ctrl_pk)
message("Dev/6k overlap: ", length(overlap_case_pk), " case(s) + ",
        length(overlap_ctrl_pk), " control(s)")

# Identify which overlap cases are insT (by matching platekeys back to positional IDs)
plate_to_positional  <- setNames(rownames(dev$psi), dev_real_ids)
overlap_case_positional <- plate_to_positional[overlap_case_pk]
overlap_insT_pk <- overlap_case_pk[overlap_case_positional %in% gel_insT_all_ids]

paired <- NULL
if (length(overlap_all_pk) > 0L) {
  s_6k <- ctrl_scores_all[match(overlap_all_pk, ctrl_scores_all$sample_id), ]
  dev_recomp_pk <- dev_recomp
  dev_recomp_pk$sample_id <- unname(positional_to_plate[dev_recomp_pk$sample_id])
  s_dev <- dev_recomp_pk[match(overlap_all_pk, dev_recomp_pk$sample_id), ]
  paired <- data.frame(
    sample_id = overlap_all_pk,
    group     = c(rep("insT case",       length(overlap_case_pk)),
                  rep("GEL dev control", length(overlap_ctrl_pk))),
    score_dev = s_dev$score,
    score_6k  = s_6k$score,
    n_dev     = s_dev$n_events,
    n_6k      = s_6k$n_events,
    stringsAsFactors = FALSE
  )
  # Override label: non-insT cases in overlap aren't the primary comparison
  paired$group[paired$sample_id %in% overlap_case_pk &
                 !paired$sample_id %in% overlap_insT_pk] <- "non-insT case"
  paired$diff <- paired$score_6k - paired$score_dev
  pair_r    <- cor(paired$score_dev, paired$score_6k, use = "complete.obs")
  pair_mad  <- mean(abs(paired$diff), na.rm = TRUE)
  pair_bias <- mean(paired$diff,      na.rm = TRUE)
  pair_sd   <- sd(  paired$diff,      na.rm = TRUE)
  loa_lo    <- pair_bias - 1.96 * pair_sd
  loa_hi    <- pair_bias + 1.96 * pair_sd
  message(sprintf("\nPaired-score: r=%.3f  MAD=%.3f  bias=%.3f  95%% LoA=[%.2f, %.2f]",
                  pair_r, pair_mad, pair_bias, loa_lo, loa_hi))
}

# --------------------------------------------------------------------------
# Dedup dev cohort from 6k matrix
# --------------------------------------------------------------------------
keep_rows    <- !rownames(ctrl_bblogit) %in% overlap_all_pk
ctrl_bblogit <- ctrl_bblogit[keep_rows, , drop = FALSE]
ctrl$samples <- ctrl$samples[keep_rows, , drop = FALSE]
ctrl_scores  <- ctrl_scores_all[keep_rows, , drop = FALSE]
ctrl_prob    <- ctrl_prob_all[keep_rows]
message("After dedup: ", nrow(ctrl_scores), " 6k samples scored")

# --------------------------------------------------------------------------
# Build dev_scores for plotting - use ins$gel_score as authoritative
# --------------------------------------------------------------------------
dev_scores <- data.frame(
  sample_id = names(ins$gel_score),
  score     = unname(ins$gel_score),
  group     = ifelse(names(ins$gel_score) %in% gel_insT_all_ids,  "insT case",
              ifelse(grepl("^gel:control:", names(ins$gel_score)), "GEL dev control",
                     "non-insT case")),
  stringsAsFactors = FALSE
)
dev_insT_mask  <- dev_scores$group == "insT case"
dev_ctrl_mask  <- dev_scores$group == "GEL dev control"
dev_insT_scores <- dev_scores$score[dev_insT_mask]
dev_ctrl_scores <- dev_scores$score[dev_ctrl_mask]

stopifnot(
  "No insT cases found in ins$gel_score"      = sum(dev_insT_mask) >= 1L,
  "No dev controls found in ins$gel_score"    = sum(dev_ctrl_mask) >= 1L
)
message("Dev insT cases: ", sum(dev_insT_mask),
        "  Dev controls: ", sum(dev_ctrl_mask),
        "  Non-insT cases: ", sum(dev_scores$group == "non-insT case"))

# --------------------------------------------------------------------------
# Glass's Delta: 6k vs dev controls as reference; insT cases as signal
# --------------------------------------------------------------------------
mu_insT <- mean(dev_insT_scores)
mu_dev  <- mean(dev_ctrl_scores)
mu_6k   <- mean(ctrl_scores$score)
sd_dev  <- sd(  dev_ctrl_scores)
sd_6k   <- sd(  ctrl_scores$score)
glass_dev <- (mu_insT - mu_dev) / sd_dev
glass_6k  <- (mu_insT - mu_6k)  / sd_6k
message("\nGlass's Delta - dev (33 ctrls):  ", round(glass_dev, 3))
message("Glass's Delta - 6k cohort:        ", round(glass_6k,  3))
message("SD ratio (6k / dev):              ", round(sd_6k / sd_dev, 3))

# --------------------------------------------------------------------------
# FP count: 6k controls scoring above min(insT case score)
# --------------------------------------------------------------------------
min_insT_score <- min(dev_insT_scores, na.rm = TRUE)

hd_platekey <- character(0)
if (file.exists(HD_CASE_FILE)) {
  hd_platekey <- read_bam_list(HD_CASE_FILE)
  hd_platekey <- hd_platekey[nzchar(hd_platekey)]
  hd_platekey <- intersect(hd_platekey, ctrl_scores$sample_id)
  message(sprintf("High-depth ReNU case(s) found in 6k pool: %d", length(hd_platekey)))
}

fp_mask <- ctrl_scores$score > min_insT_score &
           !ctrl_scores$sample_id %in% hd_platekey
fp_ids  <- ctrl_scores$sample_id[fp_mask]
n_fp    <- length(fp_ids)
message(sprintf("\nMin insT case score: %.3f", min_insT_score))
message(sprintf("6k controls > min(insT case): %d%s",
                n_fp,
                if (length(hd_platekey) > 0L)
                  sprintf(" (excluding %d HD ReNU case)", length(hd_platekey))
                else ""))

# --------------------------------------------------------------------------
# Attach manifest covariates
# --------------------------------------------------------------------------
manifest$platekey <- sub("\\.bam$", "", basename(manifest$file_path))
m <- manifest[match(ctrl_scores$sample_id, manifest$platekey), ]
ctrl_scores$participant_id   <- m$participant_id
ctrl_scores$carrier_evidence <- if ("spliceopathy_evidence" %in% names(m))
  m$spliceopathy_evidence else NA_character_
ctrl_scores$carrier_gene     <- if ("spliceopathy_gene" %in% names(m))
  m$spliceopathy_gene else NA_character_

# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11))
pdf(OUT_PDF, width = 8.5, height = 6)

# 1. Score distribution: 6k controls vs insT cases vs dev controls
# Non-insT GEL cases shown as a lighter overlay for context.
p1_df <- rbind(
  data.frame(group = "GEL 6k controls",   score = ctrl_scores$score,
             stringsAsFactors = FALSE),
  data.frame(group = dev_scores$group,     score = dev_scores$score,
             stringsAsFactors = FALSE)
)
p1_df$group <- factor(p1_df$group,
                      levels = c("GEL 6k controls", "GEL dev control",
                                 "non-insT case",   "insT case"))

p1 <- ggplot(p1_df, aes(x = score, fill = group)) +
  geom_histogram(data = subset(p1_df, group == "GEL 6k controls"),
                 bins = 80, alpha = 0.7) +
  geom_vline(data = subset(p1_df, group != "GEL 6k controls"),
             aes(xintercept = score, colour = group), linewidth = 0.4) +
  geom_vline(xintercept = min_insT_score, colour = "black",
             linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c("GEL 6k controls" = "grey60")) +
  scale_colour_manual(values = c("GEL dev control" = "#4DBBD5",
                                 "non-insT case"   = "#F39B7F",
                                 "insT case"       = "#E64B35")) +
  labs(title    = sprintf("insT-LASSO score - 6k unmatched control cohort (K=%d)", K_sel),
       subtitle = sprintf("FPs above min(insT case)=%.2f: %d  |  Glass d 6k=%.2f, dev=%.2f  |  SD ratio=%.2f",
                          min_insT_score, n_fp, glass_6k, glass_dev, sd_6k / sd_dev),
       x = "Signature score (sum of direction-oriented logit-PSI deviations)",
       y = "Count", fill = NULL, colour = NULL)
print(p1)

# 2. Carrier overlay helper
highlight_plot <- function(plate_keys, label, hl_colour, sdd_keys = NULL) {
  plate_keys <- unique(plate_keys[!is.na(plate_keys) & nzchar(plate_keys)])
  idx        <- match(plate_keys, ctrl_scores$sample_id)
  hit_mask   <- !is.na(idx)
  hit_scores <- ctrl_scores$score[idx[hit_mask]]
  n_in       <- sum(hit_mask)
  n_total    <- length(plate_keys)
  message(sprintf("  %-30s solved %d/%d in 6k scored cohort", label, n_in, n_total))

  sdd_keys <- unique(sdd_keys[!is.na(sdd_keys) & nzchar(sdd_keys)])
  sdd_keys <- setdiff(sdd_keys, plate_keys)
  sdd_idx  <- match(sdd_keys, ctrl_scores$sample_id)
  sdd_mask <- !is.na(sdd_idx)
  sdd_scores_in <- ctrl_scores$score[sdd_idx[sdd_mask]]
  sdd_n_in <- sum(sdd_mask); sdd_n_total <- length(sdd_keys)
  if (sdd_n_total > 0L) {
    message(sprintf("  %-30s sdd    %d/%d in 6k scored cohort", label, sdd_n_in, sdd_n_total))
  }
  if (n_in == 0L && sdd_n_in == 0L) return(invisible(NULL))

  hit_scores    <- hit_scores[!is.na(hit_scores)]
  sdd_scores_in <- sdd_scores_in[!is.na(sdd_scores_in)]
  if (length(hit_scores) == 0L && length(sdd_scores_in) == 0L) return(invisible(NULL))

  y_top <- max(ggplot_build(p1)$data[[1]]$count, na.rm = TRUE)
  hit_df <- data.frame(score = hit_scores,
                       y     = rep(y_top * 1.05, length(hit_scores)),
                       stringsAsFactors = FALSE)
  sdd_df <- data.frame(score = sdd_scores_in,
                       y     = rep(y_top * 0.95, length(sdd_scores_in)),
                       stringsAsFactors = FALSE)

  p <- ggplot(p1_df, aes(x = score, fill = group)) +
    geom_histogram(data = subset(p1_df, group == "GEL 6k controls"),
                   bins = 80, alpha = 0.5) +
    geom_vline(data = subset(p1_df, group != "GEL 6k controls"),
               aes(xintercept = score, colour = group),
               linewidth = 0.3, alpha = 0.4) +
    geom_vline(xintercept = min_insT_score, colour = "black",
               linetype = "dashed", linewidth = 0.3) +
    scale_fill_manual(values = c("GEL 6k controls" = "grey70")) +
    scale_colour_manual(values = c("GEL dev control" = "#4DBBD5",
                                   "non-insT case"   = "#F39B7F",
                                   "insT case"       = "#E64B35"))
  if (length(sdd_scores_in) > 0L) {
    p <- p +
      geom_vline(xintercept = sdd_scores_in, colour = hl_colour,
                 linewidth = 0.5, alpha = 0.35) +
      geom_point(data = sdd_df, aes(x = score, y = y), inherit.aes = FALSE,
                 colour = hl_colour, fill = hl_colour,
                 shape = 25, size = 2.2, alpha = 0.4)
  }
  if (length(hit_scores) > 0L) {
    p <- p +
      geom_vline(xintercept = hit_scores, colour = hl_colour, linewidth = 1.0) +
      geom_point(data = hit_df, aes(x = score, y = y), inherit.aes = FALSE,
                 colour = hl_colour, fill = hl_colour, shape = 25, size = 3.2)
  }
  subtitle <- sprintf("solved: %d of %d", n_in, n_total)
  if (sdd_n_total > 0L) subtitle <- paste0(subtitle, sprintf("  |  sdd: %d of %d", sdd_n_in, sdd_n_total))
  p <- p + labs(title    = sprintf("insT-LASSO - %s highlighted", label),
                subtitle = subtitle,
                x = "Signature score", y = "Count", fill = NULL, colour = NULL)
  print(p)
  invisible(data.frame(label = label, n_in_solved = n_in, n_total_solved = n_total,
                       n_in_sdd = sdd_n_in, n_total_sdd = sdd_n_total,
                       stringsAsFactors = FALSE))
}

message("\nCarrier highlight overlays:")
highlight_summary <- list()

if (file.exists(NDD_FILE)) {
  ndd <- read.csv(NDD_FILE, stringsAsFactors = FALSE)
  if (nrow(ndd) > 0L && all(c("hgnc_symbol", "rnaseq_platekey") %in% names(ndd))) {
    if (!"evidence" %in% names(ndd)) { warning("ndd_carriers.csv has no 'evidence' column"); ndd$evidence <- "tiered" }
    ndd_cfg <- list(
      TCF4  = list(label = "Pitt-Hopkins (TCF4)", colour = "#984EA3"),
      MECP2 = list(label = "Rett (MECP2)",        colour = "#FF7F00"),
      UBE3A = list(label = "Angelman (UBE3A)",    colour = "#33A02C")
    )
    for (g in names(ndd_cfg)) {
      sub_ndd  <- ndd[ndd$hgnc_symbol == g, , drop = FALSE]
      pk_solved <- unique(sub_ndd$rnaseq_platekey[sub_ndd$evidence == "solved"])
      pk_sdd    <- setdiff(unique(sub_ndd$rnaseq_platekey[sub_ndd$evidence == "sdd"]), pk_solved)
      if (length(pk_solved) == 0L && length(pk_sdd) == 0L) next
      res <- highlight_plot(pk_solved, ndd_cfg[[g]]$label, ndd_cfg[[g]]$colour, sdd_keys = pk_sdd)
      if (!is.null(res)) highlight_summary[[g]] <- res
    }
  }
} else {
  message("  ", NDD_FILE, " not found - skipping NDD overlay")
}

if (file.exists(BIALLELIC_FILE)) {
  bi_keys <- trimws(readLines(BIALLELIC_FILE))
  bi_keys <- bi_keys[nzchar(bi_keys) & !startsWith(bi_keys, "#")]
  res <- highlight_plot(bi_keys, "biallelic RNU4-2", "#000000")
  if (!is.null(res)) highlight_summary[["biallelic_rnu42"]] <- res
} else {
  message("  ", BIALLELIC_FILE, " not found - skipping biallelic overlay")
}

# 3. Per-event coverage in 6k cohort
# Short label: GeneSymbol (deduplicated)
short_label_full <- make.unique(sel_meta$GeneSymbol, sep = "#")
names(short_label_full) <- selected_events

ctrl_bblogit_for_cov <- ctrl$ijc[, selected_events, drop = FALSE]
ctrl_sjc_for_cov     <- ctrl$sjc[, selected_events, drop = FALSE]
ctrl_bblogit_mat     <- bblogit_at_events(ctrl_bblogit_for_cov, ctrl_sjc_for_cov)

cov_df <- data.frame(
  event_id     = selected_events,
  short_label  = short_label_full,
  coverage     = colMeans(!is.na(ctrl_bblogit_mat)),
  n_total      = nrow(ctrl_bblogit_mat),
  median_depth = apply(ctrl_bblogit_for_cov + ctrl_sjc_for_cov, 2L, median, na.rm = TRUE),
  stringsAsFactors = FALSE
)
cov_df$n_na    <- round((1 - cov_df$coverage) * cov_df$n_total)
cov_df <- cov_df[order(cov_df$coverage), ]
cov_df$short_label <- factor(cov_df$short_label, levels = cov_df$short_label)
cov_df$bar_label <- sprintf("%s NA (%s%%)", format(cov_df$n_na, big.mark = ","),
                            round(100 * (1 - cov_df$coverage), 1))

print(
  ggplot(cov_df, aes(coverage, short_label)) +
    geom_col(fill = "grey50") +
    geom_text(aes(label = bar_label), hjust = -0.05, size = 3) +
    scale_x_continuous(labels = scales::percent_format(1), limits = c(0, 1.25)) +
    labs(title = "Per-event coverage in 6k cohort (insT-LASSO selected sites)",
         x = "Fraction of samples with non-zero IJC+SJC", y = NULL)
)
print(
  ggplot(cov_df, aes(median_depth, short_label)) +
    geom_col(fill = "grey50") +
    scale_x_log10() +
    labs(title = "Per-event median depth (IJC+SJC) in 6k cohort",
         x = "Median IJC+SJC (log10 scale)", y = NULL)
)

# 4. FP oriented delta heatmap (anonymous labels)
heat_ids <- c(hd_platekey,
              fp_ids[order(-ctrl_scores$score[match(fp_ids, ctrl_scores$sample_id)])])
if (length(heat_ids) > 0L) {
  heat_ijc <- ctrl$ijc[heat_ids, selected_events, drop = FALSE]
  heat_sjc <- ctrl$sjc[heat_ids, selected_events, drop = FALSE]
  heat_bb  <- bblogit_at_events(heat_ijc, heat_sjc)
  heat_oriented <- sweep(heat_bb, 2L, mu_ctrl_sel,   "-")
  heat_oriented <- sweep(heat_oriented, 2L, direction_sel, "*")

  row_score <- ctrl_scores$score[match(heat_ids, ctrl_scores$sample_id)]
  row_label <- character(length(heat_ids))
  fp_counter <- 0L
  for (ri in seq_along(heat_ids)) {
    if (heat_ids[ri] %in% hd_platekey) {
      row_label[ri] <- sprintf("ReNU high-depth (%.2f)", row_score[ri])
    } else {
      fp_counter <- fp_counter + 1L
      row_label[ri] <- sprintf("FP-%d (%.2f)", fp_counter, row_score[ri])
    }
  }

  # Column order by dev insT-case vs ctrl effect
  dev_insT_bb <- dev_bblogit[gel_insT_all_ids, , drop = FALSE]
  dev_ctrl_bb <- dev_bblogit[grepl("^gel:control:", rownames(dev_bblogit)), , drop = FALSE]
  event_effect <- direction_sel *
    (colMeans(dev_insT_bb, na.rm = TRUE) - colMeans(dev_ctrl_bb, na.rm = TRUE))
  names(event_effect) <- selected_events
  col_order <- short_label_full[names(sort(event_effect, decreasing = TRUE))]

  heat_df <- data.frame(
    sample_id   = rep(rownames(heat_oriented), times = ncol(heat_oriented)),
    short_label = rep(short_label_full[colnames(heat_oriented)],
                      each = nrow(heat_oriented)),
    oriented    = as.vector(heat_oriented),
    stringsAsFactors = FALSE
  )
  heat_df$sample_id   <- factor(heat_df$sample_id,
                                levels = rev(heat_ids),
                                labels = rev(row_label))
  heat_df$short_label <- factor(heat_df$short_label, levels = col_order)

  print(
    ggplot(heat_df, aes(x = short_label, y = sample_id, fill = oriented)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                           midpoint = 0, na.value = "grey85",
                           name = "direction * (bblogit - mu_ctrl)") +
      labs(title    = sprintf("Oriented bb-logit deviation - %d FP control(s)%s",
                              n_fp,
                              if (length(hd_platekey) > 0L) " + HD ReNU case" else ""),
           subtitle = "Red = shifted in disease direction; rows ordered by score desc.",
           x = NULL, y = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
            axis.text.y = element_text(size = 6))
  )
} else {
  message("\nNo FPs and no HD case - skipping FP heatmap")
}

# 5. Paired-score reproducibility
if (!is.null(paired)) {
  lim <- range(c(paired$score_dev, paired$score_6k), na.rm = TRUE)
  pair_pal <- c("insT case" = "#E64B35", "non-insT case" = "#F39B7F",
                "GEL dev control" = "#4DBBD5")
  print(
    ggplot(paired, aes(x = score_dev, y = score_6k, colour = group)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(size = 2.2, alpha = 0.85) +
      scale_colour_manual(values = pair_pal) +
      coord_fixed(xlim = lim, ylim = lim) +
      labs(title = "Paired-score reproducibility (same BAM, two rMATS runs)",
           subtitle = sprintf("r=%.3f  |  MAD=%.3f  |  bias=%.3f  |  95%% LoA=[%.2f, %.2f]",
                              pair_r, pair_mad, pair_bias, loa_lo, loa_hi),
           x = "Score (dev rMATS run)", y = "Score (6k batch rMATS run)", colour = NULL)
  )
  paired$mean_score <- (paired$score_dev + paired$score_6k) / 2
  print(
    ggplot(paired, aes(x = mean_score, y = diff, colour = group)) +
      geom_hline(yintercept = 0,         linetype = "dashed", colour = "grey50") +
      geom_hline(yintercept = pair_bias, colour = "grey30") +
      geom_hline(yintercept = c(loa_lo, loa_hi), linetype = "dotted", colour = "grey30") +
      geom_point(size = 2.2, alpha = 0.85) +
      scale_colour_manual(values = pair_pal) +
      labs(title = "Bland-Altman: 6k batch score minus dev score",
           subtitle = sprintf("bias=%.3f, 95%% LoA=[%.2f, %.2f]", pair_bias, loa_lo, loa_hi),
           x = "Mean score", y = "6k - dev", colour = NULL)
  )
}

# 6. Top-N high-scoring 6k controls
top_n <- 20L
top_ctrl <- ctrl_scores[order(-ctrl_scores$score), ][seq_len(top_n), ]
top_ctrl <- top_ctrl[, c("sample_id", "score", "n_events",
                         "participant_id", "carrier_evidence", "carrier_gene")]
message("\n--- Top ", top_n, " 6k controls by insT-LASSO score ---")
print(top_ctrl, row.names = FALSE)
write.csv(top_ctrl, TOP_CSV, row.names = FALSE)
message("Wrote ", TOP_CSV)

top_str <- capture.output(print(top_ctrl, row.names = FALSE))
print(
  ggplot() +
    annotate("text", x = 0, y = 1, hjust = 0, vjust = 1,
             family = "mono", size = 3,
             label = paste(top_str, collapse = "\n")) +
    xlim(0, 1) + ylim(0, 1) + theme_void() +
    labs(title = paste("Top", top_n, "highest-scoring 6k controls (insT-LASSO)"))
)

dev.off()
message("Wrote ", OUT_PDF)

# --------------------------------------------------------------------------
# ctrl_prob_hist: same exclusions as figure_s4.R histogram (dedup + biallelic + HD)
# This is the n that figure_s4.R labels on its y-axis; figure1 panel G
# should use this field so both figures report identical n.
# --------------------------------------------------------------------------
biallelic_in_6k <- if (exists("bi_keys") && length(bi_keys) > 0L) {
  intersect(bi_keys, names(ctrl_prob))
} else {
  character(0)
}
hd_in_6k <- if (length(hd_platekey) > 0L) {
  intersect(hd_platekey, names(ctrl_prob))
} else {
  character(0)
}
hist_excl  <- unique(c(biallelic_in_6k, hd_in_6k))
ctrl_prob_hist <- ctrl_prob[!names(ctrl_prob) %in% hist_excl]
message(sprintf("ctrl_prob_hist (dedup + biallelic + HD excluded): %d samples",
                length(ctrl_prob_hist)))

# --------------------------------------------------------------------------
# Summary RDS
# --------------------------------------------------------------------------
summary_out <- list(
  n_6k                    = nrow(ctrl_scores),
  n_insT_cases            = length(dev_insT_scores),
  n_dev_ctrl              = sum(dev_ctrl_mask),
  K_selected              = K_sel,
  glass_dev               = glass_dev,
  glass_6k                = glass_6k,
  sd_dev                  = sd_dev,
  sd_6k                   = sd_6k,
  sd_ratio                = sd_6k / sd_dev,
  mu_insT                 = mu_insT,
  mu_dev                  = mu_dev,
  mu_6k                   = mu_6k,
  min_insT_score          = min_insT_score,
  n_fp_above_min_insT     = n_fp,
  per_event_coverage      = setNames(cov_df$coverage,     as.character(cov_df$event_id)),
  per_event_median_depth  = setNames(cov_df$median_depth, as.character(cov_df$event_id)),
  parity_max_abs_diff     = parity_diff,
  paired_summary          = if (!is.null(paired)) {
    list(n = nrow(paired), r = pair_r, mad = pair_mad,
         bias = pair_bias, loa_lo = loa_lo, loa_hi = loa_hi)
  } else NULL,
  paired_scores           = paired,
  top_ctrl                = top_ctrl,
  highlight_summary       = if (length(highlight_summary) > 0L)
    do.call(rbind, highlight_summary) else NULL,
  ctrl_prob               = ctrl_prob,      # LASSO P(case): dedup only
  ctrl_prob_hist          = ctrl_prob_hist  # dedup + biallelic + HD excluded (matches figure_s4.R n)
)
saveRDS(summary_out, OUT_RDS)
message("Wrote ", OUT_RDS, "\n\nDone.")
