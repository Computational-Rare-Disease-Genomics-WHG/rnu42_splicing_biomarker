#!/usr/bin/env Rscript
#
# null_lasso_within_batch.R  (Desktop, pixi)
# Robustness experiment: within a single rMATS batch of clean unmatched
# whole-blood controls, draw 11 pseudo-cases + 33 pseudo-controls, then
# after an HPC bsub'd rMATS --task stat run, fit the insT-style LASSO on
# the resulting FDR<0.05 A5SS pool and report training + holdout AUC.
#
# Two-phase run because pixi only works on Desktop and singularity only on
# HPC.  Phase is detected from what already exists in OUT_DIR:
#
#   Phase 1 (draw): draws samples + writes 0-based indices to
#                   OUT_DIR/drawn_samples.csv, prints the bsub command
#                   for the HPC rMATS step, then exits.
#
#   [User runs the printed bsub on HPC; it writes OUT_DIR/rmats_stat_od/]
#
#   Phase 2 (score): reads OUT_DIR/rmats_stat_od/A5SS.MATS.JC.txt,
#                    fits LASSO, writes null_lasso_result.rds + summary.txt.
#
# Env vars:
#   BATCH_ID   default 1     (integer batch, run check_batch_clean_availability.R first)
#   DRAW_SEED  default 1     (RNG seed; permutations use different seeds)
#   OUT_TAG    default sprintf("batch%d_seed%d", BATCH_ID, DRAW_SEED)
#
# LASSO mirrors loo_evaluate_insT_lasso.R:
#   - candidate space: A5SS FDR<0.05 pseudo-comparison, min(IJC+SJC)>=10
#     within the 44 drawn samples, MAX_PER_GENE=1 by G_AUC
#   - feature = direction * (bblogit(ijc+0.5, sjc) - mu_ctrl_pseudo)
#   - cv.glmnet(alpha=1, binomial, lower.limits=0, foldid=case-LOO,
#               standardize=TRUE, type.measure="deviance", keep=TRUE)
#   - lambda.1se primary, lambda.min fallback if 1se selects 0 sites
#   - training AUC = LOO linpred (fit.preval[, idx.1se])
#   - holdout  AUC = final-fit linpred on 3 held-out pseudo-cases vs the
#     33 pseudo-controls (mirrors figure2.R panel B negative set)

suppressPackageStartupMessages({
  library(glmnet)
})

# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
WORKSPACE     <- "/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
RMATS_TMP     <- file.path(WORKSPACE, "rmats_tmp")
RMATS_OD_ROOT <- file.path(WORKSPACE, "outputs/rmats_controls")

DEV_CASE_BAMS <- file.path(WORKSPACE, "data/rmats/inputs/renu_11_cases_paths.txt")
DEV_CTRL_BAMS <- file.path(WORKSPACE, "data/rmats/inputs/renu_11_cases_control_paths.txt")
BIALLELIC_FILE<- file.path(WORKSPACE, "data/biallelic_rnu42_platekeys.txt")
HD_CASE_FILE  <- file.path(WORKSPACE, "data/renu_highdepth_case_path.txt")

BATCH_ID  <- as.integer(Sys.getenv("BATCH_ID",  "1"))
DRAW_SEED <- as.integer(Sys.getenv("DRAW_SEED", "1"))
OUT_TAG   <- Sys.getenv("OUT_TAG",
                        sprintf("batch%d_seed%d", BATCH_ID, DRAW_SEED))

OUT_DIR   <- file.path(WORKSPACE, "outputs/null_within_batch", OUT_TAG)
STAT_OD   <- file.path(OUT_DIR, "rmats_stat_od")
DRAWN_CSV <- file.path(OUT_DIR, "drawn_samples.csv")
A5_OUT    <- file.path(STAT_OD, "A5SS.MATS.JC.txt")

# LASSO / candidate-space params
N_CASES        <- 11L
N_CTRLS        <- 33L
N_HOLDOUT      <- 3L
FDR_THRESH     <- 0.05
MIN_DEPTH      <- 10L
MAX_PER_GENE   <- 1L
ALPHA_LASSO    <- 1
LAMBDA_RULE    <- "lambda.1se"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("=== null_lasso_within_batch: BATCH_ID=%d  DRAW_SEED=%d ===\n",
            BATCH_ID, DRAW_SEED))
cat(sprintf("OUT_DIR: %s\n", OUT_DIR))

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
read_bam_list <- function(path) {
  if (!file.exists(path)) return(character(0))
  raw  <- paste(readLines(path), collapse = ",")
  bams <- trimws(unlist(strsplit(raw, ",")))
  bams <- bams[nzchar(bams)]
  sub("\\.bam$", "", basename(bams))
}
platekey_from_path <- function(path) sub("\\.bam$", "", basename(path))

bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}
bblogit_at_events <- function(ijc_mat, sjc_mat) {
  out <- bb_logit(ijc_mat, sjc_mat)
  out[(ijc_mat + sjc_mat) == 0] <- NA_real_
  out
}
parse_counts <- function(col, n) {
  t(sapply(col, function(x) {
    vals <- suppressWarnings(as.integer(strsplit(x, ",")[[1]]))
    if (length(vals) != n) rep(NA_integer_, n) else vals
  }))
}
build_coord_nostrand_a5ss <- function(df) {
  # KNOWN BUG: A5SS is chr + strand + 6 coord cols = columns 4:11 (see
  # combine_rmats_batches_cand_pool.R / export_table_s4_ipsc.R /
  # figure_s7.R for the canonical slice). This function slices 4:9,
  # truncating the last two coord cols, so the tag it emits is NOT the
  # canonical coord_nostrand. It is only written into the per-perm
  # `candidate_pool` metadata block (not used for any join or LASSO
  # feature), so the null-perm results are unaffected. Left as-is to
  # avoid rerunning the 1020-perm array; fix to 4:11 before any future
  # rerun.
  coord_key <- do.call(paste, c(df[, 4:9], sep = ":"))
  sub("^([^:]+):[^:]+:", "\\1:", coord_key)
}
auc_two_sample <- function(case_v, ctrl_v) {
  case_v <- case_v[!is.na(case_v)]; ctrl_v <- ctrl_v[!is.na(ctrl_v)]
  if (!length(case_v) || !length(ctrl_v)) return(NA_real_)
  r <- rank(c(case_v, ctrl_v))
  R1 <- sum(r[seq_along(case_v)])
  U  <- R1 - length(case_v) * (length(case_v) + 1) / 2
  U / (length(case_v) * length(ctrl_v))
}
g_auc_dir <- function(auc) {
  ifelse(is.na(auc), NA_real_, 2 * abs(auc - 0.5))
}

# --------------------------------------------------------------------------
# PHASE 1: draw samples if not done yet
# --------------------------------------------------------------------------
do_draw <- function() {
  bam_file <- file.path(RMATS_TMP, paste0("batch_", BATCH_ID),
                        "merged", "bam_list_for_post.txt")
  stopifnot("bam_list_for_post.txt missing" = file.exists(bam_file))
  batch_bams <- unlist(strsplit(readLines(bam_file), ","))
  batch_bams <- trimws(batch_bams[nzchar(trimws(batch_bams))])
  batch_platekeys <- platekey_from_path(batch_bams)
  n_batch <- length(batch_platekeys)
  cat(sprintf("\nBatch %d: %d samples in original order\n", BATCH_ID, n_batch))

  dev_case_keys  <- read_bam_list(DEV_CASE_BAMS)
  dev_ctrl_keys  <- read_bam_list(DEV_CTRL_BAMS)
  biallelic_keys <- if (file.exists(BIALLELIC_FILE))
    trimws(readLines(BIALLELIC_FILE)) else character(0)
  biallelic_keys <- biallelic_keys[nzchar(biallelic_keys)]
  hd_keys        <- read_bam_list(HD_CASE_FILE)
  excl_keys      <- unique(c(dev_case_keys, dev_ctrl_keys,
                             biallelic_keys, hd_keys))

  clean_mask       <- !batch_platekeys %in% excl_keys
  clean_platekeys  <- batch_platekeys[clean_mask]
  clean_pos_1based <- which(clean_mask)
  cat(sprintf("  clean unmatched in batch: %d/%d\n",
              length(clean_platekeys), n_batch))
  stopifnot(length(clean_platekeys) >= (N_CASES + N_CTRLS))

  set.seed(DRAW_SEED)
  drawn_local_idx <- sample.int(length(clean_platekeys),
                                N_CASES + N_CTRLS, replace = FALSE)
  case_local_idx  <- drawn_local_idx[seq_len(N_CASES)]
  ctrl_local_idx  <- drawn_local_idx[(N_CASES + 1L):(N_CASES + N_CTRLS)]

  drawn_df <- data.frame(
    role       = c(rep("case", N_CASES), rep("control", N_CTRLS)),
    platekey   = c(clean_platekeys[case_local_idx],
                   clean_platekeys[ctrl_local_idx]),
    batch_idx0 = c(clean_pos_1based[case_local_idx] - 1L,
                   clean_pos_1based[ctrl_local_idx] - 1L),
    stringsAsFactors = FALSE
  )
  write.csv(drawn_df, DRAWN_CSV, row.names = FALSE)
  cat(sprintf("  drew %d pseudo-cases + %d pseudo-controls (seed=%d)\n",
              N_CASES, N_CTRLS, DRAW_SEED))
  cat(sprintf("  wrote: %s\n", DRAWN_CSV))

  # Print the HPC command the user needs to run next
  script_sh <- file.path(WORKSPACE,
                         "scripts/gel/null_lasso_rmats_stat.sh")
  log_dir   <- file.path(WORKSPACE, "hpc_logs")
  cat("\n=========================================================\n")
  cat("PHASE 1 complete. Now run the rMATS stat step on HPC:\n\n")
  cat(sprintf(
    "  bsub -q short -P re_gecip_machine_learning -J null_%s \\\n",
    OUT_TAG))
  cat(sprintf(
    "    -o %s/null_%%J.out -e %s/null_%%J.err \\\n",
    log_dir, log_dir))
  cat(sprintf(
    "    -cwd %s -L /bin/bash \\\n", WORKSPACE))
  cat(sprintf(
    "    bash %s %d %d %s\n\n",
    script_sh, BATCH_ID, DRAW_SEED, OUT_TAG))
  cat("Then re-run this R script (same env vars) to score.\n")
  cat("=========================================================\n")
}

# --------------------------------------------------------------------------
# PHASE 2: score after rMATS stat has produced A5SS.MATS.JC.txt
# --------------------------------------------------------------------------
do_score <- function() {
  stopifnot("drawn_samples.csv missing"  = file.exists(DRAWN_CSV),
            "A5SS.MATS.JC.txt missing"   = file.exists(A5_OUT))
  drawn_df <- read.csv(DRAWN_CSV, stringsAsFactors = FALSE)
  stopifnot(nrow(drawn_df) == N_CASES + N_CTRLS)
  case_platekeys <- drawn_df$platekey[drawn_df$role == "case"]
  ctrl_platekeys <- drawn_df$platekey[drawn_df$role == "control"]

  a5 <- read.delim(A5_OUT, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("\nA5SS rows returned: %d\n", nrow(a5)))

  ijc_case <- parse_counts(a5$IJC_SAMPLE_1, N_CASES)
  sjc_case <- parse_counts(a5$SJC_SAMPLE_1, N_CASES)
  ijc_ctrl <- parse_counts(a5$IJC_SAMPLE_2, N_CTRLS)
  sjc_ctrl <- parse_counts(a5$SJC_SAMPLE_2, N_CTRLS)
  depth_case <- ijc_case + sjc_case
  depth_ctrl <- ijc_ctrl + sjc_ctrl
  min_depth_44 <- pmin(apply(depth_case, 1L, min, na.rm = TRUE),
                       apply(depth_ctrl, 1L, min, na.rm = TRUE))
  fdr <- suppressWarnings(as.numeric(a5$FDR))

  bblogit_case <- bblogit_at_events(ijc_case, sjc_case)
  bblogit_ctrl <- bblogit_at_events(ijc_ctrl, sjc_ctrl)

  auc_pseudo <- vapply(seq_len(nrow(a5)), function(j) {
    auc_two_sample(bblogit_case[j, ], bblogit_ctrl[j, ])
  }, numeric(1L))
  g_auc_pseudo <- g_auc_dir(auc_pseudo)
  direction    <- ifelse(is.na(auc_pseudo), NA_real_,
                         ifelse(auc_pseudo >= 0.5, 1, -1))

  cand_mask <- !is.na(fdr) & fdr < FDR_THRESH &
               !is.na(min_depth_44) & min_depth_44 >= MIN_DEPTH &
               !is.na(direction)

  cat("Candidate filter counts:\n")
  cat(sprintf("  FDR < %.2f                 : %d\n", FDR_THRESH,
              sum(!is.na(fdr) & fdr < FDR_THRESH)))
  cat(sprintf("  + min(IJC+SJC) >= %d       : %d\n", MIN_DEPTH,
              sum(cand_mask)))

  gene_sym <- gsub('"', '', a5$geneSymbol, fixed = TRUE)
  cand_idx  <- which(cand_mask)
  ord       <- cand_idx[order(-g_auc_pseudo[cand_idx])]
  seen_gene <- character(0)
  keep_idx  <- integer(0)
  for (i in ord) {
    g <- gene_sym[i]
    if (is.na(g) || nchar(g) == 0L) next
    if (sum(seen_gene == g) >= MAX_PER_GENE) next
    seen_gene <- c(seen_gene, g)
    keep_idx  <- c(keep_idx, i)
  }
  cat(sprintf("  + MAX_PER_GENE=%d (by G_AUC): %d\n",
              MAX_PER_GENE, length(keep_idx)))
  stopifnot("No candidate A5SS events survive filtering" = length(keep_idx) > 0L)

  cand <- data.frame(
    row_idx        = keep_idx,
    event_id       = a5$ID[keep_idx],
    GeneSymbol     = gene_sym[keep_idx],
    coord_nostrand = build_coord_nostrand_a5ss(a5)[keep_idx],
    FDR            = fdr[keep_idx],
    direction      = direction[keep_idx],
    g_auc          = g_auc_pseudo[keep_idx],
    min_depth_44   = min_depth_44[keep_idx],
    stringsAsFactors = FALSE
  )

  set.seed(DRAW_SEED)
  ho_case_idx    <- sort(sample.int(N_CASES, N_HOLDOUT, replace = FALSE))
  train_case_idx <- setdiff(seq_len(N_CASES), ho_case_idx)
  n_train_case   <- length(train_case_idx)
  cat(sprintf("\nTrain/holdout: %d train cases, %d holdout cases\n",
              n_train_case, N_HOLDOUT))

  bblogit_case_cand <- bblogit_case[keep_idx, , drop = FALSE]
  bblogit_ctrl_cand <- bblogit_ctrl[keep_idx, , drop = FALSE]
  mu_ctrl <- rowMeans(bblogit_ctrl_cand, na.rm = TRUE)

  orient_case <- sweep(sweep(bblogit_case_cand, 1L, mu_ctrl, "-"),
                       1L, cand$direction, "*")
  orient_ctrl <- sweep(sweep(bblogit_ctrl_cand, 1L, mu_ctrl, "-"),
                       1L, cand$direction, "*")

  X_train <- cbind(orient_case[, train_case_idx, drop = FALSE], orient_ctrl)
  X_train <- t(X_train)
  y_train <- c(rep(1L, n_train_case), rep(0L, N_CTRLS))
  rownames(X_train) <- c(case_platekeys[train_case_idx], ctrl_platekeys)

  if (anyNA(X_train)) {
    warning(sum(is.na(X_train)), " NA cell(s) in X_train after min-depth ",
            "filter -- imputing with column mean")
    col_mean <- colMeans(X_train, na.rm = TRUE)
    for (j in which(colSums(is.na(X_train)) > 0L)) {
      X_train[is.na(X_train[, j]), j] <- col_mean[j]
    }
  }

  case_row <- which(y_train == 1L)
  ctrl_row <- which(y_train == 0L)
  foldid   <- integer(length(y_train))
  foldid[case_row] <- seq_along(case_row)
  foldid[ctrl_row] <- ((seq_along(ctrl_row) - 1L) %% n_train_case) + 1L

  cat("\nFitting cv.glmnet ...\n")
  set.seed(DRAW_SEED)
  cv_fit <- cv.glmnet(
    x            = X_train,
    y            = y_train,
    family       = "binomial",
    alpha        = ALPHA_LASSO,
    foldid       = foldid,
    standardize  = TRUE,
    lower.limits = 0,
    type.measure = "deviance",
    keep         = TRUE
  )

  idx_1se <- which(cv_fit$lambda == cv_fit$lambda.1se)
  idx_min <- which(cv_fit$lambda == cv_fit$lambda.min)

  extract_fit <- function(cv_fit, s) {
    co <- coef(cv_fit, s = s)
    cv <- as.numeric(co); names(cv) <- rownames(co)
    intercept <- cv[1]; beta <- cv[-1]; nz <- which(beta != 0)
    list(intercept = intercept, beta = beta, n_selected = length(nz),
         event_ids = names(beta)[nz], coefs = beta[nz],
         lambda = if (s == "lambda.1se") cv_fit$lambda.1se else cv_fit$lambda.min)
  }
  fit_at_1se <- extract_fit(cv_fit, "lambda.1se")
  fit_at_min <- extract_fit(cv_fit, "lambda.min")

  chosen_label <- LAMBDA_RULE
  if (fit_at_1se$n_selected == 0L) {
    warning("lambda.1se selected 0 sites; falling back to lambda.min")
    chosen_label <- "lambda.min"
  }
  chosen_fit <- if (chosen_label == "lambda.1se") fit_at_1se else fit_at_min
  idx_chosen <- if (chosen_label == "lambda.1se") idx_1se else idx_min

  cat(sprintf("  lambda.1se=%.5g  n_selected=%d\n",
              cv_fit$lambda.1se, fit_at_1se$n_selected))
  cat(sprintf("  lambda.min=%.5g  n_selected=%d\n",
              cv_fit$lambda.min, fit_at_min$n_selected))
  cat(sprintf("  using: %s\n", chosen_label))

  preval_lp <- cv_fit$fit.preval[, idx_chosen]
  train_auc <- auc_two_sample(preval_lp[y_train == 1L],
                              preval_lp[y_train == 0L])
  cat(sprintf("\nTraining AUC (LOO linpred, %s): %.3f  (%d cases vs %d ctrls)\n",
              chosen_label, train_auc, n_train_case, N_CTRLS))

  X_holdout <- t(orient_case[, ho_case_idx, drop = FALSE])
  rownames(X_holdout) <- case_platekeys[ho_case_idx]
  if (anyNA(X_holdout)) {
    col_mean <- colMeans(X_train, na.rm = TRUE)
    for (j in which(colSums(is.na(X_holdout)) > 0L)) {
      X_holdout[is.na(X_holdout[, j]), j] <- col_mean[j]
    }
  }
  lp_holdout_case <- as.numeric(X_holdout %*% chosen_fit$beta +
                                chosen_fit$intercept)
  X_ctrl_full     <- t(orient_ctrl)
  lp_ctrl_full    <- as.numeric(X_ctrl_full %*% chosen_fit$beta +
                                chosen_fit$intercept)
  holdout_auc <- auc_two_sample(lp_holdout_case, lp_ctrl_full)
  cat(sprintf("Holdout  AUC (final-fit linpred, %s): %.3f  (%d cases vs %d ctrls)\n",
              chosen_label, holdout_auc, N_HOLDOUT, N_CTRLS))

  result <- list(
    config = list(
      batch_id = BATCH_ID, draw_seed = DRAW_SEED,
      n_cases = N_CASES, n_ctrls = N_CTRLS, n_holdout = N_HOLDOUT,
      fdr_thresh = FDR_THRESH, min_depth = MIN_DEPTH,
      max_per_gene = MAX_PER_GENE, alpha_lasso = ALPHA_LASSO,
      lambda_rule = LAMBDA_RULE, chosen_label = chosen_label
    ),
    drawn_samples  = drawn_df,
    train_case_idx = train_case_idx,
    ho_case_idx    = ho_case_idx,
    candidate_pool = cand,
    cv_fit         = cv_fit,
    fit_at_1se     = fit_at_1se,
    fit_at_min     = fit_at_min,
    chosen_fit     = chosen_fit,
    mu_ctrl        = mu_ctrl,
    preval_lp      = preval_lp,
    lp_holdout_case= setNames(lp_holdout_case, rownames(X_holdout)),
    lp_ctrl_full   = setNames(lp_ctrl_full,    ctrl_platekeys),
    auc            = c(train = train_auc, holdout = holdout_auc)
  )
  saveRDS(result, file.path(OUT_DIR, "null_lasso_result.rds"))

  writeLines(c(
    sprintf("batch_id       : %d", BATCH_ID),
    sprintf("draw_seed      : %d", DRAW_SEED),
    sprintf("candidate pool : %d events (FDR<%.2f, min_depth>=%d, MAX_PER_GENE=%d)",
            nrow(cand), FDR_THRESH, MIN_DEPTH, MAX_PER_GENE),
    sprintf("lambda.1se n_selected : %d", fit_at_1se$n_selected),
    sprintf("lambda.min n_selected : %d", fit_at_min$n_selected),
    sprintf("chosen                : %s", chosen_label),
    sprintf("training AUC (LOO)    : %.3f  (%d cases vs %d ctrls)",
            train_auc, n_train_case, N_CTRLS),
    sprintf("holdout  AUC (fit)    : %.3f  (%d cases vs %d ctrls)",
            holdout_auc, N_HOLDOUT, N_CTRLS)
  ), file.path(OUT_DIR, "summary.txt"))

  cat(sprintf("\n=== Done ===\nSaved to: %s\n", OUT_DIR))
}

# --------------------------------------------------------------------------
# Dispatch: phase 2 if rMATS output is present, phase 1 otherwise
# --------------------------------------------------------------------------
if (file.exists(A5_OUT)) {
  cat("Phase 2: A5SS.MATS.JC.txt present -- scoring.\n")
  do_score()
} else {
  cat("Phase 1: no rMATS stat output yet -- drawing samples.\n")
  do_draw()
}
