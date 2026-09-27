#!/usr/bin/env Rscript
#
# null_lasso_permutations.R  (Desktop, pixi)
# Stratified null LASSO driver: for each qualifying rMATS batch (n_clean >= 44),
# draw PER_BATCH random 11-vs-33 pseudo-comparisons, run rMATS stat per draw,
# fit the insT-style LASSO on each, aggregate training + holdout AUCs into a
# null distribution.
#
#   Phase 1 (draw_all): identify qualifying batches, write a permutation
#                       manifest (perm_id -> batch_id + within-batch seed),
#                       write per-perm drawn_samples.csv, print bsub array cmd.
#
#   [User submits the printed LSF array; each element reads the manifest for
#    its LSB_JOBINDEX and runs the rMATS stat pipeline for that draw.]
#
#   Phase 3 (score_all): for each perm load rmats_stat_od + score with
#                        the LASSO pipeline, aggregate into
#                        outputs/null_within_batch/<AGG_TAG>/null_auc_distribution.csv
#
# Env vars:
#   PER_BATCH   default 17   (draws per qualifying batch; total = PER_BATCH * n_qual)
#   N_CASES     default 11   (pseudo-cases per draw at rMATS step)
#   N_CTRLS     default 33   (pseudo-controls per draw at rMATS step)
#   N_HOLDOUT   default 3    (of the N_CASES, held out from LASSO training)
#   SAMPLE_TAG  default sprintf("n%dv%d", N_CASES, N_CTRLS)   subdir name
#   AGG_TAG     default sprintf("stratified_p%d", PER_BATCH)
#   N_CORES     default detectCores()-1
#
# All outputs go under outputs/null_within_batch/<SAMPLE_TAG>/ so that
# runs with different case/ctrl sizes never collide. Existing 11v33 output
# lives at the flat old layout and is untouched by new SAMPLE_TAG runs.
#
# Per-perm outputs:
#   outputs/null_within_batch/<SAMPLE_TAG>/batch<BATCH>_seed<SEED>/
#     drawn_samples.csv, rmats_stat_od/, null_lasso_result.rds
# Manifest + aggregate:
#   outputs/null_within_batch/<SAMPLE_TAG>/<AGG_TAG>/permutation_manifest.csv
#   outputs/null_within_batch/<SAMPLE_TAG>/<AGG_TAG>/null_auc_distribution.csv

suppressPackageStartupMessages({
  library(glmnet)
  library(data.table)   # fread is ~50x faster than read.delim on 95k rows
  library(parallel)     # mclapply for scoring 1020 perms in parallel
})

N_CORES <- as.integer(Sys.getenv(
  "N_CORES", as.character(max(1L, parallel::detectCores() - 1L))))

# --------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------
WORKSPACE     <- "/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
RMATS_TMP     <- file.path(WORKSPACE, "rmats_tmp")

DEV_CASE_BAMS <- file.path(WORKSPACE, "data/rmats/inputs/renu_11_cases_paths.txt")
DEV_CTRL_BAMS <- file.path(WORKSPACE, "data/rmats/inputs/renu_11_cases_control_paths.txt")
BIALLELIC_FILE<- file.path(WORKSPACE, "data/biallelic_rnu42_platekeys.txt")
HD_CASE_FILE  <- file.path(WORKSPACE, "data/renu_highdepth_case_path.txt")

BATCH_RANGE  <- c(1L, 61L)   # QC-pass batches
FDR_THRESH   <- 0.05
MIN_DEPTH    <- 10L
MAX_PER_GENE <- 1L
ALPHA_LASSO  <- 1
LAMBDA_RULE  <- "lambda.1se"

N_CASES   <- as.integer(Sys.getenv("N_CASES",   "11"))
N_CTRLS   <- as.integer(Sys.getenv("N_CTRLS",   "33"))
N_HOLDOUT <- as.integer(Sys.getenv("N_HOLDOUT", "3"))
stopifnot(N_HOLDOUT < N_CASES)

PER_BATCH  <- as.integer(Sys.getenv("PER_BATCH", "17"))
SAMPLE_TAG <- Sys.getenv("SAMPLE_TAG",
                          sprintf("n%dv%d", N_CASES, N_CTRLS))
AGG_TAG    <- Sys.getenv("AGG_TAG", sprintf("stratified_p%d", PER_BATCH))
OUT_ROOT   <- file.path(WORKSPACE, "outputs/null_within_batch", SAMPLE_TAG)
AGG_DIR    <- file.path(OUT_ROOT, AGG_TAG)
MANIFEST   <- file.path(AGG_DIR, "permutation_manifest.csv")
AGG_CSV    <- file.path(AGG_DIR, "null_auc_distribution.csv")

dir.create(AGG_DIR, recursive = TRUE, showWarnings = FALSE)

# out_tag stored in the manifest is a WORKSPACE-relative path fragment
# beneath outputs/null_within_batch/ so null_lasso_rmats_stat.sh (which
# composes OUT_DIR=<WORKSPACE>/outputs/null_within_batch/<OUT_TAG>)
# writes to the SAMPLE_TAG subdir automatically. For SAMPLE_TAG="" this
# collapses to the flat legacy layout.
tag_join <- function(...) {
  parts <- Filter(nzchar, c(...))
  paste(parts, collapse = "/")
}
perm_out_dir <- function(row) {
  file.path(WORKSPACE, "outputs/null_within_batch", row$out_tag)
}
perm_drawn_csv <- function(row) file.path(perm_out_dir(row), "drawn_samples.csv")
perm_a5_out    <- function(row) file.path(perm_out_dir(row),
                                          "rmats_stat_od", "A5SS.MATS.JC.txt")

cat(sprintf("=== null_lasso_permutations: PER_BATCH=%d ===\n", PER_BATCH))
cat(sprintf("AGG_DIR: %s\n", AGG_DIR))

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
  p <- (ijc + 0.5) / (ijc + sjc + 1); log(p / (1 - p))
}
bblogit_at_events <- function(ijc_mat, sjc_mat) {
  out <- bb_logit(ijc_mat, sjc_mat)
  out[(ijc_mat + sjc_mat) == 0] <- NA_real_
  out
}
parse_counts <- function(col, n) {
  # Vectorised: single strsplit call, then reshape. ~50x faster than the
  # per-row sapply variant on 95k rMATS rows. Assumes every row has exactly
  # n comma-separated values (rMATS output is uniform).
  vals <- suppressWarnings(as.integer(unlist(
    strsplit(col, ",", fixed = TRUE))))
  matrix(vals, ncol = n, byrow = TRUE)
}
auc_two_sample <- function(case_v, ctrl_v) {
  case_v <- case_v[!is.na(case_v)]; ctrl_v <- ctrl_v[!is.na(ctrl_v)]
  if (!length(case_v) || !length(ctrl_v)) return(NA_real_)
  r <- rank(c(case_v, ctrl_v))
  R1 <- sum(r[seq_along(case_v)])
  U  <- R1 - length(case_v) * (length(case_v) + 1) / 2
  U / (length(case_v) * length(ctrl_v))
}
g_auc_dir <- function(auc) ifelse(is.na(auc), NA_real_, 2 * abs(auc - 0.5))

# --------------------------------------------------------------------------
# Preload exclusion set + per-batch clean samples (shared across all perms)
# --------------------------------------------------------------------------
dev_case_keys  <- read_bam_list(DEV_CASE_BAMS)
dev_ctrl_keys  <- read_bam_list(DEV_CTRL_BAMS)
biallelic_keys <- if (file.exists(BIALLELIC_FILE))
  trimws(readLines(BIALLELIC_FILE)) else character(0)
biallelic_keys <- biallelic_keys[nzchar(biallelic_keys)]
hd_keys        <- read_bam_list(HD_CASE_FILE)
excl_keys      <- unique(c(dev_case_keys, dev_ctrl_keys, biallelic_keys, hd_keys))

load_batch_clean <- function(batch_id) {
  bam_file <- file.path(RMATS_TMP, paste0("batch_", batch_id),
                        "merged", "bam_list_for_post.txt")
  if (!file.exists(bam_file)) return(NULL)
  bams <- unlist(strsplit(readLines(bam_file), ","))
  bams <- trimws(bams[nzchar(trimws(bams))])
  platekeys <- platekey_from_path(bams)
  clean_mask <- !platekeys %in% excl_keys
  list(platekeys = platekeys[clean_mask],
       pos1based = which(clean_mask),
       n_total   = length(platekeys),
       n_clean   = sum(clean_mask))
}

cat("Scanning batches ", BATCH_RANGE[1L], "-", BATCH_RANGE[2L],
    " for qualifying (n_clean >= ", (N_CASES + N_CTRLS), ") ...\n", sep = "")
batch_ids <- BATCH_RANGE[1L]:BATCH_RANGE[2L]
batch_clean <- setNames(vector("list", length(batch_ids)),
                        as.character(batch_ids))
for (b in batch_ids) batch_clean[[as.character(b)]] <- load_batch_clean(b)

qualifying <- batch_ids[vapply(batch_ids, function(b) {
  bc <- batch_clean[[as.character(b)]]
  !is.null(bc) && bc$n_clean >= (N_CASES + N_CTRLS)
}, logical(1L))]
cat(sprintf("Qualifying batches: %d / %d\n",
            length(qualifying), length(batch_ids)))
stopifnot("No qualifying batches" = length(qualifying) > 0L)

# --------------------------------------------------------------------------
# Build (or reuse) permutation manifest
# --------------------------------------------------------------------------
if (file.exists(MANIFEST)) {
  manifest <- read.csv(MANIFEST, stringsAsFactors = FALSE)
  cat(sprintf("Reusing existing manifest: %d perms\n", nrow(manifest)))
  stopifnot(nrow(manifest) == length(qualifying) * PER_BATCH)
} else {
  manifest <- do.call(rbind, lapply(qualifying, function(b) {
    data.frame(
      perm_id           = NA_integer_,     # filled below
      batch_id          = b,
      within_batch_seed = seq_len(PER_BATCH),
      out_tag           = vapply(seq_len(PER_BATCH), function(s)
                            tag_join(SAMPLE_TAG,
                                     sprintf("batch%d_seed%d", b, s)),
                            character(1L)),
      stringsAsFactors  = FALSE
    )
  }))
  manifest$perm_id <- seq_len(nrow(manifest))
  write.csv(manifest, MANIFEST, row.names = FALSE)
  cat(sprintf("Wrote manifest: %d perms (%d batches x %d/batch)\n",
              nrow(manifest), length(qualifying), PER_BATCH))
}
N_PERMS <- nrow(manifest)

# --------------------------------------------------------------------------
# do_draw / do_score  (one perm each)
# --------------------------------------------------------------------------
do_draw_one <- function(row) {
  bc  <- batch_clean[[as.character(row$batch_id)]]
  out <- perm_out_dir(row)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  # Compose seed to be deterministic per (batch, within-batch seed).
  set.seed(row$batch_id * 100000L + row$within_batch_seed)
  drawn_local_idx <- sample.int(length(bc$platekeys),
                                N_CASES + N_CTRLS, replace = FALSE)
  case_local_idx  <- drawn_local_idx[seq_len(N_CASES)]
  ctrl_local_idx  <- drawn_local_idx[(N_CASES + 1L):(N_CASES + N_CTRLS)]
  drawn_df <- data.frame(
    role       = c(rep("case", N_CASES), rep("control", N_CTRLS)),
    platekey   = c(bc$platekeys[case_local_idx],
                   bc$platekeys[ctrl_local_idx]),
    batch_idx0 = c(bc$pos1based[case_local_idx] - 1L,
                   bc$pos1based[ctrl_local_idx] - 1L),
    stringsAsFactors = FALSE
  )
  write.csv(drawn_df, perm_drawn_csv(row), row.names = FALSE)
}

do_score_one <- function(row) {
  out_dir <- perm_out_dir(row)
  drawn_df <- read.csv(perm_drawn_csv(row), stringsAsFactors = FALSE)
  case_platekeys <- drawn_df$platekey[drawn_df$role == "case"]
  ctrl_platekeys <- drawn_df$platekey[drawn_df$role == "control"]

  # fread + column select: skip 15 unused columns, parse ~10x faster.
  a5 <- fread(perm_a5_out(row),
              select = c("ID", "geneSymbol",
                         "IJC_SAMPLE_1", "SJC_SAMPLE_1",
                         "IJC_SAMPLE_2", "SJC_SAMPLE_2", "FDR"),
              data.table = FALSE, showProgress = FALSE)
  ijc_case <- parse_counts(a5$IJC_SAMPLE_1, N_CASES)
  sjc_case <- parse_counts(a5$SJC_SAMPLE_1, N_CASES)
  ijc_ctrl <- parse_counts(a5$IJC_SAMPLE_2, N_CTRLS)
  sjc_ctrl <- parse_counts(a5$SJC_SAMPLE_2, N_CTRLS)
  depth_case <- ijc_case + sjc_case
  depth_ctrl <- ijc_ctrl + sjc_ctrl
  # Vectorised row-mins via pmin across columns (replaces apply(,1,min)
  # which is slow at 95k rows).
  row_min <- function(m) {
    m2 <- m; m2[is.na(m2)] <- .Machine$integer.max
    do.call(pmin.int, split(m2, col(m2)))
  }
  min_depth_44 <- pmin.int(row_min(depth_case), row_min(depth_ctrl))
  min_depth_44[min_depth_44 == .Machine$integer.max] <- NA_integer_
  fdr <- suppressWarnings(as.numeric(a5$FDR))
  bblogit_case <- bblogit_at_events(ijc_case, sjc_case)
  bblogit_ctrl <- bblogit_at_events(ijc_ctrl, sjc_ctrl)

  # Direction from mean bblogit difference (fast, vectorised, matches how
  # loo_evaluate_insT_lasso.R computes direction from dlogit means).
  mean_case_v <- rowMeans(bblogit_case, na.rm = TRUE)
  mean_ctrl_v <- rowMeans(bblogit_ctrl, na.rm = TRUE)
  dl          <- mean_case_v - mean_ctrl_v
  dl[!is.finite(dl)] <- NA_real_
  direction   <- ifelse(is.na(dl), NA_real_, sign(dl))
  direction[direction == 0] <- NA_real_

  cand_mask <- !is.na(fdr) & fdr < FDR_THRESH &
               !is.na(min_depth_44) & min_depth_44 >= MIN_DEPTH &
               !is.na(direction)

  # G_AUC only needed for MAX_PER_GENE ranking within the filtered set.
  # Vectorised: compute per-event Mann-Whitney U over the filtered rows
  # in one shot instead of per-event rank() calls.
  g_auc_pseudo <- rep(NA_real_, nrow(a5))
  cand_rows <- which(cand_mask)
  if (length(cand_rows)) {
    combined_c <- cbind(bblogit_case[cand_rows, , drop = FALSE],
                        bblogit_ctrl[cand_rows, , drop = FALSE])
    # Row-wise ranks. Ties get average rank; NAs propagate but are rare
    # here because the min-depth filter already guaranteed finite bblogit.
    ranks_c <- t(apply(combined_c, 1L, rank, na.last = "keep"))
    R1 <- rowSums(ranks_c[, seq_len(N_CASES), drop = FALSE], na.rm = TRUE)
    n_case_ok <- rowSums(!is.na(combined_c[, seq_len(N_CASES),
                                           drop = FALSE]))
    n_ctrl_ok <- rowSums(!is.na(combined_c[, (N_CASES + 1L):(N_CASES + N_CTRLS),
                                           drop = FALSE]))
    U <- R1 - n_case_ok * (n_case_ok + 1) / 2
    auc_c <- U / (n_case_ok * n_ctrl_ok)
    g_auc_pseudo[cand_rows] <- g_auc_dir(auc_c)
  }

  gene_sym <- gsub('"', '', a5$geneSymbol, fixed = TRUE)
  cand_idx <- which(cand_mask)
  ord      <- cand_idx[order(-g_auc_pseudo[cand_idx])]
  seen_gene <- character(0); keep_idx <- integer(0)
  for (i in ord) {
    g <- gene_sym[i]
    if (is.na(g) || nchar(g) == 0L) next
    if (sum(seen_gene == g) >= MAX_PER_GENE) next
    seen_gene <- c(seen_gene, g)
    keep_idx  <- c(keep_idx, i)
  }
  n_fdr_pass <- sum(!is.na(fdr) & fdr < FDR_THRESH)
  n_cand     <- length(keep_idx)
  base_row <- data.frame(
    perm_id      = row$perm_id, batch_id = row$batch_id,
    within_batch_seed = row$within_batch_seed,
    n_fdr_pass   = n_fdr_pass, n_candidates = n_cand,
    lambda_1se = NA_real_, K_1se = 0L,
    lambda_min = NA_real_, K_min = 0L,
    chosen = NA_character_, train_auc = NA_real_, holdout_auc = NA_real_,
    stringsAsFactors = FALSE
  )
  if (n_cand == 0L) return(base_row)

  cand <- data.frame(
    event_id   = a5$ID[keep_idx],
    GeneSymbol = gene_sym[keep_idx],
    direction  = direction[keep_idx],
    FDR        = fdr[keep_idx],
    g_auc      = g_auc_pseudo[keep_idx],
    stringsAsFactors = FALSE
  )

  perm_seed <- row$batch_id * 100000L + row$within_batch_seed
  set.seed(perm_seed)
  ho_case_idx    <- sort(sample.int(N_CASES, N_HOLDOUT, replace = FALSE))
  train_case_idx <- setdiff(seq_len(N_CASES), ho_case_idx)
  n_train_case   <- length(train_case_idx)

  bblogit_case_cand <- bblogit_case[keep_idx, , drop = FALSE]
  bblogit_ctrl_cand <- bblogit_ctrl[keep_idx, , drop = FALSE]
  mu_ctrl <- rowMeans(bblogit_ctrl_cand, na.rm = TRUE)
  orient_case <- sweep(sweep(bblogit_case_cand, 1L, mu_ctrl, "-"),
                       1L, cand$direction, "*")
  orient_ctrl <- sweep(sweep(bblogit_ctrl_cand, 1L, mu_ctrl, "-"),
                       1L, cand$direction, "*")

  X_train <- t(cbind(orient_case[, train_case_idx, drop = FALSE], orient_ctrl))
  y_train <- c(rep(1L, n_train_case), rep(0L, N_CTRLS))
  rownames(X_train) <- c(case_platekeys[train_case_idx], ctrl_platekeys)
  if (anyNA(X_train)) {
    col_mean <- colMeans(X_train, na.rm = TRUE)
    for (j in which(colSums(is.na(X_train)) > 0L)) {
      X_train[is.na(X_train[, j]), j] <- col_mean[j]
    }
  }

  case_row <- which(y_train == 1L); ctrl_row <- which(y_train == 0L)
  foldid   <- integer(length(y_train))
  foldid[case_row] <- seq_along(case_row)
  foldid[ctrl_row] <- ((seq_along(ctrl_row) - 1L) %% n_train_case) + 1L

  set.seed(perm_seed)
  cv_fit <- suppressWarnings(cv.glmnet(
    x = X_train, y = y_train, family = "binomial", alpha = ALPHA_LASSO,
    foldid = foldid, standardize = TRUE, lower.limits = 0,
    type.measure = "deviance", keep = TRUE
  ))
  idx_1se <- which(cv_fit$lambda == cv_fit$lambda.1se)
  idx_min <- which(cv_fit$lambda == cv_fit$lambda.min)

  extract_fit <- function(cv_fit, s) {
    co <- coef(cv_fit, s = s)
    cv <- as.numeric(co); names(cv) <- rownames(co)
    intercept <- cv[1]; beta <- cv[-1]; nz <- which(beta != 0)
    list(intercept = intercept, beta = beta, n_selected = length(nz),
         event_ids = names(beta)[nz], coefs = beta[nz])
  }
  fit_at_1se <- extract_fit(cv_fit, "lambda.1se")
  fit_at_min <- extract_fit(cv_fit, "lambda.min")
  chosen_label <- LAMBDA_RULE
  if (fit_at_1se$n_selected == 0L) chosen_label <- "lambda.min"
  chosen_fit <- if (chosen_label == "lambda.1se") fit_at_1se else fit_at_min
  idx_chosen <- if (chosen_label == "lambda.1se") idx_1se else idx_min

  preval_lp <- cv_fit$fit.preval[, idx_chosen]
  train_auc <- auc_two_sample(preval_lp[y_train == 1L],
                              preval_lp[y_train == 0L])

  X_holdout <- t(orient_case[, ho_case_idx, drop = FALSE])
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

  saveRDS(list(
    config = list(batch_id = row$batch_id,
                  within_batch_seed = row$within_batch_seed,
                  perm_id = row$perm_id, chosen_label = chosen_label),
    candidate_pool  = cand,
    fit_at_1se      = fit_at_1se, fit_at_min = fit_at_min,
    chosen_fit      = chosen_fit, mu_ctrl = mu_ctrl,
    preval_lp       = preval_lp,
    lp_holdout_case = lp_holdout_case,
    lp_ctrl_full    = lp_ctrl_full,
    auc = c(train = train_auc, holdout = holdout_auc)
  ), file.path(out_dir, "null_lasso_result.rds"))

  base_row$lambda_1se <- cv_fit$lambda.1se
  base_row$K_1se      <- fit_at_1se$n_selected
  base_row$lambda_min <- cv_fit$lambda.min
  base_row$K_min      <- fit_at_min$n_selected
  base_row$chosen     <- chosen_label
  base_row$train_auc  <- train_auc
  base_row$holdout_auc<- holdout_auc
  base_row
}

# --------------------------------------------------------------------------
# Dispatch: phase 3 if ALL perms have rMATS output; else phase 1
# --------------------------------------------------------------------------
has_rmats <- vapply(seq_len(N_PERMS), function(i) {
  file.exists(perm_a5_out(manifest[i, ]))
}, logical(1L))

n_present <- sum(has_rmats)
frac_present <- n_present / N_PERMS
if (n_present > 0L) {
  cat(sprintf("Phase 3: %d/%d perms have rMATS output (%.1f%%) -- scoring.\n",
              n_present, N_PERMS, 100 * frac_present))
  if (n_present < N_PERMS) {
    missing_ids <- which(!has_rmats)
    cat(sprintf("  skipping (NA row) for %d missing perm_id(s): %s\n",
                length(missing_ids),
                paste(head(missing_ids, 20L), collapse = ", ")))
  }
  na_row <- function(i) data.frame(
    perm_id = i, batch_id = manifest$batch_id[i],
    within_batch_seed = manifest$within_batch_seed[i],
    n_fdr_pass = NA_integer_, n_candidates = NA_integer_,
    lambda_1se = NA_real_, K_1se = NA_integer_,
    lambda_min = NA_real_, K_min = NA_integer_,
    chosen = NA_character_,
    train_auc = NA_real_, holdout_auc = NA_real_,
    stringsAsFactors = FALSE)
  cat(sprintf("  parallel scoring across %d cores (chunks of 50) ...\n",
              N_CORES))
  t0 <- Sys.time()
  chunk_size <- 50L
  chunks <- split(seq_len(N_PERMS),
                  ceiling(seq_len(N_PERMS) / chunk_size))
  results <- vector("list", N_PERMS)
  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    chunk_res <- mclapply(idx, function(i) {
      if (!has_rmats[i]) return(na_row(i))
      tryCatch(do_score_one(manifest[i, ]),
               error = function(e) {
                 message(sprintf("perm %d failed: %s", i, conditionMessage(e)))
                 na_row(i)
               })
    }, mc.cores = N_CORES, mc.preschedule = TRUE)
    results[idx] <- chunk_res
    elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    done <- max(idx)
    eta_min <- elapsed_min / done * (N_PERMS - done)
    cat(sprintf("  chunk %d/%d done | %d/%d perms | elapsed %.1f min | ETA %.1f min\n",
                ci, length(chunks), done, N_PERMS, elapsed_min, eta_min))
  }
  cat(sprintf("  scored %d perms in %.1f min total\n",
              N_PERMS, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  agg <- do.call(rbind, results)
  write.csv(agg, AGG_CSV, row.names = FALSE)
  cat(sprintf("\nWrote null AUC distribution: %s\n", AGG_CSV))
  cat(sprintf("Rows: %d  (failed: %d)\n",
              nrow(agg), sum(is.na(agg$train_auc))))
  cat("\nTrain AUC:\n"); print(summary(agg$train_auc))
  cat("\nHoldout AUC:\n"); print(summary(agg$holdout_auc))
  cat("\nK selected (chosen lambda):\n")
  k_chosen <- ifelse(agg$chosen == "lambda.1se", agg$K_1se, agg$K_min)
  print(summary(k_chosen))
  n_valid <- sum(!is.na(agg$train_auc) & !is.na(agg$holdout_auc))
  n_tr_1  <- sum(agg$train_auc   >= 1.0, na.rm = TRUE)
  n_ho_1  <- sum(agg$holdout_auc >= 1.0, na.rm = TRUE)
  n_both  <- sum(agg$train_auc >= 1.0 & agg$holdout_auc >= 1.0,
                 na.rm = TRUE)
  cat(sprintf("\nPerms hitting AUC = 1.0 (n_valid = %d):\n", n_valid))
  cat(sprintf("  train  only : %d (%.2f%%)\n", n_tr_1, 100 * n_tr_1 / n_valid))
  cat(sprintf("  holdout only: %d (%.2f%%)\n", n_ho_1, 100 * n_ho_1 / n_valid))
  cat(sprintf("  BOTH        : %d (%.2f%%)\n", n_both, 100 * n_both / n_valid))
} else {
  n_missing_draws <- sum(!vapply(seq_len(N_PERMS), function(i)
    file.exists(perm_drawn_csv(manifest[i, ])), logical(1L)))
  n_missing_rmats <- sum(!has_rmats)
  cat(sprintf("Phase 1: writing %d missing draws (of %d) ...\n",
              n_missing_draws, N_PERMS))
  for (i in seq_len(N_PERMS)) {
    if (!file.exists(perm_drawn_csv(manifest[i, ])))
      do_draw_one(manifest[i, ])
  }
  cat(sprintf("  all %d drawn_samples.csv now present\n", N_PERMS))

  script_sh <- file.path(WORKSPACE,
                         "scripts/gel/null_lasso_rmats_stat_array.sh")
  log_dir   <- file.path(WORKSPACE, "hpc_logs")
  cat("\n=========================================================\n")
  cat(sprintf("PHASE 2 (%d rMATS stat jobs) -- submit on HPC:\n\n",
              n_missing_rmats))
  job_tag <- if (nzchar(SAMPLE_TAG))
    sprintf("%s_%s", SAMPLE_TAG, AGG_TAG) else AGG_TAG
  cat(sprintf("  bsub -q short -P re_gecip_machine_learning \\\n"))
  cat(sprintf("    -J \"null_%s[1-%d]%%50\" \\\n", job_tag, N_PERMS))
  cat(sprintf("    -o %s/null_%s_%%J_%%I.out \\\n", log_dir, job_tag))
  cat(sprintf("    -e %s/null_%s_%%J_%%I.err \\\n", log_dir, job_tag))
  cat(sprintf("    -cwd %s -L /bin/bash \\\n", WORKSPACE))
  cat(sprintf("    bash %s %s\n\n", script_sh, MANIFEST))
  cat("(%50 caps concurrent jobs; raise/lower per HPC load.)\n")
  cat("Then re-run this R script (same env vars) to aggregate.\n")
  cat("=========================================================\n")
}
