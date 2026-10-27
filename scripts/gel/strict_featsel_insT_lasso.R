#!/usr/bin/env Rscript
#
# strict_featsel_insT_lasso.R
# Reviewer response: the insT pipeline from loo_evaluate_insT_lasso.R, with
# every feature-selection step restricted to training samples. Phase 3 of 3
# (after strict_featsel_prepare.R and strict_featsel_rmats_stat.sh).
#
# Differences from loo_evaluate_insT_lasso.R:
#   - Significance: rMATS FDR re-run on insT training cases vs training
#     controls (non-insT and holdout cases removed), per cohort.
#   - Direction concordance, depth/non-NA filters, control centering means,
#     G_AUC and LASSO: insT training cases + training controls only.
#   - Holdout = 5 insT cases + 3 GEL + 3 Nava controls
#     (data/split_strict_featsel.csv), none of which touch any step above.
# Everything else (bb-logit features, MAX_PER_GENE=1, glmnet params,
# case-LOO folds, lambda rule) is unchanged.
#
# Before fitting, checks that the all-samples stat rerun reproduces the
# original FDRs; stops if it does not.
#
# Inputs:
#   outputs/psi_matrix.rds, outputs/nava_psi_matrix.rds
#   outputs/strict_featsel/rmats_stat/{gel,nava}_{all,train}/A5SS.MATS.JC.txt
#   data/split_strict_featsel.csv
#   data/sample_genotype.csv, outputs/nava_severity_labels.csv
#   outputs/loo_results_insT_lasso.rds        (read only, for overlap)
#   outputs/gel_controls_psi_cand_pool.rds    (read only, 6k event coverage)
#
# Outputs:
#   outputs/strict_featsel/loo_results_insT_lasso_strict.rds  (GEL protected)
#   outputs/strict_featsel/strict_featsel_summary.csv         (counts + AUCs only)
#
# Run: pixi run Rscript scripts/gel/strict_featsel_insT_lasso.R

# ---- libpath -----------------------------------------------------------------
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
  if (!is.null(candidate) && dir.exists(candidate)) .libPaths(c(.libPaths(), candidate))
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
}))
suppressPackageStartupMessages({
  library(glmnet)
})

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
PSI_PATH      <- "outputs/psi_matrix.rds"
NAVA_PSI_PATH <- "outputs/nava_psi_matrix.rds"
STAT_DIR      <- "outputs/strict_featsel/rmats_stat"
SPLIT_FILE    <- "data/split_strict_featsel.csv"
GENO_FILE     <- "data/sample_genotype.csv"
NAVA_SEV_FILE <- "outputs/nava_severity_labels.csv"
CANON_RDS     <- "outputs/loo_results_insT_lasso.rds"
POOL_6K_RDS   <- "outputs/gel_controls_psi_cand_pool.rds"
OUT_DIR       <- "outputs/strict_featsel"
OUT_RDS       <- file.path(OUT_DIR, "loo_results_insT_lasso_strict.rds")
OUT_SUMMARY   <- file.path(OUT_DIR, "strict_featsel_summary.csv")

N_GEL_CASES     <- 11L
N_GEL_CONTROLS  <- 33L
N_NAVA_CASES    <- 19L
N_NAVA_CONTROLS <- 21L

INST_HGVS    <- "n.64_65insT"
EVENT_TYPE   <- "A5SS"
FDR_THRESH   <- 0.05
MIN_DEPTH    <- 10L
MAX_PER_GENE <- 1L

LAMBDA_RULE  <- "lambda.1se"
ALPHA_LASSO  <- 1
SEED         <- 2026L

stat_file <- function(run) file.path(STAT_DIR, run, "A5SS.MATS.JC.txt")
stopifnot(
  "GEL PSI matrix not found"     = file.exists(PSI_PATH),
  "Nava PSI matrix not found"    = file.exists(NAVA_PSI_PATH),
  "Strict split not found -- run strict_featsel_prepare.R" = file.exists(SPLIT_FILE),
  "Genotype file not found"      = file.exists(GENO_FILE),
  "Nava severity file not found" = file.exists(NAVA_SEV_FILE),
  "Canonical insT RDS not found" = file.exists(CANON_RDS),
  "rMATS stat outputs missing -- run strict_featsel_rmats_stat.sh" =
    all(file.exists(stat_file(c("gel_all", "nava_all", "gel_train", "nava_train"))))
)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------
# Load PSI matrices (counts are the same whichever samples rMATS tested)
# --------------------------------------------------------------------------
cat("Loading PSI matrices...\n")
gel_obj  <- readRDS(PSI_PATH)
nava_obj <- readRDS(NAVA_PSI_PATH)
gel_ijc  <- gel_obj$ijc;  gel_sjc  <- gel_obj$sjc;  gel_events  <- gel_obj$events
nava_ijc <- nava_obj$ijc; nava_sjc <- nava_obj$sjc; nava_events <- nava_obj$events

gel_case_ids  <- paste0("gel:case:",     seq_len(N_GEL_CASES))
gel_ctrl_ids  <- paste0("gel:control:",  seq_len(N_GEL_CONTROLS))
nava_case_ids <- paste0("nava:case:",    seq_len(N_NAVA_CASES))
nava_ctrl_ids <- paste0("nava:control:", seq_len(N_NAVA_CONTROLS))
stopifnot(all(c(gel_case_ids, gel_ctrl_ids)   %in% rownames(gel_ijc)),
          all(c(nava_case_ids, nava_ctrl_ids) %in% rownames(nava_ijc)))

# --------------------------------------------------------------------------
# Split + genotypes
# --------------------------------------------------------------------------
read_sep_tolerant <- function(path) {
  first <- readLines(path, n = 1L)
  sep   <- if (grepl(",", first, fixed = TRUE)) "," else "\t"
  read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
}
geno     <- read.csv(GENO_FILE, stringsAsFactors = FALSE)
nava_lbl <- read_sep_tolerant(NAVA_SEV_FILE)
nava_lbl$sample_id <- paste0("nava:case:", nava_lbl$sample_n)
hgvs <- setNames(c(geno$hgvs, nava_lbl$hgvs), c(geno$sample_id, nava_lbl$sample_id))

split_df <- read.csv(SPLIT_FILE, stringsAsFactors = FALSE)
train_ids <- split_df$sample_id[split_df$set == "train"]
ho_ids    <- split_df$sample_id[split_df$set == "holdout"]

gel_insT  <- gel_case_ids[hgvs[gel_case_ids]   == INST_HGVS]
nava_insT <- nava_case_ids[hgvs[nava_case_ids] == INST_HGVS]

gel_insT_train  <- intersect(gel_insT,  train_ids)
nava_insT_train <- intersect(nava_insT, train_ids)
gel_insT_ho     <- intersect(gel_insT,  ho_ids)
nava_insT_ho    <- intersect(nava_insT, ho_ids)
gel_ctrl_train  <- intersect(gel_ctrl_ids,  train_ids)
nava_ctrl_train <- intersect(nava_ctrl_ids, train_ids)
gel_ctrl_ho     <- intersect(gel_ctrl_ids,  ho_ids)
nava_ctrl_ho    <- intersect(nava_ctrl_ids, ho_ids)

cat(sprintf("GEL  train: %d insT vs %d ctrl | holdout: %d insT + %d ctrl\n",
            length(gel_insT_train), length(gel_ctrl_train),
            length(gel_insT_ho), length(gel_ctrl_ho)))
cat(sprintf("Nava train: %d insT vs %d ctrl | holdout: %d insT + %d ctrl\n",
            length(nava_insT_train), length(nava_ctrl_train),
            length(nava_insT_ho), length(nava_ctrl_ho)))
stopifnot(length(gel_insT_train) > 0L, length(nava_insT_train) > 0L,
          length(gel_ctrl_ho) > 0L, length(nava_ctrl_ho) > 0L,
          length(gel_insT_ho) + length(nava_insT_ho) > 0L)

gel_train_rows  <- c(gel_insT_train,  gel_ctrl_train)
nava_train_rows <- c(nava_insT_train, nava_ctrl_train)

# --------------------------------------------------------------------------
# rMATS stat outputs -> event_id (same construction as build_psi_matrix.R)
# --------------------------------------------------------------------------
read_stat <- function(run) {
  df <- read.delim(stat_file(run), stringsAsFactors = FALSE)
  coord_key      <- do.call(paste, c(df[, 4:11], sep = ":"))
  coord_nostrand <- sub("^([^:]+):[^:]+:", "\\1:", coord_key)
  ev   <- paste(EVENT_TYPE, df[[3]], coord_nostrand, sep = ":")
  keep <- !duplicated(ev)
  data.frame(event_id = ev[keep],
             FDR      = suppressWarnings(as.numeric(df$FDR[keep])),
             stringsAsFactors = FALSE)
}

# Step 1: the all-samples rerun must reproduce the original significant set
check_repro <- function(run, events) {
  st  <- read_stat(run)
  ev  <- events[events$EventType == EVENT_TYPE, c("event_id", "FDR")]
  m   <- match(ev$event_id, st$event_id)
  if (anyNA(m)) stop(run, ": ", sum(is.na(m)), " original events missing from stat rerun")
  fdr_new  <- st$FDR[m]
  max_diff <- max(abs(ev$FDR - fdr_new), na.rm = TRUE)
  sig_orig <- !is.na(ev$FDR) & ev$FDR < FDR_THRESH
  sig_new  <- !is.na(fdr_new) & fdr_new < FDR_THRESH
  cat(sprintf("  %s: %d events | sig orig=%d rerun=%d | disagree=%d | max |dFDR|=%.2g | NA mismatch=%d\n",
              run, nrow(ev), sum(sig_orig), sum(sig_new), sum(sig_orig != sig_new),
              max_diff, sum(is.na(ev$FDR) != is.na(fdr_new))))
  if (any(sig_orig != sig_new))
    stop(run, ": all-samples rerun does not reproduce the original significant set")
}
cat("\nStep 1: reproduction check (all cases vs all controls)\n")
check_repro("gel_all",  gel_events)
check_repro("nava_all", nava_events)

# --------------------------------------------------------------------------
# Candidate pool: sig in both (training-only FDR) + concordant direction
# --------------------------------------------------------------------------
cat("\nStep 2: candidate pool from training-only significance\n")
gel_tr  <- read_stat("gel_train")
nava_tr <- read_stat("nava_train")

gel_a5  <- gel_events[gel_events$EventType   == EVENT_TYPE, ]
nava_a5 <- nava_events[nava_events$EventType == EVENT_TYPE, ]
gel_a5$FDR_train  <- gel_tr$FDR[match(gel_a5$event_id,   gel_tr$event_id)]
nava_a5$FDR_train <- nava_tr$FDR[match(nava_a5$event_id, nava_tr$event_id)]

gel_sig  <- gel_a5[!is.na(gel_a5$FDR_train)   & gel_a5$FDR_train  < FDR_THRESH, ]
nava_sig <- nava_a5[!is.na(nava_a5$FDR_train) & nava_a5$FDR_train < FDR_THRESH, ]
both_keys <- intersect(gel_sig$coord_nostrand, nava_sig$coord_nostrand)
gel_m  <- gel_sig[match(both_keys,  gel_sig$coord_nostrand),  ]
nava_m <- nava_sig[match(both_keys, nava_sig$coord_nostrand), ]
cat(sprintf("  Sig (train FDR<%.2f): GEL=%d  Nava=%d  both=%d\n",
            FDR_THRESH, nrow(gel_sig), nrow(nava_sig), length(both_keys)))

bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}
gel_lp  <- bb_logit(gel_ijc[,  gel_m$event_id,  drop = FALSE],
                    gel_sjc[,  gel_m$event_id,  drop = FALSE])
nava_lp <- bb_logit(nava_ijc[, nava_m$event_id, drop = FALSE],
                    nava_sjc[, nava_m$event_id, drop = FALSE])
colnames(nava_lp) <- gel_m$event_id   # canonical names = GEL event_id

dl_gel  <- colMeans(gel_lp[gel_insT_train,    , drop = FALSE], na.rm = TRUE) -
           colMeans(gel_lp[gel_ctrl_train,    , drop = FALSE], na.rm = TRUE)
dl_nava <- colMeans(nava_lp[nava_insT_train,  , drop = FALSE], na.rm = TRUE) -
           colMeans(nava_lp[nava_ctrl_train,  , drop = FALSE], na.rm = TRUE)
dl_gel[!is.finite(dl_gel)]   <- NA_real_
dl_nava[!is.finite(dl_nava)] <- NA_real_
concordant <- !is.na(dl_gel) & !is.na(dl_nava) &
              sign(dl_gel) == sign(dl_nava) & sign(dl_gel) != 0

a5ss_df <- data.frame(
  gel_ev_id   = gel_m$event_id[concordant],
  nava_ev_id  = nava_m$event_id[concordant],
  GeneSymbol  = gel_m$GeneSymbol[concordant],
  dlogit_GEL  = unname(dl_gel[concordant]),
  dlogit_Nava = unname(dl_nava[concordant]),
  direction   = sign(unname(dl_gel[concordant])),
  stringsAsFactors = FALSE)
cat(sprintf("  Concordant A5SS pool: %d\n", nrow(a5ss_df)))
stopifnot("Empty concordant pool" = nrow(a5ss_df) > 0L)

gel_lp  <- gel_lp[,  a5ss_df$gel_ev_id, drop = FALSE]
nava_lp <- nava_lp[, a5ss_df$gel_ev_id, drop = FALSE]

# --------------------------------------------------------------------------
# Depth + non-NA filters on training samples only
# --------------------------------------------------------------------------
gel_depth  <- gel_ijc[gel_train_rows,   a5ss_df$gel_ev_id,  drop = FALSE] +
              gel_sjc[gel_train_rows,   a5ss_df$gel_ev_id,  drop = FALSE]
nava_depth <- nava_ijc[nava_train_rows, a5ss_df$nava_ev_id, drop = FALSE] +
              nava_sjc[nava_train_rows, a5ss_df$nava_ev_id, drop = FALSE]
min_gel  <- suppressWarnings(apply(gel_depth,  2, min, na.rm = TRUE))
min_nava <- suppressWarnings(apply(nava_depth, 2, min, na.rm = TRUE))
depth_pass <- is.finite(min_gel) & is.finite(min_nava) &
              min_gel >= MIN_DEPTH & min_nava >= MIN_DEPTH
nona_pass  <- colSums(is.na(gel_lp[gel_train_rows,   , drop = FALSE])) == 0L &
              colSums(is.na(nava_lp[nava_train_rows, , drop = FALSE])) == 0L

a5ss_cand <- a5ss_df[depth_pass & nona_pass, , drop = FALSE]
rownames(a5ss_cand) <- NULL
cat(sprintf("  After depth (>=%d) + non-NA filter (training samples): %d\n",
            MIN_DEPTH, nrow(a5ss_cand)))
stopifnot("No sites passed depth/NA filters" = nrow(a5ss_cand) > 0L)

gel_lp  <- gel_lp[,  a5ss_cand$gel_ev_id, drop = FALSE]
nava_lp <- nava_lp[, a5ss_cand$gel_ev_id, drop = FALSE]

# --------------------------------------------------------------------------
# G_AUC (training) -> best site per gene
# --------------------------------------------------------------------------
auc_two_sample <- function(case_scores, ctrl_scores) {
  case_scores <- case_scores[!is.na(case_scores)]
  ctrl_scores <- ctrl_scores[!is.na(ctrl_scores)]
  if (length(case_scores) == 0L || length(ctrl_scores) == 0L) return(NA_real_)
  r <- rank(c(case_scores, ctrl_scores))
  U <- sum(r[seq_along(case_scores)]) - length(case_scores) * (length(case_scores) + 1L) / 2
  U / (length(case_scores) * length(ctrl_scores))
}
auc_per_site <- function(case_mat, ctrl_mat, direction) {
  vapply(seq_len(ncol(case_mat)), function(j)
    auc_two_sample(direction[j] * case_mat[, j], direction[j] * ctrl_mat[, j]),
    numeric(1L))
}
auc_gel  <- auc_per_site(gel_lp[gel_insT_train,   , drop = FALSE],
                         gel_lp[gel_ctrl_train,   , drop = FALSE], a5ss_cand$direction)
auc_nava <- auc_per_site(nava_lp[nava_insT_train, , drop = FALSE],
                         nava_lp[nava_ctrl_train, , drop = FALSE], a5ss_cand$direction)
a1 <- auc_gel - 0.5; a2 <- auc_nava - 0.5
g_auc <- ifelse(!is.na(a1) & !is.na(a2) & a1 > 0 & a2 > 0, 0.5 + sqrt(pmax(a1, 0) * pmax(a2, 0)), NA_real_)
a5ss_cand$g_auc <- g_auc; a5ss_cand$auc_gel <- auc_gel; a5ss_cand$auc_nava <- auc_nava
a5ss_cand_pre <- a5ss_cand

eligible <- which(!is.na(g_auc))
ord      <- eligible[order(-g_auc[eligible])]
keep_idx <- ord[!duplicated(a5ss_cand$GeneSymbol[ord])]   # MAX_PER_GENE = 1
stopifnot(MAX_PER_GENE == 1L, length(keep_idx) > 0L)
a5ss_cand <- a5ss_cand[keep_idx, , drop = FALSE]
rownames(a5ss_cand) <- NULL
cat(sprintf("  After best-site-per-gene (G_AUC on training): %d\n", nrow(a5ss_cand)))

gel_lp    <- gel_lp[,  a5ss_cand$gel_ev_id, drop = FALSE]
nava_lp   <- nava_lp[, a5ss_cand$gel_ev_id, drop = FALSE]
direction <- a5ss_cand$direction

# --------------------------------------------------------------------------
# Features: centred on TRAINING control means, direction-oriented
# --------------------------------------------------------------------------
gel_ctrl_means  <- colMeans(gel_lp[gel_ctrl_train,   , drop = FALSE])
nava_ctrl_means <- colMeans(nava_lp[nava_ctrl_train, , drop = FALSE])
orient_gel  <- sweep(sweep(gel_lp,  2, gel_ctrl_means,  "-"), 2, direction, "*")
orient_nava <- sweep(sweep(nava_lp, 2, nava_ctrl_means, "-"), 2, direction, "*")

X <- rbind(orient_gel[gel_train_rows,   , drop = FALSE],
           orient_nava[nava_train_rows, , drop = FALSE])
y <- c(rep(1L, length(gel_insT_train)),  rep(0L, length(gel_ctrl_train)),
       rep(1L, length(nava_insT_train)), rep(0L, length(nava_ctrl_train)))
cohort_train <- c(rep("GEL", length(gel_train_rows)), rep("Nava", length(nava_train_rows)))
stopifnot(nrow(X) == length(y), !anyNA(X))

case_idx <- which(y == 1L); ctrl_idx <- which(y == 0L); n_cases <- length(case_idx)
foldid <- integer(length(y))
foldid[case_idx] <- seq_along(case_idx)
foldid[ctrl_idx] <- ((seq_along(ctrl_idx) - 1L) %% n_cases) + 1L

# --------------------------------------------------------------------------
# LASSO (same settings as canonical)
# --------------------------------------------------------------------------
cat(sprintf("\nStep 3: cv.glmnet on %d samples x %d sites (%d case-LOO folds)\n",
            nrow(X), ncol(X), n_cases))
set.seed(SEED)
cv_fit <- cv.glmnet(x = X, y = y, family = "binomial", alpha = ALPHA_LASSO,
                    foldid = foldid, standardize = TRUE, lower.limits = 0,
                    type.measure = "deviance", keep = TRUE)

coef_at <- function(s) { co <- coef(cv_fit, s = s); v <- as.numeric(co); names(v) <- rownames(co); v }
chosen_label <- LAMBDA_RULE
beta_all <- coef_at(chosen_label)
if (sum(beta_all[-1] != 0) == 0L) {
  warning("lambda.1se selected 0 sites; falling back to lambda.min")
  chosen_label <- "lambda.min"
  beta_all <- coef_at(chosen_label)
}
intercept <- beta_all[1]; beta <- beta_all[-1]
sel <- which(beta != 0)
K   <- length(sel)
cat(sprintf("  %s: K = %d selected sites\n", chosen_label, K))

selected_sites <- data.frame(gel_ev_id  = a5ss_cand$gel_ev_id[sel],
                             GeneSymbol = a5ss_cand$GeneSymbol[sel],
                             direction  = direction[sel],
                             g_auc      = a5ss_cand$g_auc[sel],
                             coef       = unname(beta[sel]),
                             stringsAsFactors = FALSE)
selected_sites <- selected_sites[order(-selected_sites$coef), , drop = FALSE]
print(selected_sites, row.names = FALSE)

chosen_lambda <- if (chosen_label == "lambda.1se") cv_fit$lambda.1se else cv_fit$lambda.min
preval_lp <- as.numeric(cv_fit$fit.preval[, which.min(abs(cv_fit$lambda - chosen_lambda))])

# --------------------------------------------------------------------------
# Step 4: score held-out samples (never used above)
# --------------------------------------------------------------------------
sigmoid <- function(z) 1 / (1 + exp(-z))
score_rows <- function(orient, ids, cohort, group) {
  Xh <- orient[ids, , drop = FALSE]
  if (anyNA(Xh[, sel])) stop("NA feature for held-out sample(s) in ", cohort,
                             " at selected sites -- check counts")
  lp <- as.numeric(Xh %*% beta + intercept)
  data.frame(sample_id = ids, cohort = cohort, group = group,
             linpred = lp, prob = sigmoid(lp),
             uniform_sum = rowSums(Xh[, sel, drop = FALSE]),
             stringsAsFactors = FALSE)
}
ho_scores <- rbind(
  score_rows(orient_gel,  gel_insT_ho,  "GEL",  "insT"),
  score_rows(orient_gel,  gel_ctrl_ho,  "GEL",  "control"),
  score_rows(orient_nava, nava_insT_ho, "Nava", "insT"),
  score_rows(orient_nava, nava_ctrl_ho, "Nava", "control"))

# Read depth of held-out samples at the selected sites (informational: they
# were not part of the depth filter)
ho_depth <- c(
  apply(gel_ijc[c(gel_insT_ho, gel_ctrl_ho), selected_sites$gel_ev_id, drop = FALSE] +
        gel_sjc[c(gel_insT_ho, gel_ctrl_ho), selected_sites$gel_ev_id, drop = FALSE], 1, min),
  apply(nava_ijc[c(nava_insT_ho, nava_ctrl_ho),
                 a5ss_cand$nava_ev_id[match(selected_sites$gel_ev_id, a5ss_cand$gel_ev_id)], drop = FALSE] +
        nava_sjc[c(nava_insT_ho, nava_ctrl_ho),
                 a5ss_cand$nava_ev_id[match(selected_sites$gel_ev_id, a5ss_cand$gel_ev_id)], drop = FALSE], 1, min))
ho_scores$min_depth_selected <- unname(ho_depth[ho_scores$sample_id])

cat("\nHeld-out samples:\n")
print(ho_scores[, c("cohort", "group", "prob", "uniform_sum", "min_depth_selected")],
      row.names = FALSE, digits = 3)

ho_metrics <- function(co) {
  s  <- ho_scores[ho_scores$cohort %in% co, ]
  cs <- s$prob[s$group == "insT"]; ks <- s$prob[s$group == "control"]
  cu <- s$uniform_sum[s$group == "insT"]; ku <- s$uniform_sum[s$group == "control"]
  c(n_case = length(cs), n_ctrl = length(ks),
    auc_prob = auc_two_sample(cs, ks), auc_uniform = auc_two_sample(cu, ku),
    n_case_above_all_ctrl = if (length(cs) && length(ks)) sum(cs > max(ks)) else NA)
}

# --------------------------------------------------------------------------
# Comparison with the canonical signature and 6k event coverage
# --------------------------------------------------------------------------
canon <- readRDS(CANON_RDS)
canon_events <- canon$a5ss_cand$gel_ev_id[canon$selected_idx]
canon_pool   <- canon$candidate_sites_pre$gel_ev_id
n_overlap    <- length(intersect(selected_sites$gel_ev_id, canon_events))
cat(sprintf("\nOverlap with canonical insT signature: %d of %d new sites are in the canonical %d\n",
            n_overlap, K, length(canon_events)))
cat(sprintf("Depth-filtered pool: strict=%d  canonical=%d  shared=%d\n",
            nrow(a5ss_cand_pre), length(canon_pool),
            length(intersect(a5ss_cand_pre$gel_ev_id, canon_pool))))

n_in_6k <- NA_integer_
if (file.exists(POOL_6K_RDS)) {
  pool6k  <- readRDS(POOL_6K_RDS)
  n_in_6k <- sum(selected_sites$gel_ev_id %in% colnames(pool6k$psi))
  rm(pool6k)
  cat(sprintf("Selected sites present in 6k control matrix: %d of %d%s\n", n_in_6k, K,
              if (n_in_6k < K) "  (6k scoring needs those extracted first)" else ""))
}

# --------------------------------------------------------------------------
# Summary (counts and AUCs only -- airlock-safe) + full RDS
# --------------------------------------------------------------------------
m_gel <- ho_metrics("GEL"); m_nava <- ho_metrics("Nava"); m_all <- ho_metrics(c("GEL", "Nava"))
is_case <- y == 1L
summary_df <- data.frame(
  metric = c("sig_both_train_fdr", "concordant_pool", "depth_filtered_pool",
             "best_per_gene_pool", "K_selected", "K_overlap_canonical", "K_canonical",
             "K_in_6k_matrix",
             "train_loo_auc_gel", "train_loo_auc_nava", "train_loo_auc_combined",
             paste0("holdout_gel_",  names(m_gel)),
             paste0("holdout_nava_", names(m_nava)),
             paste0("holdout_all_",  names(m_all))),
  value = c(length(both_keys), nrow(a5ss_df), nrow(a5ss_cand_pre),
            nrow(a5ss_cand), K, n_overlap, length(canon_events), n_in_6k,
            auc_two_sample(preval_lp[is_case & cohort_train == "GEL"],
                           preval_lp[!is_case & cohort_train == "GEL"]),
            auc_two_sample(preval_lp[is_case & cohort_train == "Nava"],
                           preval_lp[!is_case & cohort_train == "Nava"]),
            auc_two_sample(preval_lp[is_case], preval_lp[!is_case]),
            m_gel, m_nava, m_all),
  stringsAsFactors = FALSE)
cat("\nSummary:\n")
print(summary_df, row.names = FALSE, digits = 3)
write.csv(summary_df, OUT_SUMMARY, row.names = FALSE)

saveRDS(list(
  candidate_sites_pre = a5ss_cand_pre,
  a5ss_cand           = a5ss_cand,
  direction           = direction,
  gel_ctrl_means      = gel_ctrl_means,
  nava_ctrl_means     = nava_ctrl_means,
  X_train             = X,
  y_train             = y,
  cohort_train        = cohort_train,
  sample_ids_train    = rownames(X),
  foldid              = foldid,
  cv_fit              = cv_fit,
  chosen_lambda_label = chosen_label,
  intercept           = intercept,
  beta                = beta,
  selected_sites      = selected_sites,
  preval_lp           = setNames(preval_lp, rownames(X)),
  ho_scores           = ho_scores,
  summary             = summary_df,
  config = list(SPLIT_FILE = SPLIT_FILE, FDR_THRESH = FDR_THRESH, MIN_DEPTH = MIN_DEPTH,
                MAX_PER_GENE = MAX_PER_GENE, LAMBDA_RULE = LAMBDA_RULE,
                NOTE = "All feature selection on insT training cases + training controls only")
), OUT_RDS)
cat("\nSaved:", OUT_RDS, "\n      ", OUT_SUMMARY, "\n")
