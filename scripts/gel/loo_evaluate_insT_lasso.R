#!/usr/bin/env Rscript
#
# loo_evaluate_insT_lasso.R
# Exploratory LASSO trained on n.64_65insT carriers only vs controls.
#
# Candidate space: A5SS events from the 483 sig-in-both concordant pool
# (FDR<0.05 in both GEL and Nava all-cases rMATS, direction-concordant).
# Within that pool, A5SS events are filtered by min(IJC+SJC) >= 10 across
# insT cases and controls in both cohorts, non-NA, and MAX_PER_GENE=1
# (best G_AUC site per gene, AUC computed on insT cases vs controls).
#
# Train/holdout split: same 70/30 stratified split as v4 (split_dev_holdout_v4_mixed.csv).
# insT carriers assigned to holdout in that split remain blind; only train-set
# insT carriers enter the LASSO. Holdout AUC/Glass reported separately.
#
# Feature: direction * (bb-logit-PSI - per-cohort-ctrl-mean), per v4 formula.
# LASSO: binomial, alpha=1, lower.limits=0, foldid=case-LOO.
#
# Inputs:
#   outputs/psi_matrix.rds
#   outputs/nava_psi_matrix.rds
#   data/sample_genotype.csv
#   outputs/nava_severity_labels.csv
#
# Outputs:
#   outputs/loo_results_insT_lasso.rds
#
# Run: pixi run Rscript scripts/gel/loo_evaluate_insT_lasso.R

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
SPLIT_FILE    <- "data/split_dev_holdout_v4_mixed.csv"
GENO_FILE     <- "data/sample_genotype.csv"
NAVA_SEV_FILE <- "outputs/nava_severity_labels.csv"
OUT_RDS       <- "outputs/loo_results_insT_lasso.rds"

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

B_BOOT       <- 2000L
SEED         <- 2026L

stopifnot(
  "GEL PSI matrix not found"     = file.exists(PSI_PATH),
  "Nava PSI matrix not found"    = file.exists(NAVA_PSI_PATH),
  "Split file not found"         = file.exists(SPLIT_FILE),
  "Genotype file not found"      = file.exists(GENO_FILE),
  "Nava severity file not found" = file.exists(NAVA_SEV_FILE)
)

# --------------------------------------------------------------------------
# Load PSI matrices
# --------------------------------------------------------------------------
cat("Loading PSI matrices...\n")
gel_obj  <- readRDS(PSI_PATH)
nava_obj <- readRDS(NAVA_PSI_PATH)

stopifnot(
  all(c("ijc", "sjc") %in% names(gel_obj)),
  all(c("ijc", "sjc") %in% names(nava_obj))
)

gel_psi     <- gel_obj$psi
gel_ijc     <- gel_obj$ijc
gel_sjc     <- gel_obj$sjc
gel_events  <- gel_obj$events

nava_psi    <- nava_obj$psi
nava_ijc    <- nava_obj$ijc
nava_sjc    <- nava_obj$sjc
nava_events <- nava_obj$events

stopifnot(
  identical(dimnames(gel_ijc),  dimnames(gel_psi)),
  identical(dimnames(gel_sjc),  dimnames(gel_psi)),
  identical(dimnames(nava_ijc), dimnames(nava_psi)),
  identical(dimnames(nava_sjc), dimnames(nava_psi))
)

gel_case_ids  <- paste0("gel:case:",     seq_len(N_GEL_CASES))
gel_ctrl_ids  <- paste0("gel:control:",  seq_len(N_GEL_CONTROLS))
nava_case_ids <- paste0("nava:case:",    seq_len(N_NAVA_CASES))
nava_ctrl_ids <- paste0("nava:control:", seq_len(N_NAVA_CONTROLS))

stopifnot(
  all(gel_case_ids  %in% rownames(gel_psi)),
  all(gel_ctrl_ids  %in% rownames(gel_psi)),
  all(nava_case_ids %in% rownames(nava_psi)),
  all(nava_ctrl_ids %in% rownames(nava_psi))
)
cat(sprintf("GEL PSI:  %d samples x %d events\n", nrow(gel_psi),  ncol(gel_psi)))
cat(sprintf("Nava PSI: %d samples x %d events\n", nrow(nava_psi), ncol(nava_psi)))

# --------------------------------------------------------------------------
# Load genotype and SGE data
# --------------------------------------------------------------------------
geno <- read.csv(GENO_FILE, stringsAsFactors = FALSE)

read_sep_tolerant <- function(path) {
  first <- readLines(path, n = 1L)
  sep   <- if (grepl(",", first, fixed = TRUE)) "," else "\t"
  read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
}
nava_lbl <- read_sep_tolerant(NAVA_SEV_FILE)
nava_lbl$sample_id <- paste0("nava:case:", nava_lbl$sample_n)

all_hgvs_lookup <- setNames(
  c(geno$hgvs,      nava_lbl$hgvs),
  c(geno$sample_id, nava_lbl$sample_id)
)

gel_insT_ids   <- names(all_hgvs_lookup)[
  names(all_hgvs_lookup) %in% gel_case_ids  & all_hgvs_lookup == INST_HGVS]
nava_insT_ids  <- names(all_hgvs_lookup)[
  names(all_hgvs_lookup) %in% nava_case_ids & all_hgvs_lookup == INST_HGVS]

# Apply the same train/holdout split used for v4, restricted to insT carriers
split_df <- read.csv(SPLIT_FILE, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "cohort", "set") %in% names(split_df)),
          all(split_df$set %in% c("train", "holdout")))

train_ids_split   <- split_df$sample_id[split_df$set == "train"]
holdout_ids_split <- split_df$sample_id[split_df$set == "holdout"]

gel_insT_train_ids  <- intersect(gel_insT_ids,  train_ids_split)
gel_insT_ho_ids     <- intersect(gel_insT_ids,  holdout_ids_split)
nava_insT_train_ids <- intersect(nava_insT_ids, train_ids_split)
nava_insT_ho_ids    <- intersect(nava_insT_ids, holdout_ids_split)

cat(sprintf("GEL  insT: train=%d  holdout=%d  (of %d cases)\n",
            length(gel_insT_train_ids),  length(gel_insT_ho_ids),  N_GEL_CASES))
cat(sprintf("Nava insT: train=%d  holdout=%d  (of %d cases)\n",
            length(nava_insT_train_ids), length(nava_insT_ho_ids), N_NAVA_CASES))
stopifnot(
  "No GEL insT carriers found"       = length(gel_insT_ids)         > 0L,
  "No Nava insT carriers found"      = length(nava_insT_ids)        > 0L,
  "No GEL insT carriers in train"    = length(gel_insT_train_ids)   > 0L,
  "No Nava insT carriers in train"   = length(nava_insT_train_ids)  > 0L,
  "GEL insT IDs not in psi matrix"   = all(gel_insT_ids  %in% rownames(gel_psi)),
  "Nava insT IDs not in psi matrix"  = all(nava_insT_ids %in% rownames(nava_psi))
)

# --------------------------------------------------------------------------
# Identify candidate space: A5SS events from the 483 sig-in-both pool
# Match by EventType + coord_nostrand (verbatim from explore_dlogit_genotype_split.R).
# --------------------------------------------------------------------------
cat("\nIdentifying candidate space...\n")
gel_sig  <- gel_events[!is.na(gel_events$FDR)  & gel_events$FDR  < FDR_THRESH, ]
nava_sig <- nava_events[!is.na(nava_events$FDR) & nava_events$FDR < FDR_THRESH, ]

gel_sig$key  <- paste(gel_sig$EventType,  gel_sig$coord_nostrand,  sep = "|")
nava_sig$key <- paste(nava_sig$EventType, nava_sig$coord_nostrand, sep = "|")

both_keys <- intersect(gel_sig$key, nava_sig$key)
gel_m  <- gel_sig[match(both_keys,  gel_sig$key),  ]
nava_m <- nava_sig[match(both_keys, nava_sig$key), ]

# Concordance uses bb-logit dlogit-PSI (all cases vs all controls, per
# cohort) -- same definition as figure1.R panel D / figure_s1.R. Matches
# the transform used for LASSO features downstream, so orientation and
# candidate pool share a metric. See figure_s1.R for the rationale for
# switching from rMATS IncLevelDifference.
bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}
gel_lp_sig  <- bb_logit(gel_ijc[,  gel_m$event_id,  drop = FALSE],
                         gel_sjc[,  gel_m$event_id,  drop = FALSE])
nava_lp_sig <- bb_logit(nava_ijc[, nava_m$event_id, drop = FALSE],
                         nava_sjc[, nava_m$event_id, drop = FALSE])

dl_gel_sig  <- colMeans(gel_lp_sig[ gel_case_ids,  , drop = FALSE], na.rm = TRUE) -
               colMeans(gel_lp_sig[ gel_ctrl_ids,  , drop = FALSE], na.rm = TRUE)
dl_nava_sig <- colMeans(nava_lp_sig[nava_case_ids, , drop = FALSE], na.rm = TRUE) -
               colMeans(nava_lp_sig[nava_ctrl_ids, , drop = FALSE], na.rm = TRUE)
dl_gel_sig[ !is.finite(dl_gel_sig)]  <- NA_real_
dl_nava_sig[!is.finite(dl_nava_sig)] <- NA_real_

concordant <- !is.na(dl_gel_sig) & !is.na(dl_nava_sig) &
              sign(dl_gel_sig) == sign(dl_nava_sig) &
              sign(dl_gel_sig) != 0

sig_df <- data.frame(
  key            = both_keys[concordant],
  gel_ev_id      = gel_m$event_id[concordant],
  nava_ev_id     = nava_m$event_id[concordant],
  EventType      = gel_m$EventType[concordant],
  GeneSymbol     = gel_m$GeneSymbol[concordant],
  IncLevelDiff   = gel_m$IncLevelDiff[concordant],   # kept for reference
  dlogit_GEL     = unname(dl_gel_sig[concordant]),
  dlogit_Nava    = unname(dl_nava_sig[concordant]),
  coord_nostrand = gel_m$coord_nostrand[concordant],
  stringsAsFactors = FALSE
)
cat(sprintf("  Sig-in-both concordant (all types): %d\n", nrow(sig_df)))

# Filter to A5SS. Direction now comes from bb-logit dlogit (matches the
# feature scale used by the LASSO), not the raw rMATS IncLevelDifference.
a5ss_df <- sig_df[sig_df$EventType == EVENT_TYPE, ]
a5ss_df$direction <- sign(a5ss_df$dlogit_GEL)
cat(sprintf("  A5SS events: %d\n", nrow(a5ss_df)))
stopifnot("No A5SS events in concordant pool" = nrow(a5ss_df) > 0L)
stopifnot(!any(a5ss_df$direction == 0))

# --------------------------------------------------------------------------
# bb-logit PSI on A5SS events
# Canonical column name = gel_ev_id; nava matrix renamed to match.
# --------------------------------------------------------------------------

gel_lp_a5ss  <- bb_logit(gel_ijc[,  a5ss_df$gel_ev_id,  drop = FALSE],
                          gel_sjc[,  a5ss_df$gel_ev_id,  drop = FALSE])
nava_lp_a5ss <- bb_logit(nava_ijc[, a5ss_df$nava_ev_id, drop = FALSE],
                          nava_sjc[, a5ss_df$nava_ev_id, drop = FALSE])
colnames(nava_lp_a5ss) <- a5ss_df$gel_ev_id  # canonical names for rbind

stopifnot(identical(colnames(gel_lp_a5ss), colnames(nava_lp_a5ss)))

# --------------------------------------------------------------------------
# Depth filter: min(IJC+SJC) >= MIN_DEPTH across ALL insT cases (train +
# holdout) and controls — same scope as v4.
# --------------------------------------------------------------------------
gel_all_ids   <- c(gel_insT_ids,  gel_ctrl_ids)   # all insT train + holdout + ctrl
nava_all_ids  <- c(nava_insT_ids, nava_ctrl_ids)

gel_depth  <- gel_ijc[gel_all_ids,   a5ss_df$gel_ev_id,  drop = FALSE] +
              gel_sjc[gel_all_ids,   a5ss_df$gel_ev_id,  drop = FALSE]
nava_depth <- nava_ijc[nava_all_ids, a5ss_df$nava_ev_id, drop = FALSE] +
              nava_sjc[nava_all_ids, a5ss_df$nava_ev_id, drop = FALSE]

min_depth_gel  <- apply(gel_depth,  2, min, na.rm = TRUE)
min_depth_nava <- apply(nava_depth, 2, min, na.rm = TRUE)
depth_pass <- is.finite(min_depth_gel)  & is.finite(min_depth_nava) &
              min_depth_gel >= MIN_DEPTH & min_depth_nava >= MIN_DEPTH

# Non-NA filter across insT cases and controls
na_gel  <- colSums(is.na(gel_lp_a5ss[gel_all_ids,   , drop = FALSE])) > 0L
na_nava <- colSums(is.na(nava_lp_a5ss[nava_all_ids, , drop = FALSE])) > 0L
nona_pass <- !na_gel & !na_nava

pass <- depth_pass & nona_pass
a5ss_cand <- a5ss_df[pass, , drop = FALSE]
rownames(a5ss_cand) <- NULL

cat(sprintf("  After depth (>=%d) + non-NA filter: %d A5SS sites\n", MIN_DEPTH, nrow(a5ss_cand)))
stopifnot("No A5SS sites passed depth/NA filters" = nrow(a5ss_cand) > 0L)

gel_lp_cand  <- gel_lp_a5ss[,  a5ss_cand$gel_ev_id, drop = FALSE]
nava_lp_cand <- nava_lp_a5ss[, a5ss_cand$gel_ev_id, drop = FALSE]

# --------------------------------------------------------------------------
# Per-cohort control means
# --------------------------------------------------------------------------
gel_ctrl_means  <- colMeans(gel_lp_cand[gel_ctrl_ids,  , drop = FALSE], na.rm = TRUE)
nava_ctrl_means <- colMeans(nava_lp_cand[nava_ctrl_ids, , drop = FALSE], na.rm = TRUE)

# --------------------------------------------------------------------------
# G_AUC on insT cases for MAX_PER_GENE=1 pre-filter
# --------------------------------------------------------------------------
auc_per_site <- function(case_mat, ctrl_mat, direction) {
  n_sites <- ncol(case_mat)
  stopifnot(ncol(ctrl_mat) == n_sites, length(direction) == n_sites)
  aucs <- numeric(n_sites)
  for (j in seq_len(n_sites)) {
    d <- direction[j]
    cv <- d * case_mat[, j]; cv <- cv[!is.na(cv)]
    kv <- d * ctrl_mat[, j]; kv <- kv[!is.na(kv)]
    if (length(cv) == 0L || length(kv) == 0L) { aucs[j] <- NA_real_; next }
    r <- rank(c(cv, kv))
    U <- sum(r[seq_along(cv)]) - length(cv) * (length(cv) + 1L) / 2
    aucs[j] <- U / (length(cv) * length(kv))
  }
  aucs
}

g_auc_combine <- function(auc_gel, auc_nava) {
  a1 <- auc_gel - 0.5; a2 <- auc_nava - 0.5
  ok <- !is.na(a1) & !is.na(a2) & a1 > 0 & a2 > 0
  res <- rep(NA_real_, length(auc_gel))
  res[ok] <- 0.5 + sqrt(a1[ok] * a2[ok])
  res
}

select_sites_by_g_auc <- function(g_auc_vec, gene_sym_vec, K, max_per_gene = MAX_PER_GENE) {
  eligible <- which(!is.na(g_auc_vec))
  if (length(eligible) == 0L) return(integer(0))
  ordered_idx <- eligible[order(-g_auc_vec[eligible])]
  chosen <- integer(0); gene_count <- integer(0)
  for (idx in ordered_idx) {
    g <- gene_sym_vec[idx]; cnt <- gene_count[g]; if (is.na(cnt)) cnt <- 0L
    if (cnt < max_per_gene) {
      chosen <- c(chosen, idx); gene_count[g] <- cnt + 1L
      if (length(chosen) >= K) break
    }
  }
  chosen
}

cat("\nApplying MAX_PER_GENE=1 (G_AUC on insT training cases vs controls)...\n")
gel_ctrl_mat  <- gel_lp_cand[gel_ctrl_ids,  , drop = FALSE]
nava_ctrl_mat <- nava_lp_cand[nava_ctrl_ids, , drop = FALSE]

auc_gel  <- auc_per_site(gel_lp_cand[gel_insT_train_ids,  , drop = FALSE],
                         gel_ctrl_mat,  a5ss_cand$direction)
auc_nava <- auc_per_site(nava_lp_cand[nava_insT_train_ids, , drop = FALSE],
                         nava_ctrl_mat, a5ss_cand$direction)
g_auc    <- g_auc_combine(auc_gel, auc_nava)

keep_idx <- select_sites_by_g_auc(g_auc, a5ss_cand$GeneSymbol, K = nrow(a5ss_cand))
stopifnot(length(keep_idx) > 0L)

a5ss_cand$g_auc    <- g_auc
a5ss_cand$auc_gel  <- auc_gel
a5ss_cand$auc_nava <- auc_nava
a5ss_cand_pre      <- a5ss_cand
cat(sprintf("  Unique genes in depth pool (pre-AUC-filter): %d\n",
            length(unique(a5ss_cand_pre$GeneSymbol))))
cat(sprintf("  Sites dropped for non-positive AUC in either cohort: %d\n",
            sum(is.na(g_auc))))
a5ss_cand          <- a5ss_cand[keep_idx, , drop = FALSE]
rownames(a5ss_cand) <- NULL

gel_lp_cand     <- gel_lp_cand[,  a5ss_cand$gel_ev_id, drop = FALSE]
nava_lp_cand    <- nava_lp_cand[, a5ss_cand$gel_ev_id, drop = FALSE]
gel_ctrl_means  <- gel_ctrl_means[a5ss_cand$gel_ev_id]
nava_ctrl_means <- nava_ctrl_means[a5ss_cand$gel_ev_id]
direction       <- a5ss_cand$direction

cat(sprintf("  Sites after MAX_PER_GENE=1: %d  (down from %d)\n",
            nrow(a5ss_cand), nrow(a5ss_cand_pre)))

# --------------------------------------------------------------------------
# Direction-oriented, cohort-centered feature matrices
# --------------------------------------------------------------------------
orient_gel  <- sweep(gel_lp_cand,  2, gel_ctrl_means,  "-")
orient_nava <- sweep(nava_lp_cand, 2, nava_ctrl_means, "-")
orient_gel  <- sweep(orient_gel,   2, direction, "*")
orient_nava <- sweep(orient_nava,  2, direction, "*")

# --------------------------------------------------------------------------
# Training matrix: insT training cases + all controls (holdout excluded)
# --------------------------------------------------------------------------
train_rows_gel  <- c(gel_insT_train_ids,  gel_ctrl_ids)
train_rows_nava <- c(nava_insT_train_ids, nava_ctrl_ids)

X <- rbind(orient_gel[train_rows_gel,  , drop = FALSE],
           orient_nava[train_rows_nava, , drop = FALSE])
y <- c(rep(1L, length(gel_insT_train_ids)),  rep(0L, length(gel_ctrl_ids)),
       rep(1L, length(nava_insT_train_ids)), rep(0L, length(nava_ctrl_ids)))
sample_ids_train <- c(train_rows_gel, train_rows_nava)
cohort_train     <- c(rep("GEL",  length(train_rows_gel)),
                      rep("Nava", length(train_rows_nava)))
rownames(X)      <- sample_ids_train
stopifnot(nrow(X) == length(y), !anyNA(X))

n_cases  <- sum(y == 1L)
n_ctrls  <- sum(y == 0L)
case_idx <- which(y == 1L)
ctrl_idx <- which(y == 0L)
foldid   <- integer(length(y))
foldid[case_idx] <- seq_along(case_idx)
foldid[ctrl_idx] <- ((seq_along(ctrl_idx) - 1L) %% n_cases) + 1L
stopifnot(all(foldid > 0L), max(foldid) == n_cases)

cat(sprintf("\nTraining matrix: %d samples x %d sites\n", nrow(X), ncol(X)))
cat(sprintf("  insT train cases=%d  controls=%d  folds=%d (case-LOO)\n",
            n_cases, n_ctrls, n_cases))
cat(sprintf("  insT holdout: GEL=%d  Nava=%d\n",
            length(gel_insT_ho_ids), length(nava_insT_ho_ids)))

# --------------------------------------------------------------------------
# cv.glmnet
# --------------------------------------------------------------------------
cat("\nFitting cv.glmnet (alpha=1, binomial, lower.limits=0, keep=TRUE)...\n")
set.seed(2026L)
cv_fit <- cv.glmnet(
  x            = X,
  y            = y,
  family       = "binomial",
  alpha        = ALPHA_LASSO,
  foldid       = foldid,
  standardize  = TRUE,
  lower.limits = 0,
  type.measure = "deviance",
  keep         = TRUE
)

lambda_min <- cv_fit$lambda.min
lambda_1se <- cv_fit$lambda.1se
cat(sprintf("  lambda.min=%.5g  lambda.1se=%.5g\n", lambda_min, lambda_1se))

extract_fit <- function(cv_fit, s) {
  co <- coef(cv_fit, s = s)
  cv <- as.numeric(co); names(cv) <- rownames(co)
  intercept <- cv[1]; beta <- cv[-1]; nz <- which(beta != 0)
  list(lambda = if (s == "lambda.1se") cv_fit$lambda.1se else cv_fit$lambda.min,
       intercept = intercept, beta = beta,
       n_selected = length(nz), event_ids = names(beta)[nz], coefs = beta[nz])
}
fit_at_1se <- extract_fit(cv_fit, "lambda.1se")
fit_at_min <- extract_fit(cv_fit, "lambda.min")
cat(sprintf("  n_selected @ lambda.1se = %d   @ lambda.min = %d\n",
            fit_at_1se$n_selected, fit_at_min$n_selected))

chosen_label <- LAMBDA_RULE
if (fit_at_1se$n_selected == 0L) {
  warning("lambda.1se selected 0 sites; falling back to lambda.min")
  chosen_label <- "lambda.min"
}
chosen_fit <- if (chosen_label == "lambda.1se") fit_at_1se else fit_at_min

# --------------------------------------------------------------------------
# Selected sites
# --------------------------------------------------------------------------
beta_1se     <- fit_at_1se$beta
beta_min     <- fit_at_min$beta
selected_idx <- which((if (chosen_label == "lambda.1se") beta_1se else beta_min) != 0)
K_selected   <- length(selected_idx)
stopifnot("No sites selected" = K_selected > 0L)
cat(sprintf("K = %d selected sites\n", K_selected))

# Selected sites table (union 1se + min for display)
nz_union <- which(beta_1se != 0 | beta_min != 0)
selected_sites <- data.frame(
  gel_ev_id  = a5ss_cand$gel_ev_id[nz_union],
  GeneSymbol = a5ss_cand$GeneSymbol[nz_union],
  direction  = a5ss_cand$direction[nz_union],
  g_auc      = a5ss_cand$g_auc[nz_union],
  coef_1se   = unname(beta_1se[nz_union]),
  coef_min   = unname(beta_min[nz_union]),
  stringsAsFactors = FALSE
)
selected_sites <- selected_sites[order(-selected_sites$coef_1se,
                                       -selected_sites$coef_min), , drop = FALSE]
cat("\nSelected sites (union @1se/@min):\n")
print(selected_sites, row.names = FALSE)

# --------------------------------------------------------------------------
# LOO prevalidated linpred
# --------------------------------------------------------------------------
lambda_seq    <- cv_fit$lambda
chosen_lambda <- chosen_fit$lambda
lcol          <- which.min(abs(lambda_seq - chosen_lambda))
preval_lp     <- as.numeric(cv_fit$fit.preval[, lcol])
preval_lp_1se <- as.numeric(cv_fit$fit.preval[, which.min(abs(lambda_seq - lambda_1se))])
preval_lp_min <- as.numeric(cv_fit$fit.preval[, which.min(abs(lambda_seq - lambda_min))])

# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------
auc_two_sample <- function(case_scores, ctrl_scores) {
  case_scores <- case_scores[!is.na(case_scores)]
  ctrl_scores <- ctrl_scores[!is.na(ctrl_scores)]
  if (length(case_scores) == 0L || length(ctrl_scores) == 0L) return(NA_real_)
  r <- rank(c(case_scores, ctrl_scores))
  U <- sum(r[seq_along(case_scores)]) - length(case_scores) * (length(case_scores) + 1L) / 2
  U / (length(case_scores) * length(ctrl_scores))
}

boot_auc_ci <- function(case_scores, ctrl_scores, B = B_BOOT, seed = SEED) {
  set.seed(seed)
  vals <- replicate(B, auc_two_sample(sample(case_scores, replace = TRUE),
                                      sample(ctrl_scores, replace = TRUE)))
  c(point = auc_two_sample(case_scores, ctrl_scores),
    lo    = unname(quantile(vals, 0.025, na.rm = TRUE)),
    hi    = unname(quantile(vals, 0.975, na.rm = TRUE)))
}

is_case_train <- y == 1L
is_ctrl_train <- y == 0L
is_gel_train  <- cohort_train == "GEL"
is_nava_train <- cohort_train == "Nava"

# LOO linpred AUC with bootstrap CIs
auc_gel_loo  <- boot_auc_ci(preval_lp[is_case_train & is_gel_train],
                             preval_lp[is_ctrl_train & is_gel_train])
auc_nava_loo <- boot_auc_ci(preval_lp[is_case_train & is_nava_train],
                             preval_lp[is_ctrl_train & is_nava_train])
auc_comb_loo <- boot_auc_ci(preval_lp[is_case_train], preval_lp[is_ctrl_train])

cat("\n=== Metrics (insT, LASSO LOO linpred) ===\n")
cat(sprintf("  LOO linpred AUC: GEL=%.3f [%.3f-%.3f]  Nava=%.3f [%.3f-%.3f]  Combined=%.3f [%.3f-%.3f]  (bootstrap 95%% CI, B=%d)\n",
            auc_gel_loo["point"],  auc_gel_loo["lo"],  auc_gel_loo["hi"],
            auc_nava_loo["point"], auc_nava_loo["lo"], auc_nava_loo["hi"],
            auc_comb_loo["point"], auc_comb_loo["lo"], auc_comb_loo["hi"],
            B_BOOT))

# --------------------------------------------------------------------------
# LOO-vs-full-fit misrank diagnostic
# The final trained model separates train cases from train controls
# perfectly (AUC=1) because every case's linpred is above every control's.
# CV-AUC<1 means at least one case, when scored by the LASSO refit
# WITHOUT that case (its own fold), drops into the control-score range.
# Report per case: (a) linpred at the full fit, (b) LOO linpred, and
# (c) how many controls in the same cohort outscore it under each --
# any case with n_ctrl_above_loo > 0 is a rank-flipper contributing to
# CV-AUC < 1.
# --------------------------------------------------------------------------
linpred_full_train <- as.numeric(X %*% beta_1se + fit_at_1se$intercept)
diag_df <- data.frame(
  sample_id   = sample_ids_train,
  cohort      = cohort_train,
  is_case     = y == 1L,
  linpred_full = linpred_full_train,
  linpred_loo  = preval_lp_1se,
  stringsAsFactors = FALSE)

count_ctrls_above <- function(sub_df, scores) {
  case_scores <- scores[sub_df$is_case]
  ctrl_scores <- scores[!sub_df$is_case]
  vapply(seq_along(case_scores), function(i) sum(ctrl_scores > case_scores[i]),
         integer(1L))
}

cat("\n=== LOO vs full-fit case-rank diagnostic (linpred @ lambda.1se) ===\n")
for (co in c("GEL", "Nava")) {
  sub <- diag_df[diag_df$cohort == co, , drop = FALSE]
  case_rows <- which(sub$is_case)
  if (length(case_rows) == 0L) next
  n_full <- count_ctrls_above(sub, sub$linpred_full)
  n_loo  <- count_ctrls_above(sub, sub$linpred_loo)
  case_tbl <- data.frame(
    sample_id           = sub$sample_id[case_rows],
    linpred_full        = round(sub$linpred_full[case_rows], 4),
    linpred_loo         = round(sub$linpred_loo[case_rows],  4),
    ctrls_above_full    = n_full,
    ctrls_above_loo     = n_loo,
    stringsAsFactors    = FALSE)
  case_tbl <- case_tbl[order(-case_tbl$ctrls_above_loo,
                             case_tbl$linpred_loo), , drop = FALSE]
  n_ctrl <- sum(!sub$is_case)
  cat(sprintf("--- %s (%d cases vs %d controls) ---\n",
              co, length(case_rows), n_ctrl))
  print(case_tbl, row.names = FALSE)
  cat(sprintf("  Cases with >=1 control outscoring them at full fit : %d/%d\n",
              sum(n_full > 0L), length(n_full)))
  cat(sprintf("  Cases with >=1 control outscoring them at LOO      : %d/%d  (drives CV-AUC < 1)\n",
              sum(n_loo > 0L), length(n_loo)))
}

# --------------------------------------------------------------------------
# Lowest final-fit P(case) case (relative to all training cases). Also
# per-case sigmoid probability so it's directly comparable to Panel B.
# --------------------------------------------------------------------------
prob_full_train <- 1 / (1 + exp(-linpred_full_train))
case_prob_full <- prob_full_train[y == 1L]
case_ids_full  <- sample_ids_train[y == 1L]
min_idx <- which.min(case_prob_full)
cat(sprintf(
  "\nLowest final-fit P(case) among training cases: %s  P=%.4f  linpred=%.4f\n",
  case_ids_full[min_idx], case_prob_full[min_idx],
  linpred_full_train[y == 1L][min_idx]))

# --------------------------------------------------------------------------
# Per-fold LASSO refit: identify which sites depend on which case
# For each training case i, refit glmnet on all samples EXCEPT case i (all
# controls kept), at the same lambda.1se. Compare selected sites and
# coefficients to the full-train fit. Sites that appear only in the full
# fit and disappear in fold-i-refit are "picked because of case i"; sites
# that appear only in fold-i-refit are ones the full fit dropped when
# case i was included. Prints a wide table (rows = union of selected sites,
# cols = full + one per fold) with coefficients at lambda.1se.
# --------------------------------------------------------------------------
cat("\n=== Per-fold LASSO coefficients (drop one case, refit @ lambda.1se) ===\n")
sel_full <- which(beta_1se != 0)
site_labels <- if (exists("a5ss_cand") && "GeneSymbol" %in% names(a5ss_cand))
  a5ss_cand$GeneSymbol else colnames(X)

fold_beta_list <- vector("list", n_cases)
for (i in seq_len(n_cases)) {
  case_row <- case_idx[i]
  X_i <- X[-case_row, , drop = FALSE]
  y_i <- y[-case_row]
  fit_i <- glmnet(x = X_i, y = y_i, family = "binomial",
                  alpha = ALPHA_LASSO, standardize = TRUE,
                  lower.limits = 0)
  co_i  <- as.numeric(coef(fit_i, s = lambda_1se))[-1]
  fold_beta_list[[i]] <- co_i
}
fold_beta_mat <- do.call(cbind, fold_beta_list)
rownames(fold_beta_mat) <- site_labels
colnames(fold_beta_mat) <- sample_ids_train[case_idx]

# Union of ever-selected sites across full + all folds.
sel_any <- which(beta_1se != 0 | rowSums(fold_beta_mat != 0) > 0L)
comp <- cbind(full = beta_1se[sel_any], fold_beta_mat[sel_any, , drop = FALSE])
comp_r <- round(comp, 3)

# For each fold, list sites in full-fit that dropped out when that case
# was removed (= sites the LASSO picked BECAUSE of that case).
cat("Sites dropped from the signature when a specific case is removed:\n")
for (i in seq_len(n_cases)) {
  dropped <- which(beta_1se != 0 & fold_beta_list[[i]] == 0)
  added   <- which(beta_1se == 0 & fold_beta_list[[i]] != 0)
  if (length(dropped) == 0L && length(added) == 0L) next
  cat(sprintf("  Drop %-18s: dropped=%s ; added=%s\n",
              sample_ids_train[case_idx[i]],
              if (length(dropped)) paste(site_labels[dropped], collapse = ",") else "(none)",
              if (length(added))   paste(site_labels[added],   collapse = ",") else "(none)"))
}

cat("\nCoefficient matrix (rows = ever-selected sites; cols = full-fit + one per LOO fold):\n")
print(comp_r)

# --------------------------------------------------------------------------
# Per-case meta
# --------------------------------------------------------------------------
all_insT_ids   <- c(gel_insT_train_ids, gel_insT_ho_ids,
                    nava_insT_train_ids, nava_insT_ho_ids)
all_insT_set   <- c(rep("train",   length(gel_insT_train_ids)),
                    rep("holdout", length(gel_insT_ho_ids)),
                    rep("train",   length(nava_insT_train_ids)),
                    rep("holdout", length(nava_insT_ho_ids)))
all_insT_co    <- c(rep("GEL",  length(gel_insT_train_ids)  + length(gel_insT_ho_ids)),
                    rep("Nava", length(nava_insT_train_ids) + length(nava_insT_ho_ids)))
all_insT_hgvs  <- unname(all_hgvs_lookup[all_insT_ids])

predict_lp <- function(X_new, fit) as.numeric(X_new %*% fit$beta + fit$intercept)
sigmoid    <- function(z) 1 / (1 + exp(-z))

per_case_meta <- data.frame(
  sample_id      = all_insT_ids,
  cohort         = all_insT_co,
  set            = all_insT_set,
  hgvs           = all_insT_hgvs,
  stringsAsFactors = FALSE
)
per_case_meta$linpred_full <- NA_real_
per_case_meta$prob_full    <- NA_real_
per_case_meta$linpred_loo  <- NA_real_
for (i in seq_len(nrow(per_case_meta))) {
  sid <- per_case_meta$sample_id[i]
  co  <- per_case_meta$cohort[i]
  X_i <- if (co == "GEL") orient_gel[sid, , drop = FALSE] else orient_nava[sid, , drop = FALSE]
  per_case_meta$linpred_full[i] <- predict_lp(X_i, fit_at_1se)
  per_case_meta$prob_full[i]    <- sigmoid(per_case_meta$linpred_full[i])
  match_train <- match(sid, sample_ids_train)
  if (!is.na(match_train)) per_case_meta$linpred_loo[i] <- preval_lp[match_train]
}
per_case_meta <- per_case_meta[order(-per_case_meta$prob_full), , drop = FALSE]
cat(sprintf("\nPer-case P(insT) (n=%d, ordered by prob_full):\n", nrow(per_case_meta)))
print(data.frame(
  sample_id = per_case_meta$sample_id, cohort = per_case_meta$cohort,
  hgvs = per_case_meta$hgvs, prob = round(per_case_meta$prob_full, 3),
  row.names = NULL
))

# --------------------------------------------------------------------------
# Save RDS
# --------------------------------------------------------------------------
# Sample sizes: depth-filter scope = all insT cases (train + holdout) +
# controls; fit_scope (G_AUC ranking for best-event-per-gene, and the LASSO
# training matrix) = train cases only (holdout excluded) + controls. Both
# saved since figure_s2.R needs to annotate branch sample size and the two
# scopes genuinely differ here (unlike stem-III, which has no holdout).
sample_sizes <- list(
  depth_scope = list(
    gel_cases  = length(gel_insT_ids),       gel_ctrl  = N_GEL_CONTROLS,
    nava_cases = length(nava_insT_ids),      nava_ctrl = N_NAVA_CONTROLS
  ),
  fit_scope = list(
    gel_cases  = length(gel_insT_train_ids),  gel_ctrl  = N_GEL_CONTROLS,
    nava_cases = length(nava_insT_train_ids), nava_ctrl = N_NAVA_CONTROLS
  )
)

out <- list(
  candidate_sites_pre  = a5ss_cand_pre,
  sample_sizes         = sample_sizes,
  a5ss_cand            = a5ss_cand,
  direction            = direction,
  gel_ctrl_means       = gel_ctrl_means,
  nava_ctrl_means      = nava_ctrl_means,
  X_train              = X,
  y_train              = y,
  cohort_train         = cohort_train,
  sample_ids_train     = sample_ids_train,
  foldid               = foldid,
  cv_fit               = cv_fit,
  lambda_min           = lambda_min,
  lambda_1se           = lambda_1se,
  chosen_lambda_label  = chosen_label,
  selected_sites       = selected_sites,
  selected_idx         = selected_idx,
  K_selected           = K_selected,
  fit_at_1se           = fit_at_1se,
  fit_at_min           = fit_at_min,
  preval_lp_1se        = setNames(preval_lp_1se, sample_ids_train),
  preval_lp_min        = setNames(preval_lp_min, sample_ids_train),
  ho_scores = list(
    gel  = data.frame(sample_id = gel_insT_ho_ids,  stringsAsFactors = FALSE),
    nava = data.frame(sample_id = nava_insT_ho_ids, stringsAsFactors = FALSE)
  ),
  metrics = list(
    auc_gel_loo = auc_gel_loo, auc_nava_loo = auc_nava_loo, auc_comb_loo = auc_comb_loo
  ),
  per_case_meta = per_case_meta,
  config = list(
    LAMBDA_RULE = LAMBDA_RULE, ALPHA = ALPHA_LASSO, FAMILY = "binomial",
    LOWER_LIMITS = 0, STANDARDIZE = TRUE, FOLDS = "case-LOO",
    PSI_SCALE = "logit_bb", CENTER = "cohort_ctrl_mean",
    DIRECTION_ORIENTED = TRUE, MAX_PER_GENE = MAX_PER_GENE,
    EVENT_TYPE = EVENT_TYPE, FDR_THRESH = FDR_THRESH,
    MIN_DEPTH = MIN_DEPTH, INST_HGVS = INST_HGVS,
    CANDIDATE_POOL = "483_concordant_A5SS",
    SPLIT_FILE = SPLIT_FILE,
    NOTE = "Exploratory insT-only model; split from v4 split file"
  )
)
saveRDS(out, OUT_RDS)
cat("\nResults saved to:", OUT_RDS, "\n")
