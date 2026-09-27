#!/usr/bin/env Rscript
#
# loo_evaluate_stemIII_lasso.R
# Exploratory LASSO trained on stem-III variant carriers only vs controls.
#
# Candidate space: A5SS events from the 483 sig-in-both concordant pool
# (FDR<0.05 in both GEL and Nava all-cases rMATS, direction-concordant).
# Within that pool, A5SS events are filtered by min(IJC+SJC) >= 10 across
# stem-III cases and controls in both cohorts, non-NA, and MAX_PER_GENE=1
# (best G_AUC site per gene, AUC computed on stem-III cases vs controls).
#
# No holdout split: all stem-III carriers enter the LASSO. Performance is
# reported via case-LOO cross-validated linpred only.
#
# Feature: direction * (bb-logit-PSI - per-cohort-ctrl-mean), per v4 formula.
# LASSO: binomial, alpha=1, lower.limits=0, foldid=case-LOO.
#
# Stem-III definition: RNU4-2 nucleotide positions 75-79 (consistent with
# figure1.R STEM_III_START/END). Carriers identified by parsing the first
# integer from each HGVS string (same logic as the panel-A lollipop).
#
# Inputs:
#   outputs/psi_matrix.rds
#   outputs/nava_psi_matrix.rds
#   data/sample_genotype.csv
#   outputs/nava_severity_labels.csv
#
# Outputs:
#   outputs/loo_results_stemIII_lasso.rds
#
# Run: pixi run Rscript scripts/gel/loo_evaluate_stemIII_lasso.R

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
GENO_FILE     <- "data/sample_genotype.csv"
NAVA_SEV_FILE <- "outputs/nava_severity_labels.csv"

OUT_RDS       <- "outputs/loo_results_stemIII_lasso.rds"

N_GEL_CASES     <- 11L
N_GEL_CONTROLS  <- 33L
N_NAVA_CASES    <- 19L
N_NAVA_CONTROLS <- 21L

# Stem-III region: RNU4-2 nucleotide positions 75-79 inclusive.
# Consistent with STEM_III_START/END in figure1.R and figure_s1.R.
STEMIII_POS_MIN <- 75L
STEMIII_POS_MAX <- 79L

EVENT_TYPE   <- "A5SS"
FDR_THRESH   <- 0.05
MIN_DEPTH    <- 10L
MAX_PER_GENE <- 1L

LAMBDA_RULE  <- "lambda.1se"
ALPHA_LASSO  <- 1

stopifnot(
  "GEL PSI matrix not found"     = file.exists(PSI_PATH),
  "Nava PSI matrix not found"    = file.exists(NAVA_PSI_PATH),
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
# Load genotype data; identify stem-III carriers
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

# Identify stem-III carriers by parsing the first integer from the HGVS
# (same logic as figure1.R panel-A lollipop). Variants at positions 75-79
# are classified as stem-III.
hgvs_pos <- function(hgvs_vec) {
  as.integer(regmatches(hgvs_vec, regexpr("[0-9]+", hgvs_vec)))
}
is_stemIII <- function(hgvs_vec) {
  pos <- hgvs_pos(hgvs_vec)
  !is.na(pos) & pos >= STEMIII_POS_MIN & pos <= STEMIII_POS_MAX
}

gel_stemIII_ids  <- names(all_hgvs_lookup)[
  names(all_hgvs_lookup) %in% gel_case_ids  & is_stemIII(all_hgvs_lookup)]
nava_stemIII_ids <- names(all_hgvs_lookup)[
  names(all_hgvs_lookup) %in% nava_case_ids & is_stemIII(all_hgvs_lookup)]

# All carriers enter training (no holdout split).
cat(sprintf("GEL  stem-III carriers: %d (of %d cases)\n",
            length(gel_stemIII_ids),  N_GEL_CASES))
cat(sprintf("Nava stem-III carriers: %d (of %d cases)\n",
            length(nava_stemIII_ids), N_NAVA_CASES))

stopifnot(
  "No GEL stem-III carriers found"      = length(gel_stemIII_ids) > 0L,
  "GEL stem-III IDs not in psi matrix"  = all(gel_stemIII_ids  %in% rownames(gel_psi))
)
if (length(nava_stemIII_ids) > 0L) {
  stopifnot("Nava stem-III IDs not in psi matrix" = all(nava_stemIII_ids %in% rownames(nava_psi)))
} else {
  cat("No Nava stem-III carriers; analysis is GEL-only.\n")
}

# --------------------------------------------------------------------------
# Identify candidate space: A5SS events from the 483 sig-in-both pool
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
# cohort) -- same definition as figure1.R panel D / figure_s1.R and the
# insT LOO. Matches the transform used for LASSO features downstream, so
# orientation and candidate pool share a metric.
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
  IncLevelDiff   = gel_m$IncLevelDiff[concordant],  # kept for reference
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
# bb-logit PSI (per-event bb_logit already defined above)
# --------------------------------------------------------------------------

gel_lp_a5ss  <- bb_logit(gel_ijc[,  a5ss_df$gel_ev_id,  drop = FALSE],
                          gel_sjc[,  a5ss_df$gel_ev_id,  drop = FALSE])
nava_lp_a5ss <- bb_logit(nava_ijc[, a5ss_df$nava_ev_id, drop = FALSE],
                          nava_sjc[, a5ss_df$nava_ev_id, drop = FALSE])
colnames(nava_lp_a5ss) <- a5ss_df$gel_ev_id

stopifnot(identical(colnames(gel_lp_a5ss), colnames(nava_lp_a5ss)))

# --------------------------------------------------------------------------
# Depth filter: min(IJC+SJC) >= MIN_DEPTH across ALL stem-III cases and controls
# --------------------------------------------------------------------------
gel_all_ids  <- c(gel_stemIII_ids,  gel_ctrl_ids)
nava_all_ids <- c(nava_stemIII_ids, nava_ctrl_ids)  # may be ctrl-only if no Nava cases

gel_depth  <- gel_ijc[gel_all_ids,   a5ss_df$gel_ev_id,  drop = FALSE] +
              gel_sjc[gel_all_ids,   a5ss_df$gel_ev_id,  drop = FALSE]
nava_depth <- nava_ijc[nava_all_ids, a5ss_df$nava_ev_id, drop = FALSE] +
              nava_sjc[nava_all_ids, a5ss_df$nava_ev_id, drop = FALSE]

min_depth_gel  <- apply(gel_depth,  2, min, na.rm = TRUE)
min_depth_nava <- apply(nava_depth, 2, min, na.rm = TRUE)
depth_pass <- is.finite(min_depth_gel)  & is.finite(min_depth_nava) &
              min_depth_gel >= MIN_DEPTH & min_depth_nava >= MIN_DEPTH

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
# G_AUC on stem-III cases for MAX_PER_GENE=1 pre-filter
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

cat("\nApplying MAX_PER_GENE=1 (G_AUC on stem-III cases vs controls)...\n")
gel_ctrl_mat  <- gel_lp_cand[gel_ctrl_ids,  , drop = FALSE]
nava_ctrl_mat <- nava_lp_cand[nava_ctrl_ids, , drop = FALSE]

auc_gel <- auc_per_site(gel_lp_cand[gel_stemIII_ids, , drop = FALSE],
                        gel_ctrl_mat, a5ss_cand$direction)

# If Nava has stem-III carriers, compute G_AUC from both; else GEL G_AUC only.
if (length(nava_stemIII_ids) > 0L) {
  auc_nava <- auc_per_site(nava_lp_cand[nava_stemIII_ids, , drop = FALSE],
                           nava_ctrl_mat, a5ss_cand$direction)
  g_auc <- g_auc_combine(auc_gel, auc_nava)
} else {
  auc_nava <- rep(NA_real_, nrow(a5ss_cand))
  g_auc    <- auc_gel  # fall back to GEL AUC alone
  cat("  No Nava stem-III carriers; G_AUC = GEL AUC only.\n")
}

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
# Training matrix: ALL stem-III cases + all controls (no holdout)
# --------------------------------------------------------------------------
train_rows_gel  <- c(gel_stemIII_ids, gel_ctrl_ids)
train_rows_nava <- if (length(nava_stemIII_ids) > 0L) {
  c(nava_stemIII_ids, nava_ctrl_ids)
} else {
  nava_ctrl_ids
}

X <- rbind(orient_gel[train_rows_gel,  , drop = FALSE],
           orient_nava[train_rows_nava, , drop = FALSE])
y <- c(rep(1L, length(gel_stemIII_ids)),  rep(0L, length(gel_ctrl_ids)),
       rep(1L, length(nava_stemIII_ids)), rep(0L, length(nava_ctrl_ids)))
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
cat(sprintf("  stem-III cases=%d  controls=%d  folds=%d (case-LOO)\n",
            n_cases, n_ctrls, n_cases))
cat(sprintf("  Variant breakdown:\n"))
stemIII_hgvs_in_data <- unname(all_hgvs_lookup[c(gel_stemIII_ids, nava_stemIII_ids)])
print(table(stemIII_hgvs_in_data))

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
# Metrics (CV-only; no holdout)
# --------------------------------------------------------------------------
auc_two_sample <- function(case_scores, ctrl_scores) {
  case_scores <- case_scores[!is.na(case_scores)]
  ctrl_scores <- ctrl_scores[!is.na(ctrl_scores)]
  if (length(case_scores) == 0L || length(ctrl_scores) == 0L) return(NA_real_)
  r <- rank(c(case_scores, ctrl_scores))
  U <- sum(r[seq_along(case_scores)]) - length(case_scores) * (length(case_scores) + 1L) / 2
  U / (length(case_scores) * length(ctrl_scores))
}

is_case_train <- y == 1L
is_ctrl_train <- y == 0L
is_gel_train  <- cohort_train == "GEL"
is_nava_train <- cohort_train == "Nava"

# LOO linpred AUC (primary performance metric for this analysis)
auc_gel_loo  <- auc_two_sample(preval_lp[is_case_train & is_gel_train],
                               preval_lp[is_ctrl_train & is_gel_train])
auc_nava_loo <- if (any(is_case_train & is_nava_train))
  auc_two_sample(preval_lp[is_case_train & is_nava_train],
                 preval_lp[is_ctrl_train & is_nava_train]) else NA_real_

cat("\n=== Metrics (stem-III, LASSO LOO linpred) ===\n")
cat("  NOTE: LOO linpred AUC is primary -- no holdout; all carriers in training.\n")
cat(sprintf("  LOO linpred AUC : GEL=%.3f  Nava=%.3f\n", auc_gel_loo, auc_nava_loo))

# --------------------------------------------------------------------------
# Per-case meta
# --------------------------------------------------------------------------
all_stemIII_ids  <- c(gel_stemIII_ids, nava_stemIII_ids)
all_stemIII_co   <- c(rep("GEL",  length(gel_stemIII_ids)),
                      rep("Nava", length(nava_stemIII_ids)))
all_stemIII_hgvs <- unname(all_hgvs_lookup[all_stemIII_ids])

predict_lp <- function(X_new, fit) as.numeric(X_new %*% fit$beta + fit$intercept)
sigmoid    <- function(z) 1 / (1 + exp(-z))

per_case_meta <- data.frame(
  sample_id      = all_stemIII_ids,
  cohort         = all_stemIII_co,
  hgvs           = all_stemIII_hgvs,
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
per_case_meta$prob_loo <- sigmoid(per_case_meta$linpred_loo)
per_case_meta <- per_case_meta[order(-per_case_meta$prob_full), , drop = FALSE]
cat(sprintf("\nPer-case P(stem-III) (n=%d, ordered by prob_full):\n", nrow(per_case_meta)))
print(data.frame(
  sample_id = per_case_meta$sample_id, cohort = per_case_meta$cohort,
  hgvs      = per_case_meta$hgvs,
  prob_full = round(per_case_meta$prob_full, 3),
  prob_loo  = round(per_case_meta$prob_loo,  3),
  row.names = NULL
))

# --------------------------------------------------------------------------
# Save RDS
# --------------------------------------------------------------------------
# Sample sizes: no holdout split here (all stem-III carriers enter
# training), so depth_scope == fit_scope -- both saved anyway for a
# consistent schema with loo_evaluate_insT_lasso.R (which does differ).
sample_sizes <- list(
  depth_scope = list(
    gel_cases  = length(gel_stemIII_ids),  gel_ctrl  = N_GEL_CONTROLS,
    nava_cases = length(nava_stemIII_ids), nava_ctrl = N_NAVA_CONTROLS
  ),
  fit_scope = list(
    gel_cases  = length(gel_stemIII_ids),  gel_ctrl  = N_GEL_CONTROLS,
    nava_cases = length(nava_stemIII_ids), nava_ctrl = N_NAVA_CONTROLS
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
  metrics = list(
    auc_gel_loo = auc_gel_loo, auc_nava_loo = auc_nava_loo
  ),
  per_case_meta = per_case_meta,
  config = list(
    LAMBDA_RULE = LAMBDA_RULE, ALPHA = ALPHA_LASSO, FAMILY = "binomial",
    LOWER_LIMITS = 0, STANDARDIZE = TRUE, FOLDS = "case-LOO",
    PSI_SCALE = "logit_bb", CENTER = "cohort_ctrl_mean",
    DIRECTION_ORIENTED = TRUE, MAX_PER_GENE = MAX_PER_GENE,
    EVENT_TYPE = EVENT_TYPE, FDR_THRESH = FDR_THRESH,
    MIN_DEPTH = MIN_DEPTH,
    STEMIII_POS_MIN = STEMIII_POS_MIN, STEMIII_POS_MAX = STEMIII_POS_MAX,
    CANDIDATE_POOL = "483_concordant_A5SS",
    NOTE = "Exploratory stem-III-only model; no holdout; all carriers in training"
  )
)
saveRDS(out, OUT_RDS)
cat("\nResults saved to:", OUT_RDS, "\n")
