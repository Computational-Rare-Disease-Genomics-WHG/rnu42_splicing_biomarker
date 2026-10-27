#!/usr/bin/env Rscript
#
# strict_featsel_penalty_diagnostic.R
# Why does the strict-featsel insT LASSO select a different number of events
# from the canonical fit? Compares the two saved fits (read only):
#   1. lambda.1se / lambda.min of each fit, and K for each fit evaluated at
#      the other fit's lambda.1se.
#   2. Order in which events enter each LASSO path as lambda decreases.
#   3. Whether the strict-only events were in the canonical candidate pool.
#   4. Selection frequency: refit at each fit's own lambda.1se with each
#      training case left out in turn (same case-LOO folds as cv.glmnet).
#   5. CV deviance curves for both fits -> PDF.
#   6. Why the CV error bars differ: per-fold deviance at each fit's
#      lambda.min (which left-out case drives the spread), and how the
#      strict-featsel held-out controls scored in the canonical CV.
#   7. The hardest GEL case (the left-out case in the canonical fit's worst
#      GEL fold, from section 6): sample-level depth vs other GEL samples, and
#      per selected event its depth, feature value and contribution, next to
#      the other GEL insT training cases and the GEL training controls.
#
# Inputs:
#   outputs/loo_results_insT_lasso.rds
#   outputs/strict_featsel/loo_results_insT_lasso_strict.rds
#   data/split_strict_featsel.csv
#   outputs/psi_matrix.rds                    (depth only)
# Outputs (GEL protected - gene names):
#   outputs/strict_featsel/penalty_diagnostic_cv_curves.pdf
#   stdout
#
# Run: pixi run Rscript scripts/gel/strict_featsel_penalty_diagnostic.R

# ---- libpath -----------------------------------------------------------------
invisible(local({
  user_lib <- path.expand("~/re_gecip/machine_learning/R_packages")
  if (nzchar(Sys.getenv("PIXI_ENVIRONMENT_NAME"))) {
    keep <- !normalizePath(.libPaths(), mustWork = FALSE) %in%
            normalizePath(user_lib,    mustWork = FALSE)
    if (any(!keep)) .libPaths(.libPaths()[keep])
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
suppressPackageStartupMessages(library(glmnet))

CANON_RDS  <- "outputs/loo_results_insT_lasso.rds"
STRICT_RDS <- "outputs/strict_featsel/loo_results_insT_lasso_strict.rds"
OUT_PDF    <- "outputs/strict_featsel/penalty_diagnostic_cv_curves.pdf"
STRICT_SPLIT <- "data/split_strict_featsel.csv"
DEV_PSI    <- "outputs/psi_matrix.rds"
for (f in c(CANON_RDS, STRICT_RDS, STRICT_SPLIT, DEV_PSI)) if (!file.exists(f)) stop("Not found: ", f)

canon  <- readRDS(CANON_RDS)
strict <- readRDS(STRICT_RDS)
fits <- list(canonical = canon, strict = strict)
for (nm in names(fits))
  stopifnot(all(c("cv_fit", "X_train", "y_train", "foldid", "a5ss_cand") %in% names(fits[[nm]])))

gene_of <- function(fit, ev) fit$a5ss_cand$GeneSymbol[match(ev, fit$a5ss_cand$gel_ev_id)]
nz_at <- function(cv, s) {
  b <- as.numeric(coef(cv, s = s))[-1]
  names(b) <- rownames(coef(cv, s = s))[-1]
  names(b)[b != 0]
}

# --------------------------------------------------------------------------
# 1. Penalties and cross-evaluated K
# --------------------------------------------------------------------------
cat("=== 1. Penalty (lambda) ===\n")
for (nm in names(fits)) {
  cv <- fits[[nm]]$cv_fit
  cat(sprintf("  %-9s lambda.1se=%.4g (K=%d)  lambda.min=%.4g (K=%d)  candidates=%d\n",
              nm, cv$lambda.1se, length(nz_at(cv, "lambda.1se")),
              cv$lambda.min, length(nz_at(cv, "lambda.min")), ncol(fits[[nm]]$X_train)))
}
l_can <- canon$cv_fit$lambda.1se
l_str <- strict$cv_fit$lambda.1se
ev_str_at_can <- nz_at(strict$cv_fit, l_can)
ev_can_at_str <- nz_at(canon$cv_fit, l_str)
cat(sprintf("\n  Strict fit at canonical lambda.1se (%.4g): K=%d -> %s\n", l_can,
            length(ev_str_at_can), paste(gene_of(strict, ev_str_at_can), collapse = ", ")))
cat(sprintf("  Canonical fit at strict lambda.1se (%.4g): K=%d -> %s\n", l_str,
            length(ev_can_at_str), paste(gene_of(canon, ev_can_at_str), collapse = ", ")))

# --------------------------------------------------------------------------
# 2. Entry order along each path (largest lambda at which coef becomes non-zero)
# --------------------------------------------------------------------------
entry_order <- function(fit) {
  g  <- fit$cv_fit$glmnet.fit
  b  <- as.matrix(g$beta)
  first <- apply(b != 0, 1, function(v) if (any(v)) which(v)[1] else NA_integer_)
  d <- data.frame(gene = gene_of(fit, rownames(b)), event = rownames(b),
                  enters_at_lambda = ifelse(is.na(first), NA, g$lambda[first]),
                  stringsAsFactors = FALSE)
  d <- d[!is.na(d$enters_at_lambda), ]
  d <- d[order(-d$enters_at_lambda), ]
  d$order <- seq_len(nrow(d))
  d
}
canon_6 <- nz_at(canon$cv_fit, "lambda.1se")
cat("\n=== 2. Order events enter the model as the penalty is lowered (first 15) ===\n")
for (nm in names(fits)) {
  eo <- entry_order(fits[[nm]])
  eo$canonical_6 <- eo$event %in% canon_6
  eo$selected_1se <- eo$event %in% nz_at(fits[[nm]]$cv_fit, "lambda.1se")
  cat(sprintf("\n  %s\n", nm))
  print(head(eo[, c("order", "gene", "enters_at_lambda", "canonical_6", "selected_1se")], 15),
        row.names = FALSE, digits = 3)
}

# --------------------------------------------------------------------------
# 3. Were the strict-only events available to the canonical model?
# --------------------------------------------------------------------------
cat("\n=== 3. Strict-only events vs canonical candidate pool ===\n")
strict_sel <- nz_at(strict$cv_fit, "lambda.1se")
new_ev <- setdiff(strict_sel, canon_6)
canon_pre  <- canon$candidate_sites_pre$gel_ev_id
canon_post <- canon$a5ss_cand$gel_ev_id
print(data.frame(
  gene                    = gene_of(strict, new_ev),
  in_canonical_depth_pool = new_ev %in% canon_pre,
  in_canonical_final_pool = new_ev %in% canon_post,
  same_gene_other_event_in_canonical_final =
    gene_of(strict, new_ev) %in% canon$a5ss_cand$GeneSymbol & !new_ev %in% canon_post),
  row.names = FALSE)

# --------------------------------------------------------------------------
# 4. Selection frequency across case-LOO refits at each fit's own lambda.1se
# --------------------------------------------------------------------------
sel_freq <- function(fit) {
  X <- fit$X_train; y <- fit$y_train; fold <- fit$foldid
  lam <- fit$cv_fit$lambda.1se
  folds <- sort(unique(fold))
  hits <- setNames(numeric(ncol(X)), colnames(X))
  for (k in folds) {
    tr <- fold != k
    g <- glmnet(X[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 1,
                standardize = TRUE, lower.limits = 0, lambda = fit$cv_fit$lambda)
    b <- as.numeric(coef(g, s = lam))[-1]
    hits <- hits + (b != 0)
  }
  out <- data.frame(gene = gene_of(fit, names(hits)), n_folds = length(folds),
                    selected_in = as.integer(hits),
                    canonical_6 = names(hits) %in% canon_6,
                    stringsAsFactors = FALSE)
  out <- out[out$selected_in > 0, ]
  out[order(-out$selected_in), ]
}
cat("\n=== 4. Selection frequency, leaving each training case out in turn ===\n")
for (nm in names(fits)) {
  cat(sprintf("\n  %s (lambda.1se held fixed)\n", nm))
  print(sel_freq(fits[[nm]]), row.names = FALSE)
}

# --------------------------------------------------------------------------
# 5. CV curves
# --------------------------------------------------------------------------
pdf(OUT_PDF, width = 170 / 25.4, height = 80 / 25.4)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1), cex = 0.7)
for (nm in names(fits)) {
  plot(fits[[nm]]$cv_fit, main = "")
  title(main = nm, line = 2)
}
invisible(dev.off())
cat("\nWritten:", OUT_PDF, "\n")

# --------------------------------------------------------------------------
# 6. Per-fold deviance at lambda.min (what the CV error bars are made of)
# --------------------------------------------------------------------------
# fit.preval = out-of-fold prediction for every training sample, one column
# per lambda. glmnet versions differ on whether it is stored as probability
# or linear predictor; convert to probability either way.
preval_prob <- function(cv, s) {
  v <- cv$fit.preval[, which.min(abs(cv$lambda - s))]
  if (all(v >= 0 & v <= 1, na.rm = TRUE)) v else 1 / (1 + exp(-v))
}
fold_deviance <- function(fit) {
  cv <- fit$cv_fit
  p  <- pmin(pmax(preval_prob(cv, cv$lambda.min), 1e-15), 1 - 1e-15)
  y  <- fit$y_train
  dev_i <- -2 * (y * log(p) + (1 - y) * log(1 - p))
  ids <- rownames(fit$X_train)
  per_fold <- do.call(rbind, lapply(sort(unique(fit$foldid)), function(k) {
    r <- fit$foldid == k
    worst <- which(r)[which.max(dev_i[r])]
    data.frame(fold = k,
               left_out_case = ids[r & y == 1L][1],
               case_P = round(p[r & y == 1L][1], 3),
               max_ctrl_P = round(max(p[r & y == 0L]), 3),
               fold_deviance = round(mean(dev_i[r]), 3),
               worst_sample = ids[worst],
               worst_sample_dev = round(dev_i[worst], 3),
               stringsAsFactors = FALSE)
  }))
  per_fold[order(-per_fold$fold_deviance), ]
}
cat("\n=== 6. Per-fold CV deviance at lambda.min (sorted, worst first) ===\n")
for (nm in names(fits)) {
  pf <- fold_deviance(fits[[nm]])
  cat(sprintf("\n  %s: mean %.3f, SD across folds %.3f\n", nm,
              mean(pf$fold_deviance), sd(pf$fold_deviance)))
  print(head(pf, 6), row.names = FALSE)
}

# Strict held-out controls: were they among the hard controls in canonical CV?
split_s <- read.csv(STRICT_SPLIT, stringsAsFactors = FALSE)
ho_ctrl <- split_s$sample_id[split_s$set == "holdout" & grepl(":control:", split_s$sample_id)]
cv_c    <- canon$cv_fit
p_can   <- preval_prob(cv_c, cv_c$lambda.min)
names(p_can) <- rownames(canon$X_train)
ctrl_p  <- p_can[canon$y_train == 0L]
cat("\n  Strict held-out controls in canonical CV (out-of-fold P(case), rank among all controls):\n")
print(data.frame(sample_id = ho_ctrl,
                 canonical_P = round(ctrl_p[ho_ctrl], 4),
                 rank_of_n = paste0(rank(-ctrl_p)[ho_ctrl], "/", length(ctrl_p))),
      row.names = FALSE)

# --------------------------------------------------------------------------
# 7. The hardest GEL case: data quality first, then which events carry its signal
# --------------------------------------------------------------------------
pf_can   <- fold_deviance(canon)
FOCUS_ID <- pf_can$left_out_case[grepl("^gel:", pf_can$left_out_case)][1]
cat(sprintf("\n=== 7. Hardest GEL case in canonical CV: %s ===\n", FOCUS_ID))
dev <- readRDS(DEV_PSI)
stopifnot(FOCUS_ID %in% rownames(dev$ijc))
depth_all <- dev$ijc + dev$sjc
tot <- rowSums(depth_all, na.rm = TRUE)
gel_ids <- rownames(depth_all)
cat(sprintf("  Total junction reads over all A5SS events: %.3g (rank %d of %d GEL dev samples; median %.3g)\n",
            tot[FOCUS_ID], rank(-tot)[FOCUS_ID], length(tot), median(tot)))

event_table <- function(fit, label) {
  cv   <- fit$cv_fit
  b    <- as.numeric(coef(cv, s = "lambda.1se"))[-1]
  names(b) <- rownames(coef(cv, s = "lambda.1se"))[-1]
  sel  <- names(b)[b != 0]
  X    <- fit$X_train
  ids  <- rownames(X)
  y    <- fit$y_train
  if (!FOCUS_ID %in% ids) stop(FOCUS_ID, " not in ", label, " training set")
  is_gel  <- grepl("^gel:", ids)
  oth_cas <- ids[y == 1L & is_gel & ids != FOCUS_ID]
  ctrls   <- ids[y == 0L & is_gel]
  d <- data.frame(
    gene            = gene_of(fit, sel),
    beta            = round(b[sel], 3),
    depth           = as.integer(depth_all[FOCUS_ID, sel]),
    depth_median_gel = as.integer(apply(depth_all[gel_ids, sel, drop = FALSE], 2, median, na.rm = TRUE)),
    feature         = round(X[FOCUS_ID, sel], 2),
    other_cases_median = round(apply(X[oth_cas, sel, drop = FALSE], 2, median), 2),
    other_cases_min = round(apply(X[oth_cas, sel, drop = FALSE], 2, min), 2),
    ctrl_max        = round(apply(X[ctrls, sel, drop = FALSE], 2, max), 2),
    contribution    = round(b[sel] * X[FOCUS_ID, sel], 2),
    stringsAsFactors = FALSE)
  d <- d[order(-d$contribution), ]
  cat(sprintf("\n  %s model (K=%d). feature = direction * (logit PSI - training control mean);\n", label, length(sel)))
  cat("  positive = shifted the case way. other_cases_* = other GEL insT training cases;\n")
  cat("  ctrl_max = highest GEL training control. contribution = beta * feature.\n")
  print(d, row.names = FALSE)
  cat(sprintf("  Sum of contributions: %s %.2f | other GEL insT cases median %.2f | GEL controls max %.2f\n",
              FOCUS_ID, sum(d$contribution),
              median(as.numeric(X[oth_cas, sel, drop = FALSE] %*% b[sel])),
              max(as.numeric(X[ctrls, sel, drop = FALSE] %*% b[sel]))))
}
event_table(canon,  "canonical")
event_table(strict, "strict")
