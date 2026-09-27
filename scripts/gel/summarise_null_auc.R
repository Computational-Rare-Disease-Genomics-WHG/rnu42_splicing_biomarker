#!/usr/bin/env Rscript
#
# summarise_null_auc.R  (Desktop, pixi)
# Compact human-readable summary of a null_lasso_permutations.R aggregate:
# non-NA counts, K-selected table, chosen-lambda table, train/holdout AUC
# summary + quantiles, and the joint AUC=1 count.
#
# Env vars:
#   SAMPLE_TAG  default n11v33  (which experiment's aggregate to read)
#   AGG_TAG     default stratified_p17
#
# Reads: outputs/null_within_batch/<SAMPLE_TAG>/<AGG_TAG>/null_auc_distribution.csv
#
# Usage:
#   pixi run Rscript scripts/gel/summarise_null_auc.R
#   SAMPLE_TAG=n19v54 pixi run Rscript scripts/gel/summarise_null_auc.R
#
# Note: for the legacy 11v33 run that predates the SAMPLE_TAG subdir layout,
# the CSV lives directly under outputs/null_within_batch/stratified_p17/.
# Pass SAMPLE_TAG="" to hit that layout.

WORKSPACE  <- "/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
SAMPLE_TAG <- Sys.getenv("SAMPLE_TAG", "n11v33")
AGG_TAG    <- Sys.getenv("AGG_TAG",    "stratified_p17")

parts <- Filter(nzchar, c("outputs", "null_within_batch",
                          SAMPLE_TAG, AGG_TAG, "null_auc_distribution.csv"))
csv_path <- do.call(file.path, c(list(WORKSPACE), as.list(parts)))

stopifnot("aggregate CSV not found" = file.exists(csv_path))
d <- read.csv(csv_path, stringsAsFactors = FALSE)

cat(sprintf("=== %s / %s ===\n", SAMPLE_TAG, AGG_TAG))
cat(sprintf("CSV: %s\n", csv_path))
cat(sprintf("Rows: %d  (non-NA train: %d, holdout: %d)\n\n",
            nrow(d),
            sum(!is.na(d$train_auc)),
            sum(!is.na(d$holdout_auc))))

cat("K selected @ lambda.1se:\n")
print(table(d$K_1se, useNA = "ifany"))

cat("\nchosen lambda:\n")
print(table(d$chosen, useNA = "ifany"))

cat("\ntrain_auc:\n");   print(summary(d$train_auc))
cat("\nholdout_auc:\n"); print(summary(d$holdout_auc))

qs <- c(0.5, 0.95, 0.99, 1)
cat("\nQuantiles (50, 95, 99, 100):\n")
cat("train:   "); print(quantile(d$train_auc,   qs, na.rm = TRUE))
cat("holdout: "); print(quantile(d$holdout_auc, qs, na.rm = TRUE))

n_valid <- sum(!is.na(d$train_auc) & !is.na(d$holdout_auc))
n_tr    <- sum(d$train_auc   >= 1.0, na.rm = TRUE)
n_ho    <- sum(d$holdout_auc >= 1.0, na.rm = TRUE)
n_both  <- sum(d$train_auc >= 1.0 & d$holdout_auc >= 1.0, na.rm = TRUE)
cat(sprintf("\nPerms hitting AUC = 1.0 (n_valid = %d):\n", n_valid))
cat(sprintf("  train  only : %d (%.2f%%)\n", n_tr,   100 * n_tr   / n_valid))
cat(sprintf("  holdout only: %d (%.2f%%)\n", n_ho,   100 * n_ho   / n_valid))
cat(sprintf("  BOTH        : %d (%.2f%%)\n", n_both, 100 * n_both / n_valid))

# --------------------------------------------------------------------------
# Perm grouping by LASSO outcome:
#   1se_fit        - lambda.1se selected >=1 site (real LASSO fit)
#   min_fit        - lambda.1se selected 0 sites, fell back to lambda.min
#                    which selected >=1 site
#   intercept_only - both lambda.1se and lambda.min selected 0 sites in the
#                    full-data fit. Holdout scoring uses the full-fit -> AUC
#                    exactly 0.5. Training uses per-fold LOO fits; per-fold
#                    fits at the same lambda often still pick a handful of
#                    sites, so training AUC is NOT necessarily 0.5 -- and
#                    with direction fixed from all-data means and
#                    lower.limits=0, LOO predictions on held-out cases tend
#                    to sit systematically below controls (overfit reversal),
#                    giving training AUC well below 0.5 for this group.
# --------------------------------------------------------------------------
d$group <- ifelse(is.na(d$K_1se) | is.na(d$K_min), "no_candidates",
           ifelse(d$K_1se > 0,                     "1se_fit",
           ifelse(d$K_min > 0,                     "min_fit",
                                                    "intercept_only")))

cat("\nPerm grouping by LASSO outcome:\n")
print(table(d$group, useNA = "ifany"))

cat("\nTraining AUC by group:\n")
print(by(d$train_auc, d$group, summary))
cat("\nHoldout AUC by group:\n")
print(by(d$holdout_auc, d$group, summary))

# --------------------------------------------------------------------------
# Histogram: null train + holdout AUC distributions overlaid
# --------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggplot2))

plot_df <- rbind(
  data.frame(set = "training", auc = d$train_auc),
  data.frame(set = "holdout",  auc = d$holdout_auc)
)
plot_df <- plot_df[!is.na(plot_df$auc), ]
plot_df$set <- factor(plot_df$set, levels = c("training", "holdout"))

fill_cols <- c(training = "#2C5AA0", holdout = "#F39B7F")   # navy / peach

p <- ggplot(plot_df, aes(x = auc, fill = set)) +
  geom_histogram(binwidth = 0.02, position = "identity", alpha = 0.55,
                 colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = fill_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 1.02),
                     breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "AUC",
       y = "Number of permutations",
       title = sprintf("Null AUC distribution: %s / %s  (n = %d perms)",
                       SAMPLE_TAG, AGG_TAG, n_valid)) +
  theme_classic(base_size = 9) +
  theme(legend.position = c(0.15, 0.90),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(size = 9, face = "plain"))

pdf_path <- sub("\\.csv$", "_histogram.pdf", csv_path)
ggsave(pdf_path, p, width = 120, height = 80, units = "mm")
cat(sprintf("\nHistogram: %s\n", pdf_path))
