#!/usr/bin/env Rscript
# =============================================================================
# export_insT_signature.R
# Export the insT-only LASSO K selected A5SS splicing signature into a
# portable CSV + metadata text file.
#
# Inputs:
#   outputs/loo_results_insT_lasso.rds  (selected_idx, a5ss_cand, fit_at_1se,
#                                        gel_ctrl_means, nava_ctrl_means,
#                                        y_train, cohort_train)
#   data/rmats/A5SS.MATS.JC.txt         (ENSG + strand lookup by coordinates)
#
# Outputs:
#   outputs/insT_signature_export.csv        (K rows: coordinates, direction,
#                                             beta, gene_symbol, ENSG, rank)
#   outputs/insT_signature_export_meta.txt   (intercept, lambda.1se, K,
#                                             training counts, scoring formula)
#
# Run on GEL:
#   Rscript scripts/gel/export_insT_signature.R
#
# Scoring recipe for a new cohort:
#   1. Match rMATS A5SS events to signature sites by
#      chr + longES_0 + longEE_1 + shortES_0 + shortEE_1 (exact match).
#   2. bb-logit-PSI = logit((IJC + 0.5) / (IJC + SJC + 1)) per site.
#   3. Recompute mu_ctrl[j] on the new cohort's own controls.
#   4. oriented_contrib(i,j) = direction[j] * (bb_logit_PSI[i,j] - mu_ctrl[j])
#   5. uniform_sum_score(i)  = sum_j oriented_contrib(i,j)
#   6. lasso_P_case(i)       = sigmoid(intercept + sum_j beta[j] * oriented_contrib(i,j))
# =============================================================================

# ---- 0. libpath -------------------------------------------------------------
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

# ---- 1. Configuration -------------------------------------------------------
LOO_PATH   <- Sys.getenv("INST_LOO", "outputs/loo_results_insT_lasso.rds")
RMATS_DIR  <- Sys.getenv("GEL_RMATS_DIR", "data/rmats")
OUT_CSV    <- Sys.getenv("OUT_CSV",  "outputs/insT_signature_export.csv")
OUT_META   <- Sys.getenv("OUT_META", "outputs/insT_signature_export_meta.txt")
SCHEMA_VER <- "insT-lasso-1.0"

stopifnot("insT LOO RDS not found" = file.exists(LOO_PATH))
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

# ---- 2. Load + extract ------------------------------------------------------
cat(sprintf("Loading insT LOO RDS: %s\n", LOO_PATH))
loo <- readRDS(LOO_PATH)

required <- c("selected_idx", "a5ss_cand", "fit_at_1se",
              "gel_ctrl_means", "nava_ctrl_means")
stopifnot("insT LOO RDS missing required fields" =
            all(required %in% names(loo)))

sel_idx    <- loo$selected_idx
sel_cand   <- loo$a5ss_cand[sel_idx, , drop = FALSE]
sel_events <- sel_cand$gel_ev_id
beta_1se   <- loo$fit_at_1se$beta
beta_sel   <- unname(beta_1se[sel_idx])
intercept  <- unname(loo$fit_at_1se$intercept)

K_sel <- length(sel_idx)
cat(sprintf("K = %d selected sites at lambda.1se; intercept = %.6f\n",
            K_sel, intercept))
stopifnot("beta_sel must be >= 0 (lower.limits=0 in glmnet)" =
            all(beta_sel >= 0))

# ---- 3. Parse coordinates from event_id -------------------------------------
# event_id format: "A5SS:GeneSymbol:chr:longES_0:longEE_1:shortES_0:shortEE_1:flankES_0:flankEE_1"
parse_eid <- function(eid) {
  parts <- strsplit(eid, ":", fixed = TRUE)[[1L]]
  if (length(parts) != 9L || parts[1L] != "A5SS")
    stop("Unexpected event_id format: ", eid)
  list(
    GeneSymbol = parts[2L],
    chr        = parts[3L],
    longES_0   = as.integer(parts[4L]),
    longEE_1   = as.integer(parts[5L]),
    shortES_0  = as.integer(parts[6L]),
    shortEE_1  = as.integer(parts[7L]),
    flankES_0  = as.integer(parts[8L]),
    flankEE_1  = as.integer(parts[9L]))
}
parsed <- lapply(sel_events, parse_eid)

# Prefer values directly from a5ss_cand where available, fall back to parsed.
pick <- function(col, pkey = NULL) {
  if (col %in% names(sel_cand)) return(as.character(sel_cand[[col]]))
  if (!is.null(pkey)) return(vapply(parsed, `[[`, character(1L), pkey))
  rep(NA_character_, K_sel)
}
pick_int <- function(col, pkey = NULL) {
  if (col %in% names(sel_cand)) return(as.integer(sel_cand[[col]]))
  if (!is.null(pkey)) return(vapply(parsed, `[[`, integer(1L), pkey))
  rep(NA_integer_, K_sel)
}

gene_symbol <- pick("GeneSymbol", "GeneSymbol")
chr_v       <- pick("chr",        "chr")
longES_0    <- pick_int("longES_0",  "longES_0")
longEE_1    <- pick_int("longEE_1",  "longEE_1")
shortES_0   <- pick_int("shortES_0", "shortES_0")
shortEE_1   <- pick_int("shortEE_1", "shortEE_1")
flankES_0   <- pick_int("flankES_0", "flankES_0")
flankEE_1   <- pick_int("flankEE_1", "flankEE_1")

# ---- 4. ENSG + strand from rMATS A5SS file ----------------------------------
# Join on chr + all 6 rMATS coordinates (columns 4,6,7,8,9,10,11). A
# 4-coordinate key (long/short exon only, dropping flanking exon) can
# collide across distinct rMATS rows that reuse the same long/short exon
# boundary pair with a different flanking exon.
ensg   <- rep(NA_character_, K_sel)
strand <- rep(NA_character_, K_sel)

rmats_a5ss <- file.path(RMATS_DIR, "A5SS.MATS.JC.txt")
if (file.exists(rmats_a5ss)) {
  cat(sprintf("Parsing ENSG + strand from rMATS A5SS file: %s\n", rmats_a5ss))
  a5 <- read.delim(rmats_a5ss, stringsAsFactors = FALSE)
  a5_key  <- paste(a5[[4L]],
                   as.integer(a5[[6L]]), as.integer(a5[[7L]]),
                   as.integer(a5[[8L]]), as.integer(a5[[9L]]),
                   as.integer(a5[[10L]]), as.integer(a5[[11L]]),
                   sep = ":")
  a5_ensg <- sub("\\..*$", "", gsub('"', "", a5[[2L]]))
  a5_str  <- as.character(a5[[5L]])
  sig_key <- paste(chr_v, longES_0, longEE_1, shortES_0, shortEE_1,
                   flankES_0, flankEE_1, sep = ":")
  m       <- match(sig_key, a5_key)
  ensg    <- a5_ensg[m]
  strand  <- a5_str[m]
  cat(sprintf("  Matched ENSG  for %d / %d sites\n", sum(!is.na(ensg)),   K_sel))
  cat(sprintf("  Matched strand for %d / %d sites\n", sum(!is.na(strand)), K_sel))
} else {
  cat(sprintf("rMATS A5SS file not found at %s; ENSG + strand will be NA.\n",
              rmats_a5ss))
}

direction <- as.integer(sel_cand$direction)
stopifnot("direction must be +/- 1" = all(direction %in% c(-1L, 1L)),
          "direction length != K"   = length(direction) == K_sel)

# ---- 5. Build CSV -----------------------------------------------------------
ord <- order(-abs(beta_sel), gene_symbol)

out <- data.frame(
  rank                 = seq_along(ord),
  gene_symbol          = gene_symbol[ord],
  ensg                 = ensg[ord],
  event_id             = sel_events[ord],
  chr                  = chr_v[ord],
  strand               = strand[ord],
  longExonStart_0base  = longES_0[ord],
  longExonEnd_1based   = longEE_1[ord],
  shortES_0base        = shortES_0[ord],
  shortEE_1based       = shortEE_1[ord],
  flankingES_0base     = flankES_0[ord],
  flankingEE_1based    = flankEE_1[ord],
  direction            = direction[ord],
  beta_at_lambda1se    = beta_sel[ord],
  stringsAsFactors     = FALSE)

stopifnot("Row count != K" = nrow(out) == K_sel)
write.csv(out, OUT_CSV, row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote signature CSV: %s  (%d rows)\n", OUT_CSV, nrow(out)))

# ---- 6. Metadata text file --------------------------------------------------
ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

n_train_GEL  <- if (all(c("y_train", "cohort_train") %in% names(loo)))
  sum(loo$y_train == 1L & loo$cohort_train == "GEL")  else NA_integer_
n_train_Nava <- if (all(c("y_train", "cohort_train") %in% names(loo)))
  sum(loo$y_train == 1L & loo$cohort_train == "Nava") else NA_integer_
n_ctrl_GEL   <- if (all(c("y_train", "cohort_train") %in% names(loo)))
  sum(loo$y_train == 0L & loo$cohort_train == "GEL")  else NA_integer_
n_ctrl_Nava  <- if (all(c("y_train", "cohort_train") %in% names(loo)))
  sum(loo$y_train == 0L & loo$cohort_train == "Nava") else NA_integer_
inst_hgvs    <- if ("config" %in% names(loo)) loo$config$INST_HGVS else "n.64_65insT"

meta_lines <- c(
  "# insT-only LASSO splicing signature export",
  paste0("# schema_version  : ", SCHEMA_VER),
  paste0("# generated_at    : ", ts),
  paste0("# source_rds      : ", LOO_PATH),
  paste0("# variant_scope   : ", inst_hgvs, " carriers only (training set)"),
  "",
  paste0("K                : ", K_sel),
  paste0("intercept_at_1se : ", format(intercept,         digits = 10)),
  paste0("lambda_1se       : ", format(loo$lambda_1se,    digits = 6)),
  paste0("lambda_min       : ", format(loo$lambda_min,    digits = 6)),
  "",
  "# training cohort composition (insT cases + all controls)",
  paste0("n_train_GEL_WB   : ", n_train_GEL),
  paste0("n_train_Nava_LY  : ", n_train_Nava),
  paste0("n_ctrl_GEL_WB    : ", n_ctrl_GEL),
  paste0("n_ctrl_Nava_LY   : ", n_ctrl_Nava),
  "",
  "# scoring formula (identical to v4)",
  "# bb_logit_PSI(i,j)        = log((IJC[i,j] + 0.5) / (SJC[i,j] + 0.5))",
  "# oriented_contrib(i,j)    = direction[j] * (bb_logit_PSI[i,j] - mu_ctrl[j])",
  "# uniform_sum_score(i)     = sum over j of oriented_contrib(i,j)",
  "# lasso_linpred(i)         = intercept + sum over j of beta[j] * oriented_contrib(i,j)",
  "# lasso_P_case(i)          = 1 / (1 + exp(-lasso_linpred(i)))",
  "",
  "# IMPORTANT: recompute mu_ctrl[j] on the new cohort's own controls.",
  "# direction[j], beta[j], and intercept come from the CSV unchanged."
)
writeLines(meta_lines, OUT_META)
cat(sprintf("Wrote metadata file: %s\n", OUT_META))

# ---- 7. Console preview -----------------------------------------------------
cat("\n=== Preview (all rows) ===\n")
preview_cols <- c("rank", "gene_symbol", "ensg", "chr", "strand",
                  "longExonStart_0base", "longExonEnd_1based",
                  "shortES_0base", "shortEE_1based",
                  "direction", "beta_at_lambda1se")
print(out[, preview_cols, drop = FALSE], row.names = FALSE)

cat("\n=== Done ===\n")
cat(sprintf("Files ready for airlock:\n  %s\n  %s\n", OUT_CSV, OUT_META))
