#!/usr/bin/env Rscript
# =============================================================================
# figure_s4.R
# Supplementary Figure S2 (RNU4-2 splicing biomarker manuscript).
# NOT a plot -- Figure S2 itself is being built manually in Illustrator.
# This script just computes and records the sites-filtering funnel counts
# (both insT and Stem III) so there's a single, reproducible source of
# truth for whatever numbers go into that figure, rather than reading them
# off various script outputs by hand.
#
#   all events quantified (either cohort) -> tested in both cohorts ->
#   significant in both cohorts (FDR<0.05) -> A5SS only -> direction-
#   concordant -> [insT / Stem III fork] -> depth-filtered (that
#   subgroup's own cases+controls) -> best event per gene (that
#   subgroup's own G_AUC) -> LASSO-selected (K=6)
#
# Every count is read directly from already-computed, already-verified
# outputs -- nothing here is re-derived from raw PSI matrices, to avoid any
# risk of drifting from the canonical figure1.R / loo_evaluate_*_lasso.R
# numbers:
#   - trunk counts (all/both/sig-both/A5SS/concordant): outputs/figure1_summary.rds
#     (n_unique_combined, n_shared, n_sig_shared; summary_df's A5SS row for
#     sig_both/concordant)
#   - insT/stemIII branch counts (depth-filtered/best-per-gene/LASSO-
#     selected): outputs/loo_results_insT_lasso.rds /
#     loo_results_stemIII_lasso.rds ($candidate_sites_pre = post-depth
#     pre-MAX_PER_GENE, $a5ss_cand = post-MAX_PER_GENE, $selected_sites =
#     final K=6)
#
# insT and stemIII each do their OWN independent depth-filter and
# MAX_PER_GENE=1 selection (subgroup-specific G_AUC ranking, per
# project_rnu42_region_bounds / CLAUDE.md) -- their branch counts are NOT
# a partition of n_concordant, each subgroup narrows the same upstream
# 95-event pool separately.
#
# Inputs:
#   outputs/figure1_summary.rds
#   outputs/loo_results_insT_lasso.rds
#   outputs/loo_results_stemIII_lasso.rds
#
# Output:
#   outputs/figure_s4_summary.rds   (counts only; airlock-safe -- no gene
#                                    names, event_ids, or coefficients)
#
# DATA PROTECTION: only aggregate counts are computed/printed -- no
# event_ids, gene names, or LASSO coefficients appear anywhere here.
#
# Run: pixi run Rscript scripts/gel/figure_s4.R
# =============================================================================

# ---- 0. libpath (pixi-aware) ------------------------------------------------
invisible(local({
  user_lib <- path.expand("~/re_gecip/machine_learning/R_packages")
  if (nzchar(Sys.getenv("PIXI_ENVIRONMENT_NAME"))) {
    keep <- !normalizePath(.libPaths(), mustWork = FALSE) %in%
            normalizePath(user_lib, mustWork = FALSE)
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

# ---- 1. Config --------------------------------------------------------------
FIG1_SUMMARY_PATH <- Sys.getenv("FIG1_SUMMARY",     "outputs/figure1_summary.rds")
INST_LOO_PATH     <- Sys.getenv("INST_LOO_RDS",     "outputs/loo_results_insT_lasso.rds")
STEMIII_LOO_PATH  <- Sys.getenv("STEMIII_LOO_RDS",  "outputs/loo_results_stemIII_lasso.rds")
OUT_SUMMARY       <- "outputs/figure_s4_summary.rds"
dir.create("outputs", showWarnings = FALSE)

stopifnot(
  "outputs/figure1_summary.rds not found -- run figure1.R first" =
    file.exists(FIG1_SUMMARY_PATH),
  "insT LOO RDS not found"    = file.exists(INST_LOO_PATH),
  "stemIII LOO RDS not found" = file.exists(STEMIII_LOO_PATH)
)

# =============================================================================
# 2. Trunk counts (all event types, both cohorts) -- from figure1_summary.rds
# =============================================================================
fig1 <- readRDS(FIG1_SUMMARY_PATH)
stopifnot(
  all(c("n_gel_only", "n_nava_only", "n_shared", "n_sig_shared", "summary") %in% names(fig1))
)

a5ss_row <- fig1$summary[fig1$summary$EventType == "A5SS", , drop = FALSE]
stopifnot("No A5SS row in figure1_summary.rds$summary" = nrow(a5ss_row) == 1L)

# NOTE: fig1$n_unique_combined is NOT the union -- it's n_gel_only +
# n_nava_only, i.e. events quantified in exactly one cohort (built for the
# paragraph-1 novel-SS-rate text in figure1.R, which deliberately excludes
# shared events). The true "either cohort" union is n_gel_only + n_nava_only
# + n_shared -- computed directly here instead.
n_all        <- fig1$n_gel_only + fig1$n_nava_only + fig1$n_shared
n_both       <- fig1$n_shared
n_sigboth    <- fig1$n_sig_shared
n_a5ss       <- a5ss_row$sig_both
n_concordant <- a5ss_row$concordant

stopifnot(
  "Trunk counts must be monotonically non-increasing" =
    n_all >= n_both, n_both >= n_sigboth,
    n_sigboth >= n_a5ss, n_a5ss >= n_concordant
)

cat("=== Trunk (all event types, both cohorts) ===\n")
cat(sprintf("All events quantified (either cohort): %s\n", format(n_all, big.mark = ",")))
cat(sprintf("Tested in both cohorts:                %s\n", format(n_both, big.mark = ",")))
cat(sprintf("Significant in both (FDR<0.05):         %s\n", format(n_sigboth, big.mark = ",")))
cat(sprintf("A5SS only:                              %s\n", format(n_a5ss, big.mark = ",")))
cat(sprintf("Direction-concordant:                   %s (%d discordant dropped)\n",
            format(n_concordant, big.mark = ","), n_a5ss - n_concordant))

# =============================================================================
# 3. Branch counts -- read directly from the frozen LOO RDS objects, one
#    independent depth-filter + MAX_PER_GENE selection per subgroup.
#
# read_branch_counts() returns aggregate counts only (airlock-safe). The raw
# gel_ev_id vectors needed for the insT/stemIII overlap check in section 4
# are read separately by read_branch_ids() and never printed/saved directly
# -- only the set sizes (only-insT / shared / only-stemIII) are.
# =============================================================================
read_branch_counts <- function(path) {
  x <- readRDS(path)
  stopifnot(
    "candidate_sites_pre" %in% names(x),
    "a5ss_cand"           %in% names(x),
    "selected_idx"        %in% names(x),
    "K_selected"          %in% names(x)
  )
  cand_pre <- x$candidate_sites_pre
  stopifnot(
    "candidate_sites_pre missing GeneSymbol column" = "GeneSymbol" %in% names(cand_pre),
    "candidate_sites_pre missing g_auc column"      = "g_auc"      %in% names(cand_pre)
  )
  list(
    n_depth       = nrow(cand_pre),
    n_genes_pool  = length(unique(cand_pre$GeneSymbol)),
    n_dropped_auc = sum(is.na(cand_pre$g_auc)),
    n_bestgene    = nrow(x$a5ss_cand),
    K             = x$K_selected
  )
}

# selected_sites (nz_union of lambda.1se/lambda.min, kept for the printed
# display table) is NOT the canonical K -- use selected_idx / K_selected
# (the chosen_lambda_label selection) instead, indexed into a5ss_cand.
read_branch_ids <- function(path) {
  x <- readRDS(path)
  stopifnot("selected_idx" %in% names(x))
  list(
    depth_ids    = x$candidate_sites_pre$gel_ev_id,
    bestgene_ids = x$a5ss_cand$gel_ev_id,
    selected_ids = x$a5ss_cand$gel_ev_id[x$selected_idx]
  )
}

inst    <- read_branch_counts(INST_LOO_PATH)
stemIII <- read_branch_counts(STEMIII_LOO_PATH)

stopifnot(
  "insT branch counts must be monotonically non-increasing" =
    n_concordant >= inst$n_depth, inst$n_depth >= inst$n_bestgene,
    inst$n_bestgene >= inst$K,
  "stemIII branch counts must be monotonically non-increasing" =
    n_concordant >= stemIII$n_depth, stemIII$n_depth >= stemIII$n_bestgene,
    stemIII$n_bestgene >= stemIII$K,
  "n_genes_pool must be <= n_depth and >= n_bestgene" =
    inst$n_genes_pool    <= inst$n_depth,    inst$n_genes_pool    >= inst$n_bestgene,
    stemIII$n_genes_pool <= stemIII$n_depth, stemIII$n_genes_pool >= stemIII$n_bestgene,
  "n_dropped_auc + n_bestgene must not exceed n_depth" =
    inst$n_dropped_auc    + inst$n_bestgene    <= inst$n_depth,
    stemIII$n_dropped_auc + stemIII$n_bestgene <= stemIII$n_depth
)

cat("\n=== insT branch ===\n")
cat(sprintf("Depth-filtered (>=10):                 %d\n", inst$n_depth))
cat(sprintf("  Unique genes in pool (pre-AUC-filter): %d\n", inst$n_genes_pool))
cat(sprintf("  Dropped for non-positive AUC:          %d\n", inst$n_dropped_auc))
cat(sprintf("Best event per gene:                    %d\n", inst$n_bestgene))
cat(sprintf("LASSO-selected (K):                     %d\n", inst$K))

cat("\n=== Stem III branch ===\n")
cat(sprintf("Depth-filtered (>=10):                 %d\n", stemIII$n_depth))
cat(sprintf("  Unique genes in pool (pre-AUC-filter): %d\n", stemIII$n_genes_pool))
cat(sprintf("  Dropped for non-positive AUC:          %d\n", stemIII$n_dropped_auc))
cat(sprintf("Best event per gene:                    %d\n", stemIII$n_bestgene))
cat(sprintf("LASSO-selected (K):                     %d\n", stemIII$K))

# =============================================================================
# 4. insT/stemIII overlap at each step -- how many depth-filtered / best-per-
#    gene / LASSO-selected sites are shared between the two independently-
#    filtered tracks, vs unique to each. Matched by gel_ev_id (both tracks
#    are drawn from the same GEL rMATS event space, so ids are comparable
#    directly). Only set sizes are printed/saved -- never the ids themselves.
# =============================================================================
inst_ids    <- read_branch_ids(INST_LOO_PATH)
stemIII_ids <- read_branch_ids(STEMIII_LOO_PATH)

overlap_counts <- function(a, b) {
  list(
    only_a = length(setdiff(a, b)),
    shared = length(intersect(a, b)),
    only_b = length(setdiff(b, a))
  )
}

overlap_depth    <- overlap_counts(inst_ids$depth_ids,    stemIII_ids$depth_ids)
overlap_bestgene <- overlap_counts(inst_ids$bestgene_ids, stemIII_ids$bestgene_ids)
overlap_selected <- overlap_counts(inst_ids$selected_ids, stemIII_ids$selected_ids)

stopifnot(
  "overlap set sizes must reconcile with branch counts" =
    overlap_depth$only_a    + overlap_depth$shared    == inst$n_depth,
    overlap_depth$only_b    + overlap_depth$shared    == stemIII$n_depth,
    overlap_bestgene$only_a + overlap_bestgene$shared == inst$n_bestgene,
    overlap_bestgene$only_b + overlap_bestgene$shared == stemIII$n_bestgene,
    overlap_selected$only_a + overlap_selected$shared == inst$K,
    overlap_selected$only_b + overlap_selected$shared == stemIII$K
)

print_overlap <- function(label, ov) {
  cat(sprintf("%s -- insT only: %d, shared: %d, stemIII only: %d\n",
              label, ov$only_a, ov$shared, ov$only_b))
}

cat("\n=== insT / stemIII overlap by step ===\n")
print_overlap("Depth-filtered pool ", overlap_depth)
print_overlap("Best event per gene ", overlap_bestgene)
print_overlap("LASSO-selected      ", overlap_selected)

# =============================================================================
# 5. Sample sizes -- cases+controls each branch's pipeline was run on, for
#    annotating the "insT/Stem III Signature Development" header boxes.
#    depth_scope = cases used for the depth filter (train+holdout for insT;
#    same as fit_scope for stemIII, which has no holdout split -- see
#    loo_evaluate_stemIII_lasso.R). fit_scope = cases used for the
#    best-event-per-gene G_AUC ranking AND the LASSO training matrix (same
#    scope for both -- train-only, holdout excluded).
# =============================================================================
read_sample_sizes <- function(path) {
  x <- readRDS(path)
  stopifnot("sample_sizes" %in% names(x))
  ss <- x$sample_sizes
  stopifnot(all(c("depth_scope", "fit_scope") %in% names(ss)))
  totalize <- function(scope) with(scope, gel_cases + gel_ctrl + nava_cases + nava_ctrl)
  list(
    depth_scope = c(as.list(ss$depth_scope), total = totalize(ss$depth_scope)),
    fit_scope   = c(as.list(ss$fit_scope),   total = totalize(ss$fit_scope))
  )
}

inst_ss    <- read_sample_sizes(INST_LOO_PATH)
stemIII_ss <- read_sample_sizes(STEMIII_LOO_PATH)

cat("\n=== insT branch sample sizes ===\n")
cat(sprintf("Depth-filter scope (train+holdout cases + controls): GEL %d+%d, Nava %d+%d = %d total\n",
            inst_ss$depth_scope$gel_cases, inst_ss$depth_scope$gel_ctrl,
            inst_ss$depth_scope$nava_cases, inst_ss$depth_scope$nava_ctrl, inst_ss$depth_scope$total))
cat(sprintf("Fit scope (train-only cases + controls, used for best-event-per-gene AND LASSO): GEL %d+%d, Nava %d+%d = %d total\n",
            inst_ss$fit_scope$gel_cases, inst_ss$fit_scope$gel_ctrl,
            inst_ss$fit_scope$nava_cases, inst_ss$fit_scope$nava_ctrl, inst_ss$fit_scope$total))

cat("\n=== Stem III branch sample sizes ===\n")
cat(sprintf("Depth-filter scope (no holdout -- same as fit scope): GEL %d+%d, Nava %d+%d = %d total\n",
            stemIII_ss$depth_scope$gel_cases, stemIII_ss$depth_scope$gel_ctrl,
            stemIII_ss$depth_scope$nava_cases, stemIII_ss$depth_scope$nava_ctrl, stemIII_ss$depth_scope$total))
cat(sprintf("Fit scope (best-event-per-gene AND LASSO): GEL %d+%d, Nava %d+%d = %d total\n",
            stemIII_ss$fit_scope$gel_cases, stemIII_ss$fit_scope$gel_ctrl,
            stemIII_ss$fit_scope$nava_cases, stemIII_ss$fit_scope$nava_ctrl, stemIII_ss$fit_scope$total))

# =============================================================================
# 6. Save summary (counts only -- airlock-safe)
# =============================================================================
saveRDS(list(
  trunk   = list(n_all = n_all, n_both = n_both, n_sigboth = n_sigboth,
                 n_a5ss = n_a5ss, n_concordant = n_concordant),
  insT    = inst,
  stemIII = stemIII,
  overlap = list(depth = overlap_depth, bestgene = overlap_bestgene,
                  selected = overlap_selected),
  sample_sizes = list(insT = inst_ss, stemIII = stemIII_ss)
), OUT_SUMMARY)
cat(sprintf("\nWrote %s\n", OUT_SUMMARY))
cat("\nDone.\n")
