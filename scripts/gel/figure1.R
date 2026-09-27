#!/usr/bin/env Rscript
# Layout: patchwork only -- never gridExtra arrangeGrob/grid.arrange
# =============================================================================
# figure1.R
# Unified generator for paper Figure 1 (RNU4-2 splicing biomarker manuscript).
# Computes every number used in results paragraphs 1-3 and produces:
#
#   Figure 1   outputs/figures/figure1.pdf
#              Two-row paper figure built in memory and combined via
#              patchwork. Panels:
#                A — stacked bar plot of RNU4-2 variant carriers across
#                    cohorts. X = HGVS variant (categorical, ordered by nt
#                    position), faceted by structural region (T-loop |
#                    Other | Stem III). Y = individuals (n), fill =
#                    cohort (Whole-Blood / Lymphocytes).
#                B — proportional Euler diagram (eulerr) of significant
#                    events (FDR < 0.05 in each cohort): Whole-Blood-only,
#                    Lymphocytes-only, and sig-in-both, restricted to
#                    events quantified in BOTH cohorts. Horizontal (side-
#                    by-side) layout. The all-events overlap is reported
#                    numerically in paragraph 1, not plotted.
#                C — PCA of all samples (GEL + Nava, cases + controls) on
#                    the PSI of the n_sig_shared events significant in both
#                    cohorts (FDR < 0.05). Colour = cohort (Whole-Blood /
#                    Lymphocytes); shape = triangle for cases, circle for
#                    controls. Complete cases only (no imputation).
#                D — Delta-logit-PSI scatter (Whole-Blood vs Lymphocytes)
#                    for all cases, A5SS events only, coloured by direction
#                    concordance, with R^2 and n annotated.
#              (|dPSI| distribution numbers for paragraph 3 -- max/median/
#              90th-percentile |IncLevelDifference| among sig-in-both events
#              -- are still computed and reported to stdout/figure1_summary.rds
#              but are not plotted as a panel.)
#              The insT LASSO coefficient bar, insT P(case) by cohort, and
#              6k whole-blood control P(case) violin (formerly panels E-G)
#              have moved to figure2.R (panels A-C). Cross-cohort genotype
#              parity (TOST) and severity correlation are deferred to
#              Figure 2 and are NOT plotted here.
#
#   Numbers   outputs/figure1_summary.rds  + stdout
#              All counts, concordance/binomial p-values, Pearson/
#              Spearman/R^2 (overall + per event type), directional Storey
#              pi1, and the joined sig-in-both event table.
#
# Event matching basis: EventType + coord_nostrand (NOT raw event_id, which
# embeds GeneSymbol and can differ between cohorts) -- same join basis used
# by the insT/stemIII LASSO pipelines. Effect-size axis is Delta-logit-PSI =
# mean_case(logit) - mean_ctrl(logit) of the count-level beta-binomial-shrunk
# log-odds logit((IJC+0.5)/(IJC+SJC+1)) per cohort -- same units as the insT/
# stemIII LASSO feature.
# Novel-SS classification uses gencode v49 (Nava ran with Ensembl v106 ≈
# gencode v40; v49 is a strict superset, so applying it uniformly under-calls
# novelty in Nava — conservative for the paragraph 1 claim).
#
# -----------------------------------------------------------------------------
# RUN INSTRUCTIONS  --  via pixi (preferred)
# -----------------------------------------------------------------------------
# This script needs eulerr for the proportional Euler diagram. Run it inside
# the repo's pixi environment, which carries eulerr + gridExtra + ggplot2 and
# all their dependencies built against one consistent R -- this avoids the
# user-lib ABI mismatches that crash system R (undefined-symbol errors loading
# rlang / Rcpp when the lib was built under a different R version).
#
#   pixi install                                  # one-time, on GEL
#   pixi run Rscript scripts/gel/figure1.R
#
# The libpath block below detects pixi via PIXI_ENVIRONMENT_NAME and then uses
# only the pixi packages. Outside pixi it falls back to GEL system R + the user
# lib ~/re_gecip/machine_learning/R_packages (eulerr must then be present
# there, install.packages()'d from the GEL Desktop) -- but the pixi route is
# preferred because that fallback is fragile across R-version bumps.
#
# DATA PROTECTION: counts and figures derive from GEL protected patient data.
# Stay in workspace; airlock the PDFs + paste counts back to chat. No weights
# involved.
# =============================================================================

# ---- 0. libpath setup (pixi-aware) -----------------------------------------
# Under pixi, use only the pixi env's packages -- mixing in the GEL system
# libpath or the ~/re_gecip user lib (built against a different R) causes ABI
# mismatches: undefined-symbol crashes loading rlang / Rcpp. Outside pixi,
# fall back to GEL system R + the user lib.
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
  if (!is.null(candidate) && dir.exists(candidate)) {
    .libPaths(c(.libPaths(), candidate))
  }
  if (dir.exists(user_lib)) {
    .libPaths(c(user_lib, .libPaths()))
  }
}))

suppressPackageStartupMessages({
  library(eulerr)
  library(grid)
  library(ggplot2)
  library(patchwork)
})

# ---- 1. Configuration -------------------------------------------------------
GEL_PSI_PATH     <- Sys.getenv("GEL_PSI",      "outputs/psi_matrix.rds")
NAVA_PSI_PATH    <- Sys.getenv("NAVA_PSI",     "outputs/nava_psi_matrix.rds")
GTF_PATH         <- Sys.getenv("GTF_PATH",     "rmats_inputs/gencode.v49.annotation.gtf")
SS_CACHE_PATH    <- Sys.getenv("SS_CACHE",     "outputs/annotated_ss_v49.rds")
GENO_FILE        <- Sys.getenv("GENO_CSV",     "data/sample_genotype.csv")
NAVA_SEV_FILE    <- Sys.getenv("NAVA_SEV",     "outputs/nava_severity_labels.csv")
SPLIT_FILE       <- Sys.getenv("SPLIT_CSV",    "data/split_dev_holdout_v4_mixed.csv")

# Compute cache for the two slow derived columns (Delta-logit-PSI and the
# novel/annotated classification). Skipping these saves ~30-60s per run and
# is what makes plotting iteration fast. Invalidates when the upstream PSI
# matrices or SS cache change. DELETE THIS FILE if you change the
# event_is_novel() or compute_dlogit() logic; the cache is keyed on input
# mtimes, not on the script's hash.
EVENTS_CACHE_PATH <- Sys.getenv("EVENTS_CACHE", "outputs/figure1_events_cache.rds")

OUT_DIR_FIG     <- "outputs/figures"
OUT_SUMMARY     <- "outputs/figure1_summary.rds"
OUT_FIG_1       <- file.path(OUT_DIR_FIG, "figure1.pdf")

FDR_THRESH  <- 0.05
EVENT_TYPES <- c("SE", "A5SS", "A3SS", "MXE", "RI")

# Cohort plot labels. Internal tokens are GEL / Nava (database names); paper
# text + figures use Whole-Blood / Lymphocytes (the actual sample tissue).
# Code keeps the internal tokens; we remap only at the plot layer.
WB_LAB <- "whole blood"
LY_LAB <- "lymphocytes"
COHORT_LABELS <- c(GEL = WB_LAB, Nava = LY_LAB)

# Paper figure palette (set 2026-04-29). Colorblind-aware (Wong-style for
# upregulated/downregulated; Bang Wong CB-safe set for case-up/case-down).
PAL_WB    <- "#2C5AA0"  # Whole-Blood     (deep navy)
PAL_LY    <- "#2AA198"  # Lymphocytes     (teal / green-blue)
PAL_CONC  <- "#1F3A93"  # Concordant      (deep navy, panel C — the signal)
PAL_DISC  <- "#BBBBBB"  # Discordant      (light grey, panel C — the noise; recedes
                        #                  visually so the eye lands on concordance)
PAL_CTRL  <- "#888888"  # Ctrl / coef bar (mid grey)

stopifnot("GEL  PSI matrix not found" = file.exists(GEL_PSI_PATH),
          "Nava PSI matrix not found" = file.exists(NAVA_PSI_PATH),
          "Sample genotype CSV not found" = file.exists(GENO_FILE),
          "Nava severity labels file not found" = file.exists(NAVA_SEV_FILE),
          "Dev/holdout split CSV not found" = file.exists(SPLIT_FILE))
dir.create(OUT_DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_SUMMARY), showWarnings = FALSE, recursive = TRUE)

# ---- 2. Load PSI matrices and build shared event_key ------------------------
cat("Loading GEL  PSI matrix ...\n")
gel  <- readRDS(GEL_PSI_PATH)
cat("Loading Nava PSI matrix ...\n")
nava <- readRDS(NAVA_PSI_PATH)

stopifnot("GEL events missing coord_key/nostrand"  =
            all(c("coord_key", "coord_nostrand") %in% names(gel$events)),
          "Nava events missing coord_key/nostrand" =
            all(c("coord_key", "coord_nostrand") %in% names(nava$events)),
          "GEL  ijc/sjc/samples missing"  = all(c("ijc","sjc","samples") %in% names(gel)),
          "Nava ijc/sjc/samples missing"  = all(c("ijc","sjc","samples") %in% names(nava)),
          "mean_ctrl_depth missing GEL"   = "mean_ctrl_depth" %in% names(gel$events),
          "mean_ctrl_depth missing Nava"  = "mean_ctrl_depth" %in% names(nava$events))

gel$events$event_key  <- paste(gel$events$EventType,  gel$events$coord_nostrand, sep = "|")
nava$events$event_key <- paste(nava$events$EventType, nava$events$coord_nostrand, sep = "|")
stopifnot("GEL  event_key not unique" = !anyDuplicated(gel$events$event_key),
          "Nava event_key not unique" = !anyDuplicated(nava$events$event_key))

cat(sprintf("GEL  events: %s    Nava events: %s\n",
            format(nrow(gel$events),  big.mark = ","),
            format(nrow(nava$events), big.mark = ",")))

# ---- 3-4. Compute Delta-logit-PSI + novel/annotated classification ----------
# Both are cached together to EVENTS_CACHE_PATH because they're the only slow
# steps in the script (~30-60s combined). The cache stores per-cohort
# event_key + dlogit + is_novel and is invalidated when any of GEL_PSI,
# NAVA_PSI, or SS_CACHE_PATH is newer than the cache.
src_paths <- c(GEL_PSI_PATH, NAVA_PSI_PATH, SS_CACHE_PATH)
events_cache_valid <- file.exists(EVENTS_CACHE_PATH) &&
                      all(file.exists(src_paths)) &&
                      file.mtime(EVENTS_CACHE_PATH) > max(file.mtime(src_paths))

if (events_cache_valid) {
  cat("\nLoading cached classified events:", EVENTS_CACHE_PATH, "\n")
  cached <- readRDS(EVENTS_CACHE_PATH)
  stopifnot(
    "Events cache row count mismatch (delete EVENTS_CACHE_PATH and rerun)" =
      nrow(cached$gel)  == nrow(gel$events) &&
      nrow(cached$nava) == nrow(nava$events),
    "Events cache event_key mismatch (delete EVENTS_CACHE_PATH and rerun)" =
      identical(cached$gel$event_key,  gel$events$event_key) &&
      identical(cached$nava$event_key, nava$events$event_key))
  gel$events$dlogit    <- cached$gel$dlogit
  gel$events$is_novel  <- cached$gel$is_novel
  nava$events$dlogit   <- cached$nava$dlogit
  nava$events$is_novel <- cached$nava$is_novel
} else {
  # ---- 3. Compute per-event Delta-logit-PSI (case mean - ctrl mean) ---------
  # PSI = (IJC + 0.5) / (IJC + SJC + 1); logit(PSI) -- pseudocount-stabilized
  # (Beta(1,1) shrink). Same transform used by the insT/stemIII LASSO feature;
  # sign-equivalent to rMATS IncLevelDifference but better behaved at the rails.
  compute_dlogit <- function(obj) {
    ijc <- obj$ijc; sjc <- obj$sjc
    p <- (ijc + 0.5) / (ijc + sjc + 1)
    logit_psi <- log(p / (1 - p))
    case_idx <- which(obj$samples$group == "Case")
    ctrl_idx <- which(obj$samples$group == "Control")
    stopifnot("no Case samples"    = length(case_idx) > 0L,
              "no Control samples" = length(ctrl_idx) > 0L)
    case_mu <- colMeans(logit_psi[case_idx, , drop = FALSE], na.rm = TRUE)
    ctrl_mu <- colMeans(logit_psi[ctrl_idx, , drop = FALSE], na.rm = TRUE)
    out <- case_mu - ctrl_mu
    out[!is.finite(out)] <- NA_real_
    data.frame(event_id = colnames(ijc), dlogit = unname(out),
               stringsAsFactors = FALSE)
  }

  cat("\nComputing Delta-logit-PSI per event ...\n")
  dlogit_gel  <- compute_dlogit(gel)
  dlogit_nava <- compute_dlogit(nava)
  gel$events$dlogit  <- dlogit_gel$dlogit[ match(gel$events$event_id,  dlogit_gel$event_id) ]
  nava$events$dlogit <- dlogit_nava$dlogit[match(nava$events$event_id, dlogit_nava$event_id)]

  # ---- 4. Build / load annotated SS set from GTF --------------------------
  # An event is novel if ANY of its splice-site coords does not appear in the
  # annotated exon-boundary set on the same chromosome and strand.
  build_ss_set <- function(gtf_path) {
    cat("Parsing GTF (this may take a minute) ...\n")
    stopifnot("GTF not found" = file.exists(gtf_path))
    cmd <- sprintf("awk -F'\t' 'BEGIN{OFS=\"\\t\"} $3==\"exon\" {print $1,$4,$5,$7}' %s",
                   shQuote(gtf_path))
    exons <- read.delim(pipe(cmd), header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("chr", "start", "end", "strand"))
    cat(sprintf("  parsed %s exon rows\n", format(nrow(exons), big.mark = ",")))
    exons$start <- as.integer(exons$start)
    exons$end   <- as.integer(exons$end)
    exons$key   <- paste(exons$chr, exons$strand, sep = "|")
    pos_long <- data.frame(
      key = c(exons$key, exons$key),
      pos = c(exons$start, exons$end),
      stringsAsFactors = FALSE)
    pos_long <- pos_long[!is.na(pos_long$pos), , drop = FALSE]
    ss_set <- split(pos_long$pos, pos_long$key)
    ss_set <- lapply(ss_set, function(v) sort(unique(as.integer(v))))
    cat(sprintf("  built SS set across %d chr|strand groups\n", length(ss_set)))
    ss_set
  }

  if (file.exists(SS_CACHE_PATH) &&
      file.exists(GTF_PATH) &&
      file.mtime(SS_CACHE_PATH) > file.mtime(GTF_PATH)) {
    cat("\nLoading cached annotated-SS set: ", SS_CACHE_PATH, "\n", sep = "")
    ss_set <- readRDS(SS_CACHE_PATH)
  } else {
    ss_set <- build_ss_set(GTF_PATH)
    saveRDS(ss_set, SS_CACHE_PATH)
    cat("Cached SS set to: ", SS_CACHE_PATH, "\n", sep = "")
  }

  event_is_novel <- function(events, ss_set, event_type) {
    parts   <- strsplit(events$coord_key, ":", fixed = TRUE)
    chr     <- vapply(parts, `[`, character(1L), 1L)
    strand  <- vapply(parts, `[`, character(1L), 2L)
    n_coord <- if (event_type == "MXE") 8L else 6L
    coord_mat <- t(vapply(parts, function(p) {
      as.integer(p[3L:(2L + n_coord)])
    }, integer(n_coord)))
    start_idx <- seq(1L, n_coord, by = 2L)
    coord_mat[, start_idx] <- coord_mat[, start_idx] + 1L  # 0-based -> 1-based
    keys <- paste(chr, strand, sep = "|")
    is_novel <- logical(nrow(events))
    for (k in unique(keys)) {
      rows <- which(keys == k)
      annotated_pos <- ss_set[[k]]
      if (is.null(annotated_pos)) {
        is_novel[rows] <- TRUE
        next
      }
      sub_mat <- coord_mat[rows, , drop = FALSE]
      in_set  <- matrix(sub_mat %in% annotated_pos,
                        nrow = nrow(sub_mat), ncol = ncol(sub_mat))
      is_novel[rows] <- !apply(in_set, 1L, all)
    }
    is_novel
  }

  classify_cohort_events <- function(events, ss_set) {
    out <- logical(nrow(events))
    for (et in unique(events$EventType)) {
      sel <- events$EventType == et
      out[sel] <- event_is_novel(events[sel, , drop = FALSE], ss_set, et)
    }
    out
  }

  cat("\nClassifying GEL  events as novel/annotated ...\n")
  gel$events$is_novel  <- classify_cohort_events(gel$events,  ss_set)
  cat("Classifying Nava events as novel/annotated ...\n")
  nava$events$is_novel <- classify_cohort_events(nava$events, ss_set)

  saveRDS(list(
    gel  = data.frame(event_key = gel$events$event_key,
                      dlogit    = gel$events$dlogit,
                      is_novel  = gel$events$is_novel,
                      stringsAsFactors = FALSE),
    nava = data.frame(event_key = nava$events$event_key,
                      dlogit    = nava$events$dlogit,
                      is_novel  = nava$events$is_novel,
                      stringsAsFactors = FALSE)
  ), EVENTS_CACHE_PATH)
  cat(sprintf("Cached classified events to: %s\n", EVENTS_CACHE_PATH))
}

cat(sprintf("  GEL  dlogit non-NA: %s / %s; novel: %s\n",
            format(sum(!is.na(gel$events$dlogit)), big.mark = ","),
            format(nrow(gel$events),               big.mark = ","),
            format(sum(gel$events$is_novel),       big.mark = ",")))
cat(sprintf("  Nava dlogit non-NA: %s / %s; novel: %s\n",
            format(sum(!is.na(nava$events$dlogit)), big.mark = ","),
            format(nrow(nava$events),               big.mark = ","),
            format(sum(nava$events$is_novel),       big.mark = ",")))

# ---- 5. Event overlap (paragraph 1) -----------------------------------------
gel_keys  <- gel$events$event_key
nava_keys <- nava$events$event_key
shared    <- intersect(gel_keys, nava_keys)
gel_only  <- setdiff(gel_keys, nava_keys)
nava_only <- setdiff(nava_keys, gel_keys)

n_gel_total  <- length(gel_keys)
n_nava_total <- length(nava_keys)
n_shared     <- length(shared)
n_gel_only   <- length(gel_only)
n_nava_only  <- length(nava_only)

cat("\n=== Event totals (matched on EventType|coord_nostrand) ===\n")
cat(sprintf("GEL  events:   %s\n", format(n_gel_total,  big.mark = ",")))
cat(sprintf("Nava events:   %s\n", format(n_nava_total, big.mark = ",")))
cat(sprintf("Shared:        %s\n", format(n_shared,     big.mark = ",")))
cat(sprintf("GEL only:      %s\n", format(n_gel_only,   big.mark = ",")))
cat(sprintf("Nava only:     %s\n", format(n_nava_only,  big.mark = ",")))

gel_only_novel  <- sum(gel$events$is_novel[gel$events$event_key  %in% gel_only])
nava_only_novel <- sum(nava$events$is_novel[nava$events$event_key %in% nava_only])
n_unique_combined <- n_gel_only + n_nava_only
n_novel_combined  <- gel_only_novel + nava_only_novel
pct_novel         <- 100 * n_novel_combined / n_unique_combined

cat("\n=== Cohort-unique novel-SS summary (paragraph 1) ===\n")
cat(sprintf("Cohort-unique events (combined): %s\n",
            format(n_unique_combined, big.mark = ",")))
cat(sprintf("  involving novel SS:            %s (%.1f%%)\n",
            format(n_novel_combined, big.mark = ","), pct_novel))
cat(sprintf("    GEL-only novel:  %s / %s\n",
            format(gel_only_novel, big.mark = ","),
            format(n_gel_only,     big.mark = ",")))
cat(sprintf("    Nava-only novel: %s / %s\n",
            format(nava_only_novel, big.mark = ","),
            format(n_nava_only,     big.mark = ",")))

# ---- 6. Significance overlap among events quantified in both cohorts --------
# (paragraph 2). Restricted to the `shared` universe (events quantified in
# BOTH cohorts, section 5) so the panel-B significance Venn is restricted to
# events with a well-defined background rate in both cohorts. Events
# significant in one cohort but not quantified in the other are excluded.
sig_gel_keys  <- gel$events$event_key[!is.na(gel$events$FDR)  & gel$events$FDR  < FDR_THRESH &
                                      gel$events$event_key %in% shared]
sig_nava_keys <- nava$events$event_key[!is.na(nava$events$FDR) & nava$events$FDR < FDR_THRESH &
                                       nava$events$event_key %in% shared]
sig_shared    <- intersect(sig_gel_keys, sig_nava_keys)
sig_gel_only  <- setdiff(sig_gel_keys,  sig_nava_keys)
sig_nava_only <- setdiff(sig_nava_keys, sig_gel_keys)

n_sig_gel       <- length(sig_gel_keys)
n_sig_nava      <- length(sig_nava_keys)
n_sig_shared    <- length(sig_shared)
n_sig_gel_only  <- length(sig_gel_only)
n_sig_nava_only <- length(sig_nava_only)

cat(sprintf("\n=== Significant events among events quantified in both cohorts (FDR < %.2f) ===\n", FDR_THRESH))
cat(sprintf("Sig in GEL:        %s\n", format(n_sig_gel,       big.mark = ",")))
cat(sprintf("Sig in Nava:       %s\n", format(n_sig_nava,      big.mark = ",")))
cat(sprintf("Sig in both:       %s\n", format(n_sig_shared,    big.mark = ",")))
cat(sprintf("Sig GEL only:      %s\n", format(n_sig_gel_only,  big.mark = ",")))
cat(sprintf("Sig Nava only:     %s\n", format(n_sig_nava_only, big.mark = ",")))
if (n_sig_gel  > 0) cat(sprintf("  Sig-in-both / Sig-in-GEL  = %.1f%%\n", 100 * n_sig_shared / n_sig_gel))
if (n_sig_nava > 0) cat(sprintf("  Sig-in-both / Sig-in-Nava = %.1f%%\n", 100 * n_sig_shared / n_sig_nava))

# ---- 7. Joined event table (sig + Delta-logit + concordance) ----------------
gel_sub <- gel$events[, c("event_key", "event_id", "EventType", "GeneSymbol",
                          "FDR", "IncLevelDiff", "dlogit", "mean_ctrl_depth")]
nava_cols <- c("event_key", "event_id", "FDR", "IncLevelDiff", "dlogit", "mean_ctrl_depth")
if ("PValue" %in% names(nava$events)) nava_cols <- c(nava_cols, "PValue")
nava_sub <- nava$events[, nava_cols]

joined <- merge(gel_sub, nava_sub, by = "event_key",
                all = TRUE, suffixes = c("_GEL", "_Nava"))
stopifnot("Duplicate event_key after join" = !anyDuplicated(joined$event_key))
if ("PValue" %in% names(joined)) names(joined)[names(joined) == "PValue"] <- "PValue_Nava"

# Fill EventType for Nava-only rows (parsed from event_key prefix)
nav_only_idx <- is.na(joined$EventType)
if (any(nav_only_idx)) {
  joined$EventType[nav_only_idx] <- sub("^([^|]+)\\|.*$", "\\1",
                                        joined$event_key[nav_only_idx])
}

joined$in_GEL   <- !is.na(joined$FDR_GEL)
joined$in_Nava  <- !is.na(joined$FDR_Nava)
joined$shared   <- joined$in_GEL & joined$in_Nava
joined$sig_GEL  <- joined$in_GEL  & !is.na(joined$FDR_GEL)  & joined$FDR_GEL  < FDR_THRESH
joined$sig_Nava <- joined$in_Nava & !is.na(joined$FDR_Nava) & joined$FDR_Nava < FDR_THRESH
joined$sig_both <- joined$sig_GEL & joined$sig_Nava
joined$concordant <- joined$sig_both &
  !is.na(joined$dlogit_GEL) & !is.na(joined$dlogit_Nava) &
  sign(joined$dlogit_GEL) == sign(joined$dlogit_Nava) &
  sign(joined$dlogit_GEL) != 0

# ---- 8. Per-event-type summary table ----------------------------------------
summary_df <- do.call(rbind, lapply(EVENT_TYPES, function(et) {
  sub <- joined[!is.na(joined$EventType) & joined$EventType == et, ]
  data.frame(
    EventType     = et,
    GEL_total     = sum(sub$in_GEL),
    Nava_total    = sum(sub$in_Nava),
    shared        = sum(sub$shared),
    sig_GEL_only  = sum(sub$sig_GEL & !sub$sig_Nava),
    sig_Nava_only = sum(sub$sig_Nava & !sub$sig_GEL),
    sig_both      = sum(sub$sig_both),
    concordant    = sum(sub$concordant),
    discordant    = sum(sub$sig_both & !sub$concordant),
    stringsAsFactors = FALSE)
}))
summary_total <- data.frame(
  EventType     = "TOTAL",
  GEL_total     = sum(summary_df$GEL_total),
  Nava_total    = sum(summary_df$Nava_total),
  shared        = sum(summary_df$shared),
  sig_GEL_only  = sum(summary_df$sig_GEL_only),
  sig_Nava_only = sum(summary_df$sig_Nava_only),
  sig_both      = sum(summary_df$sig_both),
  concordant    = sum(summary_df$concordant),
  discordant    = sum(summary_df$discordant),
  stringsAsFactors = FALSE)
summary_df <- rbind(summary_df, summary_total)

cat("\n=== Per event type summary ===\n")
print(summary_df, row.names = FALSE)

concordance_rate <- if (summary_total$sig_both > 0)
  summary_total$concordant / summary_total$sig_both else NA_real_
cat(sprintf("\nOverall concordance among sig-in-both: %.3f (%d / %d)\n",
            concordance_rate, summary_total$concordant, summary_total$sig_both))

# ---- 9. Cross-cohort replication metrics (binomial concordance) ------------
n_shared_total <- sum(joined$shared)
a <- summary_total$sig_both
b <- summary_total$sig_GEL_only
c <- summary_total$sig_Nava_only
binom_both   <- binom.test(summary_total$concordant, summary_total$sig_both,
                           p = 0.5, alternative = "greater")

scatter_df <- joined[joined$shared &
                     (joined$sig_GEL | joined$sig_Nava) &
                     !is.na(joined$dlogit_GEL) &
                     !is.na(joined$dlogit_Nava), ]
scatter_df$sig_group <- with(scatter_df, ifelse(
  sig_both,
  ifelse(concordant, "Sig in both (concordant)", "Sig in both (discordant)"),
  ifelse(sig_GEL, "Sig GEL only", "Sig Nava only")))

cat("\n=== Cross-cohort replication metrics ===\n")
cat(sprintf("Shared events tested:     %d\n", n_shared_total))
cat(sprintf("Sig in GEL  (FDR<%.2f):   %d  (%.2f%% of shared)\n",
            FDR_THRESH, a + b, 100 * (a + b) / n_shared_total))
cat(sprintf("Sig in Nava (FDR<%.2f):   %d  (%.2f%% of shared)\n",
            FDR_THRESH, a + c, 100 * (a + c) / n_shared_total))
cat(sprintf("Sig in both:              %d\n", a))

# Pearson/Spearman of Delta-logit-PSI at sig-in-both events
sig_both_df <- scatter_df[scatter_df$sig_both, ]
if (nrow(sig_both_df) >= 3L) {
  pearson_r  <- cor(sig_both_df$dlogit_GEL, sig_both_df$dlogit_Nava, use = "complete.obs")
  spearman_r <- cor(sig_both_df$dlogit_GEL, sig_both_df$dlogit_Nava, use = "complete.obs",
                    method = "spearman")
} else { pearson_r <- NA_real_; spearman_r <- NA_real_ }

cat("\n--- Effect-size correlation at sig-in-both (Delta-logit-PSI) ---\n")
cat(sprintf("Pearson r:   %.3f  R^2 = %.3f  (n = %d)\n",
            pearson_r, pearson_r^2, nrow(sig_both_df)))
cat(sprintf("Spearman rho: %.3f                     (n = %d)\n",
            spearman_r, nrow(sig_both_df)))

# ---- 9b. Sig-in-both |dPSI| magnitude (paragraph 3) ------------------------
# Effect-size magnitudes for A5SS events significant in both cohorts
# (FDR < 0.05). Uses raw IncLevelDifference (PSI scale, 0-1) for reader-
# friendly units. Restricted to A5SS; not restricted to direction-concordant
# events (concordance is reported separately).
sigboth_dpsi <- joined[joined$sig_both &
                       joined$EventType == "A5SS" &
                       !is.na(joined$IncLevelDiff_GEL) &
                       !is.na(joined$IncLevelDiff_Nava), ]
n_sigboth_dpsi <- nrow(sigboth_dpsi)
abs_dpsi_gel  <- abs(sigboth_dpsi$IncLevelDiff_GEL)
abs_dpsi_nava <- abs(sigboth_dpsi$IncLevelDiff_Nava)
max_dpsi_gel  <- if (n_sigboth_dpsi > 0) max(abs_dpsi_gel,  na.rm = TRUE) else NA_real_
max_dpsi_nava <- if (n_sigboth_dpsi > 0) max(abs_dpsi_nava, na.rm = TRUE) else NA_real_
q90_dpsi_gel  <- if (n_sigboth_dpsi > 0) unname(quantile(abs_dpsi_gel,  0.90, na.rm = TRUE)) else NA_real_
q90_dpsi_nava <- if (n_sigboth_dpsi > 0) unname(quantile(abs_dpsi_nava, 0.90, na.rm = TRUE)) else NA_real_
median_dpsi_gel  <- if (n_sigboth_dpsi > 0) median(abs_dpsi_gel,  na.rm = TRUE) else NA_real_
median_dpsi_nava <- if (n_sigboth_dpsi > 0) median(abs_dpsi_nava, na.rm = TRUE) else NA_real_

cat("\n=== Sig-in-both A5SS |dPSI| magnitude (paragraph 3) ===\n")
cat(sprintf("n events:               %d\n", n_sigboth_dpsi))
cat(sprintf("max |dPSI| GEL:         %.3f\n", max_dpsi_gel))
cat(sprintf("max |dPSI| Nava:        %.3f\n", max_dpsi_nava))
cat(sprintf("median |dPSI| GEL:      %.3f\n", median_dpsi_gel))
cat(sprintf("median |dPSI| Nava:     %.3f\n", median_dpsi_nava))
cat(sprintf("90th-pct |dPSI| GEL:    %.3f\n", q90_dpsi_gel))
cat(sprintf("90th-pct |dPSI| Nava:   %.3f\n", q90_dpsi_nava))

# Identify the max-|dPSI| event in each cohort and report the paired values,
# so the prose can distinguish per-cohort maxes from a same-event paired max.
# EventType is printed so the prose can state whether the max is an A5SS.
if (n_sigboth_dpsi > 0) {
  idx_max_gel  <- which.max(abs_dpsi_gel)
  idx_max_nava <- which.max(abs_dpsi_nava)
  same_event   <- idx_max_gel == idx_max_nava

  cat("\n--- Max-|dPSI| events (paired values) ---\n")
  cat(sprintf("Max in GEL:  %s | %s | %s | GEL |dPSI|=%.3f, Nava |dPSI|=%.3f\n",
              sigboth_dpsi$GeneSymbol[idx_max_gel],
              sigboth_dpsi$EventType[idx_max_gel],
              sigboth_dpsi$event_key[idx_max_gel],
              abs_dpsi_gel[idx_max_gel],  abs_dpsi_nava[idx_max_gel]))
  cat(sprintf("Max in Nava: %s | %s | %s | GEL |dPSI|=%.3f, Nava |dPSI|=%.3f\n",
              sigboth_dpsi$GeneSymbol[idx_max_nava],
              sigboth_dpsi$EventType[idx_max_nava],
              sigboth_dpsi$event_key[idx_max_nava],
              abs_dpsi_gel[idx_max_nava], abs_dpsi_nava[idx_max_nava]))
  cat(sprintf("Same event achieves max in both cohorts? %s\n",
              if (same_event) "YES" else "NO"))

  # Paired-max alternative: largest single-event change that is large in BOTH
  # cohorts (i.e. the bigger one is small). Use min(|dPSI|_GEL, |dPSI|_Nava)
  # as the per-event "shared magnitude" and take its max.
  shared_mag <- pmin(abs_dpsi_gel, abs_dpsi_nava)
  idx_shared <- which.max(shared_mag)
  cat(sprintf("Largest paired-max event (max of min across cohorts):\n"))
  cat(sprintf("  %s | %s | %s | GEL |dPSI|=%.3f, Nava |dPSI|=%.3f\n",
              sigboth_dpsi$GeneSymbol[idx_shared],
              sigboth_dpsi$EventType[idx_shared],
              sigboth_dpsi$event_key[idx_shared],
              abs_dpsi_gel[idx_shared], abs_dpsi_nava[idx_shared]))
}

# ---- 11. Save unified summary RDS -------------------------------------------
saveRDS(list(
  fdr_threshold     = FDR_THRESH,
  # paragraph 1: overall overlap + novel SS
  n_gel_total       = n_gel_total,
  n_nava_total      = n_nava_total,
  n_shared          = n_shared,
  n_gel_only        = n_gel_only,
  n_nava_only       = n_nava_only,
  n_unique_combined = n_unique_combined,
  n_novel_combined  = n_novel_combined,
  pct_novel         = pct_novel,
  gel_only_novel    = gel_only_novel,
  nava_only_novel   = nava_only_novel,
  # paragraph 2: significance overlap
  n_sig_gel         = n_sig_gel,
  n_sig_nava        = n_sig_nava,
  n_sig_shared      = n_sig_shared,
  n_sig_gel_only    = n_sig_gel_only,
  n_sig_nava_only   = n_sig_nava_only,
  # per-type summary + scatter table for plot reuse
  summary           = summary_df,
  scatter           = scatter_df[, c("event_key", "event_id_GEL", "EventType", "GeneSymbol",
                                     "FDR_GEL", "FDR_Nava",
                                     "IncLevelDiff_GEL", "IncLevelDiff_Nava",
                                     "dlogit_GEL", "dlogit_Nava",
                                     "sig_group", "concordant")],
  # statistical tests
  concordance_rate    = concordance_rate,
  pearson_r           = pearson_r,
  pearson_r2          = pearson_r^2,
  spearman_r          = spearman_r,
  concordance_p       = binom_both$p.value,
  # paragraph 3: sig-in-both |dPSI| magnitude (reported, not plotted as a panel)
  n_sigboth_dpsi        = n_sigboth_dpsi,
  max_dpsi_gel_sigboth  = max_dpsi_gel,
  max_dpsi_nava_sigboth = max_dpsi_nava,
  median_dpsi_gel_sigboth  = median_dpsi_gel,
  median_dpsi_nava_sigboth = median_dpsi_nava,
  q90_dpsi_gel_sigboth  = q90_dpsi_gel,
  q90_dpsi_nava_sigboth = q90_dpsi_nava
), OUT_SUMMARY)
cat(sprintf("\nSummary saved to: %s\n", OUT_SUMMARY))

# ---- 12. Build panel B in memory (proportional Euler, sig events only) ------
# Single small Euler diagram of sig-in-both vs cohort-only sig events. The
# all-events overlap is reported numerically in paragraph 1 and not plotted
# (kept the figure compact). FDR threshold + total n live in the figure
# caption rather than as an in-panel title. Ellipses keep eulerr's default
# horizontal (side-by-side) layout.
quantity_label <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
cohort_fills   <- c(PAL_WB, PAL_LY)  # Whole-Blood navy, Lymphocytes teal

fit_sig <- euler(setNames(
  c(n_sig_gel_only, n_sig_nava_only, n_sig_shared),
  c(WB_LAB, LY_LAB, paste(WB_LAB, LY_LAB, sep = "&"))))

p_sig <- plot(fit_sig,
  fills      = list(fill = cohort_fills, alpha = 0.55),
  labels     = list(font = 1, cex = 0.9),
  edges      = list(col = "grey20", lwd = 1),
  quantities = list(type = "counts", cex = 0.75,
                    labels = vapply(c(n_sig_gel_only, n_sig_nava_only, n_sig_shared),
                                    quantity_label, character(1L))))

# ---- Shared publication theme (theme_classic: axis lines only, no grid) -----
pub_theme <- theme_classic(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 9),
        axis.text        = element_text(colour = "black", size = 9),
        axis.title       = element_text(size = 9),
        legend.position  = "top",
        legend.margin    = margin(0, 0, 0, 0),
        legend.box.margin= margin(-5, 0, -5, 0),
        legend.key.size  = unit(0.3, "cm"),
        legend.text      = element_text(size = 9),
        plot.margin      = margin(4, 6, 4, 4))

# HGVS + train/holdout split lookup -- authoritative, model-independent
# sources for per-case genotype identity and dev/holdout membership (per_case
# below is built from these, not from any *_loo$per_case_meta, which only
# covers that model's own training population).
geno     <- read.csv(GENO_FILE, stringsAsFactors = FALSE)
nava_lbl <- local({
  first <- readLines(NAVA_SEV_FILE, n = 1L)
  sep   <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
  read.table(NAVA_SEV_FILE, header = TRUE, sep = sep, stringsAsFactors = FALSE,
             quote = "\"", fill = TRUE)
})
nava_lbl$sample_id <- paste0("nava:case:", nava_lbl$sample_n)
hgvs_lookup <- setNames(c(geno$hgvs,      nava_lbl$hgvs),
                         c(geno$sample_id, nava_lbl$sample_id))

split_df <- read.csv(SPLIT_FILE, stringsAsFactors = FALSE)
stopifnot("Split CSV missing required fields" =
            all(c("sample_id", "cohort", "set") %in% names(split_df)),
          "Split CSV has unexpected set values" =
            all(split_df$set %in% c("train", "holdout")))
set_lookup <- setNames(split_df$set, split_df$sample_id)

gel_case_ids  <- gel$samples$sample_id[gel$samples$group  == "Case"]
nava_case_ids <- nava$samples$sample_id[nava$samples$group == "Case"]
per_case <- data.frame(
  sample_id = c(gel_case_ids, nava_case_ids),
  cohort    = c(rep("GEL",  length(gel_case_ids)),
                rep("Nava", length(nava_case_ids))),
  stringsAsFactors = FALSE)
per_case$hgvs <- unname(hgvs_lookup[per_case$sample_id])
per_case$set  <- unname(set_lookup[per_case$sample_id])
stopifnot("Some cases missing HGVS from hgvs_lookup" = !anyNA(per_case$hgvs),
          "Some cases missing train/holdout set"      = !anyNA(per_case$set))

# ---- 13b. Build panel C in memory (PCA on sig-in-both PSI) -----------------
# PCA of all samples (GEL + Nava, cases + controls) on the PSI values of the
# n_sig_shared events significant in both cohorts (FDR < 0.05). Colour =
# cohort (Whole-Blood / Lymphocytes); shape = triangle for cases, circle for
# controls. Restricted to complete cases (no imputation). Expect tissue (cohort)
# to dominate PC1/PC2; cases should cluster with their own-cohort controls.
{
  gel_col_sig  <- gel$events$event_id [match(sig_shared, gel$events$event_key)]
  nava_col_sig <- nava$events$event_id[match(sig_shared, nava$events$event_key)]
  stopifnot("sig_shared->GEL  event_id map has NAs" = !anyNA(gel_col_sig),
            "sig_shared->Nava event_id map has NAs" = !anyNA(nava_col_sig),
            "GEL  psi col missing for sig_shared"   = all(gel_col_sig  %in% colnames(gel$psi)),
            "Nava psi col missing for sig_shared"   = all(nava_col_sig %in% colnames(nava$psi)))

  gel_psi_sig  <- gel$psi [, gel_col_sig,  drop = FALSE]
  nava_psi_sig <- nava$psi[, nava_col_sig, drop = FALSE]
  colnames(gel_psi_sig) <- colnames(nava_psi_sig) <- sig_shared

  pca_sig_mat <- rbind(gel_psi_sig, nava_psi_sig)
  pca_sig_meta <- rbind(
    data.frame(sample_id = gel$samples$sample_id,  cohort = "GEL",
               group = gel$samples$group,  stringsAsFactors = FALSE),
    data.frame(sample_id = nava$samples$sample_id, cohort = "Nava",
               group = nava$samples$group, stringsAsFactors = FALSE))
  stopifnot("PCA sig matrix / meta row mismatch" = nrow(pca_sig_mat) == nrow(pca_sig_meta))

  keep_sig   <- colSums(is.na(pca_sig_mat)) == 0L
  pca_sig_cc <- pca_sig_mat[, keep_sig, drop = FALSE]
  v_sig      <- apply(pca_sig_cc, 2L, var)
  pca_sig_cc <- pca_sig_cc[, v_sig > 0, drop = FALSE]
  cat(sprintf("\nPanel C PCA: %d sig-in-both events; %d complete+nonzero-var after filtering\n",
              length(sig_shared), ncol(pca_sig_cc)))
  stopifnot("Too few complete sig-in-both events for PCA panel C" = ncol(pca_sig_cc) >= 10L)

  pca_sig      <- prcomp(pca_sig_cc, center = TRUE, scale. = FALSE)
  pct_sig      <- 100 * pca_sig$sdev^2 / sum(pca_sig$sdev^2)

  pca_sig_meta$PC1        <- pca_sig$x[, 1L]
  pca_sig_meta$PC2        <- pca_sig$x[, 2L]
  pca_sig_meta$cohort_lab <- factor(unname(COHORT_LABELS[pca_sig_meta$cohort]),
                                    levels = c(WB_LAB, LY_LAB))
  pca_sig_meta$status     <- factor(ifelse(pca_sig_meta$group == "Case",
                                           "ReNU", "Control"),
                                    levels = c("Control", "ReNU"))

  pca_sig_meta$panel_title <- sprintf("%d shared FDR < %.2f events",
                                      n_sig_shared, FDR_THRESH)

  fig1c <- ggplot(pca_sig_meta, aes(x = PC1, y = PC2,
                                     colour = cohort_lab, shape = status)) +
    geom_point(size = 1.8, alpha = 0.85, stroke = 0.5) +
    scale_colour_manual(values = setNames(c(PAL_WB, PAL_LY),
                                          c(WB_LAB, LY_LAB)),
                        name = NULL) +
    scale_shape_manual(values = c(Control = 16L, `ReNU` = 17L),
                       name = NULL) +
    facet_wrap(~ panel_title) +
    labs(x = sprintf("PC1 (%.1f%%)", pct_sig[1L]),
         y = sprintf("PC2 (%.1f%%)", pct_sig[2L])) +
    pub_theme +
    theme(legend.position = "right",
          legend.key.size = unit(3, "mm"),
          legend.text     = element_text(size = 9),
          # Bottom margin to align with panel D's bottom legend (~18pt);
          # top margin no longer needed now that C has its own strip.
          plot.margin     = margin(4, 4, 18, 4))
}

# ---- 13c. Per-event-type Delta-logit-PSI agreement (stats only, not plotted) -
# Per-event-type correlation of Delta-logit-PSI (Whole-Blood vs Lymphocytes)
# at sig-in-both events -- printed to stdout for context; the actual Panel D
# scatter below covers all cases, A5SS only.
stopifnot("Too few sig-in-both events for stats" = nrow(sig_both_df) >= 10L)
sig_both_df$EventType <- factor(sig_both_df$EventType, levels = EVENT_TYPES)

panel_cor <- do.call(rbind, lapply(EVENT_TYPES, function(et) {
  sub <- sig_both_df[sig_both_df$EventType == et, ]
  r   <- if (nrow(sub) >= 3L)
    cor(sub$dlogit_GEL, sub$dlogit_Nava, use = "complete.obs") else NA_real_
  rho <- if (nrow(sub) >= 3L)
    cor(sub$dlogit_GEL, sub$dlogit_Nava, use = "complete.obs",
        method = "spearman") else NA_real_
  data.frame(EventType = et, n = nrow(sub),
             pearson  = r, r_squared = if (is.na(r)) NA_real_ else r^2,
             spearman = rho, stringsAsFactors = FALSE)
}))
panel_cor$EventType <- factor(panel_cor$EventType, levels = EVENT_TYPES)
panel_cor$label <- ifelse(
  is.na(panel_cor$pearson),
  sprintf("n=%d", panel_cor$n),
  sprintf("R^2=%.2f\nn=%d", panel_cor$r_squared, panel_cor$n))

cat("\n--- Per-event-type Delta-logit-PSI agreement (sig-in-both) ---\n")
print(panel_cor[, c("EventType", "n", "pearson", "r_squared", "spearman")],
      row.names = FALSE)

# ---- 13f. Panel D (all-cases A5SS dlogit) prep -------------------------
# Mirrors the A5SS facet of figure_s3.R: same sig-in-both event set,
# same bb-logit PSI transform (Methods form), same all-ReNU-cases-vs-all-
# controls dlogit, same concordance rule. Restricted to A5SS.
PANEL_D_EVT_TYPES <- c("A5SS")
INST_HGVS         <- "n.64_65insT"

sigmoid <- function(z) 1 / (1 + exp(-z))

{
  # ---- Carrier counts for lollipop (panel A) --------------------------------
  d_df <- per_case[!is.na(per_case$hgvs) & nzchar(per_case$hgvs),
                   c("sample_id", "cohort", "set", "hgvs"), drop = FALSE]
  d_df$cohort_lab <- factor(unname(COHORT_LABELS[d_df$cohort]),
                            levels = c(WB_LAB, LY_LAB))
  d_df$set <- factor(d_df$set, levels = c("train", "holdout"),
                     labels = c("Training", "Holdout"))

  carrier_counts <- aggregate(sample_id ~ set + cohort_lab + hgvs,
                              data = d_df, FUN = length)
  names(carrier_counts)[ncol(carrier_counts)] <- "n"

  # ---- 13d. Panel D: delta-logit-PSI, all ReNU cases (A5SS) ------------------
  # Mirrors figure_s3.R's A5SS facet exactly:
  #   - sample set: all ReNU cases vs all controls, per cohort
  #   - event set: sig-in-both (FDR<0.05 in both), restricted to A5SS
  #   - PSI transform: bb-logit  logit((IJC+0.5)/(IJC+SJC+1))  (Methods form)
  #   - dlogit: mean(case bb-logit-PSI) - mean(ctrl bb-logit-PSI), per cohort
  #   - concordance: sign(dl_GEL) == sign(dl_Nava) & sign(dl_GEL) != 0, non-NA
  {
    gel_case_ids_pd  <- gel$samples$sample_id[ gel$samples$group  == "Case"]
    nava_case_ids_pd <- nava$samples$sample_id[nava$samples$group == "Case"]
    gel_ctrl_ids     <- gel$samples$sample_id[ gel$samples$group  == "Control"]
    nava_ctrl_ids    <- nava$samples$sample_id[nava$samples$group == "Control"]

    stopifnot("No GEL  ReNU cases"  = length(gel_case_ids_pd)  > 0L,
              "No Nava ReNU cases"  = length(nava_case_ids_pd) > 0L)
    cat(sprintf("Panel D all ReNU cases: GEL=%d  Nava=%d\n",
                length(gel_case_ids_pd), length(nava_case_ids_pd)))

    pd_sub  <- sig_both_df[!is.na(sig_both_df$EventType) &
                            sig_both_df$EventType %in% PANEL_D_EVT_TYPES, ]
    gel_ev  <- pd_sub$event_id_GEL
    nava_ev <- pd_sub$event_id_Nava

    bb_logit_pd <- function(ijc, sjc) {
      p <- (ijc + 0.5) / (ijc + sjc + 1)
      log(p / (1 - p))
    }
    gel_lp  <- bb_logit_pd(gel$ijc[,  gel_ev,  drop = FALSE],
                            gel$sjc[,  gel_ev,  drop = FALSE])
    nava_lp <- bb_logit_pd(nava$ijc[, nava_ev, drop = FALSE],
                            nava$sjc[, nava_ev, drop = FALSE])

    dl_gel  <- colMeans(gel_lp[ gel_case_ids_pd,  , drop = FALSE], na.rm = TRUE) -
               colMeans(gel_lp[ gel_ctrl_ids,     , drop = FALSE], na.rm = TRUE)
    dl_nava <- colMeans(nava_lp[nava_case_ids_pd, , drop = FALSE], na.rm = TRUE) -
               colMeans(nava_lp[nava_ctrl_ids,    , drop = FALSE], na.rm = TRUE)
    dl_gel[ !is.finite(dl_gel)]  <- NA_real_
    dl_nava[!is.finite(dl_nava)] <- NA_real_

    pd_df <- data.frame(
      EventType   = pd_sub$EventType,
      dlogit_GEL  = dl_gel,
      dlogit_Nava = dl_nava,
      concordant  = !is.na(dl_gel) & !is.na(dl_nava) &
                    sign(dl_gel) == sign(dl_nava) & sign(dl_gel) != 0L,
      stringsAsFactors = FALSE)
    pd_df <- pd_df[!is.na(pd_df$dlogit_GEL) & !is.na(pd_df$dlogit_Nava), ]
    facet_lbl <- PANEL_D_EVT_TYPES
    pd_df$facet_label <- factor(facet_lbl, levels = facet_lbl)

    pd_ann <- do.call(rbind, lapply(facet_lbl, function(fl) {
      sub <- pd_df[as.character(pd_df$facet_label) == fl, ]
      n   <- nrow(sub)
      r   <- if (n >= 3L) cor(sub$dlogit_GEL, sub$dlogit_Nava,
                              use = "complete.obs") else NA_real_
      data.frame(facet_label = fl,
                 label = if (!is.na(r)) sprintf("R^2=%.2f\nn=%d", r^2, n)
                         else           sprintf("n=%d", n),
                 stringsAsFactors = FALSE)
    }))
    pd_ann$facet_label <- factor(pd_ann$facet_label, levels = facet_lbl)

    ax_pd <- ceiling(max(abs(c(pd_df$dlogit_GEL, pd_df$dlogit_Nava)),
                         na.rm = TRUE) * 10L) / 10L

    fig1d <- ggplot(pd_df, aes(x = dlogit_GEL, y = dlogit_Nava)) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
      geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                  colour = "grey40", linewidth = 0.3) +
      geom_point(aes(colour = concordant), size = 1.0, alpha = 0.8) +
      geom_text(data = pd_ann, aes(x = -Inf, y = Inf, label = label),
                hjust = -0.1, vjust = 1.2, size = 2.4, colour = "black",
                lineheight = 0.9, inherit.aes = FALSE) +
      scale_colour_manual(
        values = c(`TRUE` = "#444444", `FALSE` = PAL_DISC),
        labels = c(`TRUE` = "Concordant", `FALSE` = "Discordant"),
        name   = NULL) +
      facet_wrap(~ facet_label, nrow = 1L) +
      coord_cartesian(xlim = c(-ax_pd, ax_pd), ylim = c(-ax_pd, ax_pd)) +
      labs(x = sprintf("delta logit-PSI (%s)", WB_LAB),
           y = sprintf("delta logit-PSI (%s)", LY_LAB)) +
      pub_theme +
      theme(legend.position  = "bottom",
            legend.key.size  = unit(3, "mm"),
            legend.text      = element_text(size = 9))
  }

  # ---- 13a. Panel A: bar plot of RNU4-2 variant carriers by HGVS ----------
  # Total carriers (train + holdout combined) per cohort per HGVS variant.
  # One stacked bar per variant (x = HGVS, categorical, ordered by nt
  # position), fill = cohort (Whole-Blood / Lymphocytes). Faceted by
  # structural region (T-loop | Other | Stem III) instead of shaded
  # background rects, so the region label comes from the facet strip and
  # gridlines stay visible across the whole panel.
  # Region coordinates consistent with figure2.R / figure_s3.R definitions.
  TLOOP_START    <- 63L; TLOOP_END    <- 70L
  STEM_III_START <- 75L; STEM_III_END <- 79L
  REGION_LEVELS  <- c("T-loop", "Other", "Stem III")

  bar_counts <- aggregate(n ~ cohort_lab + hgvs,
                          data = carrier_counts, FUN = sum)
  bar_counts$pos <- as.numeric(
    regmatches(bar_counts$hgvs,
               regexpr("[0-9]+", bar_counts$hgvs)))
  bar_counts$region <- with(bar_counts, ifelse(
    pos >= TLOOP_START & pos <= TLOOP_END, "T-loop",
    ifelse(pos >= STEM_III_START & pos <= STEM_III_END, "Stem III", "Other")))
  bar_counts$region <- factor(bar_counts$region, levels = REGION_LEVELS)

  hgvs_order <- unique(bar_counts$hgvs[order(bar_counts$pos)])
  bar_counts$hgvs <- factor(bar_counts$hgvs, levels = hgvs_order)

  # Facet strip label = region + total individuals in that region (across all
  # its variants and cohorts).
  region_totals <- aggregate(n ~ region, data = bar_counts, FUN = sum)
  region_lab_map <- setNames(
    sprintf("%s\n(n=%d)", region_totals$region, region_totals$n),
    as.character(region_totals$region))
  bar_counts$region_lab <- factor(
    unname(region_lab_map[as.character(bar_counts$region)]),
    levels = region_lab_map[REGION_LEVELS])

  # Bar-top labels: total individuals (both cohorts) per variant.
  bar_totals <- aggregate(n ~ hgvs + region_lab, data = bar_counts, FUN = sum)

  cat(sprintf("\nPanel A bar: %d variant x cohort rows; positions %d-%d\n",
              nrow(bar_counts),
              min(bar_counts$pos), max(bar_counts$pos)))
  cat("Panel A variants by region:\n")
  print(table(region = bar_counts$region[!duplicated(bar_counts$hgvs)]))

  fig1a <- ggplot(bar_counts, aes(x = hgvs, y = n, fill = cohort_lab)) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(data = bar_totals, aes(x = hgvs, y = n, label = n),
              inherit.aes = FALSE, vjust = -0.3, size = 2.8) +
    facet_grid(~ region_lab, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = setNames(c(PAL_WB, PAL_LY), c(WB_LAB, LY_LAB)),
                      name = NULL) +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "Individuals (n)",
                       breaks = scales::breaks_pretty(),
                       expand = expansion(mult = c(0, 0.12))) +
    pub_theme +
    theme(legend.position    = "bottom",
          legend.key.size    = unit(3, "mm"),
          legend.text        = element_text(size = 9),
          axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
          panel.grid.major   = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.margin        = margin(4, 6, 4, 4))

}

# ---- 14. Combine panels into Figure 1 PDF (patchwork) -----------------------
# patchwork handles axis alignment natively across panels in the same row.
# Nested row layout (row / row / row) allows different column widths per row.
# The eulerr trellis (p_sig) is wrapped with wrap_elements() so patchwork can
# place it alongside ggplots. Panel tags use labs(tag=) + plot.tag theme.
mm_to_in <- function(mm) mm / 25.4

tag_theme <- theme(
  plot.tag          = element_text(size = 14, face = "bold"),
  plot.tag.position = "topleft")

# eulerr plot() returns a gtable-like trellis; wrap_elements lets patchwork
# treat it like any other panel.
p_B <- wrap_elements(full = p_sig) + labs(tag = "B")
p_C <- fig1c + labs(tag = "C") + tag_theme

p_A <- fig1a     + labs(tag = "A") + tag_theme
p_D <- fig1d     + labs(tag = "D") + tag_theme

# Layout: 2 rows, 170 x 150 mm
#   Row 1 = A | B        (variant carrier bar plot | Euler sig-event Venn)
#   Row 2 = C | D        (PCA on sig-in-both PSI | insT dlogit scatter)
row1 <- (p_A | p_B) + plot_layout(widths = c(1, 1))
row2 <- (p_C | p_D) + plot_layout(widths = c(1, 1))

fig1_pw <- (row1 / row2) +
  plot_layout(heights = c(0.65, 0.67))

pdf(OUT_FIG_1, width = mm_to_in(170), height = mm_to_in(150))
print(fig1_pw & theme(plot.tag = element_text(size = 14, face = "bold"),
                      plot.tag.position = "topleft"))
dev.off()
cat(sprintf("Figure 1 saved to: %s\n", OUT_FIG_1))

cat("\n=== Done. Headline numbers ===\n")
cat(sprintf("Para 1: %s shared / %s GEL-only / %s Nava-only events; %.1f%% of cohort-unique involve novel SS\n",
            format(n_shared,    big.mark = ","),
            format(n_gel_only,  big.mark = ","),
            format(n_nava_only, big.mark = ","),
            pct_novel))
cat(sprintf("Para 2: %s sig in GEL, %s sig in Nava, %s sig in both; concordance %.1f%% (binom p = %.2e)\n",
            format(n_sig_gel,    big.mark = ","),
            format(n_sig_nava,   big.mark = ","),
            format(n_sig_shared, big.mark = ","),
            100 * concordance_rate, binom_both$p.value))
cat(sprintf("        Delta-logit-PSI Pearson R^2 = %.2f, Spearman rho = %.2f\n",
            pearson_r^2, spearman_r))
cat(sprintf("Para 3: sig-in-both n = %d; max |dPSI| GEL = %.3f / Nava = %.3f; 90th-pct GEL = %.3f / Nava = %.3f\n",
            n_sigboth_dpsi, max_dpsi_gel, max_dpsi_nava, q90_dpsi_gel, q90_dpsi_nava))
cat("\nNOTE: GEL-derived counts and figures stay in workspace unless airlocked.\n")
