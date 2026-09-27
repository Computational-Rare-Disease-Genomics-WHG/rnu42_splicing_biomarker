#!/usr/bin/env Rscript
#
# combine_rmats_batches_cand_pool.R
# Extract PSI + IJC + SJC at the full A5SS concordant candidate pool from all
# rMATS batch outputs and combine into sample x event matrices for the ~6k GEL
# whole-blood controls.
#
# The candidate pool = A5SS events significant (FDR < 0.05) in BOTH GEL and
# Nava rMATS and direction-concordant — currently ~483 events. This is the
# universe from which all current and future signatures are derived, so any
# signature score is a column-subset of this matrix.
#
# Pool is recomputed from psi_matrix.rds + nava_psi_matrix.rds event metadata
# (coord_nostrand, FDR, IncLevelDiff) — no LOO RDS needed.
#
# Run from workspace root: Rscript scripts/gel/combine_rmats_batches_cand_pool.R
#
# Inputs:
#   outputs/psi_matrix.rds           - dev GEL PSI (events$FDR, events$IncLevelDiff)
#   outputs/nava_psi_matrix.rds      - Nava PSI (same schema)
#   outputs/rmats_controls/batch_N/  - per-batch rMATS JC outputs
#   rmats_tmp/batch_N/merged/bam_list_for_post.txt
#   data/control_manifest.csv
#
# Output: outputs/gel_controls_psi_cand_pool.rds
#   $psi     : matrix [n_samples x n_pool_events]   raw PSI (IncLevel1)
#   $ijc     : matrix [n_samples x n_pool_events]   inclusion junction counts
#   $sjc     : matrix [n_samples x n_pool_events]   skipping junction counts
#   $samples : data.frame  sample_id (platekey), participant_id, spliceopathy*
#   $events  : data.frame  event_id, EventType, coord_nostrand, GeneSymbol,
#                          direction, gel_fdr, nava_fdr
#   $n_batches_combined : integer
#   $batch_dirs_used    : character
#
# DATA PROTECTION: event identifiers derived from GEL-protected rMATS output.
# Output stays inside the Research Environment unless approved via airlock.

local({
  v <- substr(as.character(getRversion()), 1, 3)
  candidate <- switch(
    v,
    "4.0" = "/tools/aws-workspace-apps/ce/R/4.0.2",
    "4.3" = "/tools/aws-workspace-ubuntu-apps/ce/R/4.3.3",
    "4.5" = "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2",
    NULL
  )
  if (!is.null(candidate) && dir.exists(candidate)) {
    .libPaths(c(.libPaths(), candidate))
  }
})

# Base R only - readr / tidyverse have bit64 namespace issues on GEL.

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
WORKSPACE    <- "/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
BATCH_DIR    <- file.path(WORKSPACE, "outputs/rmats_controls")
RMATS_TMP    <- file.path(WORKSPACE, "rmats_tmp")
BAM_LISTS    <- file.path(WORKSPACE, "data")
MANIFEST     <- file.path(WORKSPACE, "data/control_manifest.csv")
GEL_PSI_PATH <- file.path(WORKSPACE, "outputs/psi_matrix.rds")
NAVA_PSI_PATH<- file.path(WORKSPACE, "outputs/nava_psi_matrix.rds")
OUT_PATH     <- file.path(WORKSPACE, "outputs/gel_controls_psi_cand_pool.rds")

EVENT_TYPE  <- "A5SS"
FDR_THRESH  <- 0.05

# Batches 1-61 are the original QC-pass samples. Set to NULL to include all.
BATCH_RANGE <- c(1L, 61L)

stopifnot(
  "GEL PSI matrix not found"  = file.exists(GEL_PSI_PATH),
  "Nava PSI matrix not found" = file.exists(NAVA_PSI_PATH)
)

# --------------------------------------------------------------------------
# Reconstruct the concordant A5SS candidate pool from event metadata.
# This mirrors the candidate-space logic in loo_evaluate_insT_lasso.R /
# loo_evaluate_stemIII_lasso.R, but runs from the pre-built events frames so
# the combine script has no dependency on any LOO RDS.
# --------------------------------------------------------------------------
message("Loading PSI matrices...")
gel_full    <- readRDS(GEL_PSI_PATH)
nava_full   <- readRDS(NAVA_PSI_PATH)
gel_events  <- gel_full$events
nava_events <- nava_full$events

stopifnot(
  all(c("event_id", "EventType", "coord_nostrand",
        "FDR", "IncLevelDiff", "GeneSymbol") %in% names(gel_events)),
  all(c("event_id", "EventType", "coord_nostrand",
        "FDR", "IncLevelDiff") %in% names(nava_events)),
  all(c("ijc", "sjc", "samples") %in% names(gel_full)),
  all(c("ijc", "sjc", "samples") %in% names(nava_full))
)

gel_sig  <- gel_events[!is.na(gel_events$FDR)   & gel_events$FDR   < FDR_THRESH &
                        gel_events$EventType == EVENT_TYPE, ]
nava_sig <- nava_events[!is.na(nava_events$FDR) & nava_events$FDR  < FDR_THRESH &
                         nava_events$EventType == EVENT_TYPE, ]

both_keys <- intersect(gel_sig$coord_nostrand, nava_sig$coord_nostrand)

gel_m  <- gel_sig[ match(both_keys, gel_sig$coord_nostrand),  ]
nava_m <- nava_sig[match(both_keys, nava_sig$coord_nostrand), ]

# Concordance uses bb-logit dlogit-PSI (all cases vs all controls, per
# cohort) -- matches figure1.R panel D / figure_s1.R and the LOO scripts.
# Direction sign is taken from dlogit_GEL to keep pool orientation on the
# same metric as the LASSO feature downstream.
bb_logit_pool <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}
gel_case_all  <- gel_full$samples$sample_id[ gel_full$samples$group  == "Case"]
gel_ctrl_all  <- gel_full$samples$sample_id[ gel_full$samples$group  == "Control"]
nava_case_all <- nava_full$samples$sample_id[nava_full$samples$group == "Case"]
nava_ctrl_all <- nava_full$samples$sample_id[nava_full$samples$group == "Control"]

gel_lp_pool  <- bb_logit_pool(gel_full$ijc[,  gel_m$event_id,  drop = FALSE],
                               gel_full$sjc[,  gel_m$event_id,  drop = FALSE])
nava_lp_pool <- bb_logit_pool(nava_full$ijc[, nava_m$event_id, drop = FALSE],
                               nava_full$sjc[, nava_m$event_id, drop = FALSE])
dl_gel_pool  <- colMeans(gel_lp_pool[ gel_case_all,  , drop = FALSE], na.rm = TRUE) -
                colMeans(gel_lp_pool[ gel_ctrl_all,  , drop = FALSE], na.rm = TRUE)
dl_nava_pool <- colMeans(nava_lp_pool[nava_case_all, , drop = FALSE], na.rm = TRUE) -
                colMeans(nava_lp_pool[nava_ctrl_all, , drop = FALSE], na.rm = TRUE)
dl_gel_pool[ !is.finite(dl_gel_pool)]  <- NA_real_
dl_nava_pool[!is.finite(dl_nava_pool)] <- NA_real_

concordant <- !is.na(dl_gel_pool) & !is.na(dl_nava_pool) &
              sign(dl_gel_pool) == sign(dl_nava_pool) &
              sign(dl_gel_pool) != 0

sig <- data.frame(
  event_id       = gel_m$event_id[concordant],
  EventType      = EVENT_TYPE,
  coord_nostrand = gel_m$coord_nostrand[concordant],
  GeneSymbol     = gel_m$GeneSymbol[concordant],
  direction      = sign(unname(dl_gel_pool[concordant])),
  gel_fdr        = gel_m$FDR[concordant],
  nava_fdr       = nava_m$FDR[concordant],
  stringsAsFactors = FALSE
)

stopifnot(
  "No concordant A5SS events found" = nrow(sig) > 0L,
  "Duplicate coord_nostrand in pool" = !anyDuplicated(sig$coord_nostrand),
  "Direction has zeros" = !any(sig$direction == 0)
)
message("Concordant A5SS candidate pool: ", nrow(sig), " events")

# --------------------------------------------------------------------------
# Discover batch directories (sorted numerically)
# --------------------------------------------------------------------------
batch_dirs <- list.dirs(BATCH_DIR, recursive = FALSE, full.names = TRUE)
batch_dirs <- batch_dirs[grepl("/batch_\\d+$", batch_dirs)]
batch_nums <- as.integer(sub(".*batch_", "", batch_dirs))
batch_dirs <- batch_dirs[order(batch_nums)]
batch_nums <- sort(batch_nums)
if (!is.null(BATCH_RANGE)) {
  keep       <- batch_nums >= BATCH_RANGE[1L] & batch_nums <= BATCH_RANGE[2L]
  batch_dirs <- batch_dirs[keep]
  batch_nums <- batch_nums[keep]
  message("Restricting to batches ", BATCH_RANGE[1L], "-", BATCH_RANGE[2L],
          " (", sum(keep), " of ", length(keep), " found)")
}
stopifnot("No batch directories found under BATCH_DIR" = length(batch_dirs) > 0L)
message("\nProcessing ", length(batch_dirs), " batch directories (batches ",
        min(batch_nums), " to ", max(batch_nums), ")")

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
parse_inc <- function(col, n) {
  t(sapply(col, function(x) {
    vals <- suppressWarnings(as.numeric(strsplit(x, ",")[[1]]))
    if (length(vals) != n) rep(NA_real_, n) else vals
  }))
}
parse_counts <- function(col, n) {
  t(sapply(col, function(x) {
    vals <- suppressWarnings(as.integer(strsplit(x, ",")[[1]]))
    if (length(vals) != n) rep(NA_integer_, n) else vals
  }))
}

platekey_from_path <- function(path) sub("\\.bam$", "", basename(path))

build_coord_nostrand <- function(df, event_type) {
  n_coord   <- if (event_type == "MXE") 8L else 6L
  coord_key <- do.call(paste, c(df[, seq(4L, 5L + n_coord)], sep = ":"))
  sub("^([^:]+):[^:]+:", "\\1:", coord_key)
}

# --------------------------------------------------------------------------
# Process each batch, extracting PSI / IJC / SJC at pool events
# --------------------------------------------------------------------------
batch_psi   <- list()
batch_ijc   <- list()
batch_sjc   <- list()
batch_samps <- list()

for (i in seq_along(batch_dirs)) {
  bd   <- batch_dirs[i]
  bnum <- batch_nums[i]

  batch_tmp <- file.path(RMATS_TMP, paste0("batch_", bnum), "merged")
  bam_file  <- file.path(batch_tmp, "bam_list_for_post.txt")
  if (!file.exists(bam_file)) {
    stop("bam_list_for_post.txt not found in: ", batch_tmp,
         "\n  Re-run post step before combining.")
  }
  batch_bam_paths <- unlist(strsplit(readLines(bam_file), ","))
  batch_bam_paths <- trimws(batch_bam_paths[nzchar(trimws(batch_bam_paths))])
  n_batch         <- length(batch_bam_paths)
  batch_platekeys <- platekey_from_path(batch_bam_paths)

  orig_bam_file <- file.path(BAM_LISTS, paste0("bam_paths_batch_", bnum, ".txt"))
  if (file.exists(orig_bam_file)) {
    expected_n <- length(readLines(orig_bam_file))
    if (n_batch < expected_n) {
      stop("Batch ", bnum, ": expected ", expected_n,
           " samples but bam_list_for_post.txt has ", n_batch,
           " - ", expected_n - n_batch, " prep job(s) failed.")
    }
  }

  fname <- file.path(bd, paste0(EVENT_TYPE, ".MATS.JC.txt"))
  if (!file.exists(fname)) {
    message("  Batch ", bnum, ": ", EVENT_TYPE, " output missing - filling with NA")
    psi_mat <- matrix(NA_real_,    nrow = nrow(sig), ncol = n_batch)
    ijc_mat <- matrix(NA_integer_, nrow = nrow(sig), ncol = n_batch)
    sjc_mat <- matrix(NA_integer_, nrow = nrow(sig), ncol = n_batch)
    rownames(psi_mat) <- rownames(ijc_mat) <- rownames(sjc_mat) <- sig$event_id
    colnames(psi_mat) <- colnames(ijc_mat) <- colnames(sjc_mat) <- batch_platekeys
    batch_psi[[i]]   <- psi_mat
    batch_ijc[[i]]   <- ijc_mat
    batch_sjc[[i]]   <- sjc_mat
    batch_samps[[i]] <- batch_platekeys
    next
  }

  df <- read.delim(fname, stringsAsFactors = FALSE)

  n_in_output <- length(strsplit(df$IncLevel1[1], ",")[[1]])
  if (n_in_output != n_batch) {
    stop("Batch ", bnum, " ", EVENT_TYPE, ": IncLevel1 has ", n_in_output,
         " values but BAM list has ", n_batch, " entries - mismatch")
  }

  coord_nostrand_vec <- build_coord_nostrand(df, EVENT_TYPE)
  match_idx <- match(sig$coord_nostrand, coord_nostrand_vec)
  n_found   <- sum(!is.na(match_idx))

  psi_mat <- matrix(NA_real_,    nrow = nrow(sig), ncol = n_batch)
  ijc_mat <- matrix(NA_integer_, nrow = nrow(sig), ncol = n_batch)
  sjc_mat <- matrix(NA_integer_, nrow = nrow(sig), ncol = n_batch)

  if (n_found > 0L) {
    found_rows <- which(!is.na(match_idx))
    hits       <- match_idx[found_rows]
    psi_mat[found_rows, ] <- parse_inc(   df$IncLevel1[hits],    n_batch)
    ijc_mat[found_rows, ] <- parse_counts(df$IJC_SAMPLE_1[hits], n_batch)
    sjc_mat[found_rows, ] <- parse_counts(df$SJC_SAMPLE_1[hits], n_batch)
  }

  rownames(psi_mat) <- rownames(ijc_mat) <- rownames(sjc_mat) <- sig$event_id
  colnames(psi_mat) <- colnames(ijc_mat) <- colnames(sjc_mat) <- batch_platekeys
  batch_psi[[i]]   <- psi_mat
  batch_ijc[[i]]   <- ijc_mat
  batch_sjc[[i]]   <- sjc_mat
  batch_samps[[i]] <- batch_platekeys

  message("  Batch ", bnum, ": ", n_found, "/", nrow(sig),
          " pool events found, ", n_batch, " samples")
}

# --------------------------------------------------------------------------
# Combine batches: events x samples, then transpose to samples x events
# --------------------------------------------------------------------------
psi_es <- do.call(cbind, batch_psi)
ijc_es <- do.call(cbind, batch_ijc)
sjc_es <- do.call(cbind, batch_sjc)

stopifnot(
  identical(dim(psi_es),     dim(ijc_es)),
  identical(dim(psi_es),     dim(sjc_es)),
  identical(rownames(psi_es), rownames(ijc_es)),
  identical(colnames(psi_es), colnames(ijc_es))
)

all_sample_ids <- colnames(psi_es)
stopifnot(
  "Duplicate sample IDs in combined matrix" = !anyDuplicated(all_sample_ids),
  "event_ids in matrix don't match pool"    = setequal(rownames(psi_es), sig$event_id)
)

all_batch_samps <- unlist(batch_samps)
if (anyDuplicated(all_batch_samps)) {
  dup_ids <- all_batch_samps[duplicated(all_batch_samps)]
  stop("Duplicate sample IDs across batches: ", paste(unique(dup_ids), collapse = ", "))
}

psi_mat <- t(psi_es)
ijc_mat <- t(ijc_es)
sjc_mat <- t(sjc_es)

message("\n--- Combined matrix: ", nrow(psi_mat), " samples x ",
        ncol(psi_mat), " events ---")

# --------------------------------------------------------------------------
# Sample metadata
# --------------------------------------------------------------------------
manifest <- read.csv(MANIFEST, stringsAsFactors = FALSE)
stopifnot("file_path column missing from manifest" = "file_path" %in% names(manifest))
manifest$platekey <- platekey_from_path(manifest$file_path)

keep_cols <- intersect(
  c("platekey", "participant_id",
    "spliceopathy_flag", "spliceopathy_evidence", "spliceopathy_gene"),
  names(manifest)
)
manifest_slim <- manifest[, keep_cols, drop = FALSE]

unmatched <- setdiff(all_sample_ids, manifest_slim$platekey)
if (length(unmatched) > 0L) {
  warning(length(unmatched), " sample(s) not found in manifest: ",
          paste(head(unmatched, 5L), collapse = ", "),
          if (length(unmatched) > 5L) " ..." else "")
}

samples_meta <- data.frame(sample_id = all_sample_ids, stringsAsFactors = FALSE)
samples_meta <- merge(samples_meta, manifest_slim,
                      by.x = "sample_id", by.y = "platekey",
                      all.x = TRUE, sort = FALSE)
samples_meta <- samples_meta[match(rownames(psi_mat), samples_meta$sample_id), ]
stopifnot(identical(samples_meta$sample_id, rownames(psi_mat)))

# --------------------------------------------------------------------------
# Event metadata — preserve direction and FDRs from the pool data.frame
# --------------------------------------------------------------------------
events_meta <- sig[match(colnames(psi_mat), sig$event_id),
                   c("event_id", "EventType", "coord_nostrand",
                     "GeneSymbol", "direction", "gel_fdr", "nava_fdr"),
                   drop = FALSE]
stopifnot(identical(events_meta$event_id, colnames(psi_mat)))

# --------------------------------------------------------------------------
# Final assertions
# --------------------------------------------------------------------------
stopifnot(
  nrow(psi_mat) == nrow(samples_meta),
  ncol(psi_mat) == nrow(events_meta),
  !anyDuplicated(rownames(psi_mat)),
  !anyDuplicated(colnames(psi_mat))
)

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
out <- list(
  psi                = psi_mat,
  ijc                = ijc_mat,
  sjc                = sjc_mat,
  samples            = samples_meta,
  events             = events_meta,
  n_batches_combined = length(batch_dirs),
  batch_dirs_used    = batch_dirs
)
saveRDS(out, OUT_PATH)

message("\n=== Done ===")
message("Saved to: ", OUT_PATH)
message("Shape:    ", nrow(psi_mat), " samples x ", ncol(psi_mat), " events")
message("PSI coverage: ", round(mean(!is.na(psi_mat)) * 100L, 1L), "% non-NA")
message("IJC+SJC > 0:  ",
        round(mean((ijc_mat + sjc_mat) > 0L, na.rm = TRUE) * 100L, 1L), "% of cells")
if ("spliceopathy_flag" %in% names(samples_meta)) {
  message("\nSpliceopathy flag breakdown:")
  print(table(samples_meta$spliceopathy_flag, useNA = "ifany"))
}
message("\nPer-event coverage (frac. samples with non-NA PSI):")
print(round(summary(colMeans(!is.na(psi_mat))), 3L))
message("\nPer-event median depth (IJC+SJC):")
depth <- ijc_mat + sjc_mat
print(round(summary(apply(depth, 2L, median, na.rm = TRUE)), 1L))
