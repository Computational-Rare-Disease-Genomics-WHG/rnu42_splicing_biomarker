#!/usr/bin/env Rscript
#
# build_nava_psi_matrix.R
# Build a Nava PSI matrix from rMATS JC output, mirroring build_psi_matrix.R
# (which does the same for the GEL 11-case + 33-control rMATS run).
#
# Output schema is identical to outputs/psi_matrix.rds so downstream scripts
# (loo_evaluate_v4.R, figure-1 comparison scripts, direction concordance)
# can read both cohorts with the same code path.
#
# Inputs:
#   data/nava_rmats_tables/{SE,A5SS,A3SS,MXE,RI}.MATS.JC.txt
#
# Output: outputs/nava_psi_matrix.rds — list with:
#   $psi     : numeric matrix  [samples × events] point-estimate PSI
#                              (rows = nava:case:1..19, nava:control:1..21)
#   $ijc     : integer matrix  [samples × events] inclusion junction counts
#   $sjc     : integer matrix  [samples × events] skipping  junction counts
#   $events  : data.frame      event_id, EventType, GeneSymbol, chr, strand,
#                              coord_key, coord_nostrand, PValue, FDR,
#                              IncLevelDiff, mean_ctrl_depth, n_samples_with_psi
#   $samples : data.frame      sample_id, group (Case/Control), case_n
#
# IJC/SJC are required for count-level beta-binomial PSI shrinkage in v4
# (PSI_SCALE = "logit_bb"). For point-estimate PSI consumers, the $psi
# matrix has the same shape and semantics as before this revision.
#
# Run from workspace root: Rscript scripts/gel/build_nava_psi_matrix.R

.libPaths(c(.libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
library(tidyverse)

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
RMATS_DIR  <- "data/nava_rmats_tables"
N_CASES    <- 19L
N_CONTROLS <- 21L
OUT_PATH   <- "outputs/nava_psi_matrix.rds"

EVENT_TYPES <- c("SE", "A5SS", "A3SS", "MXE", "RI")

case_ids <- paste0("nava:case:",    seq_len(N_CASES))
ctrl_ids <- paste0("nava:control:", seq_len(N_CONTROLS))
all_ids  <- c(case_ids, ctrl_ids)

stopifnot("rMATS directory not found" = dir.exists(RMATS_DIR))
dir.create(dirname(OUT_PATH), showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------
# Helper: parse one rMATS JC file into an events × samples PSI data frame.
# Same logic as build_psi_matrix.R — MXE has 8 coord cols, others 6.
# --------------------------------------------------------------------------
parse_rmats_psi <- function(event_type, n_cases, n_controls) {
  fname <- file.path(RMATS_DIR, paste0(event_type, ".MATS.JC.txt"))
  if (!file.exists(fname)) {
    message("Skipping ", event_type, " — file not found: ", fname)
    return(NULL)
  }

  df <- read.delim(fname, stringsAsFactors = FALSE)
  message(event_type, ": ", nrow(df), " events")

  n_coord        <- if (event_type == "MXE") 8L else 6L
  coord_key      <- do.call(paste, c(df[, seq(4L, 5L + n_coord)], sep = ":"))
  coord_nostrand <- sub("^([^:]+):[^:]+:", "\\1:", coord_key)
  gene_sym       <- df[[3]]
  event_id       <- paste(event_type, gene_sym, coord_nostrand, sep = ":")

  dups <- duplicated(event_id)
  if (any(dups)) {
    message("  Removing ", sum(dups), " duplicate event_ids")
    df             <- df[!dups, ]
    coord_key      <- coord_key[!dups]
    coord_nostrand <- coord_nostrand[!dups]
    event_id       <- event_id[!dups]
    gene_sym       <- gene_sym[!dups]
  }

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

  case_mat <- parse_inc(df$IncLevel1, n_cases)
  ctrl_mat <- parse_inc(df$IncLevel2, n_controls)

  # Per-sample inclusion / skipping junction counts (cases AND controls).
  # ctrl-side aggregates (mean_ctrl_depth, min_ctrl_depth) feed v3/v4 filters;
  # the per-sample matrices are also saved for count-level beta-binomial PSI
  # shrinkage in v4 (PSI_SCALE = "logit_bb").
  ijc_case        <- parse_counts(df$IJC_SAMPLE_1, n_cases)
  sjc_case        <- parse_counts(df$SJC_SAMPLE_1, n_cases)
  ijc_ctrl        <- parse_counts(df$IJC_SAMPLE_2, n_controls)
  sjc_ctrl        <- parse_counts(df$SJC_SAMPLE_2, n_controls)
  ctrl_depth      <- ijc_ctrl + sjc_ctrl
  mean_ctrl_depth <- rowMeans(ctrl_depth, na.rm = TRUE)
  min_ctrl_depth  <- suppressWarnings(apply(ctrl_depth, 1L, min, na.rm = TRUE))
  min_ctrl_depth[!is.finite(min_ctrl_depth)] <- NA_integer_

  psi_wide <- as.data.frame(cbind(case_mat, ctrl_mat))
  colnames(psi_wide) <- all_ids

  ijc_combined <- cbind(ijc_case, ijc_ctrl)
  sjc_combined <- cbind(sjc_case, sjc_ctrl)
  colnames(ijc_combined) <- paste0("ijc__", all_ids)
  colnames(sjc_combined) <- paste0("sjc__", all_ids)
  psi_wide <- cbind(psi_wide,
                    as.data.frame(ijc_combined),
                    as.data.frame(sjc_combined))

  psi_wide$mean_ctrl_depth <- mean_ctrl_depth
  psi_wide$min_ctrl_depth  <- min_ctrl_depth
  psi_wide$event_id        <- event_id
  psi_wide$EventType       <- event_type
  psi_wide$GeneSymbol      <- gene_sym
  psi_wide$chr             <- df[[4]]
  psi_wide$strand          <- df[[5]]
  psi_wide$coord_key       <- coord_key
  psi_wide$coord_nostrand  <- coord_nostrand
  psi_wide$PValue          <- suppressWarnings(as.numeric(df$PValue))
  psi_wide$FDR             <- suppressWarnings(as.numeric(df$FDR))
  psi_wide$IncLevelDiff    <- suppressWarnings(as.numeric(df$IncLevelDifference))

  psi_wide
}

# --------------------------------------------------------------------------
# Load all event types and combine
# --------------------------------------------------------------------------
message("Building Nava PSI matrix...\n")

all_events <- map(EVENT_TYPES, parse_rmats_psi,
                  n_cases = N_CASES, n_controls = N_CONTROLS) %>%
  compact() %>%
  bind_rows()

message("\nTotal events: ", nrow(all_events))
stopifnot(
  "No events loaded — check RMATS_DIR and file names" = nrow(all_events) > 0,
  "Duplicate event_ids across event types"            = !anyDuplicated(all_events$event_id)
)

# --------------------------------------------------------------------------
# Filter: drop events with no PSI in any sample
# --------------------------------------------------------------------------
psi_only <- as.matrix(all_events[, all_ids])
ijc_only <- as.matrix(all_events[, paste0("ijc__", all_ids)])
sjc_only <- as.matrix(all_events[, paste0("sjc__", all_ids)])
has_data <- rowSums(!is.na(psi_only)) > 0

message("Events with >=1 non-NA PSI: ", sum(has_data), " / ", nrow(all_events))
message("Dropping ", sum(!has_data), " fully-missing events")

all_events <- all_events[has_data, ]
psi_only   <- psi_only[has_data, ]
ijc_only   <- ijc_only[has_data, ]
sjc_only   <- sjc_only[has_data, ]

# --------------------------------------------------------------------------
# Build output structures
# --------------------------------------------------------------------------
psi_matrix <- t(psi_only); rownames(psi_matrix) <- all_ids; colnames(psi_matrix) <- all_events$event_id
ijc_matrix <- t(ijc_only); rownames(ijc_matrix) <- all_ids; colnames(ijc_matrix) <- all_events$event_id
sjc_matrix <- t(sjc_only); rownames(sjc_matrix) <- all_ids; colnames(sjc_matrix) <- all_events$event_id

events_meta <- all_events %>%
  select(event_id, EventType, GeneSymbol, chr, strand,
         coord_key, coord_nostrand, PValue, FDR, IncLevelDiff,
         mean_ctrl_depth, min_ctrl_depth) %>%
  mutate(n_samples_with_psi = colSums(!is.na(psi_matrix)))

samples_meta <- data.frame(
  sample_id = all_ids,
  group     = c(rep("Case",    N_CASES),
                rep("Control", N_CONTROLS)),
  case_n    = c(seq_len(N_CASES), rep(NA_integer_, N_CONTROLS)),
  stringsAsFactors = FALSE
)

stopifnot(
  "PSI matrix rows != samples" = nrow(psi_matrix) == nrow(samples_meta),
  "PSI matrix cols != events"  = ncol(psi_matrix) == nrow(events_meta),
  "IJC matrix shape mismatch"  = identical(dim(ijc_matrix), dim(psi_matrix)),
  "SJC matrix shape mismatch"  = identical(dim(sjc_matrix), dim(psi_matrix)),
  "IJC dimnames mismatch"      = identical(dimnames(ijc_matrix), dimnames(psi_matrix)),
  "SJC dimnames mismatch"      = identical(dimnames(sjc_matrix), dimnames(psi_matrix))
)

out <- list(psi = psi_matrix, ijc = ijc_matrix, sjc = sjc_matrix,
            events = events_meta, samples = samples_meta)
saveRDS(out, OUT_PATH)

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
message("\n=== Nava PSI matrix built ===")
message("Saved to: ", OUT_PATH)
message("Dimensions: ", nrow(psi_matrix), " samples × ",
        ncol(psi_matrix), " events")
message("Missingness: ",
        round(mean(is.na(psi_matrix)) * 100, 1), "% NA")
message("\nEvent type breakdown:")
print(table(events_meta$EventType))
message("\nSample breakdown:")
print(table(samples_meta$group))
