#!/usr/bin/env Rscript
#
# check_batch_clean_availability.R
# For the "null LASSO within a batch" robustness experiment: report per batch
# how many samples remain after applying the figure_s4 clean-unmatched filter,
# and identify batches with >= 44 clean samples (11 pseudo-case + 33 pseudo-ctrl).
#
# Clean-unmatched pool (matches figure_s4.R lines 236-317):
#   batch samples MINUS (dev overlap ∪ biallelic RNU4-2 ∪ HD ReNU case)
#
# Run from workspace root:
#   Rscript scripts/gel/check_batch_clean_availability.R
#
# Output: prints a per-batch table and writes
#   outputs/batch_clean_availability.csv
#     batch, n_total, n_clean, has_44

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

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
WORKSPACE     <- "/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
RMATS_TMP     <- file.path(WORKSPACE, "rmats_tmp")
DEV_CASE_BAMS <- file.path(WORKSPACE, "data/rmats/inputs/renu_11_cases_paths.txt")
DEV_CTRL_BAMS <- file.path(WORKSPACE, "data/rmats/inputs/renu_11_cases_control_paths.txt")
BIALLELIC_FILE<- file.path(WORKSPACE, "data/biallelic_rnu42_platekeys.txt")
HD_CASE_FILE  <- file.path(WORKSPACE, "data/renu_highdepth_case_path.txt")
OUT_CSV       <- file.path(WORKSPACE, "outputs/batch_clean_availability.csv")

BATCH_RANGE   <- c(1L, 61L)     # QC-pass batches only
N_REQUIRED    <- 44L            # 11 pseudo-case + 33 pseudo-ctrl

# --------------------------------------------------------------------------
# Helpers (copied verbatim from figure_s4.R for identical semantics)
# --------------------------------------------------------------------------
read_bam_list <- function(path) {
  if (!file.exists(path)) return(character(0))
  raw  <- paste(readLines(path), collapse = ",")
  bams <- trimws(unlist(strsplit(raw, ",")))
  bams <- bams[nzchar(bams)]
  sub("\\.bam$", "", basename(bams))
}

platekey_from_path <- function(path) sub("\\.bam$", "", basename(path))

# --------------------------------------------------------------------------
# Build the exclusion set
# --------------------------------------------------------------------------
dev_case_keys <- read_bam_list(DEV_CASE_BAMS)
dev_ctrl_keys <- read_bam_list(DEV_CTRL_BAMS)
dev_all_keys  <- c(dev_case_keys, dev_ctrl_keys)

biallelic_keys <- if (file.exists(BIALLELIC_FILE)) {
  trimws(readLines(BIALLELIC_FILE))
} else character(0)
biallelic_keys <- biallelic_keys[nzchar(biallelic_keys)]

hd_keys <- read_bam_list(HD_CASE_FILE)

excl_keys <- unique(c(dev_all_keys, biallelic_keys, hd_keys))

cat("Exclusion set composition:\n")
cat(sprintf("  dev cases      : %d\n", length(dev_case_keys)))
cat(sprintf("  dev controls   : %d\n", length(dev_ctrl_keys)))
cat(sprintf("  biallelic RNU42: %d\n", length(biallelic_keys)))
cat(sprintf("  HD ReNU case   : %d\n", length(hd_keys)))
cat(sprintf("  TOTAL (unique) : %d\n\n", length(excl_keys)))

# --------------------------------------------------------------------------
# Discover batches
# --------------------------------------------------------------------------
batch_nums <- BATCH_RANGE[1L]:BATCH_RANGE[2L]

results <- data.frame(
  batch   = integer(0),
  n_total = integer(0),
  n_clean = integer(0),
  has_44  = logical(0),
  stringsAsFactors = FALSE
)

for (bnum in batch_nums) {
  bam_file <- file.path(RMATS_TMP, paste0("batch_", bnum),
                        "merged", "bam_list_for_post.txt")
  if (!file.exists(bam_file)) {
    warning("Batch ", bnum, ": bam_list_for_post.txt missing, skipping")
    next
  }
  bam_paths <- unlist(strsplit(readLines(bam_file), ","))
  bam_paths <- trimws(bam_paths[nzchar(trimws(bam_paths))])
  platekeys <- platekey_from_path(bam_paths)

  n_total <- length(platekeys)
  n_clean <- sum(!platekeys %in% excl_keys)

  results <- rbind(
    results,
    data.frame(batch   = bnum,
               n_total = n_total,
               n_clean = n_clean,
               has_44  = n_clean >= N_REQUIRED,
               stringsAsFactors = FALSE)
  )
}

# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
write.csv(results, OUT_CSV, row.names = FALSE)

cat("Per-batch availability:\n")
print(results, row.names = FALSE)

qualifying <- results$batch[results$has_44]
cat(sprintf("\nBatches with >= %d clean samples: %d/%d\n",
            N_REQUIRED, length(qualifying), nrow(results)))
if (length(qualifying) > 0L) {
  cat(sprintf("First qualifying batch: %d (n_clean = %d)\n",
              qualifying[1L],
              results$n_clean[results$batch == qualifying[1L]]))
  cat(sprintf("All qualifying batches: %s\n",
              paste(qualifying, collapse = ", ")))
} else {
  cat("NO batch has enough clean samples - experiment infeasible as designed.\n")
}

cat(sprintf("\nSummary of n_clean across %d batches:\n", nrow(results)))
print(summary(results$n_clean))

cat(sprintf("\nWrote: %s\n", OUT_CSV))
