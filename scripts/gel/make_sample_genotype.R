#!/usr/bin/env Rscript
#
# make_sample_genotype.R
# Build data/sample_genotype.csv by matching rMATS b1 BAM path order
# to participant genotype information.
#
# Inputs:
#   data/renu_11_cases_paths.txt   — b1 file passed to rMATS (one BAM per line,
#                                    or comma-separated; order = IncLevel1 order)
#   data/participant_genotype.tsv  — TSV with columns: participant_id, platekey,
#                                    sample_group (HGVS genotype string)
#
# Output:
#   data/sample_genotype.csv       — sample_id (gel:case:N), hgvs
#
# Run from workspace root: Rscript scripts/make_sample_genotype.R

.libPaths(c(.libPaths(), "/tools/aws-workspace-apps/ce/R/4.0.2"))
library(tidyverse)

# --------------------------------------------------------------------------
# Configuration — adjust filenames if needed
# --------------------------------------------------------------------------
PATHS_FILE    <- "data/renu_11_cases_paths.txt"
GENOTYPE_FILE <- "data/12_samples_with_renu_and_RNAseq_for_RD.tsv"
OUT_FILE      <- "data/sample_genotype.csv"
N_CASES       <- 11L

stopifnot(
  "BAM paths file not found"    = file.exists(PATHS_FILE),
  "Genotype TSV not found"      = file.exists(GENOTYPE_FILE)
)

# --------------------------------------------------------------------------
# 1. Read b1 paths file and extract platekeys in order
#
# rMATS b1 files can be one path per line OR comma-separated paths.
# We read all lines, split on commas, and flatten to get one path per element.
# --------------------------------------------------------------------------
raw_lines <- readLines(PATHS_FILE)
raw_lines <- raw_lines[nzchar(trimws(raw_lines))]  # drop blank lines

# Split on commas (handles comma-separated format; single-path lines unaffected)
all_paths <- unlist(strsplit(raw_lines, ","))
all_paths <- trimws(all_paths)
all_paths <- all_paths[nzchar(all_paths)]

cat("BAM paths found in b1 file:", length(all_paths), "\n")
stopifnot(
  "Expected 11 BAM paths in b1 file" = length(all_paths) == N_CASES
)

# Extract platekey: the folder immediately containing the BAM file
# e.g. /path/to/LP*-RNA_*/sample.bam → LP*-RNA_*
platekeys_ordered <- basename(dirname(all_paths))

cat("\nPlatekeys in rMATS b1 order (= IncLevel1 sample order):\n")
for (i in seq_along(platekeys_ordered)) {
  cat(sprintf("  gel:case:%d  %s\n", i, platekeys_ordered[i]))
}

# --------------------------------------------------------------------------
# 2. Load participant genotype TSV and join on platekey
# --------------------------------------------------------------------------
geno <- read.delim(GENOTYPE_FILE, stringsAsFactors = FALSE)

# Confirm required columns exist
required_cols <- c("platekey", "sample_group")
missing_cols  <- setdiff(required_cols, names(geno))
stopifnot(
  "Missing required columns in genotype TSV" = length(missing_cols) == 0
)

# Check all platekeys are present in the TSV
missing_keys <- setdiff(platekeys_ordered, geno$platekey)
if (length(missing_keys) > 0) {
  stop("Platekeys from b1 file not found in genotype TSV:\n  ",
       paste(missing_keys, collapse = "\n  "))
}

# Build output: one row per case, in b1 order
sample_geno <- tibble(
  sample_id = paste0("gel:case:", seq_len(N_CASES)),
  platekey  = platekeys_ordered
) %>%
  left_join(
    geno %>% select(platekey, participant_id, hgvs = sample_group),
    by      = "platekey",
    relationship = "many-to-one"
  )

stopifnot(
  "Row count changed after join"          = nrow(sample_geno) == N_CASES,
  "Some samples have missing HGVS values" = !anyNA(sample_geno$hgvs)
)

# --------------------------------------------------------------------------
# 3. Verify HGVS values look plausible (should start with "n.")
# --------------------------------------------------------------------------
unexpected_hgvs <- sample_geno$hgvs[!grepl("^n\\.", sample_geno$hgvs)]
if (length(unexpected_hgvs) > 0) {
  warning("Some HGVS values don't start with 'n.' — verify these are correct:\n  ",
          paste(unique(unexpected_hgvs), collapse = "\n  "))
}

# --------------------------------------------------------------------------
# 4. Print summary and save
# --------------------------------------------------------------------------
cat("\nSample–genotype mapping:\n")
print(sample_geno %>% select(sample_id, platekey, participant_id, hgvs))

cat("\nVariant counts:\n")
print(table(sample_geno$hgvs))

write.csv(sample_geno[, c("sample_id", "hgvs")], OUT_FILE, row.names = FALSE)
cat("\nSaved to:", OUT_FILE, "\n")
cat("Ready to join to SGE function scores on the 'hgvs' column.\n")
