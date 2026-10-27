#!/usr/bin/env Rscript
#
# strict_featsel_prepare.R
# Reviewer response: rebuild the insT signature with ALL feature selection
# (rMATS significance, direction, depth/NA filters, G_AUC, LASSO) done on
# training samples only. Phase 1 of 3.
#
#   1. Split for this analysis -> data/split_strict_featsel.csv
#      Cases: copied unchanged from split_dev_holdout_v4_mixed.csv.
#      Controls: N_HO_CTRL GEL + N_HO_CTRL Nava controls held out at random
#      (fixed seed); the rest are training. Written once, then reused.
#   2. rMATS stat input directories (A5SS only), rebuilt from the per-sample
#      counts in each cohort's A5SS.MATS.JC.txt:
#        gel_all,   nava_all   - every case vs every control. Must reproduce
#                                the original FDRs (checked in phase 3).
#        gel_train, nava_train - insT training cases vs training controls.
#
# Phase 2 (HPC):     bsub < see strict_featsel_rmats_stat.sh header
# Phase 3 (Desktop): pixi run Rscript scripts/gel/strict_featsel_insT_lasso.R
#
# Inputs:
#   data/rmats/A5SS.MATS.JC.txt
#   data/nava_rmats_tables/A5SS.MATS.JC.txt
#   data/split_dev_holdout_v4_mixed.csv
#   data/sample_genotype.csv
#   outputs/nava_severity_labels.csv
#
# Outputs:
#   data/split_strict_featsel.csv
#   outputs/strict_featsel/rmats_stat/{gel,nava}_{all,train}/  (rMATS --od)
#
# Run: pixi run Rscript scripts/gel/strict_featsel_prepare.R

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
GEL_JC        <- "data/rmats/A5SS.MATS.JC.txt"
NAVA_JC       <- "data/nava_rmats_tables/A5SS.MATS.JC.txt"
SPLIT_FILE    <- "data/split_dev_holdout_v4_mixed.csv"
GENO_FILE     <- "data/sample_genotype.csv"
NAVA_SEV_FILE <- "outputs/nava_severity_labels.csv"
NEW_SPLIT     <- "data/split_strict_featsel.csv"
STAT_DIR      <- "outputs/strict_featsel/rmats_stat"

N_GEL_CASES     <- 11L
N_GEL_CONTROLS  <- 33L
N_NAVA_CASES    <- 19L
N_NAVA_CONTROLS <- 21L

INST_HGVS <- "n.64_65insT"
N_HO_CTRL <- 3L
SEED      <- 2026L

stopifnot(
  "GEL A5SS table not found"     = file.exists(GEL_JC),
  "Nava A5SS table not found"    = file.exists(NAVA_JC),
  "Split file not found"         = file.exists(SPLIT_FILE),
  "Genotype file not found"      = file.exists(GENO_FILE),
  "Nava severity file not found" = file.exists(NAVA_SEV_FILE)
)

gel_case_ids  <- paste0("gel:case:",     seq_len(N_GEL_CASES))
gel_ctrl_ids  <- paste0("gel:control:",  seq_len(N_GEL_CONTROLS))
nava_case_ids <- paste0("nava:case:",    seq_len(N_NAVA_CASES))
nava_ctrl_ids <- paste0("nava:control:", seq_len(N_NAVA_CONTROLS))

# --------------------------------------------------------------------------
# Split: cases from the locked v4 split, plus held-out controls
# --------------------------------------------------------------------------
if (file.exists(NEW_SPLIT)) {
  cat("Reusing existing split:", NEW_SPLIT, "\n")
  split_df <- read.csv(NEW_SPLIT, stringsAsFactors = FALSE)
} else {
  case_split <- read.csv(SPLIT_FILE, stringsAsFactors = FALSE)
  stopifnot(all(c("sample_id", "cohort", "genotype", "set") %in% names(case_split)),
            all(case_split$set %in% c("train", "holdout")),
            !any(grepl(":control:", case_split$sample_id)))
  set.seed(SEED)
  gel_ho  <- sample(seq_len(N_GEL_CONTROLS),  N_HO_CTRL)
  nava_ho <- sample(seq_len(N_NAVA_CONTROLS), N_HO_CTRL)
  ctrl_split <- data.frame(
    sample_id = c(gel_ctrl_ids, nava_ctrl_ids),
    cohort    = c(rep("GEL", N_GEL_CONTROLS), rep("Nava", N_NAVA_CONTROLS)),
    genotype  = "control",
    set       = c(ifelse(seq_len(N_GEL_CONTROLS)  %in% gel_ho,  "holdout", "train"),
                  ifelse(seq_len(N_NAVA_CONTROLS) %in% nava_ho, "holdout", "train")),
    stringsAsFactors = FALSE)
  split_df <- rbind(case_split[, c("sample_id", "cohort", "genotype", "set")], ctrl_split)
  write.csv(split_df, NEW_SPLIT, row.names = FALSE)
  cat("Wrote new split:", NEW_SPLIT, "\n")
}
stopifnot(
  "Split missing some controls" = all(c(gel_ctrl_ids, nava_ctrl_ids) %in% split_df$sample_id),
  "Wrong number of held-out GEL controls" =
    sum(split_df$sample_id %in% gel_ctrl_ids  & split_df$set == "holdout") == N_HO_CTRL,
  "Wrong number of held-out Nava controls" =
    sum(split_df$sample_id %in% nava_ctrl_ids & split_df$set == "holdout") == N_HO_CTRL
)
train_ids <- split_df$sample_id[split_df$set == "train"]

# --------------------------------------------------------------------------
# insT training carriers (genotype from the same sources as the LOO script)
# --------------------------------------------------------------------------
read_sep_tolerant <- function(path) {
  first <- readLines(path, n = 1L)
  sep   <- if (grepl(",", first, fixed = TRUE)) "," else "\t"
  read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
}
geno     <- read.csv(GENO_FILE, stringsAsFactors = FALSE)
nava_lbl <- read_sep_tolerant(NAVA_SEV_FILE)
nava_lbl$sample_id <- paste0("nava:case:", nava_lbl$sample_n)
hgvs <- setNames(c(geno$hgvs, nava_lbl$hgvs), c(geno$sample_id, nava_lbl$sample_id))

gel_insT_train  <- intersect(gel_case_ids[hgvs[gel_case_ids]   == INST_HGVS], train_ids)
nava_insT_train <- intersect(nava_case_ids[hgvs[nava_case_ids] == INST_HGVS], train_ids)
gel_ctrl_train  <- intersect(gel_ctrl_ids,  train_ids)
nava_ctrl_train <- intersect(nava_ctrl_ids, train_ids)

cat(sprintf("GEL  training: %d insT cases vs %d controls\n",
            length(gel_insT_train),  length(gel_ctrl_train)))
cat(sprintf("Nava training: %d insT cases vs %d controls\n",
            length(nava_insT_train), length(nava_ctrl_train)))
stopifnot(length(gel_insT_train) > 0L, length(nava_insT_train) > 0L)

# Position of each sample within its rMATS group (SAMPLE_1 = cases, SAMPLE_2 = controls)
id_index <- function(ids) as.integer(sub(".*:", "", ids))

# --------------------------------------------------------------------------
# rMATS stat input writers
# --------------------------------------------------------------------------
FROMGTF_COLS <- c("ID", "GeneID", "geneSymbol", "chr", "strand",
                  "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
                  "flankingES", "flankingEE")
RAW_COLS <- c("ID", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
              "IncFormLen", "SkipFormLen")
# Header-only tables for the other event types: rMATS --task stat expects the
# files but tests 0 events, same trick as null_lasso_rmats_stat.sh.
OTHER_FROMGTF_COLS <- list(
  SE   = c("ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
           "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"),
  A3SS = FROMGTF_COLS,
  MXE  = c("ID", "GeneID", "geneSymbol", "chr", "strand", "1stExonStart_0base", "1stExonEnd",
           "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE",
           "downstreamES", "downstreamEE"),
  RI   = c("ID", "GeneID", "geneSymbol", "chr", "strand", "riExonStart_0base", "riExonEnd",
           "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")
)

read_jc <- function(path, n_cases, n_ctrls) {
  df <- read.delim(path, colClasses = "character", check.names = FALSE)
  stopifnot("Unexpected A5SS.MATS.JC.txt columns" =
              identical(names(df)[1:18], c(FROMGTF_COLS, RAW_COLS)))
  n1 <- lengths(strsplit(df$IJC_SAMPLE_1, ",", fixed = TRUE))
  n2 <- lengths(strsplit(df$IJC_SAMPLE_2, ",", fixed = TRUE))
  stopifnot("SAMPLE_1 count length != n cases"    = all(n1 == n_cases),
            "SAMPLE_2 count length != n controls" = all(n2 == n_ctrls))
  df
}

subset_counts <- function(x, idx) {
  vapply(strsplit(x, ",", fixed = TRUE),
         function(v) paste(v[idx], collapse = ","), character(1L))
}

write_tsv_header <- function(cols, path) writeLines(paste(cols, collapse = "\t"), path)

write_stat_od <- function(df, case_idx, ctrl_idx, od) {
  if (file.exists(file.path(od, "A5SS.MATS.JC.txt"))) {
    cat("  ", od, ": stat output already present, leaving untouched\n")
    return(invisible(NULL))
  }
  unlink(od, recursive = TRUE)
  dir.create(od, recursive = TRUE)

  # fromGTF: GeneID/geneSymbol quoted, as rMATS writes them
  gtf <- df[, 1:11]
  names(gtf) <- FROMGTF_COLS
  f_gtf <- file.path(od, "fromGTF.A5SS.txt")
  write_tsv_header(FROMGTF_COLS, f_gtf)
  write.table(gtf, f_gtf, sep = "\t", quote = c(2L, 3L), row.names = FALSE,
              col.names = FALSE, append = TRUE)

  raw <- data.frame(
    ID           = df[[12]],
    IJC_SAMPLE_1 = subset_counts(df$IJC_SAMPLE_1, case_idx),
    SJC_SAMPLE_1 = subset_counts(df$SJC_SAMPLE_1, case_idx),
    IJC_SAMPLE_2 = subset_counts(df$IJC_SAMPLE_2, ctrl_idx),
    SJC_SAMPLE_2 = subset_counts(df$SJC_SAMPLE_2, ctrl_idx),
    IncFormLen   = df$IncFormLen,
    SkipFormLen  = df$SkipFormLen,
    stringsAsFactors = FALSE)
  stopifnot(identical(raw$ID, df[[1]]))
  f_raw <- file.path(od, "JC.raw.input.A5SS.txt")
  write_tsv_header(RAW_COLS, f_raw)
  write.table(raw, f_raw, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = FALSE, append = TRUE)
  write_tsv_header(RAW_COLS, file.path(od, "JCEC.raw.input.A5SS.txt"))

  for (et in names(OTHER_FROMGTF_COLS)) {
    write_tsv_header(OTHER_FROMGTF_COLS[[et]], file.path(od, paste0("fromGTF.", et, ".txt")))
    write_tsv_header(RAW_COLS, file.path(od, paste0("JC.raw.input.",   et, ".txt")))
    write_tsv_header(RAW_COLS, file.path(od, paste0("JCEC.raw.input.", et, ".txt")))
  }
  cat(sprintf("  %s: %d events, %d vs %d samples\n",
              od, nrow(raw), length(case_idx), length(ctrl_idx)))
}

# --------------------------------------------------------------------------
# Build the four stat input directories
# --------------------------------------------------------------------------
cat("\nReading A5SS count tables...\n")
gel_jc  <- read_jc(GEL_JC,  N_GEL_CASES,  N_GEL_CONTROLS)
nava_jc <- read_jc(NAVA_JC, N_NAVA_CASES, N_NAVA_CONTROLS)

cat("Writing rMATS stat input directories...\n")
write_stat_od(gel_jc,  seq_len(N_GEL_CASES),  seq_len(N_GEL_CONTROLS),
              file.path(STAT_DIR, "gel_all"))
write_stat_od(nava_jc, seq_len(N_NAVA_CASES), seq_len(N_NAVA_CONTROLS),
              file.path(STAT_DIR, "nava_all"))
write_stat_od(gel_jc,  id_index(gel_insT_train),  id_index(gel_ctrl_train),
              file.path(STAT_DIR, "gel_train"))
write_stat_od(nava_jc, id_index(nava_insT_train), id_index(nava_ctrl_train),
              file.path(STAT_DIR, "nava_train"))

cat("\nNext (HPC): bsub -q medium -P re_gecip_machine_learning -J strict_featsel_stat",
    "-o hpc_logs/strict_featsel_%J.out -e hpc_logs/strict_featsel_%J.err",
    "bash scripts/gel/strict_featsel_rmats_stat.sh\n")
