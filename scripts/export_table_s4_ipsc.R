#!/usr/bin/env Rscript
# =============================================================================
# export_table_s4_ipsc.R
# Two-in-one export:
#
#   (1) Annotate Table 1 (outputs/candidate_table/a5ss_candidate_full_table.csv,
#       the 105-event GEL/Nava-concordant-or-discordant A5SS candidate pool)
#       against the local iPSC rMATS runs (4 DIVs, ReNU vs KOLF): quantified
#       flag, raw PValue/FDR, delta logit-PSI (ReNU - KOLF) per DIV.
#       Matching + delta-logit-PSI reuse the exact methodology figure4.R Panel
#       B uses (coord_nostrand join, bb-logit-PSI = logit((IJC+0.5)/(IJC+SJC+1))).
#       -> outputs/candidate_table/table1_ipsc_annotated.csv
#       Load-bearing: figure4.R and figure_s8.R read this file.
#
#   (2) Pivot the concordant subset into Supplementary Table S4 (long
#       format, one row per event x DIV) with novel/annotated splice-site
#       classification against gencode v49. Column names/types match
#       Table S1 so the two tables join directly on event_id.
#       -> outputs/candidate_table/table_s4_ipsc.csv
#
# Inputs:
#   outputs/candidate_table/a5ss_candidate_full_table.csv  (Table 1, 105 rows)
#   data/splicing/div{1,5,9,14}/A5SS.MATS.JC.txt            (local iPSC rMATS)
#   data/reference/gencode.v49.primary_assembly.annotation.gtf.gz
#
# Run from project root: Rscript scripts/export_table_s4_ipsc.R
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

CANDIDATE_CSV <- "outputs/candidate_table/a5ss_candidate_full_table.csv"
OUT_TABLE1    <- "outputs/candidate_table/table1_ipsc_annotated.csv"
OUT_TABLE_S4  <- "outputs/candidate_table/table_s4_ipsc.csv"
GTF_PATH      <- "data/reference/gencode.v49.primary_assembly.annotation.gtf.gz"
SS_CACHE_PATH <- "data/reference/gencode_v49_ss_set.rds"

DIVS_LC    <- c("div1", "div5", "div9", "div14")
DIVS_LABEL <- c("DIV1", "DIV5", "DIV9", "DIV14")
FDR_THRESH <- 0.05
P_THRESH   <- 0.05

stopifnot("Table 1 candidate CSV not found" = file.exists(CANDIDATE_CSV))

candidate_tbl <- read_csv(CANDIDATE_CSV, show_col_types = FALSE)
candidate_tbl$coord_nostrand <- with(candidate_tbl,
  paste(chr, longExonStart_0base, longExonEnd_1based,
        shortES_0base, shortEE_1based,
        flankingES_0base, flankingEE_1based, sep = ":"))
stopifnot("Table 1 coord_nostrand not unique" = !anyDuplicated(candidate_tbl$coord_nostrand))
cat(sprintf("Table 1: %d candidate events loaded\n", nrow(candidate_tbl)))

# Same positional key as every other local iPSC rMATS join (CLAUDE.md rMATS
# quirks): columns 4:11 = chr, strand + the 6 A5SS coordinate columns.
build_coord_nostrand <- function(df) {
  coord_key <- do.call(paste, c(df[, 4:11], sep = ":"))
  sub("^([^:]+):[+-]:(.*)", "\\1:\\2", coord_key)
}

bb_logit_mat <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}
parse_intvec <- function(s) do.call(rbind, lapply(strsplit(s, ","), as.integer))

load_div <- function(div) {
  path <- file.path("data/splicing", div, "A5SS.MATS.JC.txt")
  stopifnot("iPSC A5SS rMATS file missing" = file.exists(path))
  d <- read_tsv(path, show_col_types = FALSE,
                col_types = cols(.default = col_character()))
  d$coord_nostrand <- build_coord_nostrand(d)
  stopifnot("iPSC rMATS coord_nostrand not unique within a DIV" =
              !anyDuplicated(d$coord_nostrand))
  d$PValue <- as.numeric(d$PValue)
  d$FDR    <- as.numeric(d$FDR)
  d
}

# ---------------------------------------------------------------------------
# 1. Per-DIV annotation -> table1_ipsc_annotated.csv
# ---------------------------------------------------------------------------
out <- candidate_tbl

for (k in seq_along(DIVS_LC)) {
  div_lc  <- DIVS_LC[k]
  div_lbl <- DIVS_LABEL[k]
  d <- load_div(div_lc)
  m <- match(candidate_tbl$coord_nostrand, d$coord_nostrand)
  found <- !is.na(m)

  quantified  <- found
  pvalue_col  <- rep(NA_real_, nrow(candidate_tbl))
  fdr_col     <- rep(NA_real_, nrow(candidate_tbl))
  delta_logit <- rep(NA_real_, nrow(candidate_tbl))

  pvalue_col[found] <- d$PValue[m[found]]
  fdr_col[found]    <- d$FDR[m[found]]

  ijc1 <- parse_intvec(d$IJC_SAMPLE_1[m[found]])
  sjc1 <- parse_intvec(d$SJC_SAMPLE_1[m[found]])
  ijc2 <- parse_intvec(d$IJC_SAMPLE_2[m[found]])
  sjc2 <- parse_intvec(d$SJC_SAMPLE_2[m[found]])
  # Sample 1 = ReNU, Sample 2 = KOLF (local convention, CLAUDE.md: local
  # IncLevelDifference = ReNU - KOLF) -- mean bb-logit-PSI across replicates.
  delta_logit[found] <- rowMeans(bb_logit_mat(ijc1, sjc1), na.rm = TRUE) -
                        rowMeans(bb_logit_mat(ijc2, sjc2), na.rm = TRUE)

  out[[paste0("ipsc_quantified_", div_lbl)]]     <- quantified
  out[[paste0("ipsc_PValue_", div_lbl)]]         <- pvalue_col
  out[[paste0("ipsc_FDR_", div_lbl)]]            <- fdr_col
  out[[paste0("ipsc_FDRsig_", div_lbl)]]         <- !is.na(fdr_col) & fdr_col < FDR_THRESH
  out[[paste0("ipsc_Psig_", div_lbl)]]           <- !is.na(pvalue_col) & pvalue_col < P_THRESH
  out[[paste0("ipsc_delta_logitPSI_", div_lbl)]] <- delta_logit

  cat(sprintf("%s: quantified=%d/%d, FDR<%.2g=%d, raw P<%.2g=%d\n",
              div_lbl, sum(quantified), nrow(candidate_tbl),
              FDR_THRESH, sum(out[[paste0("ipsc_FDRsig_", div_lbl)]], na.rm = TRUE),
              P_THRESH, sum(out[[paste0("ipsc_Psig_", div_lbl)]], na.rm = TRUE)))
}

quant_cols <- paste0("ipsc_quantified_", DIVS_LABEL)
out$ipsc_quantified_any_DIV <- Reduce(`|`, out[quant_cols])
out$ipsc_quantified_all_DIV <- Reduce(`&`, out[quant_cols])

cat(sprintf("\nQuantified in >=1 DIV: %d/%d\n",
            sum(out$ipsc_quantified_any_DIV), nrow(out)))
cat(sprintf("Quantified in all 4 DIVs: %d/%d\n",
            sum(out$ipsc_quantified_all_DIV), nrow(out)))

dir.create(dirname(OUT_TABLE1), recursive = TRUE, showWarnings = FALSE)
write_csv(out, OUT_TABLE1)
cat(sprintf("\nSaved: %s (%d rows, %d cols)\n", OUT_TABLE1, nrow(out), ncol(out)))

# ---------------------------------------------------------------------------
# 2. Concordant subset -> Table S4 (long format, novel-junction annotation)
# ---------------------------------------------------------------------------
concordant_tbl <- out[out$concordant, ]
n_pool <- nrow(concordant_tbl)
cat(sprintf("\nTable S4: %d concordant candidate-pool events (Table S1 subset)\n",
            n_pool))

# Novel/annotated splice-site classification (gencode v49) -- same logic as
# figure_s8.R's check_novel(). Reuses the mtime-invalidated cache.
if (file.exists(SS_CACHE_PATH) && file.exists(GTF_PATH) &&
    file.mtime(SS_CACHE_PATH) > file.mtime(GTF_PATH)) {
  ss_set <- readRDS(SS_CACHE_PATH)
} else {
  stopifnot("gencode v49 GTF not found" = file.exists(GTF_PATH))
  cmd <- sprintf("gzcat %s | awk -F'\t' 'BEGIN{OFS=\"\t\"} $3==\"exon\" {print $1,$4,$5,$7}'",
                 shQuote(GTF_PATH))
  exons <- read.delim(pipe(cmd), header = FALSE, stringsAsFactors = FALSE,
                      col.names = c("chr", "start", "end", "strand"))
  exons$key <- paste(exons$chr, exons$strand, sep = "|")
  pos_long  <- data.frame(key = c(exons$key, exons$key),
                          pos = c(exons$start, exons$end))
  ss_set    <- split(pos_long$pos, pos_long$key)
  ss_set    <- lapply(ss_set, function(v) sort(unique(as.integer(v))))
  saveRDS(ss_set, SS_CACHE_PATH)
}
check_novel <- function(row) {
  coords1based <- c(row$longExonStart_0base + 1, row$longExonEnd_1based,
                    row$shortES_0base + 1, row$shortEE_1based,
                    row$flankingES_0base + 1, row$flankingEE_1based)
  key <- paste(row$chr, row$strand, sep = "|")
  annotated_pos <- ss_set[[key]]
  if (is.null(annotated_pos)) return(TRUE)
  !all(coords1based %in% annotated_pos)
}
concordant_tbl$is_novel_junction <- vapply(seq_len(n_pool),
  function(i) check_novel(concordant_tbl[i, ]), logical(1))
cat(sprintf("Novel junctions: %d/%d\n",
            sum(concordant_tbl$is_novel_junction), n_pool))

table_s4 <- do.call(rbind, lapply(DIVS_LABEL, function(div) {
  delta <- concordant_tbl[[paste0("ipsc_delta_logitPSI_", div)]]
  data.frame(
    event_id                  = concordant_tbl$event_id,
    gene_symbol               = concordant_tbl$gene_symbol,
    is_insT_Tloop_signature   = concordant_tbl$is_insT_Tloop_signature,
    direction                 = concordant_tbl$direction,
    is_novel_junction         = concordant_tbl$is_novel_junction,
    DIV                       = div,
    ipsc_quantified           = concordant_tbl[[paste0("ipsc_quantified_", div)]],
    ipsc_rmats_FDR            = concordant_tbl[[paste0("ipsc_FDR_", div)]],
    ipsc_delta_logitPSI       = delta,
    ipsc_direction_concordant = concordant_tbl[[paste0("ipsc_quantified_", div)]] &
      !is.na(delta) & sign(delta) == concordant_tbl$direction,
    stringsAsFactors = FALSE)
}))
table_s4$DIV <- factor(table_s4$DIV, levels = DIVS_LABEL)
table_s4 <- table_s4[order(table_s4$event_id, table_s4$DIV), ]
table_s4$DIV <- as.character(table_s4$DIV)

stopifnot("Table S4 row count != n_pool events x 4 DIVs" =
            nrow(table_s4) == n_pool * length(DIVS_LABEL))

cat("\nPer-DIV summary:\n")
print(table_s4 %>% group_by(DIV) %>%
        summarise(quantified            = sum(ipsc_quantified),
                  direction_concordant  = sum(ipsc_direction_concordant, na.rm = TRUE),
                  .groups = "drop"))

write_csv(table_s4, OUT_TABLE_S4)
cat(sprintf("\nSaved: %s (%d rows, %d cols)\n",
            OUT_TABLE_S4, nrow(table_s4), ncol(table_s4)))
