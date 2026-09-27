#!/usr/bin/env Rscript
# =============================================================================
# compare_105_vs_chen_and_nava.R
# Single consolidated comparison of our 105-event A5SS pool
# (reconstruct_105_a5ss_events.R) against two previously published
# RNU4-2-associated splicing lists:
#   - Chen et al. (FRASER, 12 events, data/Chen_et_al/Chen_table_s3.tsv)
#   - Nava et al.'s ACTUAL PUBLISHED Supplementary Table 12 ("Significant
#     alternative 5' splice sites (5'SS) events called by rMATS in patients
#     with RNU4-2 variants", 41588_2025_2184_MOESM4_ESM.xlsx) -- used
#     directly, NOT their gene-deduplicated "111" summary figure (dropped:
#     gene-level dedup hides real matches when a gene has >1 significant
#     A5SS event, e.g. MADD; the raw 148-row table is the complete,
#     unambiguous reference).
#
# Also checks the 10 distinct insT/Stem III LASSO-selected signature events
# against both lists (signature identities are no longer GEL-protected for
# this comparison -- already exported through the airlock, confirmed by
# the project owner).
#
# All matching is by EXACT coord_nostrand (chr:6-coordinates) -- gene-symbol
# overlap is never used as a match on its own; a gene can have multiple
# distinct A5SS events (see CLAUDE.md rMATS quirks / project_chen_et_al_
# comparison memory: this exact mistake overcounted an earlier version of
# this comparison by 5/30).
#
# Chen's coordinates need strand-aware anchoring + a coordinate-system
# correction to compare against rMATS's donor coordinate -- see "Chen
# donor-coordinate search" section below for the full derivation (validated
# exact, 0bp residual, across every real match in both cohorts).
#
# Inputs:
#   outputs/supplementary_tables/Table_S1_candidate_events.csv (Table S1)
#   data/Chen_et_al/Chen_table_s3.tsv
#   data/gel_rmats_tables/A5SS_summary.csv
#   data/nava_rmats_tables/A5SS.MATS.JC.txt
#   data/nava_rmats_tables/A5SS.MATS.JCEC.txt (gunzip .gz first if needed)
#   data/nava_rmats_tables/41588_2025_2184_MOESM4_ESM.xlsx
#   data/splicing_signature/insT_stemIII_export/{insT,stemIII}_signature_export.csv
#
# Outputs:
#   outputs/chen_et_al_comparison.csv               (per gene x cohort, Chen search)
#   outputs/nava_supp_table12_comparison.csv        (our 105 vs Supp Table 12, exact)
#   outputs/venn_105_vs_nava_our_only_reasons.csv   (per-event categorisation:
#                                                    first Nava row-filter that
#                                                    drops each our-only event)
#   outputs/signature_vs_chen_nava.csv              (10 signature events vs both lists)
#
# Run: Rscript scripts/compare_105_vs_chen_and_nava.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

coord_nostrand <- function(df) {
  paste(df$chr, df$longExonStart_0base, df$longExonEnd,
        df$shortES, df$shortEE, df$flankingES, df$flankingEE, sep = ":")
}

# Load the 105-event sig-in-both A5SS pool from Table S1 (single source of
# truth). Table S1's column names differ from the earlier reconstruct output
# and it doesn't carry a coord_nostrand column, so rebuild it from the 6
# rMATS coordinate columns and remap FDR names to this script's r105 schema.
s1 <- read.csv("outputs/supplementary_tables/Table_S1_candidate_events.csv",
               stringsAsFactors = FALSE)
s1_coord_nostrand <- paste(s1$chr,
                            s1$longExonStart_0base, s1$longExonEnd_1based,
                            s1$shortES_0base,       s1$shortEE_1based,
                            s1$flankingES_0base,    s1$flankingEE_1based,
                            sep = ":")
r105 <- data.frame(
  coord_nostrand = s1_coord_nostrand,
  gene_symbol    = s1$gene_symbol,
  gel_FDR        = s1$wb_rmats_FDR,
  nava_FDR       = s1$lymph_rmats_FDR,
  stringsAsFactors = FALSE)
stopifnot(nrow(r105) == 105L)

# =============================================================================
# 1. Chen et al. donor-in-pool check (12 events vs our 105 Table S1 pool)
# For each of Chen's 12 events, anchor to a single donor coordinate
# (strand-aware), then check whether that donor position exactly matches
# either the long-form or short-form donor of any of our 105 A5SS events
# on the same chr + strand. If so: Chen event is confirmed in our pool.
# =============================================================================
chen <- read.delim("data/Chen_et_al/Chen_table_s3.tsv",
                   stringsAsFactors = FALSE, check.names = FALSE)
pos <- strcapture("chr([0-9XYM]+):([0-9]+)-([0-9]+)", chen[["Genomic position"]],
                  proto = list(chr = character(), start = integer(), end = integer()))
chen$chr   <- paste0("chr", pos$chr)
chen$start <- pos$start
chen$end   <- pos$end
bad_end <- (chen$end - chen$start + 1L) != chen$width
chen$end[bad_end] <- chen$start[bad_end] + chen$width[bad_end] - 1L   # ABCA7 data-quality fix, see git history
chen$strand <- ifelse(chen$strand == "forward", "+", "-")

# Chen -> single donor coordinate, strand-aware. Empirically-derived rMATS
# anchoring (0bp residual across every real match in both cohorts):
#   + strand: donor = chen$start - 1
#   - strand: donor = chen$end (no offset)
chen$donor_ref <- ifelse(chen$strand == "+", chen$start - 1L, chen$end)

# Our 105-event A5SS pool donors, strand-aware. In Table S1 the fixed
# boundary is longExonStart_0base==shortES_0base (+strand) or
# longExonEnd_1based==shortEE_1based (-strand); the VARIABLE donor is on
# the opposite side.
r105_donors <- {
  is_plus <- s1$strand == "+"
  stopifnot(
    "+ strand: longExonStart_0base must equal shortES_0base"  =
      all(s1$longExonStart_0base[is_plus] == s1$shortES_0base[is_plus]),
    "- strand: longExonEnd_1based must equal shortEE_1based"  =
      all(s1$longExonEnd_1based[!is_plus] == s1$shortEE_1based[!is_plus])
  )
  data.frame(
    coord_nostrand = s1_coord_nostrand,
    gene_symbol    = s1$gene_symbol,
    chr            = s1$chr,
    strand         = s1$strand,
    donor_long     = ifelse(is_plus, s1$longExonEnd_1based, s1$longExonStart_0base),
    donor_short    = ifelse(is_plus, s1$shortEE_1based,     s1$shortES_0base),
    stringsAsFactors = FALSE)
}

chen_verdict <- do.call(rbind, lapply(seq_len(nrow(chen)), function(i) {
  d <- chen$donor_ref[i]; ch <- chen$chr[i]; st <- chen$strand[i]
  hit <- r105_donors[r105_donors$chr == ch & r105_donors$strand == st &
                     (r105_donors$donor_long == d | r105_donors$donor_short == d), , drop = FALSE]
  data.frame(
    gene_symbol_chen        = chen[["HGNC symbol"]][i],
    chen_donor_coord        = sprintf("%s:%d%s", ch, d, st),
    in_105_pool             = nrow(hit) > 0L,
    matched_events          = if (nrow(hit)) paste(hit$coord_nostrand, collapse = ";") else NA_character_,
    matched_gene_symbols    = if (nrow(hit)) paste(unique(hit$gene_symbol), collapse = "/") else NA_character_,
    stringsAsFactors = FALSE)
}))

cat("=== Chen et al. (n=12) vs our 105-event Table S1 pool ===\n")
cat("Match rule: Chen donor position == donor_long OR donor_short of an A5SS event on the same chr + strand.\n\n")
print(chen_verdict, row.names = FALSE)
cat(sprintf("\nMatched: %d / %d Chen events\n",
            sum(chen_verdict$in_105_pool), nrow(chen_verdict)))
write.csv(chen_verdict, "outputs/chen_et_al_comparison.csv", row.names = FALSE)

# =============================================================================
# 2. Nava Supplementary Table 12 -- direct comparison (event-level, no dedup)
# =============================================================================
t12 <- suppressMessages(read_excel(
  "data/nava_rmats_tables/41588_2025_2184_MOESM4_ESM.xlsx", sheet = "Supp. Table 12", skip = 1
))
for (col in c("longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE")) {
  t12[[col]] <- as.numeric(t12[[col]])
}
t12$key <- coord_nostrand(t12)
stopifnot("Supp Table 12 should have 148 rows" = nrow(t12) == 148L)

overlap_105_t12 <- intersect(r105$coord_nostrand, t12$key)
cat(sprintf("\n=== Our 105-event pool vs Nava Supp Table 12 (%d rows) ===\n", nrow(t12)))
cat(sprintf("Exact overlap: %d / 105\n", length(overlap_105_t12)))
cat(sprintf("Novel to our joint analysis (not in Supp Table 12): %d / 105\n", 105L - length(overlap_105_t12)))

# ---- Gap A: our 105 not in Supp Table 12 -- per-event first-failed-filter
# categorisation against Nava's own filter (from
# https://github.com/benjamin-cogne/RNU4-2_transcriptomics, tag v2.0,
# `filter_rmats(...)` in rmats-function.ipynb).
# Filter order (applied to A5SS.MATS.JCEC.txt):
#   (a) FDR                <= 0.1
#   (b) |IncLevelDifference| >= 0.05
#   (c) mean(IJC_SAMPLE_1) >= 7
#   (d) mean(SJC_SAMPLE_1) >= 7
#   (e) geneSymbol         != "NA"
#   (f) All 40 per-sample PSI values non-NA (19 case + 21 control)
# Bucket each event by the FIRST filter it fails in this order.
# na.strings = character(0): keep the literal "NA" in geneSymbol as a string so
# it can be compared against Nava's Python `geneSymbol != 'NA'` faithfully.
# quote = "": Nava's Python does field[2].replace('"','') on geneSymbol.
jcec_path <- "data/nava_rmats_tables/A5SS.MATS.JCEC.txt"
if (file.exists(jcec_path)) {
  jcec <- read.delim(jcec_path, stringsAsFactors = FALSE, check.names = FALSE,
                     quote = "", na.strings = character(0))
  jcec$geneSymbol <- gsub('"', '', jcec$geneSymbol, fixed = TRUE)
  jcec$key <- coord_nostrand(jcec)

  our_only <- setdiff(r105$coord_nostrand, overlap_105_t12)
  m <- match(our_only, jcec$key)
  stopifnot("At least one our-only event not in JCEC (coord mismatch?)" =
              all(!is.na(m)))
  jf <- jcec[m, , drop = FALSE]

  n_case <- length(strsplit(jf$IJC_SAMPLE_1[1L], ",", fixed = TRUE)[[1]])
  n_ctrl <- length(strsplit(jf$IncLevel2[1L],    ",", fixed = TRUE)[[1]])
  stopifnot("Expected 19 cases + 21 controls per Nava row" =
              n_case == 19L && n_ctrl == 21L)

  mean_num <- function(x) mean(as.numeric(x))
  n_non_na <- function(x) sum(x != "NA")
  ijc1 <- strsplit(jf$IJC_SAMPLE_1, ",", fixed = TRUE)
  sjc1 <- strsplit(jf$SJC_SAMPLE_1, ",", fixed = TRUE)
  psi1 <- strsplit(jf$IncLevel1,    ",", fixed = TRUE)
  psi2 <- strsplit(jf$IncLevel2,    ",", fixed = TRUE)
  mean_IJC1 <- vapply(ijc1, mean_num, numeric(1L))
  mean_SJC1 <- vapply(sjc1, mean_num, numeric(1L))
  n_psi_ok  <- vapply(psi1, n_non_na, integer(1L)) +
               vapply(psi2, n_non_na, integer(1L))

  fail_fdr    <- jf$FDR > 0.1
  fail_delta  <- abs(jf$IncLevelDifference) < 0.05
  fail_ijc1   <- mean_IJC1 < 7
  fail_sjc1   <- mean_SJC1 < 7
  fail_gene   <- jf$geneSymbol == "NA"
  fail_psi_na <- n_psi_ok < (n_case + n_ctrl)

  primary <- rep(NA_character_, length(our_only))
  primary[fail_fdr]                     <- "FDR > 0.1"
  primary[is.na(primary) & fail_delta]  <- "|dPSI| < 0.05"
  primary[is.na(primary) & fail_ijc1]   <- "mean(IJC case) < 7"
  primary[is.na(primary) & fail_sjc1]   <- "mean(SJC case) < 7"
  primary[is.na(primary) & fail_gene]   <- "geneSymbol == 'NA'"
  primary[is.na(primary) & fail_psi_na] <- ">=1 sample PSI is NA"
  primary[is.na(primary)]               <- "PASSES ALL NAVA FILTERS"
  stopifnot("Every our-only event should be bucketed" = !anyNA(primary))

  cat(sprintf("\n=== Our-only events (n = %d) -- primary Nava-filter that drops them ===\n",
              length(our_only)))
  print(sort(table(primary), decreasing = TRUE))

  reasons <- data.frame(
    coord_nostrand    = our_only,
    gene_symbol_jcec  = jf$geneSymbol,
    FDR_jcec          = jf$FDR,
    IncLevelDiff_jcec = jf$IncLevelDifference,
    mean_IJC_case     = mean_IJC1,
    mean_SJC_case     = mean_SJC1,
    psi_non_na_count  = n_psi_ok,
    primary_reason    = primary,
    stringsAsFactors  = FALSE)
  write.csv(reasons, "outputs/venn_105_vs_nava_our_only_reasons.csv",
            row.names = FALSE)
  cat(sprintf("Per-event reasons: outputs/venn_105_vs_nava_our_only_reasons.csv (%d rows)\n",
              nrow(reasons)))
} else {
  cat("\n(A5SS.MATS.JCEC.txt not found -- skipping gap diagnosis A; gunzip the .gz to enable)\n")
}

# ---- Gap B: Supp Table 12 events not in our 105 -- diagnose against GEL
not_in_105 <- t12[!(t12$key %in% r105$coord_nostrand), ]
gel_all <- read.csv("data/gel_rmats_tables/A5SS_summary.csv", stringsAsFactors = FALSE)
gel_all$key <- coord_nostrand(gel_all)
gm <- gel_all[match(not_in_105$key, gel_all$key), ]
reasons_b <- ifelse(is.na(gm$FDR), "not present in GEL's tested A5SS catalog",
             ifelse(gm$FDR >= 0.05, "not significant in GEL (FDR>=0.05)", "significant in GEL too (unexpected)"))
cat("\nWhy Supp Table 12 events are missing from our 105:\n")
print(table(reasons_b))

out105 <- data.frame(coord_nostrand = overlap_105_t12, stringsAsFactors = FALSE)
out105$gene_symbol <- r105$gene_symbol[match(out105$coord_nostrand, r105$coord_nostrand)]
write.csv(out105, "outputs/nava_supp_table12_comparison.csv", row.names = FALSE)

# =============================================================================
# 3. Signature (insT + Stem III, 10 distinct events) vs both lists
# =============================================================================
inst <- read.csv("data/splicing_signature/insT_stemIII_export/insT_signature_export.csv", stringsAsFactors = FALSE)
stem <- read.csv("data/splicing_signature/insT_stemIII_export/stemIII_signature_export.csv", stringsAsFactors = FALSE)
sig_key <- function(df) paste(df$chr, df$longExonStart_0base, df$longExonEnd_1based,
                               df$shortES_0base, df$shortEE_1based, df$flankingES_0base, df$flankingEE_1based, sep = ":")
inst$key <- sig_key(inst); inst$signature <- "insT"
stem$key <- sig_key(stem); stem$signature <- "stemIII"
sig <- rbind(inst[, c("gene_symbol", "key", "signature")], stem[, c("gene_symbol", "key", "signature")])
sig <- aggregate(signature ~ gene_symbol + key, sig, paste, collapse = "+")
stopifnot("Expected 10 distinct signature events (6+6, 2 shared)" = nrow(sig) == 10L)

chen_matched_events <- unique(unlist(strsplit(
  chen_verdict$matched_events[chen_verdict$in_105_pool], ";", fixed = TRUE)))
sig$in_chen <- sig$key %in% chen_matched_events
sig$in_nava_supp_table12 <- sig$key %in% t12$key

cat("\n=== Signature events (insT + Stem III, n=10) vs Chen / Nava Supp Table 12 ===\n")
print(sig[, c("gene_symbol", "signature", "in_chen", "in_nava_supp_table12")], row.names = FALSE)
write.csv(sig, "outputs/signature_vs_chen_nava.csv", row.names = FALSE)

cat("\nDone.\n")
