#!/usr/bin/env Rscript
# =============================================================================
# figure_s8.R
# Paper Supplementary Figure S8 (RNU4-2 splicing biomarker manuscript):
# genome-wide characterization and technical detail supporting Figure 4's
# cross-tissue validation story. The 3-way whole-blood/lymphocyte/iPSC Venn
# (formerly panel C) was dropped.
#
#   Panel A (50%) - FDR<0.05 (ReNU vs KOLF) iPSC event counts per
#             rMATS event type (SE/A5SS/A3SS/MXE/RI), union across the 4
#             DIVs, stacked by shared+concordant (>=2 DIVs, same direction)
#             vs shared+discordant (>=2 DIVs, direction disagrees) vs
#             singleton (1 DIV only). Ordered by total count, descending.
#   Panel B (50%) - Fisher's exact test per event type (that type vs
#             other types pooled) for the shared+concordant fraction from
#             panel A: OR + 95% CI, log scale. Ordered by OR, ascending.
#             (Moved here from figure4.R's former panel A.)
# Also writes a small missingness table (not a figure panel) for the 17/95
# candidate-pool events with >=1 unquantified DIV, categorised by likely
# cause (TPM-limited / depth-limited / unresolved / never quantified) ->
# outputs/candidate_table/table1_missingness_summary.csv.
#
# Inputs (local; no GEL data touched):
#   data/splicing/{div1,5,9,14}/{SE,A5SS,A3SS,MXE,RI}.MATS.JC.txt
#       per-DIV iPSC rMATS (ReNU vs KOLF), all 5 event types.
#   outputs/candidate_table/table1_ipsc_annotated.csv
#       Table 1 (95-event concordant candidate pool) joined against per-DIV
#       iPSC rMATS (export_table_s4_ipsc.R).
#   data/reference/gencode.v49.primary_assembly.annotation.gtf.gz
#       for the novel/annotated splice-site classification in panel B.
#
# Outputs:
#   outputs/figures/figure_s8.pdf
#   outputs/figures/figure_s8.png
#   outputs/candidate_table/table1_missingness_summary.csv
#
# Run from project root: Rscript scripts/figure_s8.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(poibin)
})

# ---- Configuration ----------------------------------------------------------
OUT_DIR_FIG <- "outputs/figures"
OUT_PDF     <- file.path(OUT_DIR_FIG, "figure_s8.pdf")
OUT_PNG     <- file.path(OUT_DIR_FIG, "figure_s8.png")

CANDIDATE_CSV <- "outputs/candidate_table/table1_ipsc_annotated.csv"
GTF_PATH      <- "data/reference/gencode.v49.primary_assembly.annotation.gtf.gz"
SS_CACHE_PATH <- "data/reference/gencode_v49_ss_set.rds"
MISSING_CSV   <- "outputs/candidate_table/table1_missingness_summary.csv"

DIVS_LC    <- c("div1", "div5", "div9", "div14")
DIVS_LABEL <- c("DIV1", "DIV5", "DIV9", "DIV14")

# Palette -- matches figure2.R / figure4.R conventions
PAL_TLOOP <- "#D55E00"   # discordant-shared / novel
PAL_DISC  <- "#BBBBBB"
PAL_CONC  <- "#1F3A93"   # deep navy

stopifnot("Table 1 iPSC-annotated table not found -- run scripts/export_table_s4_ipsc.R first" =
            file.exists(CANDIDATE_CSV))
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

pub_theme <- theme_classic(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 8, face = "bold"),
        axis.text        = element_text(size = 7),
        axis.title       = element_text(size = 8),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7),
        plot.margin      = margin(2, 2, 2, 2))

# ---- Shared helpers ----------------------------------------------------------
EVENT_TYPES_A    <- c("SE", "A5SS", "A3SS", "MXE", "RI")
FDR_THRESH_SHARE <- 0.05

build_coord_nostrand <- function(df, event_type) {
  n_coord   <- if (event_type == "MXE") 8L else 6L
  coord_key <- do.call(paste, c(df[, seq(4L, 5L + n_coord)], sep = ":"))
  sub("^([^:]+):[+-]:(.*)", "\\1:\\2", coord_key)
}

load_type_div <- function(event_type, div) {
  path <- file.path("data/splicing", div, paste0(event_type, ".MATS.JC.txt"))
  d <- read_tsv(path, show_col_types = FALSE,
               col_types = cols(.default = col_character()))
  d$coord_nostrand      <- build_coord_nostrand(d, event_type)
  d$FDR                 <- as.numeric(d$FDR)
  d$IncLevelDifference  <- as.numeric(d$IncLevelDifference)
  d
}

# ---- Panel A: FDR<0.05 event counts per type, split shared/concordant vs
# shared/discordant vs singleton ---------------------------------------------
share_rows <- lapply(EVENT_TYPES_A, function(et) {
  per_div     <- setNames(lapply(DIVS_LC, load_type_div, event_type = et),
                          DIVS_LABEL)
  all_coords  <- unique(unlist(lapply(per_div, `[[`, "coord_nostrand")))
  sig_mat <- matrix(FALSE, nrow = length(all_coords), ncol = length(DIVS_LABEL),
                    dimnames = list(all_coords, DIVS_LABEL))
  dir_mat <- matrix(NA_real_, nrow = length(all_coords), ncol = length(DIVS_LABEL),
                    dimnames = list(all_coords, DIVS_LABEL))
  fdr_mat <- matrix(NA_real_, nrow = length(all_coords), ncol = length(DIVS_LABEL),
                    dimnames = list(all_coords, DIVS_LABEL))
  for (lbl in DIVS_LABEL) {
    d   <- per_div[[lbl]]
    m   <- match(d$coord_nostrand, all_coords)
    sig <- !is.na(d$FDR) & d$FDR < FDR_THRESH_SHARE
    sig_mat[m[sig], lbl] <- TRUE
    dir_mat[m, lbl]      <- sign(d$IncLevelDifference)
    fdr_mat[m, lbl]      <- d$FDR
  }
  n_sig_div  <- rowSums(sig_mat)
  is_sig_any <- n_sig_div >= 1L
  dir_concordant <- vapply(seq_along(all_coords), function(i) {
    dv <- dir_mat[i, sig_mat[i, ], drop = TRUE]
    length(dv) > 0L && length(unique(dv)) == 1L
  }, logical(1))
  is_shared            <- is_sig_any & (n_sig_div >= 2L) & dir_concordant
  is_shared_discordant <- is_sig_any & (n_sig_div >= 2L) & !dir_concordant
  rep_direction <- vapply(seq_along(all_coords), function(i) {
    sig_idx <- which(sig_mat[i, ])
    if (length(sig_idx) == 0L) return(NA_real_)
    best <- sig_idx[which.min(fdr_mat[i, sig_idx])]
    dir_mat[i, best]
  }, numeric(1))
  data.frame(EventType         = et,
             coord_nostrand    = all_coords[is_sig_any],
             n_sig_div         = n_sig_div[is_sig_any],
             shared            = is_shared[is_sig_any],
             shared_discordant = is_shared_discordant[is_sig_any],
             direction         = rep_direction[is_sig_any],
             stringsAsFactors = FALSE)
})
share_df <- do.call(rbind, share_rows)
cat(sprintf("Panel A: FDR<%g events per type (>=1 DIV): %s\n",
            FDR_THRESH_SHARE,
            paste(sprintf("%s=%d", EVENT_TYPES_A,
                          table(factor(share_df$EventType, levels = EVENT_TYPES_A))),
                  collapse = ", ")))

count_summary <- share_df %>%
  group_by(EventType) %>%
  summarise(n_shared      = sum(shared),
            n_discordant  = sum(shared_discordant),
            n_singleton   = sum(!shared & !shared_discordant),
            n_total       = n(), .groups = "drop")
count_order <- count_summary$EventType[order(-count_summary$n_total)]
cat(sprintf("Panel A: total unique FDR<0.05 rMATS events, any DIV, all event types: %d\n",
            sum(count_summary$n_total)))
cat(sprintf("Panel A: total shared+concordant (>=2 DIVs) events, all types: %d\n",
            sum(count_summary$n_shared)))

# Is the concordance rate among shared (>=2 DIV) events higher than chance?
# Null: for an event significant in k DIVs, direction in each DIV is an
# independent fair coin flip, so P(all k agree) = 2^(1-k). Exact
# Poisson-binomial test (no normal approximation), pooling n_sig_div = 2,3,4.
conc_test_df <- share_df %>%
  filter(n_sig_div >= 2L) %>%
  mutate(p_null = 2 ^ (1 - n_sig_div))

# Plotmath expression string (geom_text(parse=TRUE)), not composed Unicode
# superscript characters -- Unicode superscript digits span two Unicode
# blocks with inconsistent glyph metrics per font (legacy Latin-1 Supplement
# U+00B9/B2/B3 for "1/2/3" vs the Superscripts-and-Subscripts block
# U+2074-2079 for "4-9"), so e.g. exponent "18" rendered as "¹⁸" mixes a
# Latin-1 "¹" with a Superscripts-block "⁸" and the two commonly differ in
# size/baseline. plotmath draws true nested superscripts from the font's
# normal digit glyphs instead, so it doesn't have this problem. Console
# cat() text (not plotted) keeps a plain-text rendering of the exponent.
format_p_sci <- function(p, lt = FALSE, plotmath = TRUE) {
  exponent <- floor(log10(p))
  mantissa <- p / 10^exponent
  rel <- if (lt) "<" else "=="
  if (plotmath) sprintf('P %s %.1f %%*%% 10^%d', rel, mantissa, exponent)
  else sprintf("P %s %.1f x 10^%d", if (lt) "<" else "=", mantissa, exponent)
}

observed <- sum(conc_test_df$shared)
expected <- sum(conc_test_df$p_null)
p_raw    <- 1 - ppoibin(observed - 1L, pp = conc_test_df$p_null, method = "RF")
eps      <- .Machine$double.eps
p_label  <- if (p_raw < eps) format_p_sci(eps, lt = TRUE, plotmath = FALSE) else format_p_sci(p_raw, plotmath = FALSE)
cat("Panel A: exact Poisson-binomial concordance-vs-chance test (shared, >=2 DIV events), overall:\n")
cat(sprintf("n_tested=%d, observed=%d, expected=%.1f, %s\n",
            nrow(conc_test_df), observed, expected, p_label))

count_long <- rbind(
  data.frame(EventType = count_summary$EventType,
            Status = "FDR < 0.05 in > 1 DIV, concordant direction of effect", Count = count_summary$n_shared),
  data.frame(EventType = count_summary$EventType,
            Status = "FDR < 0.05 in > 1 DIV, discordant direction of effect", Count = count_summary$n_discordant),
  data.frame(EventType = count_summary$EventType,
            Status = "FDR < 0.05 in 1 DIV", Count = count_summary$n_singleton)
)
count_long$EventType <- factor(count_long$EventType, levels = count_order)
count_long$Status <- factor(count_long$Status,
                            levels = c("FDR < 0.05 in 1 DIV",
                                       "FDR < 0.05 in > 1 DIV, discordant direction of effect",
                                       "FDR < 0.05 in > 1 DIV, concordant direction of effect"))

PAL_COUNT_A <- c("FDR < 0.05 in > 1 DIV, concordant direction of effect" = PAL_CONC,
                 "FDR < 0.05 in > 1 DIV, discordant direction of effect" = PAL_TLOOP,
                 "FDR < 0.05 in 1 DIV"             = PAL_DISC)

fig_a <- ggplot(count_long, aes(x = EventType, y = Count, fill = Status)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = PAL_COUNT_A, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "FDR<0.05 events\n(union across 4 DIVs)") +
  guides(fill = guide_legend(nrow = 3)) +
  pub_theme +
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"),
        legend.text     = element_text(size = 6),
        legend.box.just = "left")

# ---- Panel B: OR (event type vs other, shared in >=2 DIVs, concordant) -----
# Moved here from figure4.R -- reuses this script's own share_df (panel A's
# underlying data) directly, no extra data pull needed.
or_rows <- lapply(EVENT_TYPES_A, function(et) {
  is_et <- share_df$EventType == et
  tab <- matrix(c(sum(is_et  & share_df$shared),  sum(is_et  & !share_df$shared),
                  sum(!is_et & share_df$shared),  sum(!is_et & !share_df$shared)),
               nrow = 2, byrow = TRUE,
               dimnames = list(c(et, "Other"), c("shared", "not_shared")))
  ft <- fisher.test(tab)
  data.frame(EventType = et, OR = unname(ft$estimate),
             ci_lo = ft$conf.int[1], ci_hi = ft$conf.int[2],
             p = ft$p.value, stringsAsFactors = FALSE)
})
or_df <- do.call(rbind, or_rows)
or_df$EventType <- factor(or_df$EventType, levels = or_df$EventType[order(or_df$OR)])
or_df$p_label <- vapply(or_df$p, format_p_sci, character(1))
# A3SS and MXE (both low OR) sit close enough in y that their default-
# centred labels collide -- nudge them apart horizontally by name (order()
# above already reassigned x position by ascending OR, so this can't be a
# fixed-length vector matched to row order).
LABEL_X_NUDGE <- c(A3SS = -0.2, MXE = 0.2, SE = 0, RI = 0, A5SS = 0)
or_df$label_x_nudge <- LABEL_X_NUDGE[as.character(or_df$EventType)]
cat("Panel B: OR (event type vs other, shared in >=2 DIVs, concordant direction):\n")
print(or_df[, c("EventType", "OR", "ci_lo", "ci_hi", "p")])

fig_b_or <- ggplot(or_df, aes(x = EventType, y = OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.5, linewidth = 0.5,
                  colour = PAL_CONC) +
  geom_text(aes(x = as.numeric(EventType) + label_x_nudge, y = ci_hi, label = p_label),
           vjust = -0.6, size = 2.6, colour = PAL_CONC, parse = TRUE) +
  scale_x_discrete(expand = expansion(add = 1.0)) +
  scale_y_continuous(breaks = seq(0.6, 1.8, by = 0.2),
                     labels = scales::number_format(accuracy = 0.1),
                     expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "OR (this type vs other,\nshared in >=2 DIVs)") +
  pub_theme

# ---- Shared setup for missingness table: 95-event candidate pool ------------
candidate_tbl <- read_csv(CANDIDATE_CSV, show_col_types = FALSE)
candidate_tbl <- candidate_tbl[candidate_tbl$concordant, ]
candidate_tbl$delta_wb_ly <- rowMeans(
  cbind(candidate_tbl$wb_delta_logitPSI_Tloop, candidate_tbl$lymph_delta_logitPSI_Tloop),
  na.rm = TRUE)
n_pool <- nrow(candidate_tbl)

# Novel/annotated splice-site classification (gencode v49) -- an event is
# novel if any of its 6 A5SS coordinates isn't an actual exon boundary
# anywhere in the GTF for that chr/strand. Same logic as figure1.R's
# event_is_novel(), reimplemented locally (that version is GEL-only).
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
  pos_long  <- data.frame(key = c(exons$key, exons$key), pos = c(exons$start, exons$end))
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
candidate_tbl$is_novel <- vapply(seq_len(n_pool), function(i) check_novel(candidate_tbl[i, ]), logical(1))
cat(sprintf("Candidate pool novel events: %d/%d (not in gencode v49)\n",
            sum(candidate_tbl$is_novel), n_pool))

n_quant_per_event <- rowSums(sapply(DIVS_LABEL, function(d)
  candidate_tbl[[paste0("ipsc_quantified_", d)]]))
cat(sprintf("Quantified in all 4 DIVs: %d/%d; in >=1 DIV: %d/%d\n",
            sum(n_quant_per_event == 4L), n_pool, sum(n_quant_per_event >= 1L), n_pool))

# Restrict to events with >=1 unquantified DIV for the missingness table.
n_pool_total  <- n_pool
gap_mask      <- n_quant_per_event < 4L
candidate_tbl <- candidate_tbl[gap_mask, ]
n_quant_per_event <- n_quant_per_event[gap_mask]
n_pool        <- nrow(candidate_tbl)
cat(sprintf("Missingness table: %d/%d events with >=1 unquantified DIV\n", n_pool, n_pool_total))

# ---- Missingness table (not a figure panel) ---------------------------------
# Categorised at the gene level; genes with multiple distinct events in the
# pool may contribute rows to more than one category if their events behave
# differently.
TPM_LIMITED_GENES <- c("ALPK1", "ITGB7", "KCNK6", "RAC2", "TRAF3IP3", "WDFY4")
NEVER_QUANT_GENES <- c("DNPEP", "PGS1")

candidate_tbl$category <- with(candidate_tbl, ifelse(
  gene_symbol %in% NEVER_QUANT_GENES, "Never quantified (any DIV) despite solid TPM",
  ifelse(gene_symbol %in% TPM_LIMITED_GENES, "TPM-limited (gene expression drops)",
  ifelse(gene_symbol == "MADD", "Depth-limited, unresolved (stable per-replicate counts, cause unclear)",
         "Depth-limited (minor isoform rare, plausible sampling dropout)"))))

missing_summary <- candidate_tbl[, c("event_id", "gene_symbol", "is_novel", "category")]
missing_summary$n_DIVs_quantified <- n_quant_per_event
write_csv(missing_summary, MISSING_CSV)
cat(sprintf("\nMissingness summary (%d/%d events with >=1 unquantified DIV) written to %s:\n",
            nrow(missing_summary), n_pool_total, MISSING_CSV))
print(missing_summary %>% count(category), n = Inf)

# ---- Combine into Figure S8 with patchwork ----------------------------------
fig_s7 <- (fig_a | fig_b_or) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag          = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1))

fig_w_mm <- 170
fig_h_mm <- 75

ggsave(OUT_PDF, plot = fig_s7, width = fig_w_mm, height = fig_h_mm, units = "mm", bg = "white")
ggsave(OUT_PNG, plot = fig_s7, width = fig_w_mm, height = fig_h_mm, units = "mm", dpi = 600, bg = "white")

cat(sprintf("Figure S8 saved to: %s and %s\n", OUT_PDF, OUT_PNG))

