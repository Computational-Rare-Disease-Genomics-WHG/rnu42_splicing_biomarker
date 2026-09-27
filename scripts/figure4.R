#!/usr/bin/env Rscript
# =============================================================================
# figure4.R
# Paper Figure 4 (RNU4-2 splicing biomarker manuscript): cross-tissue
# validation story -- blood-derived ReNU splicing events recur in iPSC-
# derived neurons, the 6-event insT signature is especially reproducible,
# and the blood-trained signature classifies ReNU iPSC neurons correctly.
#
#   Panel A (row 1)      - exact hypergeometric test: is the pool's FDR<0.05
#                          rate in iPSC more than expected by chance? Null =
#                          random same-n A5SS events quantified in GEL WB +
#                          Nava + that DIV (no significance/concordance
#                          requirement). Expected vs observed % per DIV,
#                          with fold enrichment and exact p (phyper).
#   Panel B (row 2)      - 95-event candidate pool scatter, x = mean delta
#                          logit-PSI (WB+lymphocyte insT carriers - controls),
#                          y = delta logit-PSI (iPSC ReNU - KOLF), per DIV.
#                          Hollow/filled by FDR<0.05 in iPSC. R^2 (all) and
#                          R^2 (FDR<0.05 subset) annotated.
#   Panel C (row 3)      - same x/y, restricted to the 6 insT signature
#                          events only, gene-labelled, y=x line. R^2 per DIV.
#   Panel D (row 4)      - insT LASSO P(case) by cell line and DIV stage
#                          (KOLF vs ReNU, jittered points + group median
#                          line, faceted by DIV) -- the payoff: a signature
#                          trained entirely in blood classifies ReNU iPSC
#                          neurons correctly.
#
# The within-iPSC cross-DIV reproducibility OR forest (event type vs other,
# shared in >=2 DIVs) that used to open this figure as panel A moved to
# figure_s8.R (its panel B); the genome-wide event-count breakdown and
# quantifiability heatmap also live in the consolidated figure_s8.R.
#
# Inputs (local; no GEL data touched):
#   outputs/ipsc_insT_signature_scores.rds
#       built by scripts/score_ipsc_against_insT_signature.R; carries
#       per_sample (uniform_sum / lasso_prob) for panel D.
#   outputs/candidate_table/table1_ipsc_annotated.csv
#       Table 1 (GEL airlock-approved A5SS candidate pool export) joined
#       against the per-DIV iPSC rMATS (export_table_s4_ipsc.R):
#       coordinates, insT-signature flag, precomputed WB/lymphocyte delta
#       logit-PSI, and precomputed per-DIV iPSC delta logit-PSI/FDR/PValue,
#       for panels A/B/C. Regenerate that script first if Table 1 or the
#       iPSC rMATS data changes.
#
# Outputs:
#   outputs/figures/figure4.pdf
#   outputs/figures/figure4.png
#
# Run from project root:
#   Rscript scripts/score_ipsc_against_insT_signature.R   # first run only
#   Rscript scripts/export_table_s4_ipsc.R       # if Table 1 changed
#   Rscript scripts/figure4.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(patchwork)
  library(readr)
  library(dplyr)
})

# ---- Configuration ----------------------------------------------------------
SCORE_RDS      <- "outputs/ipsc_insT_signature_scores.rds"
CANDIDATE_CSV  <- "outputs/candidate_table/table1_ipsc_annotated.csv"

OUT_DIR_FIG    <- "outputs/figures"
OUT_FIG_4_PDF  <- file.path(OUT_DIR_FIG, "figure4.pdf")
OUT_FIG_4_PNG  <- file.path(OUT_DIR_FIG, "figure4.png")

DIVS_LABEL   <- c("DIV1", "DIV5", "DIV9", "DIV14")   # display

# Palette -- matches figure2.R / figure_s8.R conventions
PAL_TLOOP <- "#D55E00"   # insT signature (orange)
PAL_DISC  <- "#BBBBBB"   # KOLF / background grey
PAL_CONC  <- "#1F3A93"   # deep navy

stopifnot("iPSC insT score RDS not found" = file.exists(SCORE_RDS),
          "Table 1 iPSC-annotated table not found -- run scripts/export_table_s4_ipsc.R first" =
            file.exists(CANDIDATE_CSV))
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Shared theme ------------------------------------------------------------
pub_theme <- theme_classic(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 8, face = "bold"),
        axis.text        = element_text(size = 7),
        axis.title       = element_text(size = 8),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7),
        plot.margin      = margin(2, 2, 2, 2))

# ---- Shared setup for panels A/B/C: 95-event candidate pool -----------------
candidate_tbl <- read_csv(CANDIDATE_CSV, show_col_types = FALSE)
candidate_tbl <- candidate_tbl[candidate_tbl$concordant, ]
candidate_tbl$delta_wb_ly <- rowMeans(
  cbind(candidate_tbl$wb_delta_logitPSI_Tloop, candidate_tbl$lymph_delta_logitPSI_Tloop),
  na.rm = TRUE)
n_pool <- nrow(candidate_tbl)
cat(sprintf("Panels A/B/C: %d A5SS candidate pool events (%d insT signature)\n",
            n_pool, sum(candidate_tbl$is_insT_Tloop_signature)))

ipsc_delta_mat <- as.matrix(candidate_tbl[, paste0("ipsc_delta_logitPSI_", DIVS_LABEL)])
dimnames(ipsc_delta_mat) <- list(candidate_tbl$coord_nostrand, DIVS_LABEL)

scatter_d_base <- do.call(rbind, lapply(DIVS_LABEL, function(lbl) {
  data.frame(coord_nostrand = candidate_tbl$coord_nostrand,
             gene           = candidate_tbl$gene_symbol,
             is_insT        = candidate_tbl$is_insT_Tloop_signature,
             delta_wb_ly    = candidate_tbl$delta_wb_ly,
             delta_ip       = ipsc_delta_mat[, lbl],
             fdr_sig        = candidate_tbl[[paste0("ipsc_FDRsig_", lbl)]],
             quantified     = candidate_tbl[[paste0("ipsc_quantified_", lbl)]],
             direction      = candidate_tbl$direction,
             DIV            = factor(lbl, levels = DIVS_LABEL),
             stringsAsFactors = FALSE)
}))
scatter_d_base$concordant <- scatter_d_base$quantified & !is.na(scatter_d_base$delta_ip) &
  sign(scatter_d_base$delta_ip) == scatter_d_base$direction
scatter_d_base <- scatter_d_base[!is.na(scatter_d_base$delta_wb_ly) &
                                   is.finite(scatter_d_base$delta_ip), ]

pool_df <- scatter_d_base
pool_df$fill_col <- ifelse(pool_df$fdr_sig, PAL_DISC, NA)
inst_df <- scatter_d_base[scatter_d_base$is_insT, ]
inst_df$fill_col <- ifelse(inst_df$fdr_sig, PAL_TLOOP, NA)

ax_lim <- ceiling(max(abs(c(scatter_d_base$delta_wb_ly, scatter_d_base$delta_ip)),
                      na.rm = TRUE) * 10) / 10

# ---- Panel A: exact hypergeometric test -- is the pool's FDR<0.05 rate in
# iPSC more than expected by chance? Null: for each DIV, draw n_quant events
# (however many of the 95-event pool are QUANTIFIED in that DIV, e.g. 83 for
# DIV1) without replacement from the universe of A5SS events quantified in
# GEL WB + Nava lymphocyte + that iPSC DIV -- with NO requirement they be
# significant/concordant in blood+lymphocyte (unlike actual pool membership).
# That's exactly a hypergeometric urn draw (fixed finite population, fixed
# number of "successes" K = FDR<0.05 events in the universe), so tested via
# phyper() rather than simulation -- gives an exact p-value, and the
# expected count under the null is just n_quant * K/n_universe (so the
# expected PERCENTAGE, K/n_universe, doesn't even depend on n_quant).
#
# METHODOLOGICAL REQUIREMENT (must hold for this test to be valid): the
# hypergeometric universe must consist only of A5SS events that were
# actually ELIGIBLE to be selected into the candidate pool -- i.e. events
# quantifiable in GEL WB, Nava lymphocyte, AND that specific iPSC DIV
# simultaneously (universe_keys below, a 3-way intersection). An event never
# measurable in one of the three datasets could not have entered the
# candidate pool regardless of its true effect, so it must not be counted as
# a possible draw either -- including it would inflate n_universe with
# events that were never actually "in the running", diluting the background
# rate and overstating the enrichment. This is enforced by construction here
# (gel_wb_keys/nava_keys/ipsc$coord_nostrand are each themselves restricted
# to quantified events, not filtered by significance).
GEL_WB_DIR <- "data/gel_rmats_tables"
NAVA_DIR   <- "data/nava_rmats_tables"
DIVS_LC    <- c("div1", "div5", "div9", "div14")

build_coord_nostrand_a5ss <- function(df) {
  coord_key <- do.call(paste, c(df[, 4:11], sep = ":"))
  sub("^([^:]+):[+-]:(.*)", "\\1:\\2", coord_key)
}
load_gel_wb_a5ss <- function() {
  d <- read_csv(file.path(GEL_WB_DIR, "A5SS_summary.csv"), show_col_types = FALSE,
               col_types = cols(.default = col_character()))
  d$coord_nostrand <- build_coord_nostrand_a5ss(d)
  d
}
load_nava_a5ss <- function() {
  d <- read_tsv(file.path(NAVA_DIR, "A5SS.MATS.JC.txt"), show_col_types = FALSE,
               col_types = cols(.default = col_character()))
  d$coord_nostrand <- build_coord_nostrand_a5ss(d)
  d
}
load_ipsc_a5ss <- function(div_lc) {
  d <- read_tsv(file.path("data/splicing", div_lc, "A5SS.MATS.JC.txt"), show_col_types = FALSE,
               col_types = cols(.default = col_character()))
  d$coord_nostrand <- build_coord_nostrand_a5ss(d)
  d$FDR <- as.numeric(d$FDR)
  d
}

gel_wb_keys  <- unique(load_gel_wb_a5ss()$coord_nostrand)
nava_keys    <- unique(load_nava_a5ss()$coord_nostrand)
DIVS_LC_MAP  <- setNames(DIVS_LC, DIVS_LABEL)

hyper_result <- do.call(rbind, lapply(DIVS_LABEL, function(div_lbl) {
  ipsc <- load_ipsc_a5ss(DIVS_LC_MAP[[div_lbl]])
  universe_keys <- Reduce(intersect, list(gel_wb_keys, nava_keys, unique(ipsc$coord_nostrand)))
  m <- match(universe_keys, ipsc$coord_nostrand)
  is_sig_universe <- !is.na(ipsc$FDR[m]) & ipsc$FDR[m] < 0.05
  n_universe <- length(universe_keys)

  n_quant  <- sum(candidate_tbl[[paste0("ipsc_quantified_", div_lbl)]])
  observed <- sum(candidate_tbl[[paste0("ipsc_FDRsig_", div_lbl)]], na.rm = TRUE)
  K        <- sum(is_sig_universe)   # successes (FDR<0.05) in the universe

  p_exact <- phyper(observed - 1L, m = K, n = n_universe - K, k = n_quant, lower.tail = FALSE)
  # Exact 95% interval for the null (hypergeometric) distribution -- 2.5th/
  # 97.5th quantiles of the count, converted to %.
  ci_lo <- qhyper(0.025, m = K, n = n_universe - K, k = n_quant)
  ci_hi <- qhyper(0.975, m = K, n = n_universe - K, k = n_quant)

  data.frame(DIV = div_lbl, n_quant = n_quant, observed = observed,
             observed_pct = 100 * observed / n_quant,
             expected_pct = 100 * K / n_universe,
             expected_ci_lo = 100 * ci_lo / n_quant,
             expected_ci_hi = 100 * ci_hi / n_quant,
             K = K, n_universe = n_universe, p = p_exact, stringsAsFactors = FALSE)
}))
hyper_result$DIV <- factor(hyper_result$DIV, levels = DIVS_LABEL)
hyper_result$fold_enrichment <- hyper_result$observed_pct / hyper_result$expected_pct
hyper_result$fold_label <- sprintf("%.1fx", hyper_result$fold_enrichment)

# "P = 1.46 x 10^-16" as a plotmath expression string (geom_text(parse=TRUE))
# rather than composed Unicode superscript characters -- Unicode superscript
# digits are split across two different blocks with inconsistent glyph
# metrics per font (legacy Latin-1 Supplement U+00B9/B2/B3 for "1/2/3" vs the
# Superscripts-and-Subscripts block U+2074-2079 for "4-9"), so e.g. "16"
# rendered as "¹⁶" mixes a Latin-1 "¹" with a Superscripts-block "⁶" and the
# two glyphs commonly differ in size/baseline. plotmath draws true nested
# superscripts from the font's normal digit glyphs, so it doesn't have this
# problem.
format_p_sci_a <- function(p) {
  if (p < 1e-300) { exponent <- -300L; mantissa <- 1 } else {
    exponent <- floor(log10(p))
    mantissa <- p / 10^exponent
  }
  sprintf('P == %.2f %%*%% 10^%d', mantissa, exponent)
}
hyper_result$p_label <- vapply(hyper_result$p, format_p_sci_a, character(1))

# DIV label + n used in the hypergeometric test (= n quantified in that DIV
# among the 95-event pool), stacked on a second line via the axis text itself.
hyper_result$DIV_label <- sprintf("%s\nn = %d", hyper_result$DIV, hyper_result$n_quant)
div_label_map <- setNames(hyper_result$DIV_label, hyper_result$DIV)

cat("Panel A: exact hypergeometric test (observed vs expected FDR<0.05 % among quantified pool events):\n")
print(hyper_result[, c("DIV", "n_quant", "observed", "observed_pct", "expected_pct",
                       "expected_ci_lo", "expected_ci_hi", "fold_enrichment", "p")])

LAB_EXPECTED <- "Expected\n(hypergeometric null)"
LAB_OBSERVED <- "Candidate\npool"
div_strip_map <- setNames(sprintf("**%s**<br>n = %d", hyper_result$DIV, hyper_result$n_quant),
                          hyper_result$DIV)

# Faceted by DIV (rather than dodged along a single x axis) so its columns
# line up with the facet_wrap(~ DIV) panels used in B/C/D below. Expected vs
# observed are the two x-axis categories within each facet -- their identity
# is already conveyed by the x-axis text, so no separate colour legend.
# "n" for the Expected point is the point-estimate count implied by its own
# percentage (expected_pct% of n_quant, i.e. the hypergeometric mean
# n_quant*K/n_universe rounded) -- not itself the direct qhyper() draw size
# (that's n_quant, already shown once via the shared strip label), but the
# same kind of "how many events is this percentage actually counting"
# read-out already given for the Candidate pool point.
hyper_long <- rbind(
  data.frame(DIV = hyper_result$DIV, type = LAB_EXPECTED, pct = hyper_result$expected_pct,
             ci_lo = hyper_result$expected_ci_lo, ci_hi = hyper_result$expected_ci_hi,
             n_sig = round(hyper_result$expected_pct / 100 * hyper_result$n_quant), stringsAsFactors = FALSE),
  data.frame(DIV = hyper_result$DIV, type = LAB_OBSERVED, pct = hyper_result$observed_pct,
             ci_lo = NA_real_, ci_hi = NA_real_,
             n_sig = hyper_result$observed, stringsAsFactors = FALSE))
hyper_long$DIV  <- factor(hyper_long$DIV, levels = DIVS_LABEL)
hyper_long$type <- factor(hyper_long$type, levels = c(LAB_EXPECTED, LAB_OBSERVED))
hyper_long$label <- sprintf("%.1f%%\nn = %d", hyper_long$pct, hyper_long$n_sig)

Y_TOP <- 95  # fixed panel ceiling -- gives the fold/p label a spot clear of
             # even the tallest candidate-pool point (max observed_pct ~43%)
fig_a <- ggplot(hyper_long, aes(x = type, y = pct, colour = type)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), colour = "grey60",
               width = 0.15, linewidth = 0.4) +
  geom_point(size = 2.2) +
  geom_text(aes(y = pct + 7, label = label), vjust = 0, size = 2.6, show.legend = FALSE) +
  geom_text(data = hyper_result, aes(y = Y_TOP, label = fold_label),
            x = 1.5, vjust = 1, size = 2.6, colour = "black", inherit.aes = FALSE) +
  geom_text(data = hyper_result, aes(y = Y_TOP - 6, label = p_label),
            x = 1.5, vjust = 1, size = 2.6, colour = "black", inherit.aes = FALSE, parse = TRUE) +
  scale_colour_manual(values = setNames(c("grey50", PAL_TLOOP), c(LAB_EXPECTED, LAB_OBSERVED))) +
  facet_wrap(~ DIV, nrow = 1, labeller = as_labeller(div_strip_map)) +
  scale_y_continuous(limits = c(0, Y_TOP), expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "FDR < 0.05\nin iPSCs (%)") +
  pub_theme +
  theme(strip.text       = element_markdown(size = 8, face = "plain"),
        axis.text.x      = element_text(size = 7, lineheight = 0.85,
                                        margin = margin(t = 2)),
        legend.position  = "none",
        plot.margin      = margin(14, 2, 2, 2))

# ---- Panel B: 95-event candidate pool scatter, per DIV ----------------------
ann_b <- do.call(rbind, lapply(DIVS_LABEL, function(div) {
  pool_sub <- pool_df[pool_df$DIV == div, ]
  sig_sub  <- pool_sub[pool_sub$fdr_sig, ]
  r2 <- function(d) if (nrow(d) >= 3L)
    stats::cor(d$delta_wb_ly, d$delta_ip, use = "complete.obs")^2 else NA_real_
  data.frame(DIV = factor(div, levels = DIVS_LABEL),
             label = c(sprintf("All: R² = %.2f", r2(pool_sub)),
                      sprintf("FDR < 0.05: R² = %.2f", r2(sig_sub))),
             vjust = c(1, 2.6), stringsAsFactors = FALSE)
}))

# Direction concordance is not plotted in B (already spoken for by the R^2 /
# open-vs-filled encoding) -- reported to console only, for the manuscript text.
concordance_by_div <- do.call(rbind, lapply(DIVS_LABEL, function(div) {
  pool_sub <- pool_df[pool_df$DIV == div, ]
  data.frame(DIV = div, concordant = sum(pool_sub$concordant), n = nrow(pool_sub),
             pct = 100 * sum(pool_sub$concordant) / nrow(pool_sub))
}))
cat("\nDirection concordance by DIV (for manuscript text, not plotted in Fig 4B):\n")
print(concordance_by_div, row.names = FALSE)

cat("\nPanel B R^2 per DIV:\n"); print(ann_b[, c("DIV", "label")], row.names = FALSE)

fig_b <- ggplot(pool_df, aes(x = delta_wb_ly, y = delta_ip)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(aes(fill = fill_col), colour = PAL_DISC, shape = 21, size = 1.8,
             stroke = 0.4, alpha = 0.75) +
  scale_fill_identity() +
  geom_text(data = ann_b, aes(x = -ax_lim, y = ax_lim, label = label, vjust = vjust),
            inherit.aes = FALSE, hjust = 0, size = 2.6, colour = "grey20") +
  facet_wrap(~ DIV, nrow = 1) +
  coord_cartesian(xlim = c(-ax_lim, ax_lim), ylim = c(-ax_lim, ax_lim)) +
  labs(x = "WB + lymphocyte mean Δ logit-PSI",
       y = "Δ logit-PSI\n(iPSC ReNU - KOLF)") +
  pub_theme +
  theme(strip.text = element_text(size = 8, face = "bold"))

# ---- Panel C: insT signature (6 events) scatter, per DIV --------------------
ann_c <- do.call(rbind, lapply(DIVS_LABEL, function(div) {
  inst_sub <- inst_df[inst_df$DIV == div, ]
  r2 <- if (nrow(inst_sub) >= 3L)
    stats::cor(inst_sub$delta_wb_ly, inst_sub$delta_ip, use = "complete.obs")^2 else NA_real_
  data.frame(DIV = factor(div, levels = DIVS_LABEL), label = sprintf("R² = %.2f", r2),
             stringsAsFactors = FALSE)
}))
cat("Panel C R^2 per DIV:\n"); print(ann_c[, c("DIV", "label")], row.names = FALSE)

# Dummy off-plot points solely to generate a single shared filled/hollow
# legend for the FDR<0.05 encoding used in both B and C -- attached here so
# it renders once, directly beneath the B/C pair (not duplicated per panel).
legend_fdr <- data.frame(
  delta_wb_ly = rep(-10 * ax_lim, 2), delta_ip = rep(-10 * ax_lim, 2),
  fdr_sig     = c(TRUE, FALSE),
  DIV         = factor(DIVS_LABEL[1], levels = DIVS_LABEL))

# Gene labels sit at a FIXED (x, y) slot that is IDENTICAL across all 4 DIV
# facets, with a thin leader line to the actual per-DIV point -- rather than
# letting ggrepel re-solve independently per facet (which put the same gene's
# label in a different spot in every panel, and could still land two labels
# on top of each other e.g. RPS15A/PTBP3 in DIV1, where their points are
# almost coincident). label_y is placed relative to the y=x diagonal
# (delta_wb_ly + dy) rather than an independent rank-based slot -- the
# diagonal is the same in every facet, and alternating the offset sign for
# x-adjacent genes (CTDSP2 below vs PTBP3 above; RPS15A below vs MAP4K4
# above) keeps their fixed label slots well clear of each other even though
# their x positions (and often their real per-DIV points) sit close
# together. TOP1 and MAP4K4 also get a small dx (moved diagonally up-left)
# since a pure vertical offset still landed on top of their own point in
# some DIVs (their per-DIV points scatter close to the diagonal, i.e. close
# to whatever y a vertical-only offset would land on). Six genes, hand-tuned
# offsets -- fixed signature, safe to hardcode.
DIAG_OFFSET   <- c(TOP1 = -0.9, CTDSP2 = -1.3, PTBP3 = 1.3,
                   RPS15A = -1.3, MAP4K4 = 1.0,  STX16 = 1.6)
DIAG_X_OFFSET <- c(TOP1 = -0.2, CTDSP2 = 0,    PTBP3 = 0,
                   RPS15A = 0,    MAP4K4 = -0.7, STX16 = 0)
# TOP1 sits right at the panel's left edge, so it has no room to be centred
# on its (already leftmost) label slot without the text's right half
# reaching back onto the point -- right-align it instead so the whole word
# sits to the point's left. Every other gene keeps centred text.
LABEL_HJUST <- c(TOP1 = 1, CTDSP2 = 0.5, PTBP3 = 0.5,
                 RPS15A = 0.5, MAP4K4 = 0.5, STX16 = 0.5)
gene_label_pos <- unique(inst_df[, c("gene", "delta_wb_ly")])
gene_label_pos$label_x <- gene_label_pos$delta_wb_ly + DIAG_X_OFFSET[gene_label_pos$gene]
gene_label_pos$label_y <- gene_label_pos$delta_wb_ly + DIAG_OFFSET[gene_label_pos$gene]
gene_label_pos$hjust   <- LABEL_HJUST[gene_label_pos$gene]
inst_df <- merge(inst_df, gene_label_pos[, c("gene", "label_x", "label_y", "hjust")], by = "gene")

fig_c <- ggplot(inst_df, aes(x = delta_wb_ly, y = delta_ip)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = label_x, yend = label_y), colour = PAL_TLOOP,
              linewidth = 0.2, alpha = 0.5) +
  geom_point(aes(fill = fill_col), colour = PAL_TLOOP, shape = 21, size = 2.2, stroke = 0.6) +
  scale_fill_identity() +
  geom_point(data = legend_fdr, aes(shape = fdr_sig), colour = "grey30", size = 1.8) +
  scale_shape_manual(name = NULL, values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "iPSC FDR < 0.05", `FALSE` = "iPSC FDR ≥ 0.05")) +
  geom_text(aes(x = label_x, y = label_y, label = gene, hjust = hjust),
           colour = PAL_TLOOP, size = 2.6, fontface = "italic") +
  geom_text(data = ann_c, aes(x = -ax_lim, y = ax_lim, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.6, colour = "grey20") +
  facet_wrap(~ DIV, nrow = 1) +
  coord_cartesian(xlim = c(-ax_lim, ax_lim), ylim = c(-ax_lim, ax_lim)) +
  labs(x = "WB + lymphocyte mean Δ logit-PSI",
       y = "Δ logit-PSI\n(iPSC ReNU - KOLF)") +
  pub_theme +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom", legend.key.size = unit(2.5, "mm"),
        legend.text = element_text(size = 7))

# ---- Panel D: insT LASSO P(case) by CellLine x DIV --------------------------
cat("Loading iPSC insT scores:  ", SCORE_RDS, "\n", sep = "")
sc <- readRDS(SCORE_RDS)
required <- c("signature", "logit_psi", "per_sample", "per_event")
stopifnot("score RDS missing required fields" = all(required %in% names(sc)))
per_samp <- sc$per_sample
per_samp$DIV_label <- factor(toupper(per_samp$DIV), levels = DIVS_LABEL)
per_samp$CellLine  <- factor(per_samp$group, levels = c("KOLF", "ReNU"))

CELL_COL_D <- c(KOLF = PAL_DISC, ReNU = PAL_TLOOP)

ann_d <- per_samp %>%
  group_by(DIV_label, CellLine) %>%
  summarise(median_p = median(lasso_prob), .groups = "drop")
ann_d$label <- sprintf("median = %.2f", ann_d$median_p)
cat("Panel D median P(case) by DIV/CellLine:\n"); print(ann_d, n = Inf)

# Faceted by DIV (rather than a dodged DIV x axis) so its columns line up
# with the facet_wrap(~ DIV) panels in A/B/C -- CellLine becomes the x
# position within each facet instead of a dodge offset along one shared axis.
fig_d <- ggplot(per_samp, aes(x = CellLine, y = lasso_prob, colour = CellLine)) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.85) +
  geom_errorbar(data = ann_d, aes(x = CellLine, y = median_p, ymin = median_p, ymax = median_p,
                                  colour = CellLine),
               inherit.aes = FALSE, width = 0.35, linewidth = 0.6) +
  geom_text(data = ann_d, aes(x = CellLine, y = median_p, label = label),
            vjust = -1.4, size = 2.2, show.legend = FALSE) +
  scale_colour_manual(values = CELL_COL_D, name = NULL) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  facet_wrap(~ DIV_label, nrow = 1) +
  labs(x = NULL, y = "T-loop signature P(case)") +
  pub_theme +
  theme(strip.text       = element_text(size = 8, face = "bold"),
        axis.text.x      = element_text(size = 7),
        legend.position  = "none")

# ---- Combine into Figure 4 with patchwork ----------------------------------
fig4 <- fig_a / fig_b / fig_c / fig_d +
  plot_layout(heights = c(1.05, 1, 1, 1.25)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag          = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1))

fig_w_mm <- 170
fig_h_mm <- 190

# cairo_pdf (not the default pdf() device) -- required for the Unicode
# characters used throughout (Δ, ², ≥) to embed correctly; the base pdf()
# device silently renders them as "?". Panel A's p-value/fold-enrichment
# labels are plotmath (parse=TRUE), not Unicode, but the rest of the figure
# still needs this.
ggsave(OUT_FIG_4_PDF, plot = fig4, width = fig_w_mm, height = fig_h_mm, units = "mm",
      bg = "white", device = cairo_pdf)
ggsave(OUT_FIG_4_PNG, plot = fig4, width = fig_w_mm, height = fig_h_mm, units = "mm", dpi = 1200, bg = "white")

cat(sprintf("Figure 4 saved to: %s and %s\n", OUT_FIG_4_PDF, OUT_FIG_4_PNG))
