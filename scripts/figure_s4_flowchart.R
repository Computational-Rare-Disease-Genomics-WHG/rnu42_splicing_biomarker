#!/usr/bin/env Rscript
# =============================================================================
# figure_s4_flowchart.R
# Local, fast-iteration draft of Supplementary Figure S2 (RNU4-2 splicing
# biomarker manuscript). NOT the airlocked source of truth -- counts are
# hardcoded here (verified against a confirmed GEL Desktop run of
# scripts/gel/figure_s4.R, both chosen_lambda_label == "lambda.1se") so the
# layout/style can be iterated on locally without a GEL round trip. Figure
# S2 itself is being built manually in Illustrator -- this just gets the
# layout and numbers right first.
#
# Plain vertical flowchart (replaces an earlier proportional-sankey draft,
# abandoned -- ribbon widths spanning 5+ orders of magnitude made every node
# past "FDR<0.05" illegible, and there's no per-event transition data to
# honestly justify drawing overlap as flowing ribbons anyway). Every box is
# the same size; counts are annotated, not encoded in box size or ribbon
# width.
#
#   Shared trunk (top, single column):
#     Quantified in either cohort -> Quantified in both cohorts ->
#     Differentially spliced in both cohorts (FDR<0.05) -> A5SS events ->
#     Concordant direction
#   Split (same filtering logic, two independent genotype-specific case
#   sets -- NOT a partition of the concordant-direction event set; both
#   branches can and do draw on overlapping events, see
#   scripts/gel/figure_s4.R header):
#     n.64_65insT Signature Development     |  Stem III Signature Development
#     Depth >= 10                           |  Depth >= 10
#     Best event per gene                   |  Best event per gene
#     LASSO-selected sites                  |  LASSO-selected sites
#
# Output:
#   outputs/figures/figure_s4_flowchart.pdf / .png
#
# Run: Rscript scripts/figure_s4_flowchart.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

# ---- Verified counts (GEL Desktop run, both lambda.1se) --------------------
# "Quantified in either cohort" = 1,403,466 (n_gel_only + n_nava_only +
# n_shared -- the true union). scripts/gel/figure_s4.R previously used
# fig1$n_unique_combined for this, which is actually n_gel_only + n_nava_only
# ONLY (built in figure1.R for a different, cohort-exclusive-events
# calculation) -- that bug gave 1,154,036 and has been fixed in figure_s4.R.
TRUNK_LABS <- c("Quantified in\neither cohort", "Quantified in\nboth cohorts",
                "Differentially spliced\nin both cohorts\n(FDR < 0.05)",
                "A5SS events", "Concordant\ndirection")
TRUNK_N    <- c(1403466, 249430, 483, 105, 95)

BRANCH_LABS <- c("Depth ≥ 10", "Best event\nper gene", "LASSO-selected\nsites")
INST_HEADER <- "n.64_65insT\nSignature Development"
STEM_HEADER <- "Stem III\nSignature Development"
INST_N <- c(84, 77, 6)
STEM_N <- c(86, 52, 6)

# Sample sizes (fit scope: train-only cases + controls -- the set that
# actually produced best-event-per-gene ranking and the LASSO fit; see
# scripts/gel/figure_s4.R section 5 / loo_evaluate_*_lasso.R $sample_sizes.
# Verified GEL Desktop run: insT GEL 7+33, Nava 7+21 (14 cases, 54 controls);
# Stem III GEL 2+33, Nava 4+21 (6 cases, 54 controls; no holdout, so fit
# scope == depth-filter scope). Combined across cohorts, not split by cohort.
INST_SAMPLE_TXT <- "14 cases,\n54 controls"
STEM_SAMPLE_TXT <- "6 cases,\n54 controls"

fmt_n_compact <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

PAL_TRUNK <- "#7B94AD"   # lighter slate blue -- one colour for the whole shared funnel
PAL_INST  <- "#D55E00"   # matches figure2.R PAL_TLOOP
PAL_STEM  <- "#F4A582"   # matches figure2.R PAL_STEM3
TXT_LIGHT <- "black"
TXT_DARK  <- "black"

# =============================================================================
# Layout: 9 stacked levels (5 trunk + header + 3 branch steps), trunk
# centered at x = 0, branches at x = +-BRANCH_OFFSET below the split.
# =============================================================================
box_w <- 2.25; box_h <- 1.2; pitch <- 1.45
BRANCH_OFFSET <- 2.5
EXTRA_SPLIT_GAP <- 0.3   # small visual gap marking the split point

n_levels <- 5L + 4L
level_y  <- function(level) {
  y <- (n_levels - 1L - level) * pitch   # level 0 = top
  ifelse(level >= 5L, y - EXTRA_SPLIT_GAP, y)
}

trunk_df <- tibble(
  level = 0:4, x = 0, y = sapply(0:4, level_y),
  label = TRUNK_LABS, n = TRUNK_N, fill = PAL_TRUNK, text_col = TXT_LIGHT
)

branch_levels <- 5:8
inst_df <- tibble(
  level = branch_levels, x = -BRANCH_OFFSET, y = sapply(branch_levels, level_y),
  label = c(INST_HEADER, BRANCH_LABS), n = c(NA, INST_N),
  n_label_override = c(INST_SAMPLE_TXT, NA, NA, NA),
  fill = PAL_INST, text_col = TXT_LIGHT
)
stem_df <- tibble(
  level = branch_levels, x = BRANCH_OFFSET, y = sapply(branch_levels, level_y),
  label = c(STEM_HEADER, BRANCH_LABS), n = c(NA, STEM_N),
  n_label_override = c(STEM_SAMPLE_TXT, NA, NA, NA),
  fill = PAL_STEM, text_col = TXT_DARK
)
trunk_df <- trunk_df %>% mutate(n_label_override = NA_character_)

boxes <- bind_rows(trunk_df, inst_df, stem_df) %>%
  mutate(xmin = x - box_w / 2, xmax = x + box_w / 2,
         ymin = y - box_h / 2, ymax = y + box_h / 2,
         is_header = level == 5,
         n_label = ifelse(!is.na(n_label_override), n_label_override,
                    ifelse(is.na(n), NA_character_, sprintf("n = %s", fmt_n_compact(n)))),
         title_y  = y + ifelse(is_header, 0.3, ifelse(is.na(n_label), 0, 0.16)),
         nlabel_y = ifelse(is_header, y - 0.32, ymin + 0.2))

# ---- Straight vertical arrows within a column -------------------------------
vert_arrows <- bind_rows(
  boxes %>% filter(x == 0)              %>% arrange(level) %>%
    mutate(xend = x, yend = lead(y)) %>% filter(!is.na(yend)),
  boxes %>% filter(x == -BRANCH_OFFSET) %>% arrange(level) %>%
    mutate(xend = x, yend = lead(y)) %>% filter(!is.na(yend)),
  boxes %>% filter(x == BRANCH_OFFSET)  %>% arrange(level) %>%
    mutate(xend = x, yend = lead(y)) %>% filter(!is.na(yend))
) %>%
  transmute(x, xend, y = y - box_h / 2, yend = yend + box_h / 2)

# ---- Diagonal split arrows: concordant-direction -> each branch header -----
concordant <- boxes %>% filter(level == 4)
split_arrows <- bind_rows(
  tibble(x = concordant$x, y = concordant$y - box_h / 2,
         xend = -BRANCH_OFFSET, yend = level_y(5) + box_h / 2),
  tibble(x = concordant$x, y = concordant$y - box_h / 2,
         xend =  BRANCH_OFFSET, yend = level_y(5) + box_h / 2)
)

# =============================================================================
# Plot
# =============================================================================
p <- ggplot() +
  geom_segment(data = vert_arrows, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.09, "cm"), type = "closed"),
               linewidth = 0.4, color = "grey35") +
  geom_segment(data = split_arrows, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.09, "cm"), type = "closed"),
               linewidth = 0.4, color = "grey35") +
  geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            color = NA) +
  geom_text(data = boxes, aes(x = x, y = title_y, label = label, color = text_col),
            size = 3.0, lineheight = 0.92, fontface = "plain") +
  geom_text(data = boxes %>% filter(!is.na(n_label)),
            aes(x = x, y = nlabel_y, label = n_label, color = text_col),
            size = 3.0, fontface = "plain", lineheight = 0.92) +
  scale_fill_identity() +
  scale_color_identity() +
  coord_cartesian(xlim = c(-BRANCH_OFFSET - box_w / 2 - 0.2, BRANCH_OFFSET + box_w / 2 + 0.2),
                  ylim = c(min(boxes$ymin) - 0.2, max(boxes$ymax) + 0.2),
                  clip = "off") +
  theme_void(base_size = 9) +
  theme(plot.margin = margin(6, 6, 6, 6))

mm_to_in <- function(mm) mm / 25.4

fig_w_mm <- 170
fig_h_mm <- 165   # tightened pitch/EXTRA_SPLIT_GAP above shrank the layout;
                  # this is the smallest height that still keeps every box's
                  # text clear of its neighbours (checked by eye, print is
                  # not aspect-matched to data units -- box_h is print-only).

ggsave("outputs/figures/figure_s4_flowchart.pdf", p,
       width = mm_to_in(fig_w_mm), height = mm_to_in(fig_h_mm), units = "in")
ggsave("outputs/figures/figure_s4_flowchart.png", p,
       width = mm_to_in(fig_w_mm), height = mm_to_in(fig_h_mm), units = "in", dpi = 300)

cat("Wrote outputs/figures/figure_s4_flowchart.pdf / .png\n")
