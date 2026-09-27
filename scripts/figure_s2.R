#!/usr/bin/env Rscript
# =============================================================================
# figure_s2.R
# Paper Supplementary Figure S7 (RNU4-2 splicing biomarker manuscript):
# iPSC RNA-seq QC supporting Figure 3.
#
#   Panel A - 3x2 grid of per-sample QC metrics across the 24 iPSC samples
#             (8 cell-line/DIV combinations x 3 replicates). Violin + jitter.
#             STAR % aligned omitted (uniquely-aligned carries the same QC
#             signal).
#   Panel B - DESeq2 PCA (PC1 vs PC2 of variance-stabilised counts as written
#             by nf-core/rnaseq). Colour = DIV stage, shape = cell line.
#             Same palette + shape map as the QC report.
#
# Inputs (local; no GEL data touched):
#   data/qc/general_stats_table.tsv          - MultiQC general stats
#   data/qc/star_salmon_deseq2_pca-plot.tsv  - nf-core/rnaseq DESeq2 PCA
#
# Outputs:
#   outputs/figures/figure_s2.pdf
#   outputs/figures/figure_s2.png
#
# Run from project root: Rscript scripts/figure_s2.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# ---- Configuration ---------------------------------------------------------
QC_TSV       <- "data/qc/general_stats_table.tsv"
PCA_TSV      <- "data/qc/star_salmon_deseq2_pca-plot.tsv"

OUT_DIR_FIG  <- "outputs/figures"
OUT_PDF      <- file.path(OUT_DIR_FIG, "figure_s2.pdf")
OUT_PNG      <- file.path(OUT_DIR_FIG, "figure_s2.png")

stopifnot("QC stats TSV not found"  = file.exists(QC_TSV),
          "DESeq2 PCA TSV not found" = file.exists(PCA_TSV))
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# Palette - shared with QC report / figure3 conventions.
npg         <- pal_npg()(10)
cell_colors <- c(KOLF = npg[1], ReNU = npg[2])
div_colors  <- c(DIV1 = npg[3], DIV5 = npg[5],
                 DIV9 = npg[6], DIV14 = npg[4])
line_shapes <- c(KOLF = 16, ReNU = 17)

pub_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 8, face = "bold"),
        axis.text        = element_text(size = 7),
        axis.title       = element_text(size = 8),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7),
        plot.margin      = margin(2, 2, 2, 2))

# ---- Panel A: general QC metrics (3 x 3 grid) ------------------------------
# Duplicate column names in the MultiQC table; make.unique before reading.
raw_header      <- strsplit(readLines(QC_TSV, n = 1), "\t")[[1]]
col_names_fixed <- make.unique(raw_header, sep = "_2")

raw <- read_tsv(QC_TSV,
                col_names = col_names_fixed,
                skip      = 1,
                na        = ".",
                show_col_types = FALSE)

# Average paired _1/_2 trimming % into one per-sample value.
paired_rows <- raw %>%
  filter(str_detect(Sample, "_[12]$")) %>%
  mutate(SampleBase = str_remove(Sample, "_[12]$"))

paired_summary <- paired_rows %>%
  group_by(SampleBase) %>%
  summarise(`Trimmed bases (%)` = mean(`Trimmed bases`, na.rm = TRUE),
            .groups = "drop")

general_stats <- raw %>%
  filter(!str_detect(Sample, "_[12]$")) %>%
  left_join(paired_summary, by = c("Sample" = "SampleBase")) %>%
  rename(`STAR % uniquely aligned` = `Uniq aligned`,
         `Mapped reads (M)`        = `Reads mapped`,
         `Duplication (%)`         = `Duplication`,
         `Proper pairs (%)`        = `% Proper pairs`)

cat(sprintf("Panel A: %d samples loaded from %s\n",
            nrow(general_stats), QC_TSV))

# Helper: one violin + jitter for a single metric.
one_violin <- function(values, title, x_label,
                       fill_col = "#4C8BB5", xlim = NULL) {
  df <- tibble(x = values, y = "")
  ggplot(df, aes(x = x, y = y)) +
    geom_violin(fill = fill_col, colour = fill_col,
                alpha = 0.35, scale = "width", trim = TRUE,
                linewidth = 0.3) +
    geom_jitter(colour = fill_col, size = 1.2, alpha = 0.85,
                height = 0.08, width = 0) +
    coord_cartesian(xlim = xlim) +
    scale_x_continuous(labels = scales::comma) +
    labs(title = title, x = x_label, y = NULL) +
    pub_theme +
    theme(plot.title         = element_text(face = "bold", size = 8,
                                            margin = margin(b = 1)),
          axis.title.x       = element_text(size = 7, colour = "grey40"),
          axis.text.x        = element_text(size = 7),
          axis.text.y        = element_blank(),
          axis.ticks.y       = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92",
                                            linewidth = 0.3),
          panel.border       = element_rect(fill = NA, colour = "grey80",
                                            linewidth = 0.3),
          plot.margin        = margin(2, 5, 2, 5))
}

# NPG palette, one colour per metric for visual variety.
pal <- pal_npg()(6)

# 3x2 ordering (row-major): alignment yield -> library quality ->
# coverage/orientation -> trimming.
qc_panels <- list(
  one_violin(general_stats$`STAR % uniquely aligned`,
             "STAR % uniquely aligned", "% uniquely aligned",
             fill_col = pal[1], xlim = c(0, 100)),
  one_violin(general_stats$`Mapped reads (M)`,
             "Mapped reads (SAMtools)", "M reads mapped",
             fill_col = pal[2], xlim = c(0, NA)),
  one_violin(general_stats$`Duplication (%)`,
             "Duplication (Picard)",    "% duplicate reads",
             fill_col = pal[3], xlim = c(0, 100)),
  one_violin(general_stats$`5'-3' bias`,
             "5'-3' bias (Qualimap)",   "5'-3' bias ratio",
             fill_col = pal[4], xlim = c(0, NA)),
  one_violin(general_stats$`Proper pairs (%)`,
             "Proper pairs (SAMtools)", "% proper pairs",
             fill_col = pal[5], xlim = c(0, 100)),
  one_violin(general_stats$`Trimmed bases (%)`,
             "Trimmed bases (Cutadapt)","% bases trimmed",
             fill_col = pal[6], xlim = c(0, 100))
)

fig_a <- wrap_plots(qc_panels, ncol = 3L)

# ---- Panel B: DESeq2 PCA ---------------------------------------------------
pca_data <- read_tsv(PCA_TSV, show_col_types = FALSE) %>%
  rename(Sample = Name, PC1 = X, PC2 = Y) %>%
  mutate(CellLine = ifelse(str_detect(Sample, "KOLF"), "KOLF", "ReNU"),
         DIV      = factor(str_extract(Sample, "DIV\\d+"),
                           levels = c("DIV1", "DIV5", "DIV9", "DIV14")),
         Label    = str_remove(Sample, "^Sample-\\d+-"))

cat(sprintf("Panel B: %d samples loaded from %s\n",
            nrow(pca_data), PCA_TSV))

fig_b <- ggplot(pca_data,
                aes(x = PC1, y = PC2, colour = DIV, shape = CellLine)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_text_repel(aes(label = Label),
                  size          = 2.2,
                  max.overlaps  = Inf,
                  box.padding   = 0.35,
                  point.padding = 0.25,
                  segment.size  = 0.25,
                  segment.alpha = 0.5,
                  show.legend   = FALSE) +
  scale_colour_manual(values = div_colors,  name = "DIV stage") +
  scale_shape_manual( values = line_shapes, name = "Cell line") +
  labs(x = "PC1", y = "PC2") +
  pub_theme +
  theme(legend.position   = "right",
        legend.key.size   = unit(3, "mm"),
        panel.border      = element_rect(fill = NA, colour = "grey85",
                                         linewidth = 0.4))

# ---- Combine into Figure S7 ------------------------------------------------
# Layout: short 3x2 QC strip on top (full width), PCA below at full width,
# reduced height (no fixed aspect ratio -- fills the row's full width).
# wrap_elements() collapses the inner 3x2 patchwork so the outer A/B tagging
# skips over the inner panels.
fig_s6 <- (wrap_elements(full = fig_a) / fig_b) +
  plot_layout(heights = c(1, 1.8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag          = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1))

fig_w_mm <- 180
fig_h_mm <- 130

ggsave(OUT_PDF, plot = fig_s6,
       width = fig_w_mm, height = fig_h_mm, units = "mm", bg = "white")
ggsave(OUT_PNG, plot = fig_s6,
       width = fig_w_mm, height = fig_h_mm, units = "mm", dpi = 600,
       bg = "white")

cat(sprintf("Figure S7 saved to: %s and %s\n", OUT_PDF, OUT_PNG))
