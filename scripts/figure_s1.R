#!/usr/bin/env Rscript
# =============================================================================
# figure_s1.R
# Supplementary Figure S1: iPSC neuronal-identity marker gene expression.
# For each marker gene, log2(TPM+1) averaged across the three replicates
# within each (cell line, DIV) so there are two boxes per (gene, DIV) --
# KOLF (top) and ReNU (bottom) stacked above/below each other. Each gene
# therefore occupies two rows on the y-axis. Column axis = DIV. Gene groups
# laid out 3 rows x 2 columns via patchwork.
# Colour scale: reversed inferno, low TPM = lighter.
#
# Inputs:
#   data/gene_expression/gene_tpm_long.rds
#
# Outputs:
#   outputs/figures/figure_s1.pdf
#   outputs/figures/figure_s1.png
#
# Run from project root: Rscript scripts/figure_s1.R
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(patchwork)
})

TPM_RDS      <- "data/gene_expression/gene_tpm_long.rds"
OUT_DIR_FIG  <- "outputs/figures"
OUT_PDF      <- file.path(OUT_DIR_FIG, "figure_s1.pdf")
OUT_PNG      <- file.path(OUT_DIR_FIG, "figure_s1.png")

stopifnot("gene_tpm_long.rds not found" = file.exists(TPM_RDS))
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Marker gene groups ----------------------------------------------------
gene_groups <- tribble(
  ~gene_name, ~GeneGroup,
  "POU5F1",  "Pluripotency",
  "NANOG",   "Pluripotency",
  "SOX2",    "Pluripotency",
  "NEUROG2", "Neuronal identity",
  "TUBB3",   "Neuronal identity",
  "DCX",     "Neuronal identity",
  "MAP2",    "Neuronal identity",
  "RBFOX3",  "Neuronal identity",
  "BCL11B",  "Cortical neuronal identity",
  "SATB2",   "Cortical neuronal identity",
  "CUX1",    "Cortical neuronal identity",
  "CUX2",    "Cortical neuronal identity",
  "SLC17A7", "Glutamatergic",
  "SLC17A6", "Glutamatergic",
  "CAMK2A",  "Glutamatergic",
  "GAD1",    "GABAergic\n(negative control)",
  "GAD2",    "GABAergic\n(negative control)",
  "SLC6A1",  "GABAergic\n(negative control)",
  "SYN1",    "Synaptic maturation",
  "SNAP25",  "Synaptic maturation",
  "SYT1",    "Synaptic maturation",
  "DLG4",    "Synaptic maturation",
)
GROUP_LEVELS <- c("Pluripotency", "Neuronal identity", "Cortical neuronal identity",
                  "Glutamatergic", "GABAergic\n(negative control)", "Synaptic maturation")
gene_groups$GeneGroup <- factor(gene_groups$GeneGroup, levels = GROUP_LEVELS)

# ---- Load, subset, aggregate per (gene, DIV, CellLine) ---------------------
tpm_long <- readRDS(TPM_RDS)
marker_tpm <- tpm_long %>%
  filter(gene_name %in% gene_groups$gene_name) %>%
  inner_join(gene_groups, by = "gene_name") %>%
  mutate(log2TPM = log2(TPM + 1))
stopifnot("Some marker genes missing from gene_tpm_long.rds" =
            n_distinct(marker_tpm$gene_name) == nrow(gene_groups))

agg <- marker_tpm %>%
  group_by(GeneGroup, gene_name, DIV, CellLine) %>%
  summarise(mean_log2TPM = mean(log2TPM, na.rm = TRUE),
            n_rep        = dplyr::n(),
            .groups      = "drop")
stopifnot("Expected 3 replicates per (gene, DIV, CellLine)" =
            all(agg$n_rep == 3L))

# y-axis order: KOLF above ReNU within each gene; genes in gene_groups table
# order (top to bottom). ggplot draws first factor level at the bottom, so
# reverse the intended visual order.
row_order <- vapply(seq_len(nrow(gene_groups)), function(i) {
  c(paste(gene_groups$gene_name[i], "KOLF", sep = " | "),
    paste(gene_groups$gene_name[i], "ReNU", sep = " | "))
}, character(2L)) |> as.vector()
agg$row_lab <- paste(agg$gene_name, agg$CellLine, sep = " | ")
agg$row_lab <- factor(agg$row_lab, levels = rev(row_order))
agg$DIV     <- factor(agg$DIV, levels = c("DIV1", "DIV5", "DIV9", "DIV14"))

# ---- Theme + plot ---------------------------------------------------------
pub_theme <- theme_bw(base_size = 9) +
  theme(panel.grid       = element_blank(),
        panel.border     = element_rect(colour = "black", linewidth = 0.4),
        panel.spacing.x  = unit(1, "mm"),
        panel.spacing.y  = unit(1.5, "mm"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(size = 7, face = "bold", angle = 90),
        strip.text.x     = element_text(size = 8, face = "bold"),
        strip.placement  = "outside",
        axis.text.x      = element_text(size = 7),
        axis.text.y      = element_text(size = 6),
        axis.ticks       = element_blank(),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7),
        plot.margin      = margin(2, 2, 2, 2))

# Split into 2 columns of gene groups (3 + 3) so the figure is 3x2 rather
# than one tall single column. Shared fill scale via patchwork guides.
col1_groups <- GROUP_LEVELS[1:3]
col2_groups <- GROUP_LEVELS[4:6]

make_col <- function(agg_sub) {
  ggplot(agg_sub, aes(x = DIV, y = row_lab, fill = mean_log2TPM)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_viridis_c(name = expression(mean~log[2]*"(TPM+1)"),
                         option = "inferno", direction = -1,
                         limits = range(agg$mean_log2TPM, na.rm = TRUE)) +
    facet_grid(rows = vars(GeneGroup), scales = "free_y", space = "free_y",
               switch = "y") +
    labs(x = NULL, y = NULL) +
    pub_theme
}

fig_col1 <- make_col(agg[agg$GeneGroup %in% col1_groups, ])
fig_col2 <- make_col(agg[agg$GeneGroup %in% col2_groups, ])

fig_s6_alt <- (fig_col1 | fig_col2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

fig_w_mm <- 200
fig_h_mm <- 150

ggsave(OUT_PDF, plot = fig_s6_alt, width = fig_w_mm, height = fig_h_mm,
       units = "mm", bg = "white")
ggsave(OUT_PNG, plot = fig_s6_alt, width = fig_w_mm, height = fig_h_mm,
       units = "mm", dpi = 600, bg = "white")

cat(sprintf("Alt Figure S6 saved to: %s and %s\n", OUT_PDF, OUT_PNG))
