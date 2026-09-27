#!/usr/bin/env Rscript
# Layout: patchwork only -- never gridExtra arrangeGrob/grid.arrange
# =============================================================================
# figure_s5.R
# Supplementary Figure S3 (RNU4-2 splicing biomarker manuscript).
# Raw TPM boxplots for the insT and Stem III signature genes not
# already shown elsewhere -- same style as figure2.R's STX16 panel /
# figure3.R's (excluded-from-final-layout) CINP panel.
#
#   Panel A -- insT signature genes, excluding STX16 and RPS15A (both
#              already in figure2.R Panel D): CTDSP2, PTBP3, MAP4K4, TOP1.
#              WB control vs WB n.64_65insT.
#   Panel B -- Stem III-only signature genes (excludes MAP4K4/TOP1, shared
#              with insT and already shown in Panel A): GNB1, KMT2C, CBFB,
#              CINP. WB control vs Stem III case.
#
# Wilcoxon (Mann-Whitney U) p-values shown per gene, BH-adjusted jointly
# within each panel (4 tests in A, 4 in B) -- same convention as
# figure_s6.R's per-site tests.
#
# Inputs:
#   outputs/loo_results_insT_lasso.rds
#   outputs/loo_results_stemIII_lasso.rds
#   outputs/psi_matrix.rds
#   outputs/a5ss_candidate_full_table.csv   (ANNOT_CSV; gene_symbol -> ensg)
#   data/salmon_quant_paths.csv             (platekey -> salmon quant path)
#
# Outputs:
#   outputs/figures/figure_s5.pdf
#   outputs/figure_s5_summary.rds   (gene, group ns, MWU p-values; airlock-safe)
#
# DATA PROTECTION: TPM values are GEL-derived. Gene names here are already
# public (part of the frozen insT/stemIII signatures) but treat the PDF +
# summary RDS as RE-internal until airlocked, per project convention.
#
# Run: pixi run Rscript scripts/gel/figure_s5.R
# =============================================================================

# ---- 0. libpath (pixi-aware) ------------------------------------------------
invisible(local({
  user_lib <- path.expand("~/re_gecip/machine_learning/R_packages")
  if (nzchar(Sys.getenv("PIXI_ENVIRONMENT_NAME"))) {
    keep <- !normalizePath(.libPaths(), mustWork = FALSE) %in%
            normalizePath(user_lib, mustWork = FALSE)
    if (any(!keep)) .libPaths(.libPaths()[keep])
    return(invisible(NULL))
  }
  v <- substr(as.character(getRversion()), 1, 3)
  candidate <- switch(v,
    "4.0" = "/tools/aws-workspace-apps/ce/R/4.0.2",
    "4.3" = "/tools/aws-workspace-ubuntu-apps/ce/R/4.3.3",
    "4.5" = "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.2",
    NULL)
  if (!is.null(candidate) && dir.exists(candidate)) .libPaths(c(.libPaths(), candidate))
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
}))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- 1. Config --------------------------------------------------------------
INST_LOO_PATH    <- Sys.getenv("INST_LOO_RDS",    "outputs/loo_results_insT_lasso.rds")
STEMIII_LOO_PATH <- Sys.getenv("STEMIII_LOO_RDS", "outputs/loo_results_stemIII_lasso.rds")
GEL_PSI_PATH     <- Sys.getenv("GEL_PSI",         "outputs/psi_matrix.rds")
ANNOT_PATH       <- Sys.getenv("ANNOT_CSV",       "outputs/a5ss_candidate_full_table.csv")
QUANT_PATHS_CSV  <- Sys.getenv("QUANT_PATHS",     "data/salmon_quant_paths.csv")

OUT_DIR_FIG <- "outputs/figures"
OUT_FIG_S3  <- file.path(OUT_DIR_FIG, "figure_s5.pdf")
OUT_SUMMARY <- "outputs/figure_s5_summary.rds"
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  "insT LOO RDS not found"     = file.exists(INST_LOO_PATH),
  "stemIII LOO RDS not found"  = file.exists(STEMIII_LOO_PATH),
  "GEL PSI matrix not found"   = file.exists(GEL_PSI_PATH),
  "Candidate annotation table not found" = file.exists(ANNOT_PATH),
  "Salmon quant paths CSV not found"     = file.exists(QUANT_PATHS_CSV)
)

PAL_CTRL    <- "#888888"
PAL_TLOOP   <- "#D55E00"   # insT / T-loop colour (matches figure2.R/figure3.R)
PAL_STEMIII <- "#F4A582"   # Stem III colour (matches figure2.R/figure3.R)

pub_theme <- theme_classic(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 9, face = "bold"),
        axis.text        = element_text(size = 9),
        axis.title       = element_text(size = 9),
        plot.title       = element_text(size = 9),
        legend.title     = element_text(size = 9),
        legend.text      = element_text(size = 9),
        plot.margin      = margin(2, 2, 2, 2))

mm_to_in <- function(mm) mm / 25.4

format_p_sci <- function(p) {
  if (p >= 0.001) return(sprintf("p == %.3f", p))
  exp10  <- floor(log10(p))
  mant_r <- round(p / 10^exp10, 1)
  if (mant_r >= 10) { mant_r <- mant_r / 10; exp10 <- exp10 + 1 }
  sprintf("p == %.1f %%*%% 10^%d", mant_r, exp10)
}

# =============================================================================
# 2. Sample groups
# =============================================================================
gel <- readRDS(GEL_PSI_PATH)
gel_ctrl_ids <- gel$samples$sample_id[gel$samples$group == "Control"]

ins <- readRDS(INST_LOO_PATH)
stopifnot(all(c("sample_ids_train", "cohort_train", "y_train", "ho_scores") %in% names(ins)))
gel_insT_train_ids <- ins$sample_ids_train[ins$cohort_train == "GEL" & ins$y_train == 1L]
gel_insT_ho_ids    <- ins$ho_scores$gel$sample_id
gel_insT_all_ids   <- c(gel_insT_train_ids, gel_insT_ho_ids)

s3 <- readRDS(STEMIII_LOO_PATH)
stopifnot(all(c("sample_ids_train", "cohort_train", "y_train") %in% names(s3)))
gel_stemIII_ids <- s3$sample_ids_train[s3$cohort_train == "GEL" & s3$y_train == 1L]

cat(sprintf("GEL controls: %d | insT cases: %d (train=%d, holdout=%d) | Stem III cases: %d\n",
            length(gel_ctrl_ids), length(gel_insT_all_ids),
            length(gel_insT_train_ids), length(gel_insT_ho_ids),
            length(gel_stemIII_ids)))

# =============================================================================
# 3. Gene -> ENSG lookup (from the candidate annotation table)
# =============================================================================
annot <- read.csv(ANNOT_PATH, stringsAsFactors = FALSE)
stopifnot(all(c("gene_symbol", "ensg", "is_insT_Tloop_signature",
               "is_stemIII_signature") %in% names(annot)))

lookup_ensg <- function(gene, sig_col) {
  rows <- annot[annot[[sig_col]] & !is.na(annot$gene_symbol) &
               toupper(annot$gene_symbol) == toupper(gene), , drop = FALSE]
  if (nrow(rows) == 0L)
    stop(sprintf("Gene %s not found among %s signature events", gene, sig_col))
  rows$ensg[1L]
}

insT_genes    <- c("CTDSP2", "PTBP3", "MAP4K4", "TOP1")   # excludes STX16, RPS15A (figure2.R)
stemIII_genes <- c("GNB1", "KMT2C", "CBFB", "CINP")                # excludes MAP4K4/TOP1 (Panel A), STX16 n/a

insT_ensg    <- setNames(vapply(insT_genes,    lookup_ensg, character(1L), sig_col = "is_insT_Tloop_signature"),  insT_genes)
stemIII_ensg <- setNames(vapply(stemIII_genes, lookup_ensg, character(1L), sig_col = "is_stemIII_signature"),     stemIII_genes)

# =============================================================================
# 4. TPM extraction (mirrors figure2.R's compute_gene_tpm_f2 /
#    figure3.R's compute_single_gene_tpm_s3)
# =============================================================================
compute_gene_tpm <- function(ensg_target, sample_ids, label) {
  qp <- read.csv(QUANT_PATHS_CSV, stringsAsFactors = FALSE)
  stopifnot("quant paths CSV must have columns platekey and file_path" =
              all(c("platekey", "file_path") %in% colnames(qp)))
  qp <- qp[qp$platekey %in% sample_ids, , drop = FALSE]
  if (nrow(qp) == 0L) return(NULL)
  res <- do.call(rbind, lapply(seq_len(nrow(qp)), function(i) {
    sid <- qp$platekey[i]
    sf  <- qp$file_path[i]
    if (!file.exists(sf)) return(NULL)
    d <- read.delim(sf, stringsAsFactors = FALSE)
    ensg_base <- sub("\\.\\d+$", "", d$Name)
    hit <- d[ensg_base == ensg_target, , drop = FALSE]
    if (nrow(hit) == 0L) return(NULL)
    data.frame(sample_id = sid, tpm = sum(hit$TPM, na.rm = TRUE), stringsAsFactors = FALSE)
  }))
  if (is.null(res) || nrow(res) == 0L) return(NULL)
  res$group <- label
  res
}

# Build one small-multiple TPM boxplot per gene; returns list(plot, p_value)
make_gene_panel <- function(gene, ensg, case_ids, case_label, case_colour) {
  tpm_case <- compute_gene_tpm(ensg, case_ids,     case_label)
  tpm_ctrl <- compute_gene_tpm(ensg, gel_ctrl_ids, "WB control")
  if (is.null(tpm_case) || is.null(tpm_ctrl)) {
    cat(sprintf("  %s (%s): TPM data unavailable -- panel omitted\n", gene, ensg))
    return(list(plot = plot_spacer(), p_value = NA_real_, gene = gene,
               n_case = 0L, n_ctrl = 0L))
  }
  tpm_df <- rbind(tpm_case, tpm_ctrl)
  tpm_df$group   <- factor(tpm_df$group, levels = c("WB control", case_label))

  wt <- wilcox.test(tpm_df$tpm[tpm_df$group == case_label],
                    tpm_df$tpm[tpm_df$group == "WB control"], exact = FALSE)

  cat(sprintf("  %s (%s): n_case=%d n_ctrl=%d Mann-Whitney U p=%.3e\n",
              gene, ensg, nrow(tpm_case), nrow(tpm_ctrl), wt$p.value))

  list(tpm_df = tpm_df, p_value = wt$p.value, gene = gene,
      n_case = nrow(tpm_case), n_ctrl = nrow(tpm_ctrl), case_label = case_label,
      case_colour = case_colour)
}

plot_gene_panel <- function(res, p_display, show_stats = TRUE) {
  if (is.null(res$tpm_df)) return(plot_spacer())
  y_annot <- max(res$tpm_df$tpm, na.rm = TRUE) * 1.05
  pal <- c("WB control" = PAL_CTRL, setNames(res$case_colour, res$case_label))
  # Case group is boxplotted only when there are enough samples for IQR to
  # be meaningful (>=4); at n<=3 (e.g. Stem III cases, n=2) show points +
  # median crossbar instead. Controls always get a boxplot.
  case_lab_local <- res$case_label
  small_case <- res$n_case <= 3L
  ggplot(res$tpm_df, aes(x = group, y = tpm, colour = group)) +
    geom_boxplot(data = function(d) {
      if (small_case) d[d$group != case_lab_local, ] else d
    }, fill = NA, outlier.shape = NA, width = 0.45, linewidth = 0.35) +
    (if (small_case)
      stat_summary(data = function(d) d[d$group == case_lab_local, ],
                   fun = median, geom = "crossbar",
                   width = 0.45, fatten = 0, linewidth = 0.35)
     else NULL) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.75) +
    (if (show_stats)
      list(annotate("segment", x = 1, xend = 2, y = y_annot, yend = y_annot,
                    colour = "grey30", linewidth = 0.35),
           annotate("text", x = 1.5, y = y_annot * 1.02,
                    label = format_p_sci(p_display),
                    parse = TRUE, size = 2.4, colour = "grey20", hjust = 0.5))
     else NULL) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(x = NULL, y = "TPM",
        title = bquote(italic(.(res$gene)))) +
    pub_theme +
    theme(plot.title = element_text(size = 9, hjust = 0.5, face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# =============================================================================
# 5. Panel A -- insT signature genes (excl. STX16)
# =============================================================================
cat("=== Panel A: insT signature genes (excl. STX16, RPS15A) ===\n")
panelA_res <- lapply(insT_genes, function(g)
  make_gene_panel(g, insT_ensg[[g]], gel_insT_all_ids, "WB n.64_65insT", PAL_TLOOP))
panelA_p_raw <- vapply(panelA_res, `[[`, numeric(1L), "p_value")
panelA_p_adj <- p.adjust(panelA_p_raw, method = "BH")
cat(sprintf("Panel A BH-adjusted p: %s\n",
            paste(sprintf("%s=%.3e", insT_genes, panelA_p_adj), collapse = ", ")))
panelA_plots <- Map(plot_gene_panel, panelA_res, panelA_p_adj)
panelA <- wrap_plots(panelA_plots, nrow = 1)

# =============================================================================
# 6. Panel B -- Stem III-only signature genes
# =============================================================================
cat("\n=== Panel B: Stem III-only signature genes ===\n")
panelB_res <- lapply(stemIII_genes, function(g)
  make_gene_panel(g, stemIII_ensg[[g]], gel_stemIII_ids, "Stem III case", PAL_STEMIII))
panelB_p_raw <- vapply(panelB_res, `[[`, numeric(1L), "p_value")
panelB_p_adj <- p.adjust(panelB_p_raw, method = "BH")
cat(sprintf("Panel B BH-adjusted p: %s\n",
            paste(sprintf("%s=%.3e", stemIII_genes, panelB_p_adj), collapse = ", ")))
panelB_plots <- Map(function(res, p) plot_gene_panel(res, p, show_stats = FALSE),
                    panelB_res, panelB_p_adj)
panelB <- wrap_plots(panelB_plots, nrow = 1)

# =============================================================================
# 7. Combine + write PDF (Genome Medicine: max 170mm wide, max 225mm tall)
# =============================================================================
# wrap_elements() treats each row of small multiples as a single opaque
# grob so the tag anchors to the row's own top-left corner -- tagging the
# wrap_plots() composite directly (without wrap_elements) tags one of the
# INNER sub-plots instead (the last one added), landing the label over
# the rightmost gene rather than the row as a whole.
fig_s3 <- (patchwork::wrap_elements(panelA) + labs(tag = "A")) /
          (patchwork::wrap_elements(panelB) + labs(tag = "B"))

pdf(OUT_FIG_S3, width = mm_to_in(170), height = mm_to_in(140))
print(fig_s3 & theme(plot.tag = element_text(size = 14, face = "bold"),
                     plot.tag.position = "topleft"))
dev.off()
cat(sprintf("\nSupplementary Figure S3 saved to: %s\n", OUT_FIG_S3))

# =============================================================================
# 8. Save summary (gene, group Ns, MWU p-values -- airlock-safe: gene names
#    are already public as part of the frozen signatures, no coordinates or
#    coefficients included)
# =============================================================================
summary_df <- rbind(
  data.frame(panel = "A", gene = insT_genes,
            n_case = vapply(panelA_res, `[[`, integer(1L), "n_case"),
            n_ctrl = vapply(panelA_res, `[[`, integer(1L), "n_ctrl"),
            p_raw = panelA_p_raw, p_bh = panelA_p_adj, stringsAsFactors = FALSE),
  data.frame(panel = "B", gene = stemIII_genes,
            n_case = vapply(panelB_res, `[[`, integer(1L), "n_case"),
            n_ctrl = vapply(panelB_res, `[[`, integer(1L), "n_ctrl"),
            p_raw = panelB_p_raw, p_bh = panelB_p_adj, stringsAsFactors = FALSE)
)
saveRDS(summary_df, OUT_SUMMARY)
cat(sprintf("Wrote %s\n", OUT_SUMMARY))
cat("\nDone.\n")
