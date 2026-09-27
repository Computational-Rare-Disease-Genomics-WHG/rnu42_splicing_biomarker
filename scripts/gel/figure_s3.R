#!/usr/bin/env Rscript
# Layout: patchwork only -- never gridExtra arrangeGrob/grid.arrange
# =============================================================================
# figure_s3.R
# Supplementary Figure 1 (RNU4-2 splicing biomarker manuscript):
# Cross-cohort delta logit-PSI correlation for all sig-in-both events
# (FDR<0.05 GEL and Nava, all event types), using all ReNU cases in each
# cohort. 1 row x 5 columns (SE, A5SS, A3SS, MXE, RI). Points coloured
# by WB/Lymphocyte direction concordance.
#
# Inputs:
#   outputs/psi_matrix.rds          (GEL ijc/sjc/psi/events/samples)
#   outputs/nava_psi_matrix.rds     (Nava ijc/sjc/psi/events/samples)
#
# Outputs:
#   outputs/figures/figure_s3.pdf
#
# Run from workspace root: pixi run Rscript scripts/gel/figure_s3.R
# =============================================================================

# ---- libpath setup (pixi-aware) ---------------------------------------------
invisible(local({
  user_lib <- path.expand("~/re_gecip/machine_learning/R_packages")
  if (nzchar(Sys.getenv("PIXI_ENVIRONMENT_NAME"))) {
    keep <- !normalizePath(.libPaths(), mustWork = FALSE) %in%
            normalizePath(user_lib,    mustWork = FALSE)
    if (any(!keep)) {
      message("[pixi] removing ", sum(!keep),
              " GEL-system libpath entries to avoid ABI mismatch")
      .libPaths(.libPaths()[keep])
    }
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

suppressPackageStartupMessages(library(ggplot2))

# ---- Configuration ----------------------------------------------------------
PSI_PATH      <- Sys.getenv("GEL_PSI",  "outputs/psi_matrix.rds")
NAVA_PSI_PATH <- Sys.getenv("NAVA_PSI", "outputs/nava_psi_matrix.rds")
OUT_DIR_FIG   <- "outputs/figures"
OUT_FIG_S3    <- file.path(OUT_DIR_FIG, "figure_s3.pdf")

N_GEL_CASES     <- 11L
N_GEL_CONTROLS  <- 33L
N_NAVA_CASES    <- 19L
N_NAVA_CONTROLS <- 21L

EVENT_TYPES <- c("SE", "A5SS", "A3SS", "MXE", "RI")
FDR_THRESH  <- 0.05

# Paper display labels
WB_LAB <- "Whole-Blood"
LY_LAB <- "Lymphocytes"

# Concordant / discordant point colours (matches figure1.R panel D)
PAL_CONC <- "#444444"
PAL_DISC <- "#BBBBBB"

stopifnot(
  "GEL PSI matrix not found"  = file.exists(PSI_PATH),
  "Nava PSI matrix not found" = file.exists(NAVA_PSI_PATH)
)
dir.create(OUT_DIR_FIG, recursive = TRUE, showWarnings = FALSE)

# ---- Load inputs ------------------------------------------------------------
cat("Loading GEL PSI:   ", PSI_PATH,      "\n", sep = "")
gel_obj  <- readRDS(PSI_PATH)
cat("Loading Nava PSI:  ", NAVA_PSI_PATH, "\n", sep = "")
nava_obj <- readRDS(NAVA_PSI_PATH)

stopifnot(
  all(c("ijc", "sjc", "psi", "events", "samples") %in% names(gel_obj)),
  all(c("ijc", "sjc", "psi", "events", "samples") %in% names(nava_obj))
)

gel_ijc    <- gel_obj$ijc
gel_sjc    <- gel_obj$sjc
gel_events <- gel_obj$events

nava_ijc    <- nava_obj$ijc
nava_sjc    <- nava_obj$sjc
nava_events <- nava_obj$events

stopifnot(
  identical(rownames(gel_ijc),  rownames(gel_obj$psi)),
  identical(rownames(nava_ijc), rownames(nava_obj$psi))
)

gel_case_ids  <- paste0("gel:case:",     seq_len(N_GEL_CASES))
gel_ctrl_ids  <- paste0("gel:control:",  seq_len(N_GEL_CONTROLS))
nava_case_ids <- paste0("nava:case:",    seq_len(N_NAVA_CASES))
nava_ctrl_ids <- paste0("nava:control:", seq_len(N_NAVA_CONTROLS))

stopifnot(
  all(gel_case_ids  %in% rownames(gel_obj$psi)),
  all(gel_ctrl_ids  %in% rownames(gel_obj$psi)),
  all(nava_case_ids %in% rownames(nava_obj$psi)),
  all(nava_ctrl_ids %in% rownames(nava_obj$psi))
)
cat(sprintf("GEL PSI:  %d samples x %d events\n", nrow(gel_obj$psi), ncol(gel_obj$psi)))
cat(sprintf("Nava PSI: %d samples x %d events\n", nrow(nava_obj$psi), ncol(nava_obj$psi)))

# ---- Identify sig-in-both events (all event types) --------------------------
gel_sig  <- gel_events[!is.na(gel_events$FDR)  & gel_events$FDR  < FDR_THRESH, ]
nava_sig <- nava_events[!is.na(nava_events$FDR) & nava_events$FDR < FDR_THRESH, ]

gel_sig$key  <- paste(gel_sig$EventType,  gel_sig$coord_nostrand,  sep = "|")
nava_sig$key <- paste(nava_sig$EventType, nava_sig$coord_nostrand, sep = "|")

both_keys <- intersect(gel_sig$key, nava_sig$key)
gel_m  <- gel_sig[match(both_keys,  gel_sig$key), ]
nava_m <- nava_sig[match(both_keys, nava_sig$key), ]

sig_df <- data.frame(
  key           = both_keys,
  gel_ev_id     = gel_m$event_id,
  nava_ev_id    = nava_m$event_id,
  EventType     = gel_m$EventType,
  GeneSymbol    = gel_m$GeneSymbol,
  IncLevelDiff  = gel_m$IncLevelDiff,
  stringsAsFactors = FALSE
)
cat(sprintf("Sig-in-both events: %d\n", nrow(sig_df)))
stopifnot("No sig-in-both events found" = nrow(sig_df) > 0L)

# ---- bb-logit PSI -----------------------------------------------------------
bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}

gel_lp_sig  <- bb_logit(gel_ijc[,  sig_df$gel_ev_id,  drop = FALSE],
                         gel_sjc[,  sig_df$gel_ev_id,  drop = FALSE])
nava_lp_sig <- bb_logit(nava_ijc[, sig_df$nava_ev_id, drop = FALSE],
                         nava_sjc[, sig_df$nava_ev_id, drop = FALSE])

# ---- Per-group mean delta logit-PSI -----------------------------------------
mean_dlogit <- function(lp_mat, case_ids, ctrl_ids) {
  case_mu <- colMeans(lp_mat[case_ids, , drop = FALSE], na.rm = TRUE)
  ctrl_mu <- colMeans(lp_mat[ctrl_ids, , drop = FALSE], na.rm = TRUE)
  dl <- case_mu - ctrl_mu
  dl[!is.finite(dl)] <- NA_real_
  dl
}

# ---- Plot -------------------------------------------------------------------
mm_to_in <- function(mm) mm / 25.4

pub_theme <- theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 8),
    axis.text        = element_text(colour = "black", size = 7),
    axis.title       = element_text(size = 8),
    legend.position  = "top",
    legend.key.size  = unit(0.3, "cm"),
    legend.text      = element_text(size = 8),
    panel.grid       = element_blank(),
    plot.margin      = margin(2, 4, 2, 2)
  )

# ---- All-cases dlogit (1 row x 5 event types) ------------------------------
# Whole ReNU case cohort in each dataset (all 11 GEL + all 19 Nava cases).
dl_gel_all  <- mean_dlogit(gel_lp_sig,  gel_case_ids,  gel_ctrl_ids)
dl_nava_all <- mean_dlogit(nava_lp_sig, nava_case_ids, nava_ctrl_ids)
all_df <- data.frame(
  key         = sig_df$key,
  EventType   = factor(sig_df$EventType, levels = EVENT_TYPES),
  GeneSymbol  = sig_df$GeneSymbol,
  dlogit_GEL  = dl_gel_all,
  dlogit_Nava = dl_nava_all,
  stringsAsFactors = FALSE)
all_df$concordant <- !is.na(all_df$dlogit_GEL) & !is.na(all_df$dlogit_Nava) &
                     sign(all_df$dlogit_GEL) == sign(all_df$dlogit_Nava) &
                     sign(all_df$dlogit_GEL) != 0

all_ann <- do.call(rbind, lapply(EVENT_TYPES, function(et) {
  sub <- all_df[all_df$EventType == et &
                !is.na(all_df$dlogit_GEL) & !is.na(all_df$dlogit_Nava), ]
  n <- nrow(sub)
  n_conc <- sum(sub$concordant)
  r2_val <- if (n >= 3L) {
    r <- cor(sub$dlogit_GEL, sub$dlogit_Nava, use = "complete.obs")
    if (is.finite(r)) r^2 else NA_real_
  } else NA_real_
  data.frame(EventType = et,
             n         = n,
             n_conc    = n_conc,
             pct_conc  = if (n > 0L) 100 * n_conc / n else NA_real_,
             r_squared = r2_val,
             label     = if (!is.na(r2_val)) sprintf("R^2=%.2f\nn=%d", r2_val, n)
                         else sprintf("n=%d", n),
             stringsAsFactors = FALSE)
}))
all_ann$EventType <- factor(all_ann$EventType, levels = EVENT_TYPES)

cat("\n--- Per-event-type concordance (all cases, sig-in-both) ---\n")
print(all_ann[, c("EventType", "n", "n_conc", "pct_conc", "r_squared")],
      row.names = FALSE, digits = 3)

ax_lim <- ceiling(
  max(abs(c(all_df$dlogit_GEL, all_df$dlogit_Nava)), na.rm = TRUE) * 10) / 10

fig_s1 <- ggplot(all_df, aes(x = dlogit_GEL, y = dlogit_Nava)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey40", linewidth = 0.3) +
  geom_point(aes(colour = concordant), size = 0.8, alpha = 0.85) +
  geom_text(data = all_ann, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, size = 2.2, colour = "black",
            lineheight = 0.9, inherit.aes = FALSE) +
  scale_colour_manual(
    values = c(`TRUE` = PAL_CONC, `FALSE` = PAL_DISC),
    labels = c(`TRUE` = "Concordant", `FALSE` = "Discordant"),
    name   = NULL, drop = FALSE) +
  facet_grid(. ~ EventType) +
  coord_cartesian(xlim = c(-ax_lim, ax_lim), ylim = c(-ax_lim, ax_lim)) +
  labs(x = sprintf("delta logit-PSI (%s)", WB_LAB),
       y = sprintf("delta logit-PSI (%s)", LY_LAB)) +
  pub_theme

# 170 mm wide x 60 mm tall (1 row x 5 cols)
pdf(OUT_FIG_S3, width = mm_to_in(170), height = mm_to_in(60))
print(fig_s1)
invisible(dev.off())
cat(sprintf("Supplementary Figure S3 saved to: %s\n", OUT_FIG_S3))

cat(sprintf("\nEvent type breakdown (sig-in-both, n=%d total):\n", nrow(sig_df)))
print(table(sig_df$EventType))
