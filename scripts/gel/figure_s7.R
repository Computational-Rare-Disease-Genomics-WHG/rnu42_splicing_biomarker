#!/usr/bin/env Rscript
# Layout: patchwork only -- never gridExtra arrangeGrob/grid.arrange
# =============================================================================
# figure_s7.R
# Paper Supplementary Figure S5 (RNU4-2 splicing biomarker manuscript):
# splice-site context for the insT and Stem III signatures.
#
#   Panel A - SpliceAI donor strengths (built when SpliceAI-by-chrom data +
#             A5SS.MATS.JC.txt are available; falls back to a placeholder).
#   Panel B - 5'SS sequence logos: insT signature (top) and Stem III signature
#             (bottom); increased-usage and decreased-usage donor windows shown
#             side by side; ggseqlogo bits scale, fallback bar chart if needed.
#             Moved from figure2.R (old Panel F).
#
# Primary inputs:
#   outputs/loo_results_insT_lasso.rds
#   outputs/loo_results_stemIII_lasso.rds
#   data/rmats/A5SS.MATS.JC.txt
#   <genome FASTA> + .fai
#
# Outputs:
#   outputs/figures/figure_s7.pdf
#   outputs/figure_s7_summary.rds   (counts + stats; airlock-safe)
#   outputs/figure_s7_cache/*.rds   (workspace-only; carry gene identifiers)
#
# Run from workspace root: pixi run Rscript scripts/gel/figure_s7.R
#
# DATA PROTECTION: selected event_ids and gene symbols are GEL protected.
# Summary RDS carries only counts / aggregate stats and is airlock-safe.
# =============================================================================

# ---- libpath setup (pixi-aware) --------------------------------------------
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
  if (!is.null(candidate) && dir.exists(candidate)) {
    .libPaths(c(.libPaths(), candidate))
  }
  if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
}))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

# ---- Configuration ----------------------------------------------------------
LOO_PATH         <- Sys.getenv("LOO_RDS",       "outputs/loo_results_insT_lasso.rds")
STEMIII_LOO_PATH <- Sys.getenv("STEMIII_LOO",   "outputs/loo_results_stemIII_lasso.rds")
WB_RMATS_DIR     <- Sys.getenv("GEL_RMATS_DIR", "data/rmats")
FASTA_PATH       <- Sys.getenv("FASTA_PATH",
  "/public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa")

OUT_DIR_FIG  <- "outputs/figures"
OUT_FIG_S7   <- file.path(OUT_DIR_FIG, "figure_s7.pdf")
OUT_SUMMARY  <- "outputs/figure_s7_summary.rds"
CACHE_DIR_S5 <- "outputs/figure_s7_cache"

A5SS_WB_FILE    <- file.path(WB_RMATS_DIR, "A5SS.MATS.JC.txt")
SPLICEAI_BY_CHROM <- Sys.getenv("SPLICEAI_BY_CHROM",
                                 "data/spliceai_ref_scores/by_chrom")

PAL_INC  <- "#D55E00"
PAL_DEC  <- "#0072B2"
PAL_CTRL <- "#888888"
PAL_STEMIII <- "#F4A582"

dir.create(OUT_DIR_FIG,  recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR_S5, recursive = TRUE, showWarnings = FALSE)

pub_theme <- theme_classic(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 8, face = "bold"),
        axis.text        = element_text(size = 7),
        axis.title       = element_text(size = 8),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 7),
        plot.margin      = margin(2, 2, 2, 2))

placeholder_panel <- function(label, msg) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 3, colour = "grey40") +
    theme_void() +
    theme(panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 0.3),
          plot.margin  = margin(2, 2, 2, 2))
}

p_to_stars <- function(p) {
  if (is.na(p))         "ns"
  else if (p <= 0.0001) "****"
  else if (p <= 0.001)  "***"
  else if (p <= 0.01)   "**"
  else if (p <= 0.05)   "*"
  else                  "ns"
}

# ---- Sequence helpers -------------------------------------------------------
revcomp <- function(seq) {
  comp <- chartr("ACGTNacgtn", "TGCANtgcan", seq)
  paste(rev(strsplit(comp, NULL, fixed = TRUE)[[1]]), collapse = "")
}

.fasta_handle_s4 <- new.env(parent = emptyenv())

get_fasta_s4 <- function() {
  if (!is.null(.fasta_handle_s4$f)) return(.fasta_handle_s4$f)
  if (!requireNamespace("Rsamtools", quietly = TRUE)) return(NULL)
  if (!file.exists(FASTA_PATH) || !file.exists(paste0(FASTA_PATH, ".fai"))) return(NULL)
  .fasta_handle_s4$f <- Rsamtools::FaFile(FASTA_PATH)
  .fasta_handle_s4$f
}

fasta_fetch_s4 <- function(chr, start_1, end_1) {
  fa <- get_fasta_s4()
  if (is.null(fa)) return(NA_character_)
  gr <- tryCatch(
    GenomicRanges::GRanges(seqnames = chr,
                           ranges   = IRanges::IRanges(start = start_1, end = end_1)),
    error = function(e) NULL)
  if (is.null(gr)) return(NA_character_)
  seq <- tryCatch(toupper(unname(as.character(Rsamtools::scanFa(fa, gr)))),
                  error = function(e) NA_character_)
  if (length(seq) != 1L || is.na(seq)) return(NA_character_)
  seq
}

extract_5ss_s4 <- function(chr, donor_1, strand) {
  seq <- if (strand == "+")
    fasta_fetch_s4(chr, donor_1 - 2L, donor_1 + 6L)
  else {
    s <- fasta_fetch_s4(chr, donor_1 - 6L, donor_1 + 2L)
    if (!is.na(s)) revcomp(s) else NA_character_
  }
  if (is.na(seq) || nchar(seq) != 9L) NA_character_ else seq
}

build_a5ss_eid_lookup_s4 <- function() {
  if (!file.exists(A5SS_WB_FILE)) return(NULL)
  df  <- read.delim(A5SS_WB_FILE, stringsAsFactors = FALSE)
  ck  <- do.call(paste, c(df[, 4L:11L], sep = ":"))
  cns <- sub("^([^:]+):[^:]+:", "\\1:", ck)
  gs  <- df[[3L]]
  out <- data.frame(
    event_id       = paste("A5SS", gs, cns, sep = ":"),
    GeneSymbol     = gs,
    coord_nostrand = cns,
    coord_key      = ck,
    chr            = df[[4L]], strand = df[[5L]],
    longES_0  = as.integer(df[[6L]]), longEE_1  = as.integer(df[[7L]]),
    shortES_0 = as.integer(df[[8L]]), shortEE_1 = as.integer(df[[9L]]),
    stringsAsFactors = FALSE)
  out[!duplicated(ck), , drop = FALSE]
}

compute_donor_seqs_s4 <- function(loo_path, id_prefix, a5ss_lk) {
  loo  <- readRDS(loo_path)
  cidx <- loo$selected_idx
  cand <- loo$a5ss_cand
  ev   <- cand$gel_ev_id[cidx]
  dir  <- as.integer(cand$direction[cidx])
  m    <- match(ev, a5ss_lk$event_id)
  if (anyNA(m)) {
    cat(sprintf("compute_donor_seqs_s4 (%s): %d/%d events not in A5SS lookup\n",
                id_prefix, sum(is.na(m)), length(ev)))
    return(NULL)
  }
  hit       <- a5ss_lk[m, , drop = FALSE]
  long_pos  <- ifelse(hit$strand == "+", hit$longEE_1,  hit$longES_0 + 1L)
  short_pos <- ifelse(hit$strand == "+", hit$shortEE_1, hit$shortES_0 + 1L)
  long_seq  <- mapply(extract_5ss_s4, hit$chr, long_pos,  hit$strand)
  short_seq <- mapply(extract_5ss_s4, hit$chr, short_pos, hit$strand)
  inc_seq   <- ifelse(dir == 1L, long_seq,  short_seq)
  dec_seq   <- ifelse(dir == 1L, short_seq, long_seq)
  data.frame(event_id_short = sprintf("%s_%02d", id_prefix, seq_along(ev)),
             chr = hit$chr, strand = hit$strand, direction = dir,
             inc_donor_pos_1 = ifelse(dir == 1L, long_pos, short_pos),
             dec_donor_pos_1 = ifelse(dir == 1L, short_pos, long_pos),
             inc_seq_9mer = inc_seq, dec_seq_9mer = dec_seq,
             stringsAsFactors = FALSE)
}

# ---- Cache helper -----------------------------------------------------------
build_or_load_cache_s4 <- function(cache_path, sources, build_fn, label) {
  src <- sources[file.exists(sources)]
  fresh <- file.exists(cache_path) && length(src) > 0L &&
           file.mtime(cache_path) >= max(vapply(src, file.mtime, numeric(1L)))
  if (fresh) { cat(sprintf("[cache] %s: hit\n", label)); return(readRDS(cache_path)) }
  cat(sprintf("[cache] %s: rebuilding...\n", label))
  result <- build_fn()
  saveRDS(result, cache_path)
  result
}

# ---- Logo helpers -----------------------------------------------------------
is_clean_9mer_s4 <- function(s) !is.na(s) & nchar(s) == 9L & !grepl("[^ACGT]", s)
pos_labels_s4    <- c("-3","-2","-1","+1","+2","+3","+4","+5","+6")

ic_per_pos_s4 <- function(seqs) {
  if (!length(seqs)) return(rep(NA_real_, 9L))
  m <- do.call(rbind, strsplit(seqs, NULL, fixed = TRUE))
  apply(m, 2L, function(col) {
    p <- table(factor(col, levels = c("A","C","G","T"))) / length(col)
    p <- p[p > 0]; 2 - (-sum(p * log2(p)))
  })
}

make_logo_panel_s4 <- function(seqs_df_in, sig_title) {
  if (is.null(seqs_df_in))
    return(placeholder_panel(sig_title, paste0(sig_title, "\nlogo unavailable")))
  inc_cl <- seqs_df_in$inc_seq_9mer[is_clean_9mer_s4(seqs_df_in$inc_seq_9mer)]
  dec_cl <- seqs_df_in$dec_seq_9mer[is_clean_9mer_s4(seqs_df_in$dec_seq_9mer)]
  cat(sprintf("Panel B %s: clean 9-mers inc=%d dec=%d\n",
              sig_title, length(inc_cl), length(dec_cl)))
  if (requireNamespace("ggseqlogo", quietly = TRUE) &&
      length(inc_cl) > 0L && length(dec_cl) > 0L) {
    p <- tryCatch({
      ggseqlogo::ggseqlogo(list(`Increased usage` = inc_cl,
                                `Decreased usage` = dec_cl),
                           method = "bits", ncol = 2L) +
        scale_x_continuous(breaks = 1:9, labels = pos_labels_s4) +
        labs(x = "Position relative to 5'SS", y = "Bits", title = sig_title) +
        pub_theme +
        theme(strip.text = element_text(size = 7),
              plot.title = element_text(size = 7, hjust = 0, face = "bold"))
    }, error = function(e) NULL)
    if (!is.null(p)) return(p)
  }
  # Fallback IC bar chart
  bits_df <- data.frame(
    pos   = rep(seq_len(9L), 2L),
    bits  = c(ic_per_pos_s4(inc_cl), ic_per_pos_s4(dec_cl)),
    group = rep(c("Increased usage", "Decreased usage"), each = 9L),
    stringsAsFactors = FALSE)
  bits_df$group <- factor(bits_df$group, levels = c("Increased usage", "Decreased usage"))
  ggplot(bits_df, aes(x = pos, y = bits, fill = group)) +
    geom_col(colour = "grey20", linewidth = 0.2) +
    facet_wrap(~ group, ncol = 2L) +
    scale_x_continuous(breaks = 1:9, labels = pos_labels_s4) +
    scale_fill_manual(values = c(`Increased usage` = PAL_INC,
                                 `Decreased usage` = PAL_DEC),
                      guide = "none") +
    labs(x = "Position relative to 5'SS", y = "Bits", title = sig_title) +
    pub_theme +
    theme(plot.title = element_text(size = 7, hjust = 0, face = "bold"))
}

# ---- A5SS event-ID lookup (shared by panels A and B) -----------------------
a5ss_lk_s4 <- tryCatch(build_a5ss_eid_lookup_s4(), error = function(e) {
  cat("A5SS lookup build failed:", conditionMessage(e), "\n")
  NULL
})

# ============================================================================
# ---- Panel A helpers: SpliceAI reference donor strengths -------------------
# ============================================================================

# Donor positions for events in a LOO object (all cand, or indexed subset).
# Returns data.frame(event_id, chr, strand, direction, inc_pos, dec_pos).
get_event_positions_s4 <- function(loo, a5ss_lk, idx = NULL) {
  cand <- loo$a5ss_cand
  if (is.null(cand)) return(NULL)
  if (!is.null(idx)) cand <- cand[idx, , drop = FALSE]
  ev  <- cand$gel_ev_id
  dir <- as.integer(cand$direction)
  m   <- match(ev, a5ss_lk$event_id)
  ok  <- !is.na(m)
  if (!any(ok)) return(NULL)
  hit       <- a5ss_lk[m[ok], , drop = FALSE]
  dv        <- dir[ok]
  long_pos  <- ifelse(hit$strand == "+", hit$longEE_1,  hit$longES_0 + 1L)
  short_pos <- ifelse(hit$strand == "+", hit$shortEE_1, hit$shortES_0 + 1L)
  data.frame(
    event_id = ev[ok],
    chr      = hit$chr,
    strand   = hit$strand,
    direction = dv,
    inc_pos  = ifelse(dv == 1L, long_pos,  short_pos),
    dec_pos  = ifelse(dv == 1L, short_pos, long_pos),
    stringsAsFactors = FALSE
  )
}

# Find zip file(s) for a given chromosome from the by_chrom directory.
# For split chromosomes (chr4/chr5), returns the top-level part zips
# (chr4_part1, chr4_part2) and excludes sub-split (part1a, part1b).
get_chr_zips_s4 <- function(chr_name, by_chrom_dir) {
  chr_num <- sub("^chr", "", chr_name)
  direct  <- file.path(by_chrom_dir,
                       sprintf("MANE_44_spliceAI_chr%s.zip", chr_num))
  if (file.exists(direct)) return(direct)
  # Split chromosomes: prefer partN over partNa/partNb (no trailing letter)
  parts <- sort(list.files(by_chrom_dir,
    pattern = sprintf("^MANE_44_spliceAI_chr%s_part[0-9]+\\.zip$", chr_num),
    full.names = TRUE))
  if (length(parts) > 0L) return(parts)
  character(0L)
}

# Scan by_chrom zip files to find donor_prob at each (chr, pos) in target_df.
# target_df must have chr_str (e.g. "chr12") and pos_int (integer).
# Returns named numeric vector "chr_str:pos_int" -> donor_prob.
build_spliceai_pos_lookup_s4 <- function(target_df, by_chrom_dir) {
  ud     <- unique(target_df[c("chr_str","pos_int")])
  keys   <- paste(ud$chr_str, ud$pos_int, sep = ":")
  result <- setNames(rep(NA_real_, length(keys)), keys)

  for (chr_name in unique(ud$chr_str)) {
    zips   <- get_chr_zips_s4(chr_name, by_chrom_dir)
    if (length(zips) == 0L) { cat(sprintf("  SpliceAI: no zip for %s\n", chr_name)); next }
    t_pos  <- ud$pos_int[ud$chr_str == chr_name]
    t_max  <- max(t_pos)
    t_keys <- paste(chr_name, t_pos, sep = ":")
    found  <- setNames(rep(FALSE, length(t_keys)), t_keys)

    for (zip_path in zips) {
      file_list <- unzip(zip_path, list = TRUE)$Name
      for (ef in file_list) {
        if (all(found)) break
        tryCatch({
          # Peek: read header + first data row to get min position
          con1     <- gzcon(unz(zip_path, ef))
          peek_raw <- readLines(con1, n = 2L)
          close(con1)
          if (length(peek_raw) < 2L) next
          peek_fields <- strsplit(peek_raw[2L], "\t", fixed = TRUE)[[1L]]
          if (length(peek_fields) < 2L) next
          min_pos <- suppressWarnings(as.integer(peek_fields[2L]))
          if (is.na(min_pos) || min_pos > t_max) next

          # Full read (readLines needed: gzcon is binary, read.table needs text)
          con2  <- gzcon(unz(zip_path, ef))
          lines <- readLines(con2)
          close(con2)
          dat   <- read.table(
            text      = lines[!startsWith(lines, "#")],
            header    = TRUE,
            col.names = c("chr","pos","ref_nt","acceptor_prob","donor_prob"),
            colClasses = c("character","integer","character","numeric","numeric"))
          hits <- dat[dat$pos %in% t_pos, , drop = FALSE]
          for (hi in seq_len(nrow(hits))) {
            k <- paste(chr_name, hits$pos[hi], sep = ":")
            if (!is.na(result[k])) next
            result[k] <- hits$donor_prob[hi]
            found[k]  <- TRUE
          }
        }, error = function(e) NULL)
      }
    }
    cat(sprintf("  SpliceAI %s: %d/%d positions found\n",
                chr_name, sum(found), length(t_pos)))
  }
  result
}

# ============================================================================
# ---- Panel A: SpliceAI reference donor strengths ---------------------------
# Two groups: full insT signature (n=6) and full Stem III signature (n=6);
# shared events (MAP4K4, TOP1) appear in BOTH bars.
# For each event: donor_prob at decreased 5'SS (left) and increased 5'SS (right).
# ============================================================================
fig_s4_panelA <- placeholder_panel("A", "SpliceAI donor strengths\nloading...")
panel_a_wide  <- NULL   # per-event wide table: group, event_id, inc_prob, dec_prob

if (!is.null(a5ss_lk_s4) && dir.exists(SPLICEAI_BY_CHROM) &&
    file.exists(LOO_PATH)) {
  tryCatch({
    SPLICEAI_CACHE_A <- file.path(CACHE_DIR_S5, "spliceai_panel_a_v3.rds")
    src_a <- c(LOO_PATH, STEMIII_LOO_PATH, A5SS_WB_FILE)
    src_a <- src_a[file.exists(src_a)]

    panel_a_wide <- if (
      file.exists(SPLICEAI_CACHE_A) && length(src_a) > 0L &&
      file.mtime(SPLICEAI_CACHE_A) >= max(vapply(src_a, file.mtime, numeric(1L)))) {
      cat("[cache] panel A SpliceAI: hit\n")
      readRDS(SPLICEAI_CACHE_A)
    } else {
      cat("[cache] panel A SpliceAI: building (scanning zip files)...\n")
      loo_a    <- readRDS(LOO_PATH)
      pos_inst <- get_event_positions_s4(loo_a, a5ss_lk_s4, loo_a$selected_idx)
      pos_siii <- if (file.exists(STEMIII_LOO_PATH)) {
        loo_siii <- readRDS(STEMIII_LOO_PATH)
        get_event_positions_s4(loo_siii, a5ss_lk_s4, loo_siii$selected_idx)
      } else NULL

      # Both groups keep their full signature membership; shared events
      # (MAP4K4, TOP1) appear in BOTH insT and Stem III bars.

      # Bulk position lookup across both groups
      mk_tgt <- function(pd) {
        if (is.null(pd)) return(NULL)
        data.frame(chr_str = pd$chr,
                   pos_int = c(pd$inc_pos, pd$dec_pos),
                   stringsAsFactors = FALSE)
      }
      all_tgt <- do.call(rbind, list(mk_tgt(pos_inst), mk_tgt(pos_siii)))
      all_tgt <- all_tgt[!is.na(all_tgt$pos_int), ]
      lkp     <- build_spliceai_pos_lookup_s4(all_tgt, SPLICEAI_BY_CHROM)

      mk_wide <- function(pd, grp) {
        if (is.null(pd) || nrow(pd) == 0L) return(NULL)
        data.frame(
          group    = grp,
          event_id = pd$event_id,
          inc_prob = lkp[paste(pd$chr, pd$inc_pos, sep = ":")],
          dec_prob = lkp[paste(pd$chr, pd$dec_pos, sep = ":")],
          stringsAsFactors = FALSE
        )
      }
      wide <- rbind(
        mk_wide(pos_inst, "insT"),
        mk_wide(pos_siii, "Stem III")
      )
      saveRDS(wide, SPLICEAI_CACHE_A)
      wide
    }

    if (!is.null(panel_a_wide) && nrow(panel_a_wide) > 0L) {
      sig_levels <- c("insT", "Stem III")
      sig_levels <- sig_levels[sig_levels %in% panel_a_wide$group]
      panel_a_wide$group <- factor(panel_a_wide$group, levels = sig_levels)

      # Long format: x = donor type (Decreased / Increased), fill = signature group
      plot_long <- rbind(
        data.frame(sig_group  = panel_a_wide$group,
                   event_id   = panel_a_wide$event_id,
                   donor_type = "Decreased 5'SS",
                   prob       = panel_a_wide$dec_prob,
                   stringsAsFactors = FALSE),
        data.frame(sig_group  = panel_a_wide$group,
                   event_id   = panel_a_wide$event_id,
                   donor_type = "Increased 5'SS",
                   prob       = panel_a_wide$inc_prob,
                   stringsAsFactors = FALSE)
      )
      plot_long <- plot_long[!is.na(plot_long$prob), ]
      plot_long$donor_type <- factor(plot_long$donor_type,
                                     levels = c("Decreased 5'SS", "Increased 5'SS"))
      plot_long$sig_group  <- factor(plot_long$sig_group, levels = sig_levels)

      pal_pa <- c("insT" = PAL_INC, "Stem III" = PAL_STEMIII)

      fig_s4_panelA <- ggplot(plot_long,
                               aes(x = donor_type, y = prob, fill = sig_group)) +
        geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35, alpha = 0.7,
                     position = position_dodge(0.7)) +
        geom_jitter(aes(colour = sig_group),
                    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
                    size = 1.0, alpha = 0.7, show.legend = FALSE) +
        scale_fill_manual(values = pal_pa, name = NULL) +
        scale_colour_manual(values = pal_pa, guide = "none") +
        scale_y_continuous(limits = c(0, 1.08),
                           breaks = c(0, 0.25, 0.5, 0.75, 1),
                           name = "Raw SpliceAI score") +
        labs(x = NULL) +
        pub_theme +
        theme(legend.position   = "top",
              legend.key.size   = unit(3, "mm"),
              legend.text       = element_text(size = 6.5),
              legend.margin     = margin(0, 0, 0, 0),
              legend.box.margin = margin(-4, 0, -4, 0),
              axis.text.x       = element_text(size = 7))

      n_inst <- sum(panel_a_wide$group == "insT")
      n_siii <- sum(panel_a_wide$group == "Stem III")
      cat(sprintf("Panel A: insT n=%d, Stem III n=%d\n", n_inst, n_siii))
    }
  }, error = function(e) {
    cat("Panel A error:", conditionMessage(e), "\n")
  })
}

# ============================================================================
# ---- Panel B: 5'SS sequence logos ------------------------------------------
# Sequence logos for insT and Stem III signatures (stacked 2 rows).
# Moved from figure2.R (old Panel F).
# ============================================================================
fig_s4_panelB <- patchwork::plot_spacer()

if (!is.null(a5ss_lk_s4) && file.exists(FASTA_PATH) &&
    file.exists(paste0(FASTA_PATH, ".fai"))) {
  seqs_insT_s4 <- tryCatch(
    build_or_load_cache_s4(
      file.path(CACHE_DIR_S5, "donor_seqs_insT.rds"),
      sources  = c(LOO_PATH, A5SS_WB_FILE),
      build_fn = function() compute_donor_seqs_s4(LOO_PATH, "insT", a5ss_lk_s4),
      label    = "panel B insT seqs"),
    error = function(e) { cat("Panel B insT:", conditionMessage(e), "\n"); NULL })

  seqs_sIII_s4 <- tryCatch(
    build_or_load_cache_s4(
      file.path(CACHE_DIR_S5, "donor_seqs_stemIII.rds"),
      sources  = c(STEMIII_LOO_PATH, A5SS_WB_FILE),
      build_fn = function() compute_donor_seqs_s4(STEMIII_LOO_PATH, "stemIII", a5ss_lk_s4),
      label    = "panel B stemIII seqs"),
    error = function(e) { cat("Panel B stemIII:", conditionMessage(e), "\n"); NULL })

  logo_insT_s4    <- make_logo_panel_s4(seqs_insT_s4,  "insT signature")
  logo_stemIII_s4 <- make_logo_panel_s4(seqs_sIII_s4,  "Stem III signature")
  fig_s4_panelB   <- patchwork::wrap_elements(full = logo_insT_s4 / logo_stemIII_s4)
} else {
  cat("Panel B: FASTA or A5SS rMATS file missing -- panel shows spacer.\n")
  seqs_insT_s4  <- NULL
  seqs_sIII_s4  <- NULL
}

# ============================================================================
# ---- Combine into Figure S5 PDF --------------------------------------------
# Two-panel layout: A (SpliceAI boxplots) | B (5'SS logos), 50/50.
# ============================================================================
mm_to_in <- function(mm) mm / 25.4

fig_s4_combined <- (fig_s4_panelA + labs(tag = "A")) |
                   (fig_s4_panelB + labs(tag = "B"))
fig_s4_combined <- fig_s4_combined + plot_layout(widths = c(0.35, 0.65))

pdf(OUT_FIG_S7, width = mm_to_in(170), height = mm_to_in(96))
print(fig_s4_combined & theme(plot.tag = element_text(size = 14, face = "bold"),
                              plot.tag.position = "topleft"))
dev.off()
cat(sprintf("Figure S7 saved to: %s\n", OUT_FIG_S7))

# ---- Summary RDS (airlock-safe) --------------------------------------------
clean_inc_insT  <- if (!is.null(seqs_insT_s4))
  seqs_insT_s4$inc_seq_9mer[is_clean_9mer_s4(seqs_insT_s4$inc_seq_9mer)] else character(0L)
clean_dec_insT  <- if (!is.null(seqs_insT_s4))
  seqs_insT_s4$dec_seq_9mer[is_clean_9mer_s4(seqs_insT_s4$dec_seq_9mer)] else character(0L)
clean_inc_sIII  <- if (!is.null(seqs_sIII_s4))
  seqs_sIII_s4$inc_seq_9mer[is_clean_9mer_s4(seqs_sIII_s4$inc_seq_9mer)] else character(0L)
clean_dec_sIII  <- if (!is.null(seqs_sIII_s4))
  seqs_sIII_s4$dec_seq_9mer[is_clean_9mer_s4(seqs_sIII_s4$dec_seq_9mer)] else character(0L)

summary_out <- list(
  # Panel A: SpliceAI counts per group
  panel_a_n_insT_only = if (!is.null(panel_a_wide))
    sum(!is.na(panel_a_wide$inc_prob) & !is.na(panel_a_wide$dec_prob) &
        panel_a_wide$group == "insT") else 0L,
  panel_a_n_stemIII   = if (!is.null(panel_a_wide))
    sum(!is.na(panel_a_wide$inc_prob) & !is.na(panel_a_wide$dec_prob) &
        panel_a_wide$group == "Stem III") else 0L,
  # Panel B: logo 9-mer counts
  panel_b_n_clean_inc_insT    = length(clean_inc_insT),
  panel_b_n_clean_dec_insT    = length(clean_dec_insT),
  panel_b_n_clean_inc_stemIII = length(clean_inc_sIII),
  panel_b_n_clean_dec_stemIII = length(clean_dec_sIII))
saveRDS(summary_out, OUT_SUMMARY)
cat(sprintf("Summary saved to: %s\n", OUT_SUMMARY))

cat("\n=== Done. Figure S5 ===\n")
cat(sprintf("Panel A: insT n=%d, Stem III n=%d events with SpliceAI scores\n",
            summary_out$panel_a_n_insT_only, summary_out$panel_a_n_stemIII))
cat(sprintf("Panel B logos: insT inc=%d dec=%d; stemIII inc=%d dec=%d clean 9-mers\n",
            summary_out$panel_b_n_clean_inc_insT, summary_out$panel_b_n_clean_dec_insT,
            summary_out$panel_b_n_clean_inc_stemIII, summary_out$panel_b_n_clean_dec_stemIII))
cat("\nNOTE: GEL-derived figure stays in workspace unless airlocked.\n")
