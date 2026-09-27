#!/usr/bin/env Rscript
# =============================================================================
# export_a5ss_candidate_full_table.R
# Build one master table with every A5SS event significant in both cohorts
# (FDR<0.05 in both GEL whole-blood and Nava lymphocyte rMATS -- the ~105-
# event shared pool). NOT restricted to direction-concordant events -- a
# `concordant` column flags agreement so discordant events stay visible
# rather than being silently dropped. NOT the 483-event candidate pool that
# feeds the LASSO models (that pool additionally applies a depth filter and
# is built by combine_rmats_batches_cand_pool.R). One row per event.
#
# Simplified column set: identity (event_id, gene, chr, strand, EventType),
# competing 5'SS coordinates + donor distance, direction/concordance,
# GTF/CDS annotation (towards/away transcript come from match_full_event(),
# the adjacency-checked method -- see below), NDD panel flag, insT/stemIII
# signature membership, and per-cohort (whole-blood / lymphocyte) carrier
# counts, delta logit-PSI, FDR, P-value for T-loop and Stem III carriers.
# Blood-side Stem III stats (n_wb_StemIII, wb_delta_logitPSI_StemIII) are
# computed here for the first time -- GEL does have 2 Stem III carriers
# among the 11 ReNU cases (see project_rnu42_region_bounds / stemIII-LASSO
# training counts), this just wasn't wired up previously.
#
# Two transcript-matching strategies feed the GTF/CDS annotation columns:
#   - match_tx(): per-coordinate only. Tests each of the 6 rMATS coordinates
#     independently against any exon boundary of any transcript for the
#     gene -- doesn't require the flanking exon and the alt exon to belong
#     to the same transcript. Feeds long_donor_tx/short_donor_tx (reported
#     for comparison only -- no longer drives CDS/protein-consequence calls).
#   - match_full_event(): joint check. Requires the flanking exon and the
#     alt exon (long- or short-form) to be the SAME transcript's consecutive
#     exons (adjacent exon_number) -- i.e. that transcript actually
#     reproduces the whole flanking-intron-altExon splicing pattern. Feeds
#     long_form_best_transcript/short_form_best_transcript (+ full
#     transcript lists/counts), which are what towards_tx/away_tx/
#     classify_cds/protein_consequence are computed from.
#
# "T-loop carriers" / "Stem III carriers" = all case samples whose HGVS
# falls in nt 63-70 / 75-79 respectively (T-loop includes insT plus other
# T-loop variants), NOT just the 6-event insT-only or stemIII-only LASSO
# models' training subsets. See project_rnu42_region_bounds memory.
#
# Inputs:
#   outputs/psi_matrix.rds             (GEL blood PSI/IJC/SJC + events + samples)
#   outputs/nava_psi_matrix.rds        (Nava lymphocyte, same schema + PValue)
#   data/rmats/A5SS.MATS.JC.txt        (GEL PValue + ENSG + strand lookup)
#   outputs/loo_results_insT_lasso.rds     (insT/"T-loop" signature membership)
#   outputs/loo_results_stemIII_lasso.rds  (Stem III signature membership)
#   data/sample_genotype.csv           (GEL case HGVS)
#   outputs/nava_severity_labels.csv   (Nava case HGVS)
#   rmats_inputs/gencode.v49.annotation.gtf (GTF_PATH env)
#   data/panelapp/panelappAus_ID_green_ensg_grch38.txt
#
# Output: outputs/a5ss_candidate_full_table.csv (one row per candidate event)
#
# GTF caches (built on first run):
#   outputs/gtf_cache_tx_meta.rds, outputs/gtf_cache_cds.rds, outputs/gtf_cache_exons.rds
#
# DATA PROTECTION: event_ids, gene symbols, and per-event GEL case-derived
# stats (dPSI, FDR, P-values, PSI in T-loop carriers) are GEL-protected
# patient-derived data. RE-internal only -- do not airlock without review.
#
# Run: pixi run Rscript scripts/gel/export_a5ss_candidate_full_table.R
#   or Rscript scripts/gel/export_a5ss_candidate_full_table.R
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

# ---- 1. Config --------------------------------------------------------------
GEL_PSI_PATH     <- Sys.getenv("GEL_PSI",         "outputs/psi_matrix.rds")
NAVA_PSI_PATH    <- Sys.getenv("NAVA_PSI",        "outputs/nava_psi_matrix.rds")
RMATS_DIR        <- Sys.getenv("GEL_RMATS_DIR",   "data/rmats")
GTF_PATH         <- Sys.getenv("GTF_PATH",        "rmats_inputs/gencode.v49.annotation.gtf")
PANEL_PATH       <- Sys.getenv("PANEL_NDD",       "data/panelapp/panelappAus_ID_green_ensg_grch38.txt")
INST_LOO_PATH    <- Sys.getenv("INST_LOO_RDS",    "outputs/loo_results_insT_lasso.rds")
STEMIII_LOO_PATH <- Sys.getenv("STEMIII_LOO_RDS", "outputs/loo_results_stemIII_lasso.rds")
GENO_FILE        <- Sys.getenv("GENO_CSV",        "data/sample_genotype.csv")
NAVA_SEV_FILE    <- Sys.getenv("NAVA_SEV",        "outputs/nava_severity_labels.csv")

META_CACHE <- "outputs/gtf_cache_tx_meta.rds"
CDS_CACHE  <- "outputs/gtf_cache_cds.rds"
EXON_CACHE <- "outputs/gtf_cache_exons.rds"

OUT_CSV <- Sys.getenv("OUT_CSV", "outputs/a5ss_candidate_full_table.csv")

EVENT_TYPE  <- "A5SS"
FDR_THRESH  <- 0.05
RNU4_TLOOP  <- c(63L, 70L)
RNU4_STEM3  <- c(75L, 79L)

stopifnot(
  "GEL PSI matrix not found"    = file.exists(GEL_PSI_PATH),
  "Nava PSI matrix not found"   = file.exists(NAVA_PSI_PATH),
  "rMATS A5SS file not found"   = file.exists(file.path(RMATS_DIR, "A5SS.MATS.JC.txt")),
  "GTF not found"               = file.exists(GTF_PATH),
  "insT LOO RDS not found"      = file.exists(INST_LOO_PATH),
  "stemIII LOO RDS not found"   = file.exists(STEMIII_LOO_PATH)
)
dir.create("outputs", showWarnings = FALSE)

# =============================================================================
# 2. Reconstruct the full shared (sig-in-both) A5SS pool -- concordant AND
#    discordant events kept, flagged via `concordant`.
# =============================================================================
cat("Loading PSI matrix event metadata...\n")
gel  <- readRDS(GEL_PSI_PATH)
nava <- readRDS(NAVA_PSI_PATH)

stopifnot(
  all(c("event_id","EventType","coord_nostrand","FDR","IncLevelDiff","GeneSymbol",
        "chr","strand") %in% names(gel$events)),
  all(c("event_id","EventType","coord_nostrand","FDR","IncLevelDiff","PValue") %in%
        names(nava$events))
)

gel_sig  <- gel$events[ !is.na(gel$events$FDR)  & gel$events$FDR  < FDR_THRESH &
                         gel$events$EventType == EVENT_TYPE, ]
nava_sig <- nava$events[!is.na(nava$events$FDR) & nava$events$FDR < FDR_THRESH &
                         nava$events$EventType == EVENT_TYPE, ]

both_keys <- intersect(gel_sig$coord_nostrand, nava_sig$coord_nostrand)
gel_m  <- gel_sig[ match(both_keys, gel_sig$coord_nostrand),  ]
nava_m <- nava_sig[match(both_keys, nava_sig$coord_nostrand), ]

# Concordance uses bb-logit dlogit-PSI (all cases vs all controls, per
# cohort) -- matches figure1.R panel D / figure_s3.R and the LOO scripts.
# Site direction is sign(dlogit_GEL), on the same metric as the LASSO
# feature -- not rMATS's raw IncLevelDifference.
gel_case_all  <- gel$samples$sample_id[ gel$samples$group  == "Case"]
gel_ctrl_all  <- gel$samples$sample_id[ gel$samples$group  == "Control"]
nava_case_all <- nava$samples$sample_id[nava$samples$group == "Case"]
nava_ctrl_all <- nava$samples$sample_id[nava$samples$group == "Control"]

bb_logit_ecc <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}
gel_lp_sig  <- bb_logit_ecc(gel$ijc[,  gel_m$event_id,  drop = FALSE],
                             gel$sjc[,  gel_m$event_id,  drop = FALSE])
nava_lp_sig <- bb_logit_ecc(nava$ijc[, nava_m$event_id, drop = FALSE],
                             nava$sjc[, nava_m$event_id, drop = FALSE])
dl_gel_sig  <- colMeans(gel_lp_sig[ gel_case_all,  , drop = FALSE], na.rm = TRUE) -
               colMeans(gel_lp_sig[ gel_ctrl_all,  , drop = FALSE], na.rm = TRUE)
dl_nava_sig <- colMeans(nava_lp_sig[nava_case_all, , drop = FALSE], na.rm = TRUE) -
               colMeans(nava_lp_sig[nava_ctrl_all, , drop = FALSE], na.rm = TRUE)
dl_gel_sig[ !is.finite(dl_gel_sig)]  <- NA_real_
dl_nava_sig[!is.finite(dl_nava_sig)] <- NA_real_

concordant <- !is.na(dl_gel_sig) & !is.na(dl_nava_sig) &
              sign(dl_gel_sig) == sign(dl_nava_sig) &
              sign(dl_gel_sig) != 0

sig <- data.frame(
  event_id       = gel_m$event_id,
  nava_event_id  = nava_m$event_id,
  gene_symbol    = gel_m$GeneSymbol,
  chr            = gel_m$chr,
  strand         = gel_m$strand,
  coord_nostrand = gel_m$coord_nostrand,
  direction      = sign(unname(dl_gel_sig)),
  concordant     = concordant,
  dPSI_GEL       = gel_m$IncLevelDiff,         # raw rMATS (kept for reference)
  dPSI_Nava      = nava_m$IncLevelDiff,
  dlogit_GEL     = unname(dl_gel_sig),
  dlogit_Nava    = unname(dl_nava_sig),
  FDR_GEL        = gel_m$FDR,
  FDR_Nava       = nava_m$FDR,
  PValue_Nava    = nava_m$PValue,
  stringsAsFactors = FALSE
)
n_cand <- nrow(sig)
stopifnot(
  "No shared A5SS events found" = n_cand > 0L,
  !anyDuplicated(sig$event_id),
  # nava_event_id need not be globally unique here (Nava GeneSymbol
  # annotation can differ from GEL's for the same coord_nostrand), but the
  # (event, nava_event_id) pairing itself must be 1:1 for matrix indexing.
  !anyDuplicated(sig$nava_event_id),
  "direction has zeros" = !any(sig$direction == 0)
)
cat(sprintf("Shared (sig-in-both) A5SS pool: %d events (%d concordant, %d discordant)\n",
            n_cand, sum(sig$concordant), sum(!sig$concordant)))

# ---- 2b. Parse the 6 rMATS coordinates from coord_nostrand ------------------
# coord_nostrand format: chr:longES_0:longEE_1:shortES_0:shortEE_1:flankES_0:flankEE_1
ck_parts <- strsplit(sig$coord_nostrand, ":", fixed = TRUE)
stopifnot(all(lengths(ck_parts) == 7L))
sig$longExonStart_0base <- as.integer(vapply(ck_parts, `[`, character(1L), 2L))
sig$longExonEnd_1based  <- as.integer(vapply(ck_parts, `[`, character(1L), 3L))
sig$shortES_0base       <- as.integer(vapply(ck_parts, `[`, character(1L), 4L))
sig$shortEE_1based      <- as.integer(vapply(ck_parts, `[`, character(1L), 5L))
sig$flankingES_0base    <- as.integer(vapply(ck_parts, `[`, character(1L), 6L))
sig$flankingEE_1based   <- as.integer(vapply(ck_parts, `[`, character(1L), 7L))

# Competing 5'SS genomic positions (1-based). Donor distance is NOT a
# separate column -- it's mathematically identical to alt_len_nt (section 3
# below): the alt region IS the sequence between the two competing donor
# sites, so its length always equals their distance. alt_len_nt is the one
# carried through to protein_delta/cds_class/the figures.
sig$long_5ss_pos  <- ifelse(sig$strand == "+", sig$longExonEnd_1based, sig$longExonStart_0base + 1L)
sig$short_5ss_pos <- ifelse(sig$strand == "+", sig$shortEE_1based,     sig$shortES_0base + 1L)

# ---- 2c. ENSG + PValue (GEL) from the raw rMATS A5SS file -------------------
cat("Parsing rMATS A5SS file for ENSG + GEL PValue...\n")
a5 <- read.delim(file.path(RMATS_DIR, "A5SS.MATS.JC.txt"), stringsAsFactors = FALSE, check.names = FALSE)
a5_ensg      <- sub("\\..*$", "", gsub('"', "", a5[[2L]]))
# Full 6-coordinate key (chr + long/short exon + flanking exon) -- a 4-
# coordinate key (dropping flanking exon) can collide across distinct rMATS
# rows that reuse the same long/short exon boundary pair with a different
# flanking exon, silently pulling the wrong row's PValue via match().
a5_joinkey_ns <- paste(a5[[4L]], as.integer(a5[[6L]]), as.integer(a5[[7L]]),
                       as.integer(a5[[8L]]), as.integer(a5[[9L]]),
                       as.integer(a5[[10L]]), as.integer(a5[[11L]]), sep = ":")
sig_joinkey <- paste(sig$chr, sig$longExonStart_0base, sig$longExonEnd_1based,
                     sig$shortES_0base, sig$shortEE_1based,
                     sig$flankingES_0base, sig$flankingEE_1based, sep = ":")
hit_idx <- match(sig_joinkey, a5_joinkey_ns)
sig$ensg        <- a5_ensg[hit_idx]
sig$PValue_GEL  <- suppressWarnings(as.numeric(a5$PValue))[hit_idx]
cat(sprintf("  Matched ENSG for %d / %d candidates; matched GEL PValue for %d / %d\n",
            sum(!is.na(sig$ensg)), n_cand, sum(!is.na(sig$PValue_GEL)), n_cand))

# =============================================================================
# 3. NDD panel + GTF-derived transcript/CDS annotation
#    (operating on `sig`)
# =============================================================================
ndd_ensg <- if (file.exists(PANEL_PATH)) {
  v <- unique(trimws(readLines(PANEL_PATH))); v[nzchar(v)]
} else {
  warning("NDD panel file not found: ", PANEL_PATH, " -- in_ndd_panel will be NA")
  character(0L)
}

priority_label <- c(
  "0" = "MANE_Select",
  "1" = "appris_principal_1", "2" = "appris_principal_2",
  "3" = "appris_principal_3", "4" = "appris_principal_4", "5" = "appris_principal_5",
  "6" = "appris_alternative_1", "7" = "appris_alternative_2",
  "8" = "Ensembl_canonical", "9" = "basic", "99" = "none")

rank_of <- function(tags) {
  if ("MANE_Select" %in% tags) return(0L)
  pr <- intersect(paste0("appris_principal_", 1:5), tags)
  if (length(pr)) return(as.integer(sub(".*_", "", pr[1L])))
  if ("appris_alternative_1" %in% tags) return(6L)
  if ("appris_alternative_2" %in% tags) return(7L)
  if ("Ensembl_canonical"    %in% tags) return(8L)
  if ("basic"                %in% tags) return(9L)
  99L
}

if (!file.exists(META_CACHE) || file.mtime(GTF_PATH) > file.mtime(META_CACHE)) {
  message("Building transcript metadata cache from ", GTF_PATH, " ...")
  tx_raw <- read.delim(
    pipe(sprintf("awk -F'\\t' '$3==\"transcript\"' %s", shQuote(GTF_PATH))),
    header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE,
    col.names = c("seqnames","source","type","start","end","score","strand","frame","attr"))
  attr_str <- tx_raw$attr
  tx_id    <- sub('^.*transcript_id "([^"]+)".*$', '\\1', attr_str)
  biotype  <- sub('^.*transcript_type "([^"]+)".*$', '\\1', attr_str)
  biotype[!grepl('transcript_type "', attr_str, fixed = TRUE)] <- NA_character_
  ccdsid   <- sub('^.*ccdsid "([^"]+)".*$', '\\1', attr_str)
  ccdsid[!grepl('ccdsid "', attr_str, fixed = TRUE)] <- NA_character_
  tags_list <- regmatches(attr_str, gregexpr('(?<=tag ")[^"]+', attr_str, perl = TRUE))
  tx_priority <- data.frame(
    transcript_id = tx_id, priority = vapply(tags_list, rank_of, integer(1L)),
    transcript_type = biotype, ccdsid = ccdsid, stringsAsFactors = FALSE)
  saveRDS(tx_priority, META_CACHE)
  message(sprintf("  Cached %d transcript records -> %s", nrow(tx_priority), META_CACHE))
} else {
  tx_priority <- readRDS(META_CACHE)
  message("Loaded transcript metadata from cache.")
}

if (!file.exists(CDS_CACHE) || file.mtime(GTF_PATH) > file.mtime(CDS_CACHE)) {
  message("Building CDS cache from ", GTF_PATH, " ...")
  cds_raw <- read.delim(
    pipe(sprintf("awk -F'\\t' '$3==\"CDS\"' %s", shQuote(GTF_PATH))),
    header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE,
    col.names = c("seqnames","source","type","start","end","score","strand","frame","attr"))
  cds_df <- data.frame(
    transcript_id = sub('^.*transcript_id "([^"]+)".*$', '\\1', cds_raw$attr),
    seqnames = cds_raw$seqnames, start = as.integer(cds_raw$start), end = as.integer(cds_raw$end),
    stringsAsFactors = FALSE)
  saveRDS(cds_df, CDS_CACHE)
  message(sprintf("  Cached %d CDS segments -> %s", nrow(cds_df), CDS_CACHE))
} else {
  cds_df <- readRDS(CDS_CACHE)
  message("Loaded CDS cache.")
}

# exon_number retained for match_full_event() below (the flanking-exon/
# alt-exon adjacency check) -- GENCODE GTF exon_number always runs 5'->3'
# along the transcript, so exons N and N+1 are guaranteed adjacent
# regardless of strand. NB exon_number is an unquoted GTF attribute
# ('exon_number 1;'), unlike the quoted string attributes above. Cache is
# GENCODE GTF exon_number always runs 5-prime to 3-prime along the transcript.
build_exon_cache <- function() {
  message("Building exon boundary cache from ", GTF_PATH, " ...")
  exon_raw <- read.delim(
    pipe(sprintf("awk -F'\\t' '$3==\"exon\"' %s", shQuote(GTF_PATH))),
    header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE,
    col.names = c("seqnames","source","type","start","end","score","strand","frame","attr"))
  exon_attr <- exon_raw$attr
  exons_df <- data.frame(
    seqnames = exon_raw$seqnames, start = as.integer(exon_raw$start), end = as.integer(exon_raw$end),
    strand = exon_raw$strand, transcript_id = sub('^.*transcript_id "([^"]+)".*$', '\\1', exon_attr),
    gene_id_base = sub("\\.\\d+$", "", sub('^.*gene_id "([^"]+)".*$', '\\1', exon_attr)),
    exon_number = as.integer(sub('^.*exon_number ([0-9]+);.*$', '\\1', exon_attr)),
    stringsAsFactors = FALSE)
  saveRDS(exons_df, EXON_CACHE)
  message(sprintf("  Cached %d exon records -> %s", nrow(exons_df), EXON_CACHE))
  exons_df
}

if (!file.exists(EXON_CACHE) || file.mtime(GTF_PATH) > file.mtime(EXON_CACHE)) {
  exons_df <- build_exon_cache()
} else {
  exons_df <- readRDS(EXON_CACHE)
  if (!"exon_number" %in% names(exons_df)) {
    message("Exon cache predates exon_number -- rebuilding.")
    exons_df <- build_exon_cache()
  } else {
    message("Loaded exon cache.")
  }
}

target_ensgs <- unique(sig$ensg[!is.na(sig$ensg)])
exons_scoped <- exons_df[exons_df$gene_id_base %in% target_ensgs, ]
cat(sprintf("Exons scoped to %d target genes: %d rows\n", length(target_ensgs), nrow(exons_scoped)))

coord_order <- c("longExonStart_0base","longExonEnd_1based","shortES_0base","shortEE_1based",
                 "flankingES_0base","flankingEE_1based")
coord_boundary <- c(
  longExonStart_0base = "start", longExonEnd_1based = "end",
  shortES_0base = "start", shortEE_1based = "end",
  flankingES_0base = "start", flankingEE_1based = "end")

long_list <- lapply(coord_order, function(cn) {
  d <- sig[, c("gene_symbol","ensg","event_id","chr","strand","direction")]
  d$coord_name <- cn
  d$rmats_pos  <- sig[[cn]]
  d$boundary   <- coord_boundary[[cn]]
  d$gtf_pos    <- if (coord_boundary[[cn]] == "start") d$rmats_pos + 1L else d$rmats_pos
  d
})
long <- do.call(rbind, long_list)
rownames(long) <- NULL

match_tx <- function(chr, strand_, ensg, boundary, gtf_pos) {
  if (is.na(ensg))
    return(list(transcripts = NA_character_, best_transcript = NA_character_, best_tier = NA_character_))
  rows <- exons_scoped[
    exons_scoped$gene_id_base == ensg & exons_scoped$seqnames == chr & exons_scoped$strand == strand_ &
    (if (boundary == "start") exons_scoped$start == gtf_pos else exons_scoped$end == gtf_pos), ]
  if (nrow(rows) == 0L) {
    list(transcripts = NA_character_, best_transcript = NA_character_, best_tier = NA_character_)
  } else {
    tx <- sort(unique(rows$transcript_id))
    pr <- tx_priority$priority[match(tx, tx_priority$transcript_id)]
    pr[is.na(pr)] <- 99L
    ord <- order(pr, tx)
    list(transcripts = paste(tx[ord], collapse = ";"), best_transcript = tx[ord[1L]],
         best_tier = unname(priority_label[as.character(pr[ord[1L]])]))
  }
}

cat(sprintf("Matching transcript boundaries (%d coordinate x event pairs)...\n", nrow(long)))
m_list <- vector("list", nrow(long))
for (i in seq_len(nrow(long)))
  m_list[[i]] <- match_tx(long$chr[i], long$strand[i], long$ensg[i], long$boundary[i], long$gtf_pos[i])
long$transcripts     <- vapply(m_list, `[[`, character(1L), "transcripts")
long$best_transcript <- vapply(m_list, `[[`, character(1L), "best_transcript")
long$best_tier       <- vapply(m_list, `[[`, character(1L), "best_tier")

# Stronger check than match_tx():
# match_tx() tests whether a single coordinate (e.g. the alt-exon donor)
# exists as *some* exon boundary in *some* transcript of the gene --
# independently for each of the 6 rMATS coordinates, with no requirement
# they belong to the same transcript. match_full_event() instead requires
# that the flanking exon and the alt exon (long- or short-form) are the
# SAME transcript's consecutive exons (exon_number differs by exactly 1) --
# i.e. that transcript actually reproduces the full flanking-intron-altExon
# splicing pattern being tested, not just an isolated boundary coincidence.
match_full_event <- function(chr, strand_, ensg, flank_start, flank_end, alt_start, alt_end) {
  empty <- list(transcripts = NA_character_, best_transcript = NA_character_, best_tier = NA_character_)
  if (is.na(ensg)) return(empty)
  tx_exons <- exons_scoped[
    exons_scoped$gene_id_base == ensg & exons_scoped$seqnames == chr & exons_scoped$strand == strand_, ]
  if (nrow(tx_exons) == 0) return(empty)
  flank_rows <- tx_exons[tx_exons$start == flank_start & tx_exons$end == flank_end, ]
  if (nrow(flank_rows) == 0) return(empty)
  matched_tx <- character(0)
  for (i in seq_len(nrow(flank_rows))) {
    tid   <- flank_rows$transcript_id[i]
    f_num <- flank_rows$exon_number[i]
    # either neighbour -- flanking exon can sit on either side of the alt
    # exon depending on strand, so don't hardcode which
    nbrs <- tx_exons[tx_exons$transcript_id == tid &
                     tx_exons$exon_number %in% c(f_num - 1L, f_num + 1L), ]
    if (any(nbrs$start == alt_start & nbrs$end == alt_end)) matched_tx <- c(matched_tx, tid)
  }
  matched_tx <- sort(unique(matched_tx))
  if (length(matched_tx) == 0) return(empty)
  pr  <- tx_priority$priority[match(matched_tx, tx_priority$transcript_id)]
  pr[is.na(pr)] <- 99L
  ord <- order(pr, matched_tx)
  list(transcripts     = paste(matched_tx[ord], collapse = ";"),
       best_transcript = matched_tx[ord[1L]],
       best_tier       = unname(priority_label[as.character(pr[ord[1L]])]))
}

cat("Matching full events (flanking + alt exon adjacency, long and short form)...\n")
flank_start_v <- sig$flankingES_0base + 1L
flank_end_v   <- sig$flankingEE_1based
long_start_v  <- sig$longExonStart_0base + 1L
long_end_v    <- sig$longExonEnd_1based
short_start_v <- sig$shortES_0base + 1L
short_end_v   <- sig$shortEE_1based

long_full_list  <- vector("list", n_cand)
short_full_list <- vector("list", n_cand)
for (i in seq_len(n_cand)) {
  long_full_list[[i]]  <- match_full_event(sig$chr[i], sig$strand[i], sig$ensg[i],
                                           flank_start_v[i], flank_end_v[i], long_start_v[i], long_end_v[i])
  short_full_list[[i]] <- match_full_event(sig$chr[i], sig$strand[i], sig$ensg[i],
                                           flank_start_v[i], flank_end_v[i], short_start_v[i], short_end_v[i])
}
long_form_transcripts      <- vapply(long_full_list,  `[[`, character(1L), "transcripts")
long_form_best_transcript  <- vapply(long_full_list,  `[[`, character(1L), "best_transcript")
long_form_best_tier        <- vapply(long_full_list,  `[[`, character(1L), "best_tier")
n_long_form_transcripts    <- ifelse(is.na(long_form_transcripts), 0L,
                                     lengths(strsplit(long_form_transcripts, ";", fixed = TRUE)))
short_form_transcripts     <- vapply(short_full_list, `[[`, character(1L), "transcripts")
short_form_best_transcript <- vapply(short_full_list, `[[`, character(1L), "best_transcript")
short_form_best_tier       <- vapply(short_full_list, `[[`, character(1L), "best_tier")
n_short_form_transcripts   <- ifelse(is.na(short_form_transcripts), 0L,
                                     lengths(strsplit(short_form_transcripts, ";", fixed = TRUE)))
cat(sprintf(
  "Full-event match: long form in >=1 transcript for %d/%d events; short form for %d/%d events\n",
  sum(n_long_form_transcripts  > 0L), n_cand, sum(n_short_form_transcripts > 0L), n_cand))

classify_cds <- function(tx, chr, alt_start, alt_end, strand) {
  alt_len <- alt_end - alt_start + 1L
  if (is.na(tx)) return(list(class = "unannotated", cds_nt = NA_integer_))
  rows <- cds_df[cds_df$transcript_id == tx & cds_df$seqnames == chr, ]
  if (nrow(rows) == 0L) return(list(class = "non_coding", cds_nt = 0L))
  ov_nt <- as.integer(sum(pmax(0L, pmin(rows$end, alt_end) - pmax(rows$start, alt_start) + 1L)))
  cls <- if (ov_nt > 0L && ov_nt == alt_len) {
    "CDS"
  } else if (ov_nt > 0L) {
    "partial_CDS"
  } else {
    cds_min <- min(rows$start); cds_max <- max(rows$end)
    is_five_prime <- if (strand == "+") alt_end < cds_min else alt_start > cds_max
    if (is_five_prime) "5UTR" else "3UTR"
  }
  list(class = cls, cds_nt = ov_nt)
}

frame_text <- function(alt_len, direction, cds_class, cds_nt) {
  sign_ <- if (direction > 0) "+" else "-"
  switch(cds_class,
    "CDS" = if (alt_len %% 3L == 0L)
              sprintf("in-frame %s%d AA (%s%d nt CDS)", sign_, alt_len %/% 3L, sign_, alt_len)
            else sprintf("frameshift (%s%d nt CDS)", sign_, alt_len),
    "5UTR"        = sprintf("5'UTR-only, no protein change (%s%d nt)", sign_, alt_len),
    "3UTR"        = sprintf("3'UTR-only, no protein change (%s%d nt)", sign_, alt_len),
    "partial_CDS" = sprintf("spans CDS/UTR boundary: %s%d nt CDS + %s%d nt UTR", sign_, cds_nt, sign_, alt_len - cds_nt),
    "non_coding"  = sprintf("non-coding transcript (%s%d nt)", sign_, alt_len),
    "unannotated" = sprintf("unannotated splice junction, protein effect uncertain (%s%d nt)", sign_, alt_len))
}

lu_key  <- paste(long$event_id, long$coord_name, sep = "|||")
lu_tx   <- setNames(long$best_transcript, lu_key)
lu_tier <- setNames(long$best_tier,       lu_key)

long_donor_coord  <- ifelse(sig$strand == "+", "longExonEnd_1based",  "longExonStart_0base")
short_donor_coord <- ifelse(sig$strand == "+", "shortEE_1based",      "shortES_0base")
ld_key <- paste(sig$event_id, long_donor_coord,  sep = "|||")
sd_key <- paste(sig$event_id, short_donor_coord, sep = "|||")

long_donor_tx    <- unname(lu_tx[ld_key])
long_donor_tier  <- unname(lu_tier[ld_key])
short_donor_tx   <- unname(lu_tx[sd_key])
short_donor_tier <- unname(lu_tier[sd_key])

# towards/away transcript basis: the full-event (adjacency-checked) match,
# not the single-boundary long_donor_tx/short_donor_tx -- a transcript
# should only inform CDS/protein-consequence calls if it actually
# reproduces the flanking+alt exon splicing pattern being tested, kept
# alongside the raw single-boundary columns for comparison.
towards_tx   <- ifelse(sig$direction > 0, long_form_best_transcript,  short_form_best_transcript)
towards_tier <- ifelse(sig$direction > 0, long_form_best_tier,        short_form_best_tier)
away_tx      <- ifelse(sig$direction > 0, short_form_best_transcript, long_form_best_transcript)
away_tier    <- ifelse(sig$direction > 0, short_form_best_tier,       long_form_best_tier)

tx_lu_bio  <- setNames(tx_priority$transcript_type, tx_priority$transcript_id)
tx_lu_ccds <- setNames(tx_priority$ccdsid,          tx_priority$transcript_id)
towards_biotype <- unname(tx_lu_bio[towards_tx]); towards_ccds <- unname(tx_lu_ccds[towards_tx])
away_biotype    <- unname(tx_lu_bio[away_tx]);    away_ccds    <- unname(tx_lu_ccds[away_tx])

alt_start  <- ifelse(sig$strand == "+", sig$shortEE_1based + 1L, sig$longExonStart_0base + 1L)
alt_end    <- ifelse(sig$strand == "+", sig$longExonEnd_1based,  sig$shortES_0base)
alt_len_nt <- alt_end - alt_start + 1L

cat("Classifying CDS overlap and protein consequence...\n")
# alt region exists as exonic sequence ONLY in the long-form transcript, so
# we test overlap against long_form_best_transcript -- the full-event
# (flanking + alt exon adjacency-checked) match, not the single-boundary
# long_donor_tx.
cds_results <- vector("list", n_cand)
for (i in seq_len(n_cand))
  cds_results[[i]] <- classify_cds(long_form_best_transcript[i], sig$chr[i], alt_start[i], alt_end[i], sig$strand[i])
cds_class      <- vapply(cds_results, `[[`, character(1L), "class")
cds_overlap_nt <- vapply(cds_results, function(x) { v <- x$cds_nt; if (is.null(v) || is.na(v)) NA_integer_ else as.integer(v) }, integer(1L))

protein_delta <- vapply(seq_len(n_cand), function(i)
  frame_text(alt_len_nt[i], sig$direction[i], cds_class[i], cds_overlap_nt[i]), character(1L))

# Simplified to a binary NMD call -- the only thing figure2.R/figure3.R's
# event schematic actually reads from this column (is_nmd <- grepl("NMD",
# ...)). NMD predicted only when there IS a towards-transcript with a known
# biotype flagged nonsense_mediated_decay / protein_coding_CDS_not_defined;
# UTR-only/non-coding/unannotated events are not NMD calls at all, so they
# fall to "no NMD predicted" same as ordinary protein-coding events.
protein_consequence <- vapply(seq_len(n_cand), function(i) {
  cc <- cds_class[i]; tb <- towards_biotype[i]
  is_nmd <- !is.na(towards_tx[i]) &&
            !cc %in% c("5UTR", "3UTR", "non_coding", "unannotated") &&
            !is.na(tb) && tb %in% c("nonsense_mediated_decay", "protein_coding_CDS_not_defined")
  if (is_nmd) "NMD predicted" else "no NMD predicted"
}, character(1L))

tier_label_for <- function(tier) if (is.na(tier)) "an unannotated splice site" else tier
annotation <- vapply(seq_len(n_cand), function(i) {
  td <- if (is.na(towards_tx[i])) "" else
    paste0(" (", towards_tx[i], ", ", if (is.na(towards_biotype[i])) "biotype unknown" else towards_biotype[i],
           if (!is.na(towards_ccds[i])) paste0(", ", towards_ccds[i]) else "", ")")
  ad <- if (is.na(away_tx[i])) "" else
    paste0(" (", away_tx[i], ", ", if (is.na(away_biotype[i])) "biotype unknown" else away_biotype[i],
           if (!is.na(away_ccds[i])) paste0(", ", away_ccds[i]) else "", ")")
  sprintf("ReNU shifts splicing towards %s%s; away from %s%s. Naive frame call: %s. Refined: %s.",
          tier_label_for(towards_tier[i]), td, tier_label_for(away_tier[i]), ad,
          protein_delta[i], protein_consequence[i])
}, character(1L))

sig$in_ndd_panel          <- if (length(ndd_ensg) > 0L) sig$ensg %in% ndd_ensg else NA
sig$favoured_form_in_renu <- ifelse(sig$direction > 0, "long", "short")
sig$alt_len_nt          <- alt_len_nt
sig$cds_class           <- cds_class
sig$cds_overlap_nt      <- cds_overlap_nt
sig$protein_delta       <- protein_delta
sig$protein_consequence <- protein_consequence
sig$towards_tier        <- towards_tier
sig$towards_tx          <- towards_tx
sig$towards_biotype     <- towards_biotype
sig$towards_ccds        <- towards_ccds
sig$away_tier           <- away_tier
sig$away_tx             <- away_tx
sig$away_biotype        <- away_biotype
sig$away_ccds           <- away_ccds
sig$long_donor_tier     <- long_donor_tier
sig$long_donor_tx       <- long_donor_tx
sig$short_donor_tier    <- short_donor_tier
sig$short_donor_tx      <- short_donor_tx
sig$long_form_best_tier        <- long_form_best_tier
sig$long_form_best_transcript  <- long_form_best_transcript
sig$n_long_form_transcripts    <- n_long_form_transcripts
sig$long_form_transcripts      <- long_form_transcripts
sig$short_form_best_tier       <- short_form_best_tier
sig$short_form_best_transcript <- short_form_best_transcript
sig$n_short_form_transcripts   <- n_short_form_transcripts
sig$short_form_transcripts     <- short_form_transcripts
sig$annotation          <- annotation

# =============================================================================
# 4. insT ("T-loop") and Stem III signature membership
# =============================================================================
insT_loo    <- readRDS(INST_LOO_PATH)
stemIII_loo <- readRDS(STEMIII_LOO_PATH)
stopifnot(all(c("a5ss_cand","selected_idx") %in% names(insT_loo)),
          all(c("a5ss_cand","selected_idx") %in% names(stemIII_loo)))
insT_sel_events    <- insT_loo$a5ss_cand$gel_ev_id[insT_loo$selected_idx]
stemIII_sel_events <- stemIII_loo$a5ss_cand$gel_ev_id[stemIII_loo$selected_idx]

sig$is_insT_Tloop_signature <- sig$event_id %in% insT_sel_events
sig$is_stemIII_signature    <- sig$event_id %in% stemIII_sel_events
cat(sprintf("insT/T-loop signature events found in pool: %d / %d\n",
            sum(sig$is_insT_Tloop_signature), length(insT_sel_events)))
cat(sprintf("Stem III signature events found in pool: %d / %d\n",
            sum(sig$is_stemIII_signature), length(stemIII_sel_events)))

# =============================================================================
# 5. Per-cohort PSI stats: controls vs T-loop carriers vs Stem III carriers
#    (blood + lymphocyte)
# =============================================================================
classify_rnu4_region <- function(hgvs) {
  out <- rep(NA_character_, length(hgvs))
  for (i in seq_along(hgvs)) {
    h <- hgvs[i]
    if (is.na(h) || !nzchar(h)) next
    coord <- regmatches(h, regexpr("n\\.[0-9]+(_[0-9]+)?", h))
    if (length(coord) == 0L) next
    pos <- as.integer(regmatches(coord, gregexpr("[0-9]+", coord))[[1]])
    if (length(pos) == 0L) next
    if (any(pos >= RNU4_TLOOP[1] & pos <= RNU4_TLOOP[2])) out[i] <- "T-loop"
    else if (any(pos >= RNU4_STEM3[1] & pos <= RNU4_STEM3[2])) out[i] <- "Stem III"
  }
  out
}

bb_logit <- function(ijc, sjc) {
  p <- (ijc + 0.5) / (ijc + sjc + 1)
  log(p / (1 - p))
}

hgvs_gel <- setNames(character(0), character(0))
if (file.exists(GENO_FILE)) {
  geno <- read.csv(GENO_FILE, stringsAsFactors = FALSE)
  stopifnot(all(c("sample_id","hgvs") %in% names(geno)))
  hgvs_gel <- setNames(geno$hgvs, geno$sample_id)
} else {
  warning("GENO_FILE not found: ", GENO_FILE, " -- GEL T-loop carrier group will be empty")
}

hgvs_nava <- setNames(character(0), character(0))
if (file.exists(NAVA_SEV_FILE)) {
  nava_lbl <- local({
    first <- readLines(NAVA_SEV_FILE, n = 1L)
    sep   <- if (grepl("\t", first, fixed = TRUE)) "\t" else ","
    read.table(NAVA_SEV_FILE, header = TRUE, sep = sep, stringsAsFactors = FALSE,
               quote = "\"", fill = TRUE)
  })
  nava_lbl$sample_id <- paste0("nava:case:", nava_lbl$sample_n)
  hgvs_nava <- setNames(nava_lbl$hgvs, nava_lbl$sample_id)
} else {
  warning("NAVA_SEV_FILE not found: ", NAVA_SEV_FILE, " -- Nava T-loop carrier group will be empty")
}

gel_ctrl_ids     <- gel$samples$sample_id[gel$samples$group == "Control"]
gel_case_ids     <- gel$samples$sample_id[gel$samples$group == "Case"]
gel_region       <- classify_rnu4_region(hgvs_gel[gel_case_ids])
gel_tloop_ids    <- gel_case_ids[!is.na(gel_region) & gel_region == "T-loop"]
gel_stemIII_ids  <- gel_case_ids[!is.na(gel_region) & gel_region == "Stem III"]

nava_ctrl_ids    <- nava$samples$sample_id[nava$samples$group == "Control"]
nava_case_ids    <- nava$samples$sample_id[nava$samples$group == "Case"]
nava_region      <- classify_rnu4_region(hgvs_nava[nava_case_ids])
nava_tloop_ids   <- nava_case_ids[!is.na(nava_region) & nava_region == "T-loop"]
nava_stemIII_ids <- nava_case_ids[!is.na(nava_region) & nava_region == "Stem III"]

cat(sprintf("GEL (blood): %d controls, %d T-loop carriers, %d Stem III carriers (of %d cases)\n",
            length(gel_ctrl_ids), length(gel_tloop_ids), length(gel_stemIII_ids), length(gel_case_ids)))
cat(sprintf("Nava (lymphocyte): %d controls, %d T-loop carriers, %d Stem III carriers (of %d cases)\n",
            length(nava_ctrl_ids), length(nava_tloop_ids), length(nava_stemIII_ids), length(nava_case_ids)))
stopifnot("No GEL T-loop carriers found"    = length(gel_tloop_ids)    > 0L,
          "No GEL Stem III carriers found"  = length(gel_stemIII_ids)  > 0L,
          "No Nava T-loop carriers found"   = length(nava_tloop_ids)   > 0L,
          "No Nava Stem III carriers found" = length(nava_stemIII_ids) > 0L)

ev      <- sig$event_id
nava_ev <- sig$nava_event_id  # Nava's own event_id space -- differs from
                              # GEL's when GeneSymbol annotation differs for
                              # the same coord_nostrand; same row order as sig.
gel_logit_psi  <- bb_logit(gel$ijc[,  ev,      drop = FALSE], gel$sjc[,  ev,      drop = FALSE])
nava_logit_psi <- bb_logit(nava$ijc[, nava_ev, drop = FALSE], nava$sjc[, nava_ev, drop = FALSE])

col_mean <- function(mat, ids) {
  ids <- intersect(ids, rownames(mat))
  if (length(ids) == 0L) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[ids, , drop = FALSE], na.rm = TRUE)
}

gel_psi_sub <- gel$psi[, ev, drop = FALSE]
sig$n_gel_ctrl              <- length(gel_ctrl_ids)
sig$n_gel_Tloop_carrier     <- length(gel_tloop_ids)
sig$gel_mean_PSI_ctrl       <- col_mean(gel_psi_sub, gel_ctrl_ids)
sig$gel_mean_PSI_Tloop      <- col_mean(gel_psi_sub, gel_tloop_ids)
sig$gel_delta_PSI           <- sig$gel_mean_PSI_Tloop - sig$gel_mean_PSI_ctrl
sig$gel_mean_logitPSI_ctrl  <- col_mean(gel_logit_psi, gel_ctrl_ids)
sig$gel_mean_logitPSI_Tloop <- col_mean(gel_logit_psi, gel_tloop_ids)
sig$gel_delta_logitPSI      <- sig$gel_mean_logitPSI_Tloop - sig$gel_mean_logitPSI_ctrl

sig$n_gel_stemIII_carrier      <- length(gel_stemIII_ids)
sig$gel_mean_PSI_stemIII       <- col_mean(gel_psi_sub, gel_stemIII_ids)
sig$gel_delta_PSI_stemIII      <- sig$gel_mean_PSI_stemIII - sig$gel_mean_PSI_ctrl
sig$gel_mean_logitPSI_stemIII  <- col_mean(gel_logit_psi, gel_stemIII_ids)
sig$gel_delta_logitPSI_stemIII <- sig$gel_mean_logitPSI_stemIII - sig$gel_mean_logitPSI_ctrl

sig$n_nava_ctrl              <- length(nava_ctrl_ids)
sig$n_nava_Tloop_carrier      <- length(nava_tloop_ids)
nava_psi_sub <- nava$psi[, nava_ev, drop = FALSE]
sig$nava_mean_PSI_ctrl        <- col_mean(nava_psi_sub, nava_ctrl_ids)
sig$nava_mean_PSI_Tloop       <- col_mean(nava_psi_sub, nava_tloop_ids)
sig$nava_delta_PSI            <- sig$nava_mean_PSI_Tloop - sig$nava_mean_PSI_ctrl
sig$nava_mean_logitPSI_ctrl   <- col_mean(nava_logit_psi, nava_ctrl_ids)
sig$nava_mean_logitPSI_Tloop  <- col_mean(nava_logit_psi, nava_tloop_ids)
sig$nava_delta_logitPSI       <- sig$nava_mean_logitPSI_Tloop - sig$nava_mean_logitPSI_ctrl

sig$n_nava_stemIII_carrier         <- length(nava_stemIII_ids)
sig$nava_mean_PSI_stemIII          <- col_mean(nava_psi_sub, nava_stemIII_ids)
sig$nava_delta_PSI_stemIII         <- sig$nava_mean_PSI_stemIII - sig$nava_mean_PSI_ctrl
sig$nava_mean_logitPSI_stemIII     <- col_mean(nava_logit_psi, nava_stemIII_ids)
sig$nava_delta_logitPSI_stemIII    <- sig$nava_mean_logitPSI_stemIII - sig$nava_mean_logitPSI_ctrl

# =============================================================================
# 6. Assemble + write final table
# =============================================================================
# Simplified export -- one row per shared (sig-in-both) A5SS event,
# concordant and discordant alike (see `concordant` column). Dropped
# relative to the full internal `sig` table: raw PSI means (kept as
# delta-logit-PSI only), CDS overlap size, per-boundary (non-adjacency-
# checked) donor tx/tier columns, transcript ID lists/counts, and the
# free-text `annotation` sentence -- all still computed above and
# available in `sig` if needed.
# alt_len_nt, protein_consequence, towards_tier, away_tier are kept
# (despite being droppable-looking) because figure2.R/figure3.R's
# make_event_schematic() and CINP lookup read this file directly --
# see project note on retiring annotate_candidate_pool.R in favour of
# this table as the single source for those panels.
out <- data.frame(
  event_id                  = sig$event_id,
  gene_symbol                = sig$gene_symbol,
  ensg                       = sig$ensg,
  chr                        = sig$chr,
  strand                     = sig$strand,
  EventType                  = EVENT_TYPE,
  longExonStart_0base        = sig$longExonStart_0base,
  longExonEnd_1based         = sig$longExonEnd_1based,
  shortES_0base               = sig$shortES_0base,
  shortEE_1based               = sig$shortEE_1based,
  flankingES_0base             = sig$flankingES_0base,
  flankingEE_1based            = sig$flankingEE_1based,
  long_5ss_pos                 = sig$long_5ss_pos,
  short_5ss_pos                 = sig$short_5ss_pos,
  direction                    = sig$direction,
  concordant                   = sig$concordant,
  favoured_form_in_renu        = sig$favoured_form_in_renu,
  in_ndd_panel                  = sig$in_ndd_panel,
  alt_len_nt                    = sig$alt_len_nt,
  cds_class                     = sig$cds_class,
  protein_delta                 = sig$protein_delta,
  protein_consequence           = sig$protein_consequence,
  towards_tier                  = sig$towards_tier,
  towards_transcript            = sig$towards_tx,
  towards_biotype                = sig$towards_biotype,
  away_tier                       = sig$away_tier,
  away_transcript                 = sig$away_tx,
  away_biotype                     = sig$away_biotype,
  is_insT_Tloop_signature           = sig$is_insT_Tloop_signature,
  is_stemIII_signature               = sig$is_stemIII_signature,
  n_wb_ctrl                          = sig$n_gel_ctrl,
  n_wb_Tloop                         = sig$n_gel_Tloop_carrier,
  n_wb_StemIII                        = sig$n_gel_stemIII_carrier,
  wb_delta_logitPSI_Tloop              = sig$gel_delta_logitPSI,
  wb_delta_logitPSI_StemIII             = sig$gel_delta_logitPSI_stemIII,
  wb_rmats_FDR                          = sig$FDR_GEL,
  wb_rmats_PValue                       = sig$PValue_GEL,
  n_lymph_ctrl                          = sig$n_nava_ctrl,
  n_lymph_Tloop                         = sig$n_nava_Tloop_carrier,
  n_lymph_StemIII                       = sig$n_nava_stemIII_carrier,
  lymph_delta_logitPSI_Tloop            = sig$nava_delta_logitPSI,
  lymph_delta_logitPSI_StemIII          = sig$nava_delta_logitPSI_stemIII,
  lymph_rmats_FDR                       = sig$FDR_Nava,
  lymph_rmats_PValue                    = sig$PValue_Nava,
  stringsAsFactors                      = FALSE
)

stopifnot(nrow(out) == n_cand, !anyDuplicated(out$event_id))
write.csv(out, OUT_CSV, row.names = FALSE)

cat(sprintf("\nWrote %s (%d events, %d columns)\n", OUT_CSV, nrow(out), ncol(out)))
cat("\ncds_class summary:\n"); print(sort(table(out$cds_class), decreasing = TRUE))
cat("\nconcordant summary:\n"); print(table(out$concordant, useNA = "ifany"))
cat("\nin_ndd_panel summary:\n"); print(table(out$in_ndd_panel, useNA = "ifany"))
cat(sprintf("\nT-loop signature events: %d; Stem III signature events: %d\n",
            sum(out$is_insT_Tloop_signature), sum(out$is_stemIII_signature)))
cat("\nNOTE: event_ids, gene names, and per-event GEL patient-derived stats ",
    "are GEL-protected. RE-internal only.\n", sep = "")
