#!/usr/bin/env Rscript
# =============================================================================
# export_rmats_summary_allTypes.R
# Strip the GEL whole-blood dev-cohort rMATS output (11 ReNU cases vs 33 GEL
# controls, all 5 event types) down to group-level summary stats, for airlock
# export. Drops every per-sample column (IJC/SJC per-sample junction counts,
# IncLevel1/2 per-sample PSI lists) -- keeps only event identity, coordinates,
# and the aggregate group-comparison stats (PValue, FDR, IncLevelDifference).
#
# Purpose: extends figure4.R panel D (cross-tissue direction-concordance
# replication, currently iPSC vs Nava lymphocyte only) to also compare
# against GEL whole blood. Panel D needs unfiltered per-event FDR + direction
# for ALL events, all 5 types -- not a pre-selected candidate pool -- so it
# can test arbitrary iPSC-significant coordinates against GEL significance
# without circularity.
#
# Columns dropped (per event type): IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2,
# SJC_SAMPLE_2, IncLevel1, IncLevel2. Everything else retained verbatim,
# including column order and the duplicate "ID" column, so
# figure4.R::build_coord_nostrand() (positional column indexing) works on
# the exported file unmodified.
#
# Input:  data/rmats/{SE,A5SS,A3SS,MXE,RI}.MATS.JC.txt
#         (whole-blood dev-cohort rMATS: 11 ReNU cases vs 33 GEL controls)
# Output: outputs/rmats_export/{SE,A5SS,A3SS,MXE,RI}_summary.csv
#
# DATA PROTECTION: gene names, coordinates, and per-event GEL FDR/dPSI are
# GEL patient-derived data -- RE-internal only until airlocked.
#
# Run: pixi run Rscript scripts/gel/export_rmats_summary_allTypes.R
#   or Rscript scripts/gel/export_rmats_summary_allTypes.R
# =============================================================================

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

RMATS_DIR   <- Sys.getenv("GEL_RMATS_DIR", "data/rmats")
OUT_DIR     <- Sys.getenv("OUT_DIR", "outputs/rmats_export")
EVENT_TYPES <- c("SE", "A5SS", "A3SS", "MXE", "RI")

DROP_COLS <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
               "IncLevel1", "IncLevel2")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (et in EVENT_TYPES) {
  in_path  <- file.path(RMATS_DIR, paste0(et, ".MATS.JC.txt"))
  out_path <- file.path(OUT_DIR, paste0(et, "_summary.csv"))
  stopifnot("rMATS input file not found" = file.exists(in_path))

  d <- read.delim(in_path, stringsAsFactors = FALSE, check.names = FALSE)
  present <- intersect(DROP_COLS, names(d))
  stopifnot("Expected per-sample columns missing -- check rMATS header" =
              length(present) == length(DROP_COLS))

  d_out <- d[, !(names(d) %in% DROP_COLS), drop = FALSE]

  stopifnot(
    "Per-sample columns leaked into export" = !any(DROP_COLS %in% names(d_out)),
    "Row count changed during strip"        = nrow(d_out) == nrow(d)
  )

  write.csv(d_out, out_path, row.names = FALSE)
  cat(sprintf("%s: %d events, %d -> %d columns -> %s\n",
              et, nrow(d_out), ncol(d), ncol(d_out), out_path))
}

cat("\nNOTE: gene names, coordinates, and per-event GEL FDR/dPSI are GEL",
    "patient-derived data. RE-internal only -- do not airlock without review.\n")
