#!/usr/bin/env Rscript
# =============================================================================
# export_supplementary_tables_xlsx.R
# Bundles the four paper-ready supplementary table CSVs
# (outputs/supplementary_tables/Table_S{1..4}_*.csv) into one workbook,
# one sheet per table. Run this after regenerating any of the four CSVs.
#
# Inputs:
#   outputs/supplementary_tables/Table_S1_candidate_events.csv
#   outputs/supplementary_tables/Table_S2_insT_signature.csv
#   outputs/supplementary_tables/Table_S3_stemIII_signature.csv
#   outputs/supplementary_tables/Table_S4_ipsc_validation.csv
#
# Output:
#   outputs/supplementary_tables/all_supplementary_tables.xlsx
#
# Run from project root: Rscript scripts/export_supplementary_tables_xlsx.R
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(readr)
})

TABLE_DIR <- "outputs/supplementary_tables"
OUT_XLSX  <- file.path(TABLE_DIR, "supplementary_tables.xlsx")

tables <- c(
  Table_S1 = "Table_S1_candidate_events.csv",
  Table_S2 = "Table_S2_insT_signature.csv",
  Table_S3 = "Table_S3_stemIII_signature.csv",
  Table_S4 = "Table_S4_ipsc_validation.csv"
)
paths <- file.path(TABLE_DIR, tables)
stopifnot("Missing supplementary table CSV(s)" = all(file.exists(paths)))

wb <- createWorkbook()
for (i in seq_along(tables)) {
  sheet <- names(tables)[i]
  df <- read_csv(paths[i], show_col_types = FALSE)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df, headerStyle = createStyle(textDecoration = "bold"))
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = "auto")
  cat(sprintf("%s: %d rows, %d cols (%s)\n", sheet, nrow(df), ncol(df), tables[i]))
}

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat(sprintf("\nSaved: %s\n", OUT_XLSX))
