#!/usr/bin/env bash
# build_salmon_quant_paths.sh
# Build data/salmon_quant_paths.csv (platekey,file_path) for analysis C of
# scripts/gel/v4_signature_extra_analyses.R.
#
# Reads:
#   data/renu_11_cases_paths.txt          (single comma-separated line of BAM paths)
#   data/renu_11_cases_control_paths.txt  (likewise)
# Writes:
#   data/salmon_quant_paths.csv           (45 lines: header + 11 cases + 33 ctrls)
#
# The platekey column uses rMATS sample IDs (gel:case:N, gel:ctrl:N) so the
# downstream R script's join against psi_matrix$samples$sample_id works.
# The file_path column swaps `.bam` -> `.quant.genes.sf` so we point at GEL's
# pre-aggregated gene-level salmon output.
#
# Run from workspace root:
#   bash scripts/gel/build_salmon_quant_paths.sh
#
# Why a script file instead of an inline one-liner: pasting multi-line
# commands into VS Code's terminal can capture shell-integration escape
# sequences (OSC 633 markers) into redirected output, corrupting the CSV.

set -euo pipefail

CASES_IN=data/renu_11_cases_paths.txt
CTRLS_IN=data/renu_11_cases_control_paths.txt
OUT=data/salmon_quant_paths.csv

for f in "$CASES_IN" "$CTRLS_IN"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

{
    echo "platekey,file_path"
    tr ',' '\n' < "$CASES_IN" | grep -E '\.bam$' | \
        awk '{i++; sub(/\.bam$/, ".quant.genes.sf"); print "gel:case:" i "," $0}'
    tr ',' '\n' < "$CTRLS_IN" | grep -E '\.bam$' | \
        awk '{i++; sub(/\.bam$/, ".quant.genes.sf"); print "gel:control:" i "," $0}'
} > "$OUT"

n=$(wc -l < "$OUT")
echo "Wrote $OUT (${n} lines, expected 45)"
head -3 "$OUT"
