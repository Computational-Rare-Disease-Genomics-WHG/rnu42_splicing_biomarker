#!/usr/bin/env bash
#
# null_lasso_rmats_stat.sh  (HPC, submitted via bsub)
# Middle step of the null_lasso_within_batch experiment: takes the drawn
# 0-based indices written by null_lasso_within_batch.R phase 1, runs
# rMATS prepare_stat_inputs.py + rmats.py --task stat inside the rMATS
# singularity container to produce case-vs-control LRT p-values (identical
# to what rMATS would emit end-to-end).
#
# Args (positional; env is wiped by -L /bin/bash on GEL LSF):
#   $1  BATCH_ID   integer batch number (must match phase 1)
#   $2  DRAW_SEED  integer seed (used only for OUT_TAG default)
#   $3  OUT_TAG    directory tag; default batch${BATCH_ID}_seed${DRAW_SEED}
#
# Reads: <WORKSPACE>/outputs/null_within_batch/<OUT_TAG>/drawn_samples.csv
# Writes:<WORKSPACE>/outputs/null_within_batch/<OUT_TAG>/rmats_stat_od/
#
# Submit example (printed by phase 1):
#   bsub -q short -J null_batch1_seed1 \
#     -o hpc_logs/null_%J.out -e hpc_logs/null_%J.err \
#     -env "all,BATCH_ID=1,DRAW_SEED=1,OUT_TAG=batch1_seed1" \
#     scripts/gel/null_lasso_rmats_stat.sh
#
# NOTE: LSF -env may strip vars when combined with -L /bin/bash; hardcode
# WORKSPACE as default so the script is self-contained.

set -euo pipefail

module load singularity/4.1.1

WORKSPACE="/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"

if [[ $# -lt 1 ]]; then
    echo "usage: $0 BATCH_ID [DRAW_SEED] [OUT_TAG]" >&2
    exit 1
fi
BATCH_ID="$1"
DRAW_SEED="${2:-1}"
OUT_TAG="${3:-batch${BATCH_ID}_seed${DRAW_SEED}}"

OUT_DIR="${WORKSPACE}/outputs/null_within_batch/${OUT_TAG}"
STAT_OD="${OUT_DIR}/rmats_stat_od"
STAT_TMP="${OUT_DIR}/rmats_stat_tmp"
DRAWN_CSV="${OUT_DIR}/drawn_samples.csv"
OLD_OD="${WORKSPACE}/outputs/rmats_controls/batch_${BATCH_ID}"
SIF="${WORKSPACE}/rmats_inputs/rmats.sif"

if [[ ! -f "$DRAWN_CSV" ]]; then
    echo "ERROR: drawn_samples.csv missing at $DRAWN_CSV" >&2
    echo "       Run phase 1 (Rscript scripts/gel/null_lasso_within_batch.R) on Desktop first." >&2
    exit 1
fi
if [[ ! -d "$OLD_OD" ]]; then
    echo "ERROR: original rMATS --od missing at $OLD_OD" >&2
    exit 1
fi
if [[ ! -f "$SIF" ]]; then
    echo "ERROR: rMATS singularity image missing at $SIF" >&2
    exit 1
fi

# --new-output-dir must NOT pre-exist (prepare_stat_inputs.py mkdirs it,
# and rmats.py --task stat re-uses it).
rm -rf "$STAT_OD" "$STAT_TMP"
mkdir -p "$STAT_TMP"

# Extract 0-based indices for cases and controls from drawn_samples.csv.
# CSV columns: role,platekey,batch_idx0
# R's write.csv wraps character fields in double quotes -- strip them
# before matching.
CASE_IDX=$(awk -F, 'NR>1 {gsub(/"/,"",$1); if ($1=="case")    print $3}' "$DRAWN_CSV" | paste -sd,)
CTRL_IDX=$(awk -F, 'NR>1 {gsub(/"/,"",$1); if ($1=="control") print $3}' "$DRAWN_CSV" | paste -sd,)

if [[ -z "$CASE_IDX" || -z "$CTRL_IDX" ]]; then
    echo "ERROR: could not parse case/control indices from $DRAWN_CSV" >&2
    exit 1
fi

N_CASE=$(awk -F, 'NR>1 {gsub(/"/,"",$1); if ($1=="case")    c++} END{print c+0}' "$DRAWN_CSV")
N_CTRL=$(awk -F, 'NR>1 {gsub(/"/,"",$1); if ($1=="control") c++} END{print c+0}' "$DRAWN_CSV")

echo "=== null_lasso_rmats_stat ==="
echo "BATCH_ID   : $BATCH_ID"
echo "OUT_TAG    : $OUT_TAG"
echo "OLD_OD     : $OLD_OD"
echo "STAT_OD    : $STAT_OD"
echo "cases (${N_CASE}) idx: $CASE_IDX"
echo "ctrls (${N_CTRL}) idx: $CTRL_IDX"
echo "============================"

echo "Running prepare_stat_inputs.py ..."
singularity exec \
    --bind "${WORKSPACE}:${WORKSPACE}" \
    "$SIF" \
    python /rmats/rMATS_P/prepare_stat_inputs.py \
        --old-output-dir   "$OLD_OD" \
        --new-output-dir   "$STAT_OD" \
        --group-1-indices  "$CASE_IDX" \
        --group-2-indices  "$CTRL_IDX"

# We only score A5SS downstream. rMATS --task stat has no per-event-type
# switch, so gut the other four event types down to header-only files
# before the stat call -- 0 events => LRT skips, cuts wall time from ~1h
# to seconds.
echo "Truncating non-A5SS event tables to header-only ..."
for et in SE A3SS RI MXE; do
    for f in "fromGTF.${et}.txt" "JC.raw.input.${et}.txt" "JCEC.raw.input.${et}.txt"; do
        if [[ -f "${STAT_OD}/${f}" ]]; then
            head -n1 "${STAT_OD}/${f}" > "${STAT_OD}/${f}.tmp" && \
                mv "${STAT_OD}/${f}.tmp" "${STAT_OD}/${f}"
        fi
    done
done

echo "Running rmats.py --task stat ..."
singularity exec \
    --bind "${WORKSPACE}:${WORKSPACE}" \
    "$SIF" \
    python /rmats/rmats.py \
        --od   "$STAT_OD" \
        --tmp  "$STAT_TMP" \
        --task stat

echo "rMATS stat complete. Output in: $STAT_OD"
ls -1 "$STAT_OD"/*.MATS.JC.txt 2>/dev/null | sed 's/^/  /'

echo "Now re-run phase 2 on Desktop:"
echo "  BATCH_ID=$BATCH_ID DRAW_SEED=$DRAW_SEED pixi run Rscript scripts/gel/null_lasso_within_batch.R"
