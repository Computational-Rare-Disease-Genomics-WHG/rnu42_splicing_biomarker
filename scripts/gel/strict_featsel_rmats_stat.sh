#!/usr/bin/env bash
#
# strict_featsel_rmats_stat.sh  (HPC, submitted via bsub)
# Phase 2 of the strict feature-selection rerun: runs rmats.py --task stat in
# each input directory written by strict_featsel_prepare.R. Runs that already
# have an A5SS.MATS.JC.txt are skipped.
#
# Submit from the workspace root:
#   bsub -q medium -P re_gecip_machine_learning -J strict_featsel_stat \
#     -o hpc_logs/strict_featsel_%J.out -e hpc_logs/strict_featsel_%J.err \
#     bash scripts/gel/strict_featsel_rmats_stat.sh
#
# Paths are hardcoded: LSF -env vars are reset by the -L /bin/bash login shell.

set -eo pipefail

module load singularity/4.1.1

WORKSPACE="/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
STAT_DIR="${WORKSPACE}/outputs/strict_featsel/rmats_stat"
SIF="${WORKSPACE}/rmats_inputs/rmats.sif"
RUNS="gel_all nava_all gel_train nava_train"

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: rMATS singularity image missing at $SIF" >&2
    exit 1
fi

for run in $RUNS; do
    OD="${STAT_DIR}/${run}"
    TMP="${STAT_DIR}/${run}_tmp"
    if [[ ! -f "${OD}/JC.raw.input.A5SS.txt" ]]; then
        echo "ERROR: ${OD} has no inputs -- run strict_featsel_prepare.R first" >&2
        exit 1
    fi
    if [[ -f "${OD}/A5SS.MATS.JC.txt" ]]; then
        echo "=== ${run}: output already present, skipping ==="
        continue
    fi
    echo "=== ${run}: rmats.py --task stat ==="
    rm -rf "$TMP"
    mkdir -p "$TMP"
    singularity exec \
        --bind "${WORKSPACE}:${WORKSPACE}" \
        "$SIF" \
        python /rmats/rmats.py \
            --od   "$OD" \
            --tmp  "$TMP" \
            --task stat
    if [[ ! -f "${OD}/A5SS.MATS.JC.txt" ]]; then
        echo "ERROR: ${run} finished without writing A5SS.MATS.JC.txt" >&2
        exit 1
    fi
    echo "  wrote ${OD}/A5SS.MATS.JC.txt"
done

echo "All runs done. Next (Desktop): pixi run Rscript scripts/gel/strict_featsel_insT_lasso.R"
