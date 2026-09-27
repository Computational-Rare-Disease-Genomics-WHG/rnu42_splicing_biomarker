#!/usr/bin/env bash
# resubmit_post_rerun_batches.sh
# Re-submit rMATS post for the QC-fail rerun batches after duplicate .rmats
# cleanup. For each batch:
#   1. Delete merged/ (recreated fresh from the now-deduplicated prep dirs)
#   2. Delete the existing --od dir (rMATS post fails with FileExistsError if present)
#   3. Submit post using the existing _rerun wrapper
#
# Run from workspace root on GEL HPC after running cleanup_dup_rmats_prep.sh:
#   bash scripts/gel/resubmit_post_rerun_batches.sh [batch_n ...]
#   DRY_RUN=1 bash scripts/gel/resubmit_post_rerun_batches.sh

set -eo pipefail

WORKSPACE="/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
WRAPPER_DIR="${WORKSPACE}/rmats_inputs/batch_wrappers"
DRY_RUN="${DRY_RUN:-0}"

DEFAULT_BATCHES=(62 63 64 65 66 67 68 70 71 72 73 74 75 76 77 79 80)
if [[ $# -gt 0 ]]; then
    BATCHES=("$@")
else
    BATCHES=("${DEFAULT_BATCHES[@]}")
fi

echo "Workspace: $WORKSPACE"
echo "Dry run:   $DRY_RUN"
echo "Batches:   ${BATCHES[*]}"
echo ""

for BATCH in "${BATCHES[@]}"; do
    TMP_DIR="${WORKSPACE}/rmats_tmp/batch_${BATCH}"
    MERGED_DIR="${TMP_DIR}/merged"
    OD_DIR="${WORKSPACE}/outputs/rmats_controls/batch_${BATCH}"
    WRAPPER="${WRAPPER_DIR}/post_batch_${BATCH}_rerun.sh"
    JOB_NAME="rmats_post_b${BATCH}_rerun2"

    if [[ ! -f "$WRAPPER" ]]; then
        echo "Batch ${BATCH}: wrapper not found at ${WRAPPER} — skipping" >&2
        continue
    fi

    # 1. Remove merged/ so it is rebuilt from the deduplicated prep dirs
    if [[ -d "$MERGED_DIR" ]]; then
        echo "Batch ${BATCH}: rm -rf ${MERGED_DIR}"
        [[ "$DRY_RUN" -eq 0 ]] && rm -rf "$MERGED_DIR"
    else
        echo "Batch ${BATCH}: no merged/ dir to remove"
    fi

    # 2. Remove --od dir so rMATS can create it (Python 2 makedirs no exist_ok)
    if [[ -d "$OD_DIR" ]]; then
        echo "Batch ${BATCH}: rm -rf ${OD_DIR}"
        [[ "$DRY_RUN" -eq 0 ]] && rm -rf "$OD_DIR"
    else
        echo "Batch ${BATCH}: no od dir to remove"
    fi

    # 3. Submit post — no prep dependency, prep is already done
    CMD=(bsub
        -J "$JOB_NAME"
        -q short
        -P re_gecip_machine_learning
        -o "${WORKSPACE}/hpc_logs/${JOB_NAME}.%J.out"
        -e "${WORKSPACE}/hpc_logs/${JOB_NAME}.%J.err"
        -cwd "$WORKSPACE"
        -L /bin/bash
        -R "rusage[mem=12000]" -M 12000
        bash "$WRAPPER")

    echo "Batch ${BATCH}: submitting ${JOB_NAME}"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        "${CMD[@]}"
    else
        echo "  [DRY RUN] ${CMD[*]}"
    fi
    echo ""
done

echo "Done. Monitor with:"
echo "  bjobs -J 'rmats_post_b*_rerun2'"
echo "  grep -c 'rMATS post complete' hpc_logs/rmats_post_b*_rerun2.*.out 2>/dev/null"
