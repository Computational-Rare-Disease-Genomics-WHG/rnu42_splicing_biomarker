#!/usr/bin/env bash
# rerun_failed_qcfail_batches.sh
# For each batch in FAILED_BATCHES:
#   1. Identify which prep array indices failed (missing or empty _od dir)
#   2. Remove partial _od dirs so rMATS won't crash with FileExistsError
#   3. Write a prep rerun wrapper and submit it for the failed indices only
#   4. Submit post with -w "ended(rerun_prep_job)" so it runs once reruns finish
#
# Run from workspace root on GEL HPC:
#   bash scripts/gel/rerun_failed_qcfail_batches.sh [batch_n ...]
#
# Pass batch numbers as arguments to override the default list:
#   bash scripts/gel/rerun_failed_qcfail_batches.sh 62
#   bash scripts/gel/rerun_failed_qcfail_batches.sh 62 63 64
#
# Dry-run mode (print what would happen, submit nothing):
#   DRY_RUN=1 bash scripts/gel/rerun_failed_qcfail_batches.sh 62

set -eo pipefail

WORKSPACE="/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
WRAPPER_DIR="${WORKSPACE}/rmats_inputs/batch_wrappers"
DRY_RUN="${DRY_RUN:-0}"

# Batches with one or more failed prep indices — override by passing args
DEFAULT_BATCHES=(62 63 64 65 66 67 68 70 71 72 73 74 75 76 77 79 80)
if [[ $# -gt 0 ]]; then
    FAILED_BATCHES=("$@")
else
    FAILED_BATCHES=("${DEFAULT_BATCHES[@]}")
fi

echo "Workspace:  $WORKSPACE"
echo "Wrappers:   $WRAPPER_DIR"
echo "Dry run:    $DRY_RUN"
echo "Batches:    ${FAILED_BATCHES[*]}"
echo ""

for BATCH in "${FAILED_BATCHES[@]}"; do
    BAM_LIST="${WORKSPACE}/data/bam_paths_batch_${BATCH}.txt"
    TMP_DIR="${WORKSPACE}/rmats_tmp/batch_${BATCH}"
    OD_DIR="${WORKSPACE}/outputs/rmats_controls/batch_${BATCH}"

    if [[ ! -f "$BAM_LIST" ]]; then
        echo "Batch ${BATCH}: BAM list not found at ${BAM_LIST} — skipping" >&2
        continue
    fi

    N=$(wc -l < "$BAM_LIST")

    # ------------------------------------------------------------------
    # 1. Identify failed prep indices: tmp dir missing or contains no
    #    .rmats files. rMATS prep writes .rmats to --tmp (${i}/), NOT
    #    to --od (${i}_od/). Checking _od was a bug that caused all
    #    indices to be re-submitted even when prep had succeeded.
    # ------------------------------------------------------------------
    FAILED_INDICES=()
    for i in $(seq 1 "$N"); do
        tmp="${TMP_DIR}/${i}"
        # Mirror rmats_post.sh: only count .rmats files that are non-empty (-s).
        # A 0-byte .rmats is skipped by post's [[ -s "$f" ]] check and causes
        # the "missing N sample(s)" error even though the file exists.
        if [[ ! -d "$tmp" ]] || \
           [[ -z "$(find "$tmp" -maxdepth 1 -name "*.rmats" -size +0c 2>/dev/null)" ]]; then
            FAILED_INDICES+=("$i")
        fi
    done

    if [[ ${#FAILED_INDICES[@]} -eq 0 ]]; then
        echo "Batch ${BATCH}: all ${N} prep jobs succeeded — skipping"
        continue
    fi

    INDICES_CSV=$(IFS=','; echo "${FAILED_INDICES[*]}")
    echo "Batch ${BATCH}: ${#FAILED_INDICES[@]}/${N} failed — indices: ${INDICES_CSV}"

    # ------------------------------------------------------------------
    # 2. Remove partial _od dirs so rMATS can recreate them cleanly
    # ------------------------------------------------------------------
    for i in "${FAILED_INDICES[@]}"; do
        od="${TMP_DIR}/${i}_od"
        if [[ -d "$od" ]]; then
            echo "  rm -rf ${od}"
            [[ "$DRY_RUN" -eq 0 ]] && rm -rf "$od"
        fi
    done

    # ------------------------------------------------------------------
    # 3. Write prep rerun wrapper (same env vars as original, new name)
    # ------------------------------------------------------------------
    RERUN_WRAPPER="${WRAPPER_DIR}/prep_batch_${BATCH}_rerun.sh"
    echo "  Writing ${RERUN_WRAPPER}"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        cat > "$RERUN_WRAPPER" <<WRAPPER
#!/usr/bin/env bash
export BAM_LIST="${BAM_LIST}"
export TMP_DIR="${TMP_DIR}"
export OD_DIR="${OD_DIR}"
exec bash "${WORKSPACE}/scripts/gel/rmats_prep_array.sh"
WRAPPER
        chmod +x "$RERUN_WRAPPER"
    fi

    # ------------------------------------------------------------------
    # 4. Submit prep rerun array for failed indices only
    # ------------------------------------------------------------------
    RERUN_JOB_NAME="rmats_prep_b${BATCH}_rerun"
    POST_JOB_NAME="rmats_post_b${BATCH}_rerun"
    POST_WRAPPER="${WRAPPER_DIR}/post_batch_${BATCH}.sh"

    PREP_CMD=(bsub
        -J "${RERUN_JOB_NAME}[${INDICES_CSV}]%200"
        -q short
        -P re_gecip_machine_learning
        -o "${WORKSPACE}/hpc_logs/rmats_prep_b${BATCH}_rerun.%J.%I.out"
        -e "${WORKSPACE}/hpc_logs/rmats_prep_b${BATCH}_rerun.%J.%I.err"
        -cwd "$WORKSPACE"
        -L /bin/bash
        -R "rusage[mem=8000]" -M 8000
        bash "$RERUN_WRAPPER")

    echo "  Submitting prep rerun: ${RERUN_JOB_NAME}[${INDICES_CSV}]"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        "${PREP_CMD[@]}"
    else
        echo "  [DRY RUN] ${PREP_CMD[*]}"
    fi

    # ------------------------------------------------------------------
    # 5. Write post rerun wrapper + submit — depends on the rerun prep
    #    job name so it only starts once all rerun preps have ended.
    #    New wrapper so logs are clearly labelled _rerun in bjobs/hpc_logs.
    # ------------------------------------------------------------------
    ORIG_POST_WRAPPER="${WRAPPER_DIR}/post_batch_${BATCH}.sh"
    if [[ ! -f "$ORIG_POST_WRAPPER" ]]; then
        echo "  WARNING: original post wrapper not found at ${ORIG_POST_WRAPPER} — skipping post resubmit" >&2
        continue
    fi

    # Read EXPECTED_N from the original post wrapper so the rerun carries
    # the same value (avoids re-deriving it from the BAM list).
    EXPECTED_N=$(grep 'export EXPECTED_N=' "$ORIG_POST_WRAPPER" | \
                 sed 's/.*export EXPECTED_N="\([0-9]*\)".*/\1/')
    if [[ -z "$EXPECTED_N" ]]; then
        echo "  WARNING: could not parse EXPECTED_N from ${ORIG_POST_WRAPPER} — skipping post resubmit" >&2
        continue
    fi

    RERUN_POST_WRAPPER="${WRAPPER_DIR}/post_batch_${BATCH}_rerun.sh"
    echo "  Writing ${RERUN_POST_WRAPPER}"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        cat > "$RERUN_POST_WRAPPER" <<WRAPPER
#!/usr/bin/env bash
export EXPECTED_N="${EXPECTED_N}"
export TMP_DIR="${TMP_DIR}"
export OD_DIR="${OD_DIR}"
# Clear stale merged/ and output dir so post rebuilds from current prep files
rm -rf "${TMP_DIR}/merged"
rm -rf "${OD_DIR}"
exec bash "${WORKSPACE}/scripts/gel/rmats_post.sh"
WRAPPER
        chmod +x "$RERUN_POST_WRAPPER"
    fi

    POST_CMD=(bsub
        -J "$POST_JOB_NAME"
        -q short
        -P re_gecip_machine_learning
        -o "${WORKSPACE}/hpc_logs/rmats_post_b${BATCH}_rerun.%J.out"
        -e "${WORKSPACE}/hpc_logs/rmats_post_b${BATCH}_rerun.%J.err"
        -cwd "$WORKSPACE"
        -L /bin/bash
        -R "rusage[mem=12000]" -M 12000
        -w "ended(\"${RERUN_JOB_NAME}\")"
        bash "$RERUN_POST_WRAPPER")

    echo "  Submitting post rerun: ${POST_JOB_NAME} (waits on ${RERUN_JOB_NAME})"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        "${POST_CMD[@]}"
    else
        echo "  [DRY RUN] ${POST_CMD[*]}"
    fi

    echo ""
done

echo "Done. Monitor with:"
echo "  bjobs -J 'rmats_prep_b*_rerun'"
echo "  bjobs -J 'rmats_post_b*_rerun'"
