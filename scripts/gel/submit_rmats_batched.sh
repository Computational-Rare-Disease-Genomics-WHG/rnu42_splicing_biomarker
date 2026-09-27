#!/usr/bin/env bash
# submit_rmats_batched.sh
# Submit rMATS prep+post as batched LSF jobs (BATCH_SIZE BAMs per batch).
#
# Why batched: rMATS post memory scales superlinearly with --novelSS.
# 100-sample batches keep post memory ~5-10 GB (manageable on GEL HPC).
# Each batch has its own tmp and output directory to prevent cross-contamination.
#
# Run from workspace root on GEL HPC:
#   bash scripts/gel/submit_rmats_batched.sh
#
# Before running:
#   - data/bam_paths.txt must exist (run submit_rmats_controls.sh or generate manually)
#   - Run a single batch first to verify: set N_BATCHES_MAX=1 below
#
# Env-var overrides (used to add a second wave of batches for the ~1000
# QC-failing samples without disturbing the existing batch_1..N outputs):
#   MANIFEST_CSV   csv to derive BAM_LIST from (default: data/control_manifest.csv)
#   BAM_LIST       BAM paths file              (default: data/bam_paths.txt)
#   BATCH_OFFSET   integer added to each batch number for output dir + job names
#                  (default: 0). Set to the highest existing batch_N when adding
#                  more samples so new batches land at batch_(N+1), batch_(N+2), ...
#                  next to the originals; combine_rmats_batches*.R picks them up
#                  automatically.
#
# Example — submit QC-fail samples as continuation batches:
#   python -u scripts/gel/get_control_manifest.py --mode fail
#   MANIFEST_CSV=data/qcfail_manifest.csv \
#   BAM_LIST=data/bam_paths_qcfail.txt \
#   BATCH_OFFSET=$(ls -d outputs/rmats_controls/batch_* | sed 's|.*batch_||' | sort -n | tail -1) \
#   bash scripts/gel/submit_rmats_batched.sh
#
# After all batches complete, combine outputs with scripts/gel/combine_rmats_batches_cand_pool.R

set -eo pipefail

WORKSPACE="/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
MANIFEST_CSV="${MANIFEST_CSV:-${WORKSPACE}/data/control_manifest.csv}"
BAM_LIST="${BAM_LIST:-${WORKSPACE}/data/bam_paths.txt}"
BATCH_OFFSET="${BATCH_OFFSET:-0}"
BATCH_SIZE=100

# Resolve relative paths against the workspace
[[ "$MANIFEST_CSV" != /* ]] && MANIFEST_CSV="${WORKSPACE}/${MANIFEST_CSV}"
[[ "$BAM_LIST"     != /* ]] && BAM_LIST="${WORKSPACE}/${BAM_LIST}"

# Generate BAM list from manifest if not already present
if [[ ! -f "$BAM_LIST" ]]; then
    echo "Generating $BAM_LIST from $(basename "$MANIFEST_CSV")..."
    /resources/conda/miniconda3/envs/py3pypirev2/bin/python -c "
import pandas as pd
df = pd.read_csv('${MANIFEST_CSV}')
df['file_path'].to_csv('${BAM_LIST}', index=False, header=False)
print(f'  {len(df)} BAMs written')
"
fi

# Set START_BATCH to resume from a specific batch (e.g. 3 to skip already-completed 1+2)
START_BATCH=1

# Set to a small number to submit only the first N batches (e.g. 1 for a verification run)
# Set to 9999 to submit all batches
N_BATCHES_MAX=9999

WRAPPER_DIR="${WORKSPACE}/rmats_inputs/batch_wrappers"
mkdir -p "${WRAPPER_DIR}" "${WORKSPACE}/hpc_logs"

N_TOTAL=$(wc -l < "$BAM_LIST")
N_BATCHES=$(( (N_TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))
N_SUBMIT=$(( N_BATCHES < N_BATCHES_MAX ? N_BATCHES : N_BATCHES_MAX ))
ACTUAL_SUBMIT=$(( N_SUBMIT - START_BATCH + 1 ))

echo "Manifest CSV:   $MANIFEST_CSV"
echo "BAM list:       $BAM_LIST"
echo "Total BAMs:     $N_TOTAL"
echo "Batch size:     $BATCH_SIZE"
echo "Batch offset:   $BATCH_OFFSET"
echo "Total batches:  $N_BATCHES"
echo "Start batch:    $START_BATCH"
echo "Submitting:     $ACTUAL_SUBMIT (local ${START_BATCH}–${N_SUBMIT}, absolute $((START_BATCH+BATCH_OFFSET))–$((N_SUBMIT+BATCH_OFFSET)))"
echo ""

for BATCH in $(seq "$START_BATCH" "$N_SUBMIT"); do
    START=$(( (BATCH - 1) * BATCH_SIZE + 1 ))
    END=$(( BATCH * BATCH_SIZE ))
    [[ $END -gt $N_TOTAL ]] && END=$N_TOTAL
    N_BATCH=$(( END - START + 1 ))

    ABS_BATCH=$(( BATCH + BATCH_OFFSET ))
    BATCH_BAM_LIST="${WORKSPACE}/data/bam_paths_batch_${ABS_BATCH}.txt"
    BATCH_TMP="${WORKSPACE}/rmats_tmp/batch_${ABS_BATCH}"
    BATCH_OD="${WORKSPACE}/outputs/rmats_controls/batch_${ABS_BATCH}"

    if [[ -d "$BATCH_OD" ]]; then
        echo "Batch ${ABS_BATCH}: output dir already exists ($BATCH_OD) — skipping" >&2
        continue
    fi

    # Write per-batch BAM list (absolute paths only, no header)
    sed -n "${START},${END}p" "$BAM_LIST" > "$BATCH_BAM_LIST"

    mkdir -p "$BATCH_TMP"
    # Do NOT pre-create BATCH_OD — rMATS post calls makedirs() without exist_ok
    # and will crash with FileExistsError if the output directory already exists

    # Write prep wrapper — exports env vars before calling the prep array worker.
    # This avoids the -env flag which is incompatible with -L /bin/bash login shells.
    PREP_WRAPPER="${WRAPPER_DIR}/prep_batch_${ABS_BATCH}.sh"
    cat > "$PREP_WRAPPER" <<WRAPPER
#!/usr/bin/env bash
export BAM_LIST="${BATCH_BAM_LIST}"
export TMP_DIR="${BATCH_TMP}"
export OD_DIR="${BATCH_OD}"
exec bash "${WORKSPACE}/scripts/gel/rmats_prep_array.sh"
WRAPPER
    chmod +x "$PREP_WRAPPER"

    # Write post wrapper
    POST_WRAPPER="${WRAPPER_DIR}/post_batch_${ABS_BATCH}.sh"
    cat > "$POST_WRAPPER" <<WRAPPER
#!/usr/bin/env bash
export EXPECTED_N="${N_BATCH}"
export TMP_DIR="${BATCH_TMP}"
export OD_DIR="${BATCH_OD}"
exec bash "${WORKSPACE}/scripts/gel/rmats_post.sh"
WRAPPER
    chmod +x "$POST_WRAPPER"

    PREP_JOB_NAME="rmats_prep_b${ABS_BATCH}"
    POST_JOB_NAME="rmats_post_b${ABS_BATCH}"

    # Submit prep array — cap at 200 concurrent across all arrays
    bsub \
        -J "${PREP_JOB_NAME}[1-${N_BATCH}]%200" \
        -q short \
        -P re_gecip_machine_learning \
        -o "${WORKSPACE}/hpc_logs/rmats_prep_b${ABS_BATCH}.%J.%I.out" \
        -e "${WORKSPACE}/hpc_logs/rmats_prep_b${ABS_BATCH}.%J.%I.err" \
        -cwd "$WORKSPACE" \
        -L /bin/bash \
        -R "rusage[mem=8000]" -M 8000 \
        bash "$PREP_WRAPPER"

    # Submit post job — waits for its prep array, not other batches
    bsub \
        -J "$POST_JOB_NAME" \
        -q short \
        -P re_gecip_machine_learning \
        -o "${WORKSPACE}/hpc_logs/rmats_post_b${ABS_BATCH}.%J.out" \
        -e "${WORKSPACE}/hpc_logs/rmats_post_b${ABS_BATCH}.%J.err" \
        -cwd "$WORKSPACE" \
        -L /bin/bash \
        -R "rusage[mem=12000]" -M 12000 \
        -w "ended(\"${PREP_JOB_NAME}\")" \
        bash "$POST_WRAPPER"

    echo "Batch ${ABS_BATCH} (local ${BATCH}/${N_SUBMIT}): ${N_BATCH} BAMs (lines ${START}-${END})"
done

echo ""
echo "All $N_SUBMIT batch(es) submitted."
echo ""
echo "Monitor prep:  bjobs -J 'rmats_prep_b*'"
echo "Monitor post:  bjobs -J 'rmats_post_b*'"
echo "Failed jobs:   bjobs -J 'rmats_prep_b*' | grep EXIT"
echo ""
echo "After all post jobs complete, combine outputs with:"
echo "  Rscript scripts/gel/combine_rmats_batches_cand_pool.R"
