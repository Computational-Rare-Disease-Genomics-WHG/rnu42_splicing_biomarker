#!/usr/bin/env bash
# rmats_post.sh
# rMATS post step — merges all prep outputs into final PSI tables.
# Submitted via submit_rmats_batched.sh. Do not submit directly.
#
# Environment variables (set by submission script):
#   GTF_PATH  — path to GTF file
#   TMP_DIR   — rMATS tmp directory   (default: ${WORKSPACE}/rmats_tmp)
#   OD_DIR    — rMATS output directory (default: ${WORKSPACE}/outputs/rmats_controls)
#   WORKSPACE — absolute path to workspace root
#
# BAM list is derived from the .rmats files in TMP_DIR rather than passed
# externally — this guarantees exact matching and avoids ordering mismatches.

module load singularity/4.1.1

WORKSPACE="${WORKSPACE:-/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature}"
GTF_PATH="${GTF_PATH:-${WORKSPACE}/rmats_inputs/gencode.v49.annotation.gtf}"
TMP_DIR="${TMP_DIR:-${WORKSPACE}/rmats_tmp}"
OD_DIR="${OD_DIR:-${WORKSPACE}/outputs/rmats_controls}"

# Do NOT mkdir OD_DIR here — rMATS calls makedirs() without exist_ok=True
# (Python 2 compat) and will crash with FileExistsError if --od already exists.
# rMATS creates --od itself during the post step.

# Merge per-job .rmats files into a shared directory for post.
# Each prep job writes to its own --tmp subdirectory (TMP_DIR/INDEX/) to prevent
# filename conflicts when concurrent jobs run at the same second. Here we copy
# all .rmats files into a flat merged/ dir, prefixed with the job index so
# filenames are guaranteed unique.
MERGED_TMP="${TMP_DIR}/merged"
mkdir -p "$MERGED_TMP"

for job_dir in "${TMP_DIR}"/[0-9]*/; do
    job_idx=$(basename "$job_dir")
    for f in "${job_dir}"*.rmats; do
        [[ -s "$f" ]] && cp "$f" "${MERGED_TMP}/${job_idx}_$(basename "$f")"
    done
done

# Check that all prep jobs produced a non-empty .rmats file
actual_rmats=$(find "${MERGED_TMP}" -name "*.rmats" -size +0c 2>/dev/null | wc -l)
if [[ -n "$EXPECTED_N" && "$EXPECTED_N" -gt 0 ]]; then
    if [[ "$actual_rmats" -lt "$EXPECTED_N" ]]; then
        echo "ERROR: prep incomplete — expected ${EXPECTED_N} .rmats files, found ${actual_rmats} in ${MERGED_TMP}" >&2
        echo "       Missing: $(( EXPECTED_N - actual_rmats )) sample(s)" >&2
        echo "       Check prep logs in hpc_logs/ to identify failed array indices" >&2
        exit 1
    fi
    echo "Prep check passed: ${actual_rmats}/${EXPECTED_N} .rmats files in ${MERGED_TMP}"
else
    echo "WARNING: EXPECTED_N not set — skipping prep completeness check (found ${actual_rmats} .rmats files)"
fi

# Build comma-separated BAM list from the original per-batch BAM list.
# rMATS --b1 requires a single line, comma-separated.
# Using bam_paths_batch_N.txt guarantees paths match exactly what was used in prep.
BAM_LIST="${MERGED_TMP}/bam_list_for_post.txt"
BATCH_NUM=$(basename "$TMP_DIR" | grep -oP '\d+$')
ORIG_BAM_LIST="${WORKSPACE}/data/bam_paths_batch_${BATCH_NUM}.txt"
if [[ ! -f "$ORIG_BAM_LIST" ]]; then
    echo "ERROR: original BAM list not found: $ORIG_BAM_LIST" >&2
    exit 1
fi
paste -sd ',' "$ORIG_BAM_LIST" > "$BAM_LIST"

echo "=== rmats_post diagnostics ==="
echo "WORKSPACE:  $WORKSPACE"
echo "GTF:        $GTF_PATH"
echo "BAM list:   $BAM_LIST ($(wc -l < "$BAM_LIST") BAMs)"
echo "TMP:        $MERGED_TMP (merged from per-job subdirs)"
echo "Output:     $OD_DIR"
echo "SIF:        ${WORKSPACE}/rmats_inputs/rmats.sif"
echo "--- file existence checks ---"
[[ -f "$GTF_PATH" ]] && echo "GTF:     EXISTS" || echo "GTF:     MISSING"
[[ -f "${WORKSPACE}/rmats_inputs/rmats.sif" ]] && echo "SIF:     EXISTS" || echo "SIF:     MISSING"
[[ -d "$TMP_DIR" ]] && echo "TMP_DIR: EXISTS ($(ls "$TMP_DIR" | wc -l) files)" || echo "TMP_DIR: MISSING"
echo "==============================="

singularity exec \
    --bind /nas/weka.gel.zone/pgen_int_data_resources:/nas/weka.gel.zone/pgen_int_data_resources \
    --bind /gel_data_resources:/gel_data_resources \
    --bind "${WORKSPACE}:${WORKSPACE}" \
    "${WORKSPACE}/rmats_inputs/rmats.sif" \
    python /rmats/rmats.py \
        --b1 "$BAM_LIST" \
        --gtf "$GTF_PATH" \
        -t paired \
        --anchorLength 1 \
        --libType fr-firststrand \
        --novelSS \
        --readLength 100 \
        --variable-read-length \
        --allow-clipping \
        --statoff \
        --task post \
        --od "$OD_DIR" \
        --tmp "$MERGED_TMP"

echo "rMATS post complete. Output in: $OD_DIR"
