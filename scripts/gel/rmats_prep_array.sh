#!/usr/bin/env bash
# rmats_prep_array.sh
# rMATS prep step — one job per BAM, run as LSF job array.
# Submitted by submit_rmats_controls.sh or submit_rmats_trial.sh.
# Do not submit directly.
#
# Environment variables (set by submission script):
#   WORKSPACE — absolute path to workspace root (set by submission script)
#   GTF_PATH  — path to GTF file
#   BAM_LIST  — path to BAM paths file (default: ${WORKSPACE}/data/bam_paths.txt)
#   TMP_DIR   — rMATS tmp directory   (default: ${WORKSPACE}/rmats_tmp)
#   OD_DIR    — rMATS output directory (default: ${WORKSPACE}/outputs/rmats_controls)

module load singularity/4.1.1

WORKSPACE="${WORKSPACE:-/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature}"
GTF_PATH="${GTF_PATH:-${WORKSPACE}/rmats_inputs/gencode.v49.annotation.gtf}"
BAM_LIST="${BAM_LIST:-${WORKSPACE}/data/bam_paths.txt}"
TMP_DIR="${TMP_DIR:-${WORKSPACE}/rmats_tmp}"
OD_DIR="${OD_DIR:-${WORKSPACE}/outputs/rmats_controls}"

# Get BAM path for this array index
BAM_PATH=$(sed -n "${LSB_JOBINDEX}p" "$BAM_LIST")

if [[ -z "$BAM_PATH" ]]; then
    echo "ERROR: no BAM path found for index $LSB_JOBINDEX in $BAM_LIST" >&2
    exit 1
fi

# Each job uses its own --tmp and --od subdirectories.
# --tmp: prevents .rmats filename conflicts (rMATS uses timestamp-based names;
#        concurrent jobs at the same second silently overwrite each other).
# --od:  rMATS calls os.makedirs() without exist_ok=True (Python 2 compat);
#        if two jobs share --od, they race to makedirs(--od/tmp) and the loser
#        crashes with FileExistsError. Per-job --od dirs are never pre-created
#        so rMATS creates them fresh with no conflict.
JOB_TMP="${TMP_DIR}/${LSB_JOBINDEX}"
JOB_OD="${TMP_DIR}/${LSB_JOBINDEX}_od"
B1_FILE="${JOB_TMP}/b1_${LSB_JOBINDEX}.txt"
mkdir -p "$JOB_TMP"   # must pre-exist for --tmp; rMATS does not makedirs --tmp
# Do NOT mkdir JOB_OD — rMATS must create it (makedirs fails on pre-existing dirs)
echo "$BAM_PATH" > "$B1_FILE"

echo "=== rmats_prep_array diagnostics ==="
echo "Index:      $LSB_JOBINDEX"
echo "WORKSPACE:  $WORKSPACE"
echo "BAM:        $BAM_PATH"
echo "GTF:        $GTF_PATH"
echo "BAM_LIST:   $BAM_LIST"
echo "TMP:        $TMP_DIR"
echo "OD:         $OD_DIR"
echo "SIF:        ${WORKSPACE}/rmats_inputs/rmats.sif"
echo "--- file existence checks ---"
[[ -f "$GTF_PATH" ]]  && echo "GTF:     EXISTS" || echo "GTF:     MISSING"
[[ -f "${WORKSPACE}/rmats_inputs/rmats.sif" ]] && echo "SIF:     EXISTS" || echo "SIF:     MISSING"
[[ -f "$BAM_PATH" ]]  && echo "BAM:     EXISTS" || echo "BAM:     MISSING (may be remote)"
echo "--- singularity version ---"
singularity --version
echo "====================================="

singularity exec \
    --bind /nas/weka.gel.zone/pgen_int_data_resources:/nas/weka.gel.zone/pgen_int_data_resources \
    --bind /gel_data_resources:/gel_data_resources \
    --bind "${WORKSPACE}:${WORKSPACE}" \
    "${WORKSPACE}/rmats_inputs/rmats.sif" \
    python /rmats/rmats.py \
        --b1 "$B1_FILE" \
        --gtf "$GTF_PATH" \
        -t paired \
        --anchorLength 1 \
        --libType fr-firststrand \
        --novelSS \
        --readLength 100 \
        --variable-read-length \
        --allow-clipping \
        --statoff \
        --task prep \
        --od "$JOB_OD" \
        --tmp "$JOB_TMP"
