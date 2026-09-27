#!/usr/bin/env bash
#
# null_lasso_rmats_stat_array.sh  (HPC, submitted as LSF array)
# One array element = one permutation. Reads the perm's batch_id +
# out_tag from the manifest by LSB_JOBINDEX and delegates to
# null_lasso_rmats_stat.sh.
#
# Args (positional; env is wiped by -L /bin/bash on GEL LSF):
#   $1  MANIFEST  path to permutation_manifest.csv (columns:
#                 perm_id,batch_id,within_batch_seed,out_tag)
#
# Submit example (printed by null_lasso_permutations.R phase 1):
#   bsub -q short -P re_gecip_machine_learning \
#     -J "null_stratified_p17[1-1020]%50" \
#     -o hpc_logs/null_stratified_p17_%J_%I.out \
#     -e hpc_logs/null_stratified_p17_%J_%I.err \
#     -cwd <WORKSPACE> -L /bin/bash \
#     bash scripts/gel/null_lasso_rmats_stat_array.sh \
#          <WORKSPACE>/outputs/null_within_batch/stratified_p17/permutation_manifest.csv

set -euo pipefail

WORKSPACE="/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"

if [[ $# -lt 1 ]]; then
    echo "usage: $0 MANIFEST_CSV" >&2
    exit 1
fi
MANIFEST="$1"

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: manifest not found at $MANIFEST" >&2
    exit 1
fi
if [[ -z "${LSB_JOBINDEX:-}" ]]; then
    echo "ERROR: LSB_JOBINDEX not set -- must be run as an LSF array job." >&2
    exit 1
fi

# Read manifest row matching perm_id == LSB_JOBINDEX.
# CSV columns: perm_id,batch_id,within_batch_seed,out_tag
# R writes character cols with double quotes; awk strips them.
row=$(awk -F, -v pid="$LSB_JOBINDEX" '
    NR==1 { next }
    { gsub(/"/, ""); if ($1 == pid) { print; exit } }
' "$MANIFEST")

if [[ -z "$row" ]]; then
    echo "ERROR: no manifest row for perm_id=$LSB_JOBINDEX" >&2
    exit 1
fi

IFS=, read -r PERM_ID BATCH_ID DRAW_SEED OUT_TAG <<< "$row"

echo "=== array element ==="
echo "PERM_ID    : $PERM_ID"
echo "BATCH_ID   : $BATCH_ID"
echo "DRAW_SEED  : $DRAW_SEED"
echo "OUT_TAG    : $OUT_TAG"

A5_OUT="${WORKSPACE}/outputs/null_within_batch/${OUT_TAG}/rmats_stat_od/A5SS.MATS.JC.txt"
if [[ -f "$A5_OUT" ]]; then
    echo "A5SS.MATS.JC.txt already present -- skipping."
    exit 0
fi

# -l = login shell, so /etc/profile is sourced and `module` is defined
# for the `module load singularity/4.1.1` line inside the inner script.
bash -l "${WORKSPACE}/scripts/gel/null_lasso_rmats_stat.sh" \
    "$BATCH_ID" "$DRAW_SEED" "$OUT_TAG"
