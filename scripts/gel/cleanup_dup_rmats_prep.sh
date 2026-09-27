#!/bin/bash
# cleanup_dup_rmats_prep.sh
# For the rerun QC-fail batches: remove older .rmats and read_outcomes files
# when a prep tmp dir has >1 from a previous run.
# Then verifies each batch has 100 prep dirs each with exactly 1 .rmats file.
#
# Note: only .rmats duplicates cause rMATS post errors; read_outcomes duplicates
# are harmless but cleaned for tidiness.
#
# Usage (dry run first!):
#   bash scripts/gel/cleanup_dup_rmats_prep.sh --dry-run
#   bash scripts/gel/cleanup_dup_rmats_prep.sh

set -euo pipefail

RERUN_BATCHES=(62 63 64 65 66 67 68 70 71 72 73 74 75 76 77 79 80)
DRY_RUN=false
BASE="rmats_tmp"

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=true ;;
    *) BASE="$arg" ;;
  esac
done

if [[ "$DRY_RUN" == "true" ]]; then
  echo "[DRY-RUN] No files will be deleted"
else
  echo "WILL DELETE older prep files in rerun batches"
fi
echo ""

# --- cleanup -----------------------------------------------------------------

ndirs=0
ndeleted=0

dedup_pattern() {
  local dir="$1" pattern="$2"
  mapfile -t files < <(
    find "$dir" -maxdepth 1 -name "$pattern" -printf "%T@ %p\n" \
      | sort -rn | awk '{print $2}'
  )
  local n=${#files[@]}
  if (( n > 1 )); then
    echo "    [$pattern] keep:   ${files[0]##*/}"
    for f in "${files[@]:1}"; do
      echo "    [$pattern] remove: ${f##*/}"
      [[ "$DRY_RUN" == "false" ]] && rm -f "$f"
      ndeleted=$(( ndeleted + 1 ))
    done
    return 0   # had duplicates
  fi
  return 1     # no duplicates
}

for batch in "${RERUN_BATCHES[@]}"; do
  batch_dir="$BASE/batch_${batch}"
  if [[ ! -d "$batch_dir" ]]; then
    echo "WARNING: $batch_dir does not exist — skipping"
    continue
  fi
  while IFS= read -r dir; do
    had_dup=false
    dedup_pattern "$dir" "*.rmats"         && had_dup=true
    dedup_pattern "$dir" "*read_outcomes*" && had_dup=true
    if [[ "$had_dup" == "true" ]]; then
      echo "  ^ $dir"
      ndirs=$(( ndirs + 1 ))
    fi
  done < <(find "$batch_dir" -mindepth 1 -maxdepth 1 -type d ! -name "*_od" ! -name "merged" | sort -V)
done

echo ""
echo "Cleanup complete."
echo "  Dirs with duplicates: $ndirs"
if [[ "$DRY_RUN" == "true" ]]; then
  echo "  Files would be removed: $ndeleted"
else
  echo "  Files removed: $ndeleted"
fi

# --- verification ------------------------------------------------------------

echo ""
echo "=========================================="
echo "Verification (expecting 100 prep dirs / batch, each with 1 .rmats)"
echo "=========================================="

all_ok=true

for batch in "${RERUN_BATCHES[@]}"; do
  batch_dir="$BASE/batch_${batch}"
  if [[ ! -d "$batch_dir" ]]; then
    echo "MISSING  batch_${batch}: directory not found"
    all_ok=false
    continue
  fi

  # expected = number of BAMs in the batch list (last batch may be < 100)
  bam_list="${WORKSPACE:-/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature}/data/bam_paths_batch_${batch}.txt"
  if [[ -f "$bam_list" ]]; then
    expected=$(wc -l < "$bam_list")
  else
    expected=100
  fi

  mapfile -t prep_dirs < <(
    find "$batch_dir" -mindepth 1 -maxdepth 1 -type d ! -name "*_od" ! -name "merged" | sort -V
  )
  total=${#prep_dirs[@]}

  n_ok=0; n_zero=0; n_multi=0
  zero_list=(); multi_list=()

  for dir in "${prep_dirs[@]}"; do
    mapfile -t rf < <(find "$dir" -maxdepth 1 -name "*.rmats")
    cnt=${#rf[@]}
    if   (( cnt == 1 )); then n_ok=$(( n_ok + 1 ))
    elif (( cnt == 0 )); then n_zero=$(( n_zero + 1 )); zero_list+=("${dir##*/}")
    else                      n_multi=$(( n_multi + 1 )); multi_list+=("${dir##*/}")
    fi
  done

  if (( total == expected && n_ok == expected )); then
    echo "OK       batch_${batch}: ${expected}/${expected} prep dirs have 1 .rmats"
  else
    echo "PROBLEM  batch_${batch}: ${total}/${expected} prep dirs total | ok=${n_ok} | missing_rmats=${n_zero} | still_dup=${n_multi}"
    [[ ${#zero_list[@]}  -gt 0 ]] && echo "         missing .rmats in indices: ${zero_list[*]}"
    [[ ${#multi_list[@]} -gt 0 ]] && echo "         still >1 .rmats in indices: ${multi_list[*]}"
    all_ok=false
  fi
done

echo ""
if [[ "$all_ok" == "true" ]]; then
  echo "All batches OK — safe to re-submit post jobs."
else
  echo "Some batches need attention (see PROBLEM lines above)."
fi
