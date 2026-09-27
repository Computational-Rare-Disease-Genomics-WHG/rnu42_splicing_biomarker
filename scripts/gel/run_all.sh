#!/usr/bin/env bash
# =============================================================================
# run_all.sh  (GEL workspace)
# Orchestrator: runs the scripts in this folder in dependency order.
# Grouped by phase; each phase is idempotent (skips a step whose expected
# output file already exists).
#
# Usage:
#   bash scripts/gel/run_all.sh [phase ...]
#
# phase = one of:
#   data       PSI matrices, per-sample genotype, salmon paths, batch pool
#   signature  LOO fits + score 6k WB controls
#   export     airlock exports (signature CSVs, table S1, summary tables)
#   figures    all paper figure scripts
#   null       null-LASSO robustness driver (phase 1 draw; array must be
#              submitted on HPC separately, see the printed bsub)
#   all        every phase above, in order (default)
#
# Prerequisites (checked, not run by this script):
#   - rMATS batches for the 6k WB controls have finished (see §1 of
#     scripts/gel/README.md; submit_rmats_batched.sh spawns LSF jobs that
#     take days).
#   - Pixi env active (`pixi shell`) or scripts invoked via `pixi run ...`.
#     R scripts here need pixi's r-base + package set.
#
# Notes:
#   - Scripts marked WORKSPACE-only run under the pixi env on GEL Desktop.
#   - Scripts marked HPC-only are submitted via bsub inside their own
#     wrappers and this driver does not block on them.
#   - To rerun a step whose output already exists, delete that output first.
# =============================================================================

set -euo pipefail

WORKSPACE="/re_gecip/machine_learning/CRDG_Shared/renu_splicing_signature"
cd "$WORKSPACE"

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
say() { printf "\n\033[1;36m=== %s ===\033[0m\n" "$*"; }
step() { printf "  \033[0;36m> %s\033[0m\n" "$*"; }
skip() { printf "  \033[0;33m- skip: %s (output exists)\033[0m\n" "$*"; }

run_if_missing() {
  # run_if_missing <output_path> <human_label> <command...>
  local out="$1"; shift
  local label="$1"; shift
  if [[ -e "$out" ]]; then
    skip "$label"
    return 0
  fi
  step "$label"
  "$@"
}

# --------------------------------------------------------------------------
# Phase: data processing (PSI matrices, sample genotype, salmon paths,
#        candidate-pool matrix from batch outputs)
# --------------------------------------------------------------------------
phase_data() {
  say "PHASE: data processing"

  if [[ ! -d outputs/rmats_controls ]]; then
    echo "ERROR: outputs/rmats_controls/ is missing." >&2
    echo "  Run §1 (submit_rmats_batched.sh) first -- see scripts/gel/README.md." >&2
    exit 1
  fi

  run_if_missing "data/control_manifest.csv" \
    "control manifest" \
    bash scripts/gel/run_get_control_manifest.sh

  run_if_missing "outputs/psi_matrix.rds" \
    "GEL dev-cohort PSI matrix" \
    pixi run Rscript scripts/gel/build_psi_matrix.R

  run_if_missing "outputs/nava_psi_matrix.rds" \
    "Nava PSI matrix" \
    pixi run Rscript scripts/gel/build_nava_psi_matrix.R

  run_if_missing "data/sample_genotype.csv" \
    "sample <-> HGVS genotype table" \
    pixi run Rscript scripts/gel/make_sample_genotype.R

  run_if_missing "data/salmon_quant_paths.csv" \
    "salmon quant paths for dev cohort" \
    bash scripts/gel/build_salmon_quant_paths.sh

  run_if_missing "outputs/gel_controls_psi_cand_pool.rds" \
    "6k WB candidate-pool PSI matrix" \
    pixi run Rscript scripts/gel/combine_rmats_batches_cand_pool.R
}

# --------------------------------------------------------------------------
# Phase: signature creation (LOO fits + score 6k controls)
# --------------------------------------------------------------------------
phase_signature() {
  say "PHASE: signature creation"

  run_if_missing "outputs/loo_results_insT_lasso.rds" \
    "insT LASSO LOO fit" \
    pixi run Rscript scripts/gel/loo_evaluate_insT_lasso.R

  run_if_missing "outputs/loo_results_stemIII_lasso.rds" \
    "stemIII LASSO LOO fit" \
    pixi run Rscript scripts/gel/loo_evaluate_stemIII_lasso.R

  run_if_missing "outputs/gel_controls_insT_robustness.rds" \
    "score 6k WB controls against insT signature" \
    pixi run Rscript scripts/gel/score_gel_controls_insT.R

  run_if_missing "outputs/gel_controls_stemIII_robustness.rds" \
    "score 6k WB controls against stemIII signature" \
    pixi run Rscript scripts/gel/score_gel_controls_stemIII.R
}

# --------------------------------------------------------------------------
# Phase: airlock exports
# --------------------------------------------------------------------------
phase_export() {
  say "PHASE: airlock exports"

  run_if_missing "data/splicing_signature/insT_stemIII_export/insT_signature_export.csv" \
    "insT signature export" \
    pixi run Rscript scripts/gel/export_insT_signature.R

  run_if_missing "data/splicing_signature/insT_stemIII_export/stemIII_signature_export.csv" \
    "stemIII signature export" \
    pixi run Rscript scripts/gel/export_stemIII_signature.R

  run_if_missing "outputs/rmats_export/A5SS_summary.csv" \
    "dev-cohort rMATS summary (all event types)" \
    pixi run Rscript scripts/gel/export_rmats_summary_allTypes.R

  run_if_missing "outputs/a5ss_candidate_full_table.csv" \
    "A5SS 105-event full annotated table (Table S1)" \
    pixi run Rscript scripts/gel/export_a5ss_candidate_full_table.R

  # Airlock scope = every overlay cohort included (dev + biallelic + HD + 6k).
  # Emits both v4_signature_participants.csv and airlock_participant_ids.txt.
  run_if_missing "outputs/airlock_participant_ids.txt" \
    "participant IDs + airlock ID list" \
    python scripts/gel/export_rnaseq_participants.py \
        --include-biallelic --include-hd --include-6k
}

# --------------------------------------------------------------------------
# Phase: paper figures
# --------------------------------------------------------------------------
phase_figures() {
  say "PHASE: paper figures"

  run_if_missing "outputs/figures/figure1.pdf" \
    "Figure 1"       pixi run Rscript scripts/gel/figure1.R
  run_if_missing "outputs/figures/figure2.pdf" \
    "Figure 2"       pixi run Rscript scripts/gel/figure2.R
  run_if_missing "outputs/figures/figure3.pdf" \
    "Figure 3"       pixi run Rscript scripts/gel/figure3.R
  run_if_missing "outputs/figures/figure_s3.pdf" \
    "Figure S3"      pixi run Rscript scripts/gel/figure_s3.R
  run_if_missing "outputs/figure_s4_summary.rds" \
    "Figure S4 funnel-count summary (no plot; figure built manually)" \
    pixi run Rscript scripts/gel/figure_s4.R
  run_if_missing "outputs/figures/figure_s5.pdf" \
    "Figure S5"      pixi run Rscript scripts/gel/figure_s5.R
  run_if_missing "outputs/figures/figure_s6.pdf" \
    "Figure S6"      pixi run Rscript scripts/gel/figure_s6.R
  run_if_missing "outputs/figures/figure_s7.pdf" \
    "Figure S7"      pixi run Rscript scripts/gel/figure_s7.R
}

# --------------------------------------------------------------------------
# Phase: null-LASSO robustness experiment
# --------------------------------------------------------------------------
# Only the Desktop-side draw + aggregate phases are run here. The rMATS
# stat step is an LSF array that must be submitted separately -- the
# permutations driver prints the exact bsub command when phase 1 finishes.
phase_null() {
  say "PHASE: null-LASSO robustness (Desktop-side only; HPC array is manual)"

  step "batch-availability check"
  pixi run Rscript scripts/gel/check_batch_clean_availability.R

  step "permutation driver (phase 1 draw OR phase 3 aggregate, depending on state)"
  pixi run Rscript scripts/gel/null_lasso_permutations.R

  echo "  If the driver just wrote drawn_samples.csv files and printed a"
  echo "  bsub command, submit it on HPC now, wait for completion, then"
  echo "  rerun this phase to aggregate."
}

# --------------------------------------------------------------------------
# Dispatch
# --------------------------------------------------------------------------
PHASES=("$@")
if [[ ${#PHASES[@]} -eq 0 ]] || [[ "${PHASES[0]}" == "all" ]]; then
  PHASES=(data signature export figures)
fi

for p in "${PHASES[@]}"; do
  case "$p" in
    data)      phase_data      ;;
    signature) phase_signature ;;
    export)    phase_export    ;;
    figures)   phase_figures   ;;
    null)      phase_null      ;;
    *) echo "unknown phase: $p" >&2; exit 1 ;;
  esac
done

say "run_all.sh complete"
