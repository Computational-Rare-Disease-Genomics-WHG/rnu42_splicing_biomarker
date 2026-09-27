# scripts/gel/ — GEL research workspace scripts

Everything here runs on the GEL research workspace only (never locally — see
project CLAUDE.md). Flat folder by design; this file groups scripts by role
and notes what has to run before what. Run from workspace root.

**End-to-end orchestrator:** `bash scripts/gel/run_all.sh [phase ...]` runs
the scripts below in dependency order, idempotently (skips a step whose
expected output already exists). Phases: `data | signature | figures | null
| export`. Default `all` runs data → signature → export → figures (dependency
order); pass any subset to run only those.

## 1. rMATS batch orchestration (6k whole-blood control cohort)

Run in this order:

1. **`get_control_manifest.py`** (wrapper: `run_get_control_manifest.sh`) —
   builds `data/control_manifest.csv` (+ `data/qcfail_manifest.csv` in
   `--mode fail`). Must run before batch submission — `submit_rmats_batched.sh`
   derives its BAM list from this manifest.
2. **`submit_rmats_batched.sh`** — submits rMATS prep+post as batched LSF jobs
   (100 BAMs/batch). Internally drives:
   - `rmats_prep_array.sh` (one job per BAM array index — never submit directly)
   - `rmats_post.sh` (merges prep outputs per batch — never submit directly)
3. **`combine_rmats_batches_cand_pool.R`** (see §2) — run directly (no LSF
   wrapper) after all `rmats_post_b*` jobs finish.

**Exception-handling scripts** (only needed if a batch fails, not part of the
normal path):
- `cleanup_dup_rmats_prep.sh` — removes duplicate `.rmats`/read_outcomes files
  in a prep tmp dir after a rerun.
- `rerun_failed_qcfail_batches.sh` — identifies and resubmits failed prep
  indices for QC-fail continuation batches.
- `resubmit_post_rerun_batches.sh` — resubmits rMATS post after prep dedup.

## 2. Data processing (PSI matrices, salmon quant, candidate pool)

- **`build_psi_matrix.R`** — builds `outputs/psi_matrix.rds` from the GEL
  11-case + 33-control rMATS run. Foundational input to almost everything
  below and to `figure1.R`/`figure2.R`/the LOO scripts.
- **`build_nava_psi_matrix.R`** — same, mirrored for Nava, → `outputs/nava_psi_matrix.rds`.
- **`combine_rmats_batches_cand_pool.R`** — extracts PSI+IJC+SJC at the full
  483-event A5SS candidate pool from the 6k-cohort rMATS batches →
  `outputs/gel_controls_psi_cand_pool.rds`. Requires §1 (batches) to be done.
  This is the current input for `score_gel_controls_insT.R`/`_stemIII.R` and
  `figure_s6.R`.
- **`make_sample_genotype.R`** — builds `data/sample_genotype.csv` (sample_id
  → HGVS) by matching the rMATS `b1` BAM order to participant genotype info.
  Needed before the LOO scripts (§3) and any HGVS-based carrier classification.
- **`build_salmon_quant_paths.sh`** — builds `data/salmon_quant_paths.csv`
  (platekey → salmon quant path) for the 44-sample dev cohort, from local
  BAM path lists. Needed by `figure2.R` Panel F (STX16 TPM) and
  `figure3.R` Panel C (CINP TPM) — its `platekey` column uses rMATS
  sample IDs (`gel:case:N`/`gel:control:N`), which is what those scripts
  filter on.
## 3. Signature creation (LOO fitting + scoring)

Run after §2's PSI matrices, `sample_genotype.csv`, and the locked split CSV
(`data/split_dev_holdout_v4_mixed.csv`) exist:

- **`loo_evaluate_insT_lasso.R`** → `outputs/loo_results_insT_lasso.rds`. Canonical.
- **`loo_evaluate_stemIII_lasso.R`** → `outputs/loo_results_stemIII_lasso.rds`. Canonical.

Then, after §2's `combine_rmats_batches_cand_pool.R` output exists:

- **`score_gel_controls_insT.R`** → scores the ~6k whole-blood controls
  against the insT signature (`outputs/gel_controls_insT_robustness.rds`).
- **`score_gel_controls_stemIII.R`** → same, for stemIII.

## 4. Paper figures

Each pulls together outputs from §2, §3, and §6. See each script's own header
for its exact input list; summary:

- **`figure1.R`** — needs `psi_matrix.rds`/`nava_psi_matrix.rds` (§2),
  `sample_genotype.csv`/`nava_severity_labels.csv`/split CSV (§2). Panels A-D
  only (variant carrier bar plot, sig-event Venn, PCA, insT dlogit scatter).
- **`figure2.R`** — needs `loo_results_insT_lasso.rds` (§3),
  `gel_controls_insT_robustness.rds` (§3), both PSI matrices (§2),
  `gel_controls_psi_cand_pool.rds` (§2), `a5ss_candidate_full_table.csv` (§6),
  `salmon_quant_paths.csv` (§2), `control_manifest.csv` (§1). Panels A-C
  (insT LASSO coef bar / P(case) by cohort / unmatched-control violin, moved
  from old `figure1.R` E-G) plus D-F (P(case) by variant group / event
  annotation table / STX16 TPM).
- **`figure3.R`** — needs both LOO results (§3), both PSI matrices (§2),
  `gel_controls_psi_cand_pool.rds` (§2), `a5ss_candidate_full_table.csv`
  (§6), `salmon_quant_paths.csv` (§2). Panels A-C (Stem III classifier +
  event table + CINP TPM) plus D (insT + stemIII heatmap).
- **`figure_s3.R`** — GEL delta-logit-PSI scatter (whole-blood vs lymphocytes)
  for the sig-in-both event set. Needs both PSI matrices (§2),
  `sample_genotype.csv` (§2).
- **`figure_s4.R`** — NOT a plot (Figure S4 itself is built manually in
  Illustrator). Needs `figure1_summary.rds` (produced by `figure1.R`),
  both LOO results (§3). Computes and prints/saves the sites-filtering
  funnel counts (all quantified events down to the insT and Stem III K=6
  signatures) as a single reproducible record for whatever goes into that
  figure. All counts read from already-computed outputs, nothing
  re-derived from raw PSI matrices.
- **`figure_s5.R`** — log2(TPM+1) boxplots (WB control vs carrier) for
  insT/stemIII signature genes not already shown in
  `figure2.R`/`figure3.R`'s own TPM panels. Needs both LOO results (§3),
  `psi_matrix.rds` (§2), `a5ss_candidate_full_table.csv` (§6),
  `salmon_quant_paths.csv` (§2).
- **`figure_s6.R`** — validation and specificity of the ReNU splicing
  signature in the ~5,984 whole-blood control cohort. Needs
  `loo_results_insT_lasso.rds` (§3), `gel_controls_psi_cand_pool.rds` (§2),
  `psi_matrix.rds` (§2).
- **`figure_s7.R`** — sequence-motif panels for the signature events.
  Needs both LOO results (§3), `data/rmats/A5SS.MATS.JC.txt`.

## 5. Robustness experiment — null LASSO within batch

See paper methods for the design. Runs in three phases: draw (Desktop, pixi)
→ rMATS statistical test (HPC, LSF array) → LASSO scoring + aggregate
(Desktop, pixi).

- **`check_batch_clean_availability.R`** — per-batch count of samples
  surviving the clean-unmatched filter (6k minus dev-overlap ∪ biallelic ∪
  HD case). Reports which batches have ≥ N clean samples for a given draw
  size. Run first to confirm feasibility.
- **`null_lasso_within_batch.R`** — single-draw driver. Env vars `BATCH_ID`,
  `DRAW_SEED`. Phase 1 draws 11 pseudo-cases + 33 pseudo-controls within a
  batch and prints a bsub for phase 2; phase 3 re-runs the same script to
  score once rMATS output is present.
- **`null_lasso_permutations.R`** — batch-of-draws driver (stratified across
  qualifying batches). Env vars `PER_BATCH`, `N_CASES`, `N_CTRLS`,
  `N_HOLDOUT`, `SAMPLE_TAG`. Phase 1 writes a permutation manifest +
  drawn_samples.csv per perm, prints an LSF array bsub; phase 3 scores each
  perm and aggregates AUCs into
  `outputs/null_within_batch/<SAMPLE_TAG>/<AGG_TAG>/null_auc_distribution.csv`.
- **`null_lasso_rmats_stat.sh`** — HPC worker for one draw: runs
  `prepare_stat_inputs.py` + `rmats.py --task stat` inside the rMATS
  singularity container to produce case-vs-control LRT p-values identical
  to a native rMATS run. Truncates non-A5SS event tables to header-only
  before `--task stat` to keep runtime in seconds.
- **`null_lasso_rmats_stat_array.sh`** — LSF-array wrapper: reads
  `permutation_manifest.csv`, resolves its `LSB_JOBINDEX` row, delegates to
  `null_lasso_rmats_stat.sh`. Idempotent (skips a perm whose A5SS output
  already exists).
- **`summarise_null_auc.R`** — reads a `null_auc_distribution.csv`, prints
  the K-selected/chosen-lambda tables, per-group (1se_fit / min_fit /
  intercept_only) AUC summaries, joint AUC=1 counts, and writes
  `null_auc_distribution_histogram.pdf` alongside the CSV.

## 6. Airlock exports

Run after §3 (need LOO results) unless noted. Note: `export_a5ss_candidate_full_table.R`
must actually run before §4 figures — its output `a5ss_candidate_full_table.csv`
feeds `figure2.R`/`figure3.R`/`figure_s5.R`. The other exports are airlock-
only and can run at any point after §3.

- **`export_insT_signature.R`** / **`export_stemIII_signature.R`** — the
  canonical signature exports (event coords, direction, beta) →
  `data/splicing_signature/insT_stemIII_export/`. GEL-protected until airlocked.
- **`export_rmats_summary_allTypes.R`** — strips the dev-cohort rMATS output
  (all 5 event types) down to group-level stats only, for airlock.
- **`export_a5ss_candidate_full_table.R`** — full ~105-event shared
  (sig-in-both) A5SS pool, concordant and discordant, with GTF protein
  consequence/transcript tier/NDD panel annotation → `outputs/a5ss_candidate_full_table.csv`.
  Needs §3 (both LOO results, for signature membership flags) to exist first.
  Feeds `figure2.R`/`figure3.R`'s event-schematic and CINP/STX16 TPM panels.
- **`export_rnaseq_participants.py`** — participant ID / cohort-role
  exports for airlock requests. Writes `outputs/v4_signature_participants.csv`
  (structured, one row per platekey) and `outputs/airlock_participant_ids.txt`
  (de-duplicated participant_ids, one per line — the airlock payload).
  Overlay cohorts opt-in via `--include-biallelic --include-hd --include-6k`;
  pass all three for the full Figure 1-2 + S3/S6 export scope. Independent
  of the signature LOO results.
