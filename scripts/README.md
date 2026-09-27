# scripts/ — local (non-GEL) scripts

Everything here runs locally (not on the GEL research workspace — see
`scripts/gel/README.md` for that side). Flat folder by design; this file
groups scripts by role and notes what has to run before what. Run from
project root.

## 1. iPSC signature scoring

- **`score_ipsc_against_insT_signature.R`** — applies the insT-only LASSO
  signature (K=6 A5SS events) to the iPSC rMATS output (`data/splicing/div{1,5,9,14}/A5SS.MATS.JC.txt`)
  → `outputs/ipsc_insT_signature_scores.rds`. Needs
  `data/splicing_signature/insT_stemIII_export/insT_signature_export.csv`
  (+ `_meta.txt`), airlocked from GEL.
- **`score_ipsc_against_stemIII_signature.R`** — same, for the stemIII
  signature → `outputs/ipsc_stemIII_signature_scores.rds`.

## 2. Paper figures (iPSC)

- **`export_table_s4_ipsc.R`** — two-in-one export.
  (1) Joins Table S1 (`outputs/candidate_table/a5ss_candidate_full_table.csv`,
  GEL airlock export) against the per-DIV iPSC rMATS
  (`data/splicing/div{1,5,9,14}/A5SS.MATS.JC.txt`): per-event, per-DIV
  quantified flag, FDR/PValue, and delta logit-PSI (ReNU − KOLF) →
  `outputs/candidate_table/table1_ipsc_annotated.csv` (consumed by
  `figure4.R` and `figure_s8.R`).
  (2) Pivots the concordant subset into Supplementary Table S4 (long
  format, one row per event × DIV) with novel/annotated splice-site
  classification from the gencode v49 GTF cache; column names/types match
  Table S1 for a direct `event_id` join →
  `outputs/candidate_table/table_s4_ipsc.csv`.
  Rerun if Table S1 or the iPSC rMATS data changes.
- **`figure4.R`** — Figure 4: cross-tissue validation of the 96-event
  (direction-concordant) A5SS candidate pool in iPSC neuronal
  differentiation. Panel A: exact hypergeometric test — is the pool's
  FDR<0.05 rate in each iPSC DIV higher than a same-n random draw from the
  3-way GEL WB / Nava / iPSC-quantified universe (`phyper`/`qhyper`, own
  rMATS-loading code, not read from `figure_s8.R`). Panel B: candidate-pool
  delta-logit-PSI scatter (WB+lymphocyte vs iPSC, per DIV), hollow/filled by
  iPSC FDR. Panel C: same scatter restricted to the 6 insT signature events,
  gene-labelled at a fixed cross-DIV position with a leader line to each
  DIV's actual point. Panel D: insT LASSO P(case) by cell line/DIV (points +
  group median). Needs §1's `ipsc_insT_signature_scores.rds`, Nava + GEL
  rMATS summary tables, iPSC per-DIV rMATS tables, and
  `table1_ipsc_annotated.csv` (above).
- **`figure_s8.R`** — Supplementary
  Figure S8: Panel A, genome-wide FDR<0.05 iPSC event counts by rMATS type
  (shared+concordant/discordant/singleton across DIVs). Panel B, per-type
  Fisher OR for the shared+concordant fraction. Panel C, the 17/96 candidate-pool events NOT quantified
  in every DIV: gene name / annotated-novel-junction marker / per-DIV
  Reads+TPM cell pair (two independent colour scales via `ggnewscale`). Also
  writes `outputs/candidate_table/table1_missingness_summary.csv`. Needs
  `table1_ipsc_annotated.csv` (above) and the gencode v49 GTF/cache under
  `data/reference/`; otherwise self-contained (own copies of the
  rMATS-loading helpers, no shared state with `figure4.R`).
- **`export_supplementary_tables_xlsx.R`** — bundles the four paper-ready
  `outputs/supplementary_tables/Table_S{1..4}_*.csv` into one workbook
  (`supplementary_tables.xlsx`, one sheet per table). Run after touching any
  of the four CSVs — see `outputs/supplementary_tables/README.md`.
- **`figure_s1.R`** — Supplementary Figure S1: iPSC cell-identity marker
  expression, log2(TPM+1) averaged across the three replicates within each
  (cell line, DIV) → two stacked rows per gene (KOLF top / ReNU bottom),
  DIVs across columns. Six gene groups (pluripotency-loss / neuronal /
  cortical / glutamatergic / GABAergic negative-control / synaptic-
  maturation) laid out via patchwork. Needs
  `data/gene_expression/gene_tpm_long.rds` only — independent of §1.
- **`figure_s2.R`** — Supplementary Figure S2: iPSC RNA-seq QC (per-sample MultiQC metrics,
  DESeq2 PCA) supporting Figure 4. Needs `data/qc/general_stats_table.tsv`
  and `data/qc/star_salmon_deseq2_pca-plot.tsv` only — independent of §1/§2.
- **`figure_s4_flowchart.R`** — Supplementary Figure S4: flowchart of the
  insT / stemIII development pipeline (candidate-pool filtering funnel).
  Local, fast-iteration draft; the airlocked counts come from
  `scripts/gel/figure_s4.R`. Independent of §1/§2.

## `gel/` subfolder

Everything under `scripts/gel/` runs on the GEL research workspace only —
see `scripts/gel/README.md`.
