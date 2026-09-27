# RNU4-2 Splicing Biomarker

Analysis pipeline for a treatment-response splicing biomarker for RNU4-2-related
neurodevelopmental syndrome (ReNU syndrome). The biomarker is a LASSO-based
case-vs-control classifier (P(case)) over a panel of splice events, designed
for monitoring antisense oligonucleotide (ASO) therapy response.

## Datasets

| Dataset | Tissue | Samples | Purpose |
|---------|--------|---------|---------|
| iPSC neurons | iPSC-derived neurons | 24 (3 KOLF + 3 ReNU × 4 DIV stages) | Biomarker application across neuronal differentiation |
| GEL whole blood | Whole blood | 11 ReNU cases + ~6k controls | Biomarker development (LOO), robustness testing |
| Nava lymphocytes | Lymphocyte culture | 19 cases + 21 controls | Co-development cohort (joint LOO with GEL) |

## Structure

```
scripts/                        Local (non-GEL) scripts — see scripts/README.md for
                                 full script-by-script breakdown (iPSC signature
                                 scoring, paper figures)
  gel/                          GEL Research Environment scripts (run on GEL only)
                                 — see scripts/gel/README.md

03-analysis/                    rMATS run configs (iPSC)

data/                           Local data (gitignored)
outputs/                        Analysis outputs (gitignored)
```

## Splicing signature

Two genotype-restricted LASSO signatures, each selected from a shared 483-event
A5SS candidate pool (FDR<0.05 in both GEL and Nava, direction-concordant) by
`cv.glmnet` (alpha=1, lower.limits=0, λ.1se) on a case-LOO fold structure:

- **insT** (K=6) — fit on carriers of the recurrent founder variant (n.64_65insT,
  also the iPSC ReNU line genotype). Uses the locked 70/30 dev/holdout split;
  holdout stays blind.
- **stemIII** (K=6) — fit on carriers in the Stem III structural region
  (HGVS nt 75-79). Cohort too small for a holdout split; case-LOO linpred AUC
  is the primary metric.

Each signature produces a single score, **LASSO P(case)** — a case-vs-control
classifier:
```
P(case) = sigmoid(intercept + Sigma_j  beta_j * direction_j * (bb_logit_PSI_sample - mean_bb_logit_PSI_control))
```

## Pipeline dependencies

**Alignment:** nf-core/rnaseq 3.22.2 (STAR + Salmon), run via Seqera Platform

**Splicing:** rMATS-turbo v4.3.0, run via Seqera/Nextflow

**R packages:** tidyverse, patchwork, ggsci, ggrepel, pheatmap, RColorBrewer,
UpSetR, rtracklayer, DESeq2

**GEL workspace** (`scripts/gel/`, never run locally): pixi env (`pixi.toml`,
workspace root) — r-ggplot2, r-patchwork, r-biocmanager, bioconductor-rsamtools,
etc. See CLAUDE.md "GEL pixi" for setup/gotchas.
