#!/usr/bin/env python3
"""
get_control_manifest.py
Query GEL LabKey to build a QC-filtered, spliceopathy-annotated manifest of
whole-blood RNA-seq BAMs for robustness testing of the splicing signature
(Phase 2 / Q3 stress test).

Modes (--mode):
  pass (default)  Samples meeting all QC thresholds → data/control_manifest.csv.
  fail            Samples failing at least one QC threshold →
                  data/qcfail_manifest.csv. Used to add rMATS batches for the
                  ~1000 previously-excluded samples so downstream analysis can
                  filter on qc_pass + the rnaseqc_* columns.

Outputs:
  data/control_manifest.csv       — one row per sample; spliceopathy_evidence
                                    is an ordered label (none < tiered < sdd <
                                    solved) giving the strongest source of
                                    evidence that the participant has a
                                    spliceopathy-gene variant. 'solved' =
                                    gmc_exit_questionnaire case_solved_family
                                    IN ('yes','partially'); 'sdd' = a row in
                                    submitted_diagnostic_discovery; 'tiered' =
                                    only a TIER1/TIER2 tiering variant (a
                                    candidate, not a confirmed diagnosis).
                                    spliceopathy_flag is a back-compat boolean
                                    TRUE whenever evidence != 'none'. All
                                    samples retained.
  data/spliceopathy_diagnoses.csv — long format; one row per variant per
                                    source (exit_questionnaire, SDD, tiering)
                                    for any participant with RNA-seq and at
                                    least one hit; source column lets you
                                    filter to solved-only.

Authentication: reads ~/.netrc automatically (machine labkey-embassy.gel.zone).

IMPORTANT: Run only inside the GEL Research Environment.
Outputs must not leave GEL without airlock approval.

Run from workspace root: python -u scripts/get_control_manifest.py
(-u disables output buffering so progress prints appear immediately)
"""

import argparse
import sys
import os
import gzip
import warnings
import labkey
import pandas as pd

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

LABKEY_HOST   = "labkey-embassy.gel.zone"
LABKEY_FOLDER = "/main-programme/main-programme_v19_2024-10-31"
LABKEY_SCHEMA = "lists"

# QC thresholds — aligned with Lord et al. 2026 (medRxiv 10.64898/2026.03.19.26348811)
# RIN is NOT used: unreliable for PAXgene blood RNA; not output by RNASeQC2.
# pct_correct_strand is NOT used: no matching column in rnaseq_qc_metrics.
#
# mapped_unique_reads is DERIVED (not a LabKey column):
#   rnaseqc_total_reads * rnaseqc_mapping_rate * (1 - rnaseqc_duplicate_rate_of_mapped)
MIN_TOTAL_READS   = 120e6   # rnaseqc_total_reads >= 120M (Lord et al.)
MIN_MAPPED_UNIQUE = 30e6    # derived unique mapped reads >= 30M (Lord et al.)
MIN_MAPPING_RATE  = 0.95    # rnaseqc_mapping_rate >= 0.95 (unmapped <= 5%)
MAX_RRNA_RATE     = 0.12    # rnaseqc_rrna_rate <= 12% (Lord et al.)

BATCH_SIZE = 200            # genes per IN clause
BAM_BATCH  = 500            # platekeys per IN clause (large lists need batching)

GENES_FILE    = "data/spliceopathy_genes/pathogenic_variants.vcf.gz"

# Output paths depend on --mode (set in main()).
OUT_PATHS = {
    "pass": ("data/control_manifest.csv",
             "data/spliceopathy_diagnoses.csv"),
    "fail": ("data/qcfail_manifest.csv",
             "data/spliceopathy_diagnoses_qcfail.csv"),
}

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def labkey_to_df(sql_query, maxrows=10_000_000):
    """Version-aware LabKey SQL query → DataFrame (matches GEL RE working pattern)."""
    ver = labkey.__version__
    if ver < "2.4.0":
        server_context = labkey.utils.create_server_context(
            domain         = LABKEY_HOST,
            container_path = LABKEY_FOLDER,
            context_path   = "labkey",
            use_ssl        = True,
        )
        results = labkey.query.execute_sql(
            server_context,
            schema_name = LABKEY_SCHEMA,
            sql         = sql_query,
            max_rows    = maxrows,
        )
    else:
        from labkey.api_wrapper import APIWrapper
        api = APIWrapper(
            LABKEY_HOST, LABKEY_FOLDER, "labkey", use_ssl=True
        )
        results = api.query.execute_sql(
            sql            = sql_query,
            schema_name    = LABKEY_SCHEMA,
            max_rows       = maxrows,
            waf_encode_sql = False,
        )
    if results is None or "rows" not in results:
        return pd.DataFrame()
    return pd.DataFrame(results["rows"])


def execute_sql(sql, label="query"):
    """Wrapper around labkey_to_df with error handling."""
    try:
        return labkey_to_df(sql)
    except Exception as e:
        warnings.warn(f"{label} failed: {e}")
        return pd.DataFrame()


def batched(items, size):
    """Yield successive chunks of `size` from `items`."""
    for i in range(0, len(items), size):
        yield items[i : i + size]


def sql_in(values):
    """Return a SQL IN list string: 'val1','val2',..."""
    escaped = [str(v).replace("'", "''") for v in values]
    return ",".join(f"'{v}'" for v in escaped)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=["pass", "fail"], default="pass",
                        help=("'pass' (default): keep samples meeting all QC "
                              "thresholds → data/control_manifest.csv. "
                              "'fail': keep samples that fail at least one QC "
                              "threshold → data/qcfail_manifest.csv. Use 'fail' "
                              "to add rMATS batches for the previously-excluded "
                              "samples so downstream can include/exclude on "
                              "qc_pass + the rnaseqc_* columns."))
    args = parser.parse_args()
    mode = args.mode
    out_manifest, out_diagnoses = OUT_PATHS[mode]

    print(f"LabKey host:   {LABKEY_HOST}")
    print(f"Folder:        {LABKEY_FOLDER}")
    print(f"Mode:          {mode}\n")

    # ----------------------------------------------------------------------
    # Step 1: RNA-seq QC metrics
    #
    # Table: rnaseq_qc_metrics (confirmed 2026-03-25)
    # All records are whole blood — no tissue_type filter needed.
    # ----------------------------------------------------------------------
    print("Querying rnaseq_qc_metrics...")

    qc_sql = """
        SELECT
            participant_id,
            rnaseq_platekey,
            rnaseqc_total_reads,
            rnaseqc_mapping_rate,
            rnaseqc_duplicate_rate_of_mapped,
            rnaseqc_rrna_rate
        FROM rnaseq_qc_metrics
    """
    qc_raw = execute_sql(qc_sql, "rnaseq_qc_metrics")
    if qc_raw.empty:
        sys.exit("ERROR: rnaseq_qc_metrics returned no rows — check folder path and table name.")
    print(f"  Rows returned: {len(qc_raw)}")

    # Numeric coercion (LabKey may return strings)
    for col in ["rnaseqc_total_reads", "rnaseqc_mapping_rate",
                "rnaseqc_duplicate_rate_of_mapped", "rnaseqc_rrna_rate"]:
        qc_raw[col] = pd.to_numeric(qc_raw[col], errors="coerce")

    # Derive mapped_unique_reads (Lord et al. secondary filter)
    qc_raw["mapped_unique_reads"] = (
        qc_raw["rnaseqc_total_reads"]
        * qc_raw["rnaseqc_mapping_rate"]
        * (1 - qc_raw["rnaseqc_duplicate_rate_of_mapped"])
    )

    # Apply QC filters
    qc_pass_mask = (
        qc_raw["rnaseqc_total_reads"].notna()  & (qc_raw["rnaseqc_total_reads"]  >= MIN_TOTAL_READS)   &
        qc_raw["mapped_unique_reads"].notna()  & (qc_raw["mapped_unique_reads"]  >= MIN_MAPPED_UNIQUE)  &
        qc_raw["rnaseqc_mapping_rate"].notna() & (qc_raw["rnaseqc_mapping_rate"] >= MIN_MAPPING_RATE)   &
        qc_raw["rnaseqc_rrna_rate"].notna()    & (qc_raw["rnaseqc_rrna_rate"]    <= MAX_RRNA_RATE)
    )
    qc_raw["qc_pass"] = qc_pass_mask

    if mode == "pass":
        qc_keep = qc_raw[qc_pass_mask].copy()
        print(f"  After QC filter (pass): {len(qc_keep)} samples\n")
        assert len(qc_keep) > 0, "No samples passed QC filter — check column names"
    else:
        qc_keep = qc_raw[~qc_pass_mask].copy()
        print(f"  After QC filter (fail): {len(qc_keep)} samples\n")
        assert len(qc_keep) > 0, "No samples failed QC filter — check column names"

    # participant_id → rnaseq_platekey lookup (used to enrich diagnoses table)
    pid_to_platekey = qc_keep.set_index("participant_id")["rnaseq_platekey"].to_dict()

    # ----------------------------------------------------------------------
    # Step 2: Spliceopathy gene list from pathogenic_variants.vcf.gz
    # SYMBOL INFO field contains HGNC gene symbols (comma-separated).
    # ----------------------------------------------------------------------
    assert os.path.exists(GENES_FILE), f"Spliceopathy genes file not found: {GENES_FILE}"

    splice_genes_set = set()
    with gzip.open(GENES_FILE, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            info = line.split("\t")[7]
            for field in info.split(";"):
                if field.startswith("SYMBOL="):
                    for symbol in field[7:].split(","):
                        s = symbol.strip()
                        if s:
                            splice_genes_set.add(s)
                    break
    splice_genes = sorted(splice_genes_set)
    print(f"Spliceopathy gene list: {len(splice_genes)} unique genes")

    # ----------------------------------------------------------------------
    # Step 3: Query diagnosis evidence in three layers
    #
    # PRIMARY   gmc_exit_questionnaire        — case_solved_family IN
    #                                           ('yes','partially'); this is
    #                                           the only layer that means
    #                                           "solved" in GEL's own records.
    # SECONDARY submitted_diagnostic_discovery — a variant was submitted as
    #                                           diagnostic but may not yet be
    #                                           confirmed as solving the case.
    # TERTIARY  tiering_data TIER1/TIER2      — pipeline-ranked candidates;
    #                                           most never become diagnoses.
    #
    # gmc_exit_questionnaire columns (v19, confirmed 2026-04-15):
    #   participant_id, case_solved_family, gene_name, chromosome, position,
    #   reference, alternate, acmg_classification, assembly
    #   NB: uses 'assembly', not 'genome_build', and 'gene_name', not
    #   'genomic_feature_hgnc'.
    # ----------------------------------------------------------------------

    # -- 3a: gmc_exit_questionnaire (solved cases) -------------------------
    print(f"Querying gmc_exit_questionnaire (batches of {BATCH_SIZE})...")
    eq_frames = []
    for batch in batched(splice_genes, BATCH_SIZE):
        sql = f"""
            SELECT DISTINCT
                participant_id,
                gene_name,
                case_solved_family,
                chromosome,
                position,
                reference,
                alternate,
                acmg_classification,
                assembly
            FROM gmc_exit_questionnaire
            WHERE gene_name IN ({sql_in(batch)})
              AND case_solved_family IN ('yes','partially')
        """
        df = execute_sql(sql, "gmc_exit_questionnaire")
        if not df.empty:
            eq_frames.append(df)

    if eq_frames:
        eq_hits = pd.concat(eq_frames, ignore_index=True).drop_duplicates()
        print(f"  Exit-questionnaire solved hits: {len(eq_hits)} participant-variant rows")
    else:
        warnings.warn("gmc_exit_questionnaire returned no results — check table/column names")
        eq_hits = pd.DataFrame()

    # -- 3b: tiering_data TIER1/TIER2 (candidates, not confirmed diagnoses)
    print(f"Querying tiering_data TIER1/TIER2 (batches of {BATCH_SIZE})...")
    tier12_frames = []
    for batch in batched(splice_genes, BATCH_SIZE):
        sql = f"""
            SELECT DISTINCT
                participant_id,
                genomic_feature_hgnc,
                tier,
                chromosome,
                position,
                reference,
                alternate,
                assembly,
                genotype,
                mode_of_inheritance,
                phenotype
            FROM tiering_data
            WHERE genomic_feature_hgnc IN ({sql_in(batch)})
              AND tier IN ('TIER1','TIER2')
        """
        df = execute_sql(sql, "tiering_data TIER1/2")
        if not df.empty:
            tier12_frames.append(df)

    if tier12_frames:
        tier12_hits = pd.concat(tier12_frames, ignore_index=True).drop_duplicates()
        print(f"  Tier 1/2 hits: {len(tier12_hits)} participant-variant rows")
    else:
        warnings.warn("tiering_data TIER1/TIER2 returned no results — check table/column names")
        tier12_hits = pd.DataFrame()

    # -- 3c: submitted_diagnostic_discovery --------------------------------
    print("Querying submitted_diagnostic_discovery...")
    sdd_frames = []
    for batch in batched(splice_genes, BATCH_SIZE):
        sql = f"""
            SELECT DISTINCT
                participant_id,
                gene_of_interest,
                chromosome,
                position,
                reference,
                alternate,
                genome_build,
                variant_of_interest,
                zygosity
            FROM submitted_diagnostic_discovery
            WHERE gene_of_interest IN ({sql_in(batch)})
        """
        df = execute_sql(sql, "submitted_diagnostic_discovery")
        if not df.empty:
            sdd_frames.append(df)

    if sdd_frames:
        sdd_hits = pd.concat(sdd_frames, ignore_index=True).drop_duplicates()
        print(f"  submitted_diagnostic_discovery hits: {len(sdd_hits)} participant-variant rows")
    else:
        warnings.warn("submitted_diagnostic_discovery returned no results — check table/column names")
        sdd_hits = pd.DataFrame()

    # ----------------------------------------------------------------------
    # Step 4: Build diagnoses table (long format, one row per hit per source)
    #
    # Restrict to QC-passing RNA-seq samples; add rnaseq_platekey.
    # Harmonise column names across sources so 'hgnc_symbol' exists everywhere.
    # ----------------------------------------------------------------------
    diag_parts = []

    if not eq_hits.empty:
        eq = eq_hits.rename(columns={"gene_name": "hgnc_symbol",
                                     "assembly": "genome_build"}).copy()
        eq["source"] = "gmc_exit_questionnaire"
        eq = eq[eq["participant_id"].isin(pid_to_platekey)]
        eq["rnaseq_platekey"] = eq["participant_id"].map(pid_to_platekey)
        diag_parts.append(eq)

    if not sdd_hits.empty:
        sd = sdd_hits.rename(columns={"gene_of_interest": "hgnc_symbol"}).copy()
        sd["source"] = "submitted_diagnostic_discovery"
        sd = sd[sd["participant_id"].isin(pid_to_platekey)]
        sd["rnaseq_platekey"] = sd["participant_id"].map(pid_to_platekey)
        diag_parts.append(sd)

    if not tier12_hits.empty:
        td = tier12_hits.rename(columns={"genomic_feature_hgnc": "hgnc_symbol",
                                         "assembly": "genome_build"}).copy()
        td["source"] = "tiering_data"
        td = td[td["participant_id"].isin(pid_to_platekey)]
        td["rnaseq_platekey"] = td["participant_id"].map(pid_to_platekey)
        diag_parts.append(td)

    if diag_parts:
        diagnoses = pd.concat(diag_parts, ignore_index=True, sort=False)
        # Consistent column order
        col_order = [
            "participant_id", "rnaseq_platekey", "hgnc_symbol",
            "source", "case_solved_family", "acmg_classification", "tier",
            "chromosome", "position", "reference", "alternate", "genome_build",
            "genotype", "mode_of_inheritance", "phenotype",
            "variant_of_interest", "zygosity",
        ]
        diagnoses = diagnoses.reindex(columns=[c for c in col_order if c in diagnoses.columns])
        n_carriers = diagnoses["participant_id"].nunique()
        print(f"\nDiagnoses table: {len(diagnoses)} rows, {n_carriers} unique participants with RNA-seq")
    else:
        warnings.warn("No diagnosis data retrieved — diagnoses table will be empty")
        diagnoses = pd.DataFrame()
        n_carriers = 0

    # ----------------------------------------------------------------------
    # Step 5: Collapse to per-participant evidence flag.
    #
    # spliceopathy_evidence is ordered: 'solved' > 'sdd' > 'tiered' > 'none'.
    # Each participant gets their strongest-available source. The same source
    # ranking is used downstream (score_gel_controls.R) to decide whether to
    # tint/label a dot as a confirmed carrier or a weaker candidate.
    # ----------------------------------------------------------------------
    SOURCE_TO_LEVEL = {
        "gmc_exit_questionnaire":         "solved",
        "submitted_diagnostic_discovery": "sdd",
        "tiering_data":                   "tiered",
    }
    LEVEL_RANK = {"none": 0, "tiered": 1, "sdd": 2, "solved": 3}

    if not diagnoses.empty:
        d = diagnoses.copy()
        d["_level"] = d["source"].map(SOURCE_TO_LEVEL).fillna("tiered")
        d["_rank"]  = d["_level"].map(LEVEL_RANK)

        # Strongest level per participant
        best = (
            d.sort_values("_rank", ascending=False)
             .drop_duplicates(subset="participant_id", keep="first")
             [["participant_id", "_level"]]
             .rename(columns={"_level": "spliceopathy_evidence"})
        )

        # Gene list: collapse all genes seen at the participant's STRONGEST level
        # (so a 'solved' participant's gene column only shows solved genes, not
        # every candidate they were ever tiered for).
        strongest = d.merge(best, on="participant_id")
        strongest = strongest[strongest["_level"] == strongest["spliceopathy_evidence"]]
        gene_df = (
            strongest.groupby("participant_id")["hgnc_symbol"]
            .apply(lambda g: ";".join(sorted(set(g.dropna()))))
            .reset_index()
            .rename(columns={"hgnc_symbol": "spliceopathy_gene"})
        )
        splice_flag_df = best.merge(gene_df, on="participant_id", how="left")
    else:
        splice_flag_df = pd.DataFrame(
            columns=["participant_id", "spliceopathy_evidence", "spliceopathy_gene"]
        )

    # ----------------------------------------------------------------------
    # Step 6: Join QC-passing samples with spliceopathy evidence
    # ----------------------------------------------------------------------
    eligible = qc_keep.merge(splice_flag_df, on="participant_id", how="left")
    eligible["spliceopathy_evidence"] = eligible["spliceopathy_evidence"].fillna("none")
    eligible["spliceopathy_gene"]     = eligible["spliceopathy_gene"].fillna("")
    # Back-compat boolean: TRUE whenever there's ANY evidence. Downstream
    # code that only cares about "any hit" keeps working; code that wants to
    # be strict should filter on spliceopathy_evidence == 'solved'.
    eligible["spliceopathy_flag"] = eligible["spliceopathy_evidence"] != "none"

    assert len(eligible) == len(qc_keep), "Row count changed after flag join"

    print(f"\n{mode.title()}-mode whole-blood samples: {len(eligible)}")
    print(f"  of which confirmed spliceopathy-flagged: {eligible['spliceopathy_flag'].sum()}")

    # ----------------------------------------------------------------------
    # Step 7: Fetch BAM file paths (batched — platekey list may be large)
    #
    # Table: transcriptome_file_paths_and_types (confirmed 2026-03-25)
    # BAMs: file_type = 'Assembly', file_sub_type = 'BAM'
    # ----------------------------------------------------------------------
    print(f"\nQuerying transcriptome_file_paths_and_types for BAM paths (batches of {BAM_BATCH})...")

    platekeys = eligible["rnaseq_platekey"].dropna().tolist()
    file_frames = []
    for batch in batched(platekeys, BAM_BATCH):
        sql = f"""
            SELECT participant_id, platekey, file_path, delivery_date
            FROM transcriptome_file_paths_and_types
            WHERE file_type     = 'Assembly'
              AND file_sub_type = 'BAM'
              AND platekey IN ({sql_in(batch)})
        """
        df = execute_sql(sql, "transcriptome_file_paths_and_types")
        if not df.empty:
            file_frames.append(df)

    if not file_frames:
        sys.exit("ERROR: No BAM paths found — check file_type/file_sub_type values and platekey column")

    files_raw = pd.concat(file_frames, ignore_index=True).drop_duplicates(subset="platekey")
    print(f"  BAM paths returned: {len(files_raw)}")
    assert len(files_raw) > 0, "No BAM paths found"

    # ----------------------------------------------------------------------
    # Step 7b: Participant demographics (sex, ethnicity, year of birth)
    #
    # Table: participant (confirmed 2026-03-28)
    # Columns: participant_phenotypic_sex, participant_ethnic_category, year_of_birth
    # Age at RNA-seq derived as: delivery_date.year - year_of_birth
    # ----------------------------------------------------------------------
    print(f"\nQuerying participant demographics (batches of {BAM_BATCH})...")
    participant_ids = qc_keep["participant_id"].dropna().tolist()
    demo_frames = []
    for batch in batched(participant_ids, BAM_BATCH):
        sql = f"""
            SELECT
                participant_id,
                participant_phenotypic_sex,
                participant_ethnic_category,
                year_of_birth
            FROM participant
            WHERE participant_id IN ({sql_in(batch)})
        """
        df = execute_sql(sql, "participant demographics")
        if not df.empty:
            demo_frames.append(df)

    if demo_frames:
        demographics = pd.concat(demo_frames, ignore_index=True).drop_duplicates(subset="participant_id")
        print(f"  Demographics returned: {len(demographics)} participants")
    else:
        warnings.warn("participant table returned no demographics — sex/ethnicity/age will be missing")
        demographics = pd.DataFrame(columns=["participant_id", "participant_phenotypic_sex",
                                              "participant_ethnic_category", "year_of_birth"])

    # ----------------------------------------------------------------------
    # Step 8: Join QC metrics + BAM paths + demographics → final manifest
    # ----------------------------------------------------------------------
    manifest = eligible.merge(
        files_raw[["platekey", "file_path", "delivery_date"]],
        left_on  = "rnaseq_platekey",
        right_on = "platekey",
        how      = "inner",
    ).drop(columns=["platekey"])

    manifest = manifest.merge(demographics, on="participant_id", how="left")

    # Age at RNA-seq = year of delivery - year of birth
    manifest["delivery_year"] = pd.to_datetime(manifest["delivery_date"], errors="coerce").dt.year
    manifest["year_of_birth"]  = pd.to_numeric(manifest["year_of_birth"], errors="coerce")
    manifest["age_at_rnaseq"]  = manifest["delivery_year"] - manifest["year_of_birth"]
    manifest = manifest.drop(columns=["delivery_year"])

    print(f"Final manifest: {len(manifest)} samples with BAM paths")
    assert len(manifest) > 0.8 * len(eligible), \
        f"Lost >20% of samples in BAM path join ({len(manifest)} / {len(eligible)})"

    # ----------------------------------------------------------------------
    # Step 9: Write outputs
    # ----------------------------------------------------------------------
    manifest.to_csv(out_manifest, index=False)
    print(f"\nManifest written to:   {out_manifest}")

    if not diagnoses.empty:
        diagnoses.to_csv(out_diagnoses, index=False)
        print(f"Diagnoses written to:  {out_diagnoses}")

    print("\n=== Summary ===")
    print(f"  Mode:                                {mode}")
    print(f"  Whole-blood samples in scope:        {len(qc_keep)}")
    print(f"  With BAM paths (total manifest):     {len(manifest)}")
    print(f"  Diagnoses table rows:                {len(diagnoses)}")
    print(f"  Unique participants with a hit:      {n_carriers}")
    print("\nSpliceopathy evidence breakdown (per sample, strongest source):")
    for lvl in ("solved", "sdd", "tiered", "none"):
        n = (manifest["spliceopathy_evidence"] == lvl).sum()
        print(f"  {lvl:>7}: {n}")
    print(f"\nManifest columns: {', '.join(manifest.columns)}")
    print("\nNOTE: flagged samples are RETAINED — apply signature to all, then")
    print("      compare 'solved' vs 'none' to test specificity. 'sdd' / 'tiered'")
    print("      are weaker evidence (candidate variants, not confirmed diagnoses).")
    print("\nNOTE: outputs must not leave GEL without airlock approval.")


if __name__ == "__main__":
    main()
