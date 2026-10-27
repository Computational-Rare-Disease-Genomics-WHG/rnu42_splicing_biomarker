#!/usr/bin/env python3
"""
donor_variant_get_vcf_paths.py
Step 1 of the donor-site variant check: find the germline small-variant VCF
for each of the 44 GEL dev samples (11 cases + 33 matched controls).

  dev positional IDs (gel:case:N / gel:control:N)
    -> RNA platekey (BAM basename, rMATS b1/b2 order)
    -> participant_id (control_manifest.csv, qcfail_manifest.csv, case TSV)
    -> germline VCF path (LabKey genome_file_paths_and_types)

The LabKey table is queried with SELECT * so no column names are assumed in
SQL. The script prints the table's columns and the distinct values of its
type columns, then keeps one germline "Standard VCF" per participant,
preferring GRCh38 and then the latest delivery. GRCh37-only genomes are kept
(donor_variant_check.R maps coordinates for them). Stops if any participant
has no germline standard VCF.

Run on the HPC (same python env as get_control_manifest.py):
  /resources/conda/miniconda3/envs/py3pypirev2/bin/python -u scripts/gel/donor_variant_get_vcf_paths.py

Output (GEL protected - participant identifiers):
  data/dev_germline_vcf_paths.csv
    sample_id, rna_platekey, participant_id, dna_platekey, genome_build, vcf_path
"""

import os
import sys
import pandas as pd
import labkey

LABKEY_HOST   = "labkey-embassy.gel.zone"
LABKEY_FOLDER = "/main-programme/main-programme_v19_2024-10-31"
LABKEY_SCHEMA = "lists"
FILE_TABLE    = "genome_file_paths_and_types"

DEV_CASE_BAMS = "data/rmats/inputs/renu_11_cases_paths.txt"
DEV_CTRL_BAMS = "data/rmats/inputs/renu_11_cases_control_paths.txt"
MANIFESTS     = ["data/control_manifest.csv", "data/qcfail_manifest.csv"]
CASE_TSV      = "data/12_samples_with_renu_and_RNAseq_for_RD.tsv"
OUT_CSV       = "data/dev_germline_vcf_paths.csv"

N_CASES, N_CONTROLS = 11, 33


def labkey_to_df(sql_query, maxrows=10_000_000):
    """Version-aware LabKey SQL query -> DataFrame (as get_control_manifest.py)."""
    if labkey.__version__ < "2.4.0":
        ctx = labkey.utils.create_server_context(
            domain=LABKEY_HOST, container_path=LABKEY_FOLDER,
            context_path="labkey", use_ssl=True)
        res = labkey.query.execute_sql(ctx, schema_name=LABKEY_SCHEMA,
                                       sql=sql_query, max_rows=maxrows)
    else:
        from labkey.api_wrapper import APIWrapper
        api = APIWrapper(LABKEY_HOST, LABKEY_FOLDER, "labkey", use_ssl=True)
        res = api.query.execute_sql(sql=sql_query, schema_name=LABKEY_SCHEMA,
                                    max_rows=maxrows, waf_encode_sql=False)
    if res is None or "rows" not in res:
        return pd.DataFrame()
    return pd.DataFrame(res["rows"])


def strip_bam(p):
    b = os.path.basename(p)
    return b[:-4] if b.endswith(".bam") else b


def read_bam_list(path):
    with open(path) as fh:
        raw = ",".join(fh.read().splitlines())
    bams = [b.strip() for b in raw.split(",") if b.strip()]
    return [strip_bam(b) for b in bams]


def main():
    for f in [DEV_CASE_BAMS, DEV_CTRL_BAMS, CASE_TSV] + MANIFESTS:
        if not os.path.exists(f):
            sys.exit(f"ERROR: input not found: {f}")

    # ---- dev positional IDs -> RNA platekey ---------------------------------
    case_keys, ctrl_keys = read_bam_list(DEV_CASE_BAMS), read_bam_list(DEV_CTRL_BAMS)
    if (len(case_keys), len(ctrl_keys)) != (N_CASES, N_CONTROLS):
        sys.exit(f"ERROR: expected {N_CASES}+{N_CONTROLS} BAMs, got "
                 f"{len(case_keys)}+{len(ctrl_keys)}")
    dev = pd.DataFrame({
        "sample_id": [f"gel:case:{i+1}" for i in range(N_CASES)] +
                     [f"gel:control:{i+1}" for i in range(N_CONTROLS)],
        "rna_platekey": case_keys + ctrl_keys})

    # ---- RNA platekey -> participant_id -------------------------------------
    rows = []
    for f in MANIFESTS:
        m = pd.read_csv(f, dtype=str)
        for col in ["participant_id", "rnaseq_platekey", "file_path"]:
            if col not in m.columns:
                sys.exit(f"ERROR: {f} has no column {col}")
        bam_key = m["file_path"].map(strip_bam)
        rows.append(pd.DataFrame({"key": m["rnaseq_platekey"], "participant_id": m["participant_id"]}))
        rows.append(pd.DataFrame({"key": bam_key, "participant_id": m["participant_id"]}))
    cg = pd.read_csv(CASE_TSV, sep="\t", dtype=str)
    rows.append(pd.DataFrame({"key": cg["platekey"], "participant_id": cg["participant_id"]}))
    lookup = pd.concat(rows).dropna().drop_duplicates()
    clash = lookup[lookup.duplicated("key", keep=False)]
    if not clash.empty:
        sys.exit(f"ERROR: platekey maps to >1 participant_id:\n{clash}")
    dev = dev.merge(lookup, left_on="rna_platekey", right_on="key", how="left").drop(columns="key")
    if dev["participant_id"].isna().any():
        sys.exit("ERROR: no participant_id for:\n"
                 f"{dev[dev['participant_id'].isna()].to_string(index=False)}")

    # ---- LabKey file table --------------------------------------------------
    pids = ",".join(f"'{p}'" for p in dev["participant_id"].unique())
    files = labkey_to_df(f"SELECT * FROM {FILE_TABLE} WHERE participant_id IN ({pids})")
    if files.empty:
        sys.exit(f"ERROR: {FILE_TABLE} returned no rows for the dev participants")
    files = files.astype(str)
    print(f"{FILE_TABLE}: {len(files)} rows for {files['participant_id'].nunique()} participants")
    print("Columns:", ", ".join(files.columns))
    type_cols = [c for c in files.columns
                 if any(t in c.lower() for t in ["type", "build", "assembly"])]
    for c in type_cols:
        print(f"  distinct {c}: {sorted(files[c].unique())}")

    if "file_path" not in files.columns:
        sys.exit("ERROR: no file_path column - see column list above")

    # Identify columns by their values (names vary between releases)
    def col_with(pattern):
        hits = [c for c in files.columns
                if files[c].str.match(f"(?:{pattern})$").mean() > 0.5]
        if len(hits) != 1:
            sys.exit(f"ERROR: expected one column matching {pattern!r}, found {hits}")
        return hits[0]
    sub_col   = col_with(r"Standard VCF|Genomic VCF|Structural VCF|BAM|Repeat VCF|"
                         r"Array Genotyping VCF|.*VCF.*|BAM")
    build_col = col_with(r"GRCh3[78]")
    ver_col   = col_with(r"V\d+")
    type_col  = [c for c in files.columns if files[c].str.contains("germline").mean() > 0.5]
    print(f"Using columns: file type={sub_col}, build={build_col}, delivery={ver_col}, "
          f"germline={type_col}")

    cand = files[(files[sub_col] == "Standard VCF") &
                 files["file_path"].str.endswith(".vcf.gz")]
    if type_col:
        cand = cand[cand[type_col[0]].str.contains("germline")]
    cand = cand.assign(
        build_rank=(cand[build_col] == "GRCh38").astype(int),
        ver_num=pd.to_numeric(cand[ver_col].str.extract(r"(\d+)", expand=False),
                              errors="coerce").fillna(0))
    # Prefer GRCh38, then the latest delivery
    cand = (cand.sort_values(["participant_id", "build_rank", "ver_num"],
                             ascending=[True, False, False])
                .drop_duplicates("participant_id"))
    print(f"\nSelected one germline standard VCF per participant: {len(cand)}")
    print(cand[build_col].value_counts().to_string())

    missing = sorted(set(dev["participant_id"]) - set(cand["participant_id"]))
    if missing:
        print(files[files["participant_id"].isin(missing)]
              [["participant_id", "file_path", sub_col, build_col, ver_col]].to_string(index=False))
        sys.exit(f"ERROR: {len(missing)} participant(s) with no germline standard VCF")

    cand = cand.rename(columns={"platekey": "dna_platekey", build_col: "genome_build"}) \
        if "platekey" in cand.columns else \
        cand.rename(columns={build_col: "genome_build"}).assign(dna_platekey="")
    out = dev.merge(cand[["participant_id", "dna_platekey", "genome_build", "file_path"]],
                    on="participant_id", how="left").rename(columns={"file_path": "vcf_path"})
    absent = out[~out["vcf_path"].map(os.path.exists)]
    if not absent.empty:
        sys.exit(f"ERROR: VCF path not readable from here:\n{absent.to_string(index=False)}")
    out.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {OUT_CSV}: {len(out)} samples")


if __name__ == "__main__":
    main()
