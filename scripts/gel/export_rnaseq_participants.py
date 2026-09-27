#!/usr/bin/env python
"""
Export the list of GEL participants (and their RNA-seq platekeys) used in
the v4-LASSO splicing biomarker manuscript: dev cases, dev controls, the
6k whole-blood control cohort used in phase-2 robustness, biallelic
RNU4-2 carriers, and the high-depth ReNU case.

Inputs
------
  data/rmats/inputs/renu_11_cases_paths.txt          (11 dev case BAMs)
  data/rmats/inputs/renu_11_cases_control_paths.txt  (33 dev ctrl BAMs)
  data/biallelic_rnu42_platekeys.txt                 (biallelic carriers)
  data/renu_highdepth_case_path.txt                  (HD ReNU case BAM)
  data/control_manifest.csv                          (6k WB cohort - already
                                                      carries participant_id)

LabKey
------
  Authentication via ~/.netrc (machine labkey-embassy.gel.zone), same
  setup as get_control_manifest.py. The
  transcriptome_file_paths_and_types table maps platekey to
  participant_id (one platekey may map to multiple file rows; we
  deduplicate to one (platekey, participant_id) pair).

Outputs
-------
  outputs/v4_signature_participants.csv
    columns: cohort_role, group, rnaseq_platekey, participant_id
      cohort_role : "dev_case" / "dev_control" / "wb_6k_control" /
                    "biallelic_RNU4_2" / "hd_renu_case"
      group       : "case" / "control" / "carrier" (for biallelic)
      rnaseq_platekey, participant_id
  outputs/airlock_participant_ids.txt
    de-duplicated participant_ids, one per line, sorted. Airlock payload
    for the Figure 1-2 + S3/S6 export. Written iff --participant-txt is
    set (default: outputs/airlock_participant_ids.txt; pass NONE to skip).

Run (from workspace root, on GEL Desktop):
  # Dev cohort only (default):
  python -u scripts/gel/export_rnaseq_participants.py
  # Airlock scope (all overlay cohorts):
  python -u scripts/gel/export_rnaseq_participants.py \
      --include-biallelic --include-hd --include-6k

DATA PROTECTION: participant_ids are GEL-protected patient identifiers.
Both outputs are GEL-side only; airlock is required if they need to
leave the RE.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import warnings
from pathlib import Path
from typing import Iterable

import labkey
import pandas as pd

# ----------------------------------------------------------------------------
# Paths + LabKey config
# ----------------------------------------------------------------------------
DEV_CASE_BAMS  = Path("data/rmats/inputs/renu_11_cases_paths.txt")
DEV_CTRL_BAMS  = Path("data/rmats/inputs/renu_11_cases_control_paths.txt")
BIALLELIC_FILE = Path("data/biallelic_rnu42_platekeys.txt")
HD_CASE_FILE   = Path("data/renu_highdepth_case_path.txt")
WB_6K_MANIFEST = Path("data/control_manifest.csv")
OUT_CSV        = Path("outputs/v4_signature_participants.csv")
OUT_TXT        = Path("outputs/airlock_participant_ids.txt")

LABKEY_HOST    = "labkey-embassy.gel.zone"
LABKEY_FOLDER  = "/main-programme/main-programme_v19_2024-10-31"
LABKEY_SCHEMA  = "lists"
BAM_BATCH      = 500   # platekeys per SQL IN clause


def labkey_to_df(sql_query: str, maxrows: int = 10_000_000) -> pd.DataFrame:
    """Version-aware LabKey SQL → DataFrame."""
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
        api = APIWrapper(LABKEY_HOST, LABKEY_FOLDER, "labkey", use_ssl=True)
        results = api.query.execute_sql(
            sql            = sql_query,
            schema_name    = LABKEY_SCHEMA,
            max_rows       = maxrows,
        )
    return pd.DataFrame(results["rows"])


def sql_in(values: Iterable[str]) -> str:
    """Build a SQL IN-list literal from an iterable of strings."""
    quoted = [f"'{v}'" for v in values]
    return ",".join(quoted)


def batched(seq: list, n: int):
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


def read_bam_list(path: Path) -> list[str]:
    """Read a comma-separated BAM-list file, return platekeys (basename - .bam)."""
    if not path.exists():
        warnings.warn(f"BAM list missing: {path}")
        return []
    raw = ",".join(line.strip() for line in path.read_text().splitlines())
    bams = [b.strip() for b in raw.split(",") if b.strip()]
    return [Path(b).stem for b in bams]


def read_one_per_line(path: Path) -> list[str]:
    if not path.exists():
        warnings.warn(f"Platekey list missing: {path}")
        return []
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip()
    ]


def lookup_participants(platekeys: list[str]) -> dict[str, str]:
    """Query LabKey transcriptome_file_paths_and_types to map platekey ->
    participant_id. One platekey can have multiple file rows but should map
    to a single participant_id; we keep the first seen and warn on
    contradictions.
    """
    if not platekeys:
        return {}
    rows: list[pd.DataFrame] = []
    for chunk in batched(sorted(set(platekeys)), BAM_BATCH):
        sql = f"""
            SELECT DISTINCT participant_id, platekey
            FROM transcriptome_file_paths_and_types
            WHERE platekey IN ({sql_in(chunk)})
        """
        try:
            rows.append(labkey_to_df(sql))
        except Exception as e:
            warnings.warn(f"LabKey lookup failed for chunk: {e}")
    if not rows:
        return {}
    df = pd.concat(rows, ignore_index = True)
    df = df.dropna(subset = ["platekey", "participant_id"])
    # Coerce participant_id to str so the CSV doesn't carry numeric formatting.
    df["participant_id"] = df["participant_id"].astype(str)
    df = df.drop_duplicates(subset = ["platekey"])
    mapping = dict(zip(df["platekey"], df["participant_id"]))
    return mapping


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default = str(OUT_CSV),
                        help = "Output CSV path")
    parser.add_argument("--include-biallelic", action = "store_true",
                        help = ("Add biallelic RNU4-2 carriers (n ~3) "
                                "to the output. Default: skip "
                                "(they're validation, not score "
                                "development)."))
    parser.add_argument("--include-hd", action = "store_true",
                        help = ("Add the high-depth ReNU case (n=1) to "
                                "the output. Default: skip "
                                "(validation overlay, not training)."))
    parser.add_argument("--include-6k", action = "store_true",
                        help = ("Include the ~6k WB control cohort from "
                                "control_manifest.csv. Default: skip "
                                "(phase-2 robustness only)."))
    parser.add_argument("--participant-txt", default = str(OUT_TXT),
                        help = ("Also write a de-duplicated participant_ids "
                                "list (one per line) to this path. Pass "
                                "'NONE' to skip. Default: %(default)s."))
    args = parser.parse_args()

    # ---- 1. Collect platekeys per cohort role -------------------------------
    # Default scope = the 44 score-development samples only (11 dev cases
    # + 33 dev controls). Biallelic / HD / 6k are validation overlays
    # opted in via CLI flags above.
    cohort_keys: dict[str, list[str]] = {
        "dev_case":    read_bam_list(DEV_CASE_BAMS),
        "dev_control": read_bam_list(DEV_CTRL_BAMS),
    }
    if args.include_biallelic:
        cohort_keys["biallelic_RNU4_2"] = read_one_per_line(BIALLELIC_FILE)
    if args.include_hd:
        cohort_keys["hd_renu_case"]     = read_bam_list(HD_CASE_FILE)
    for role, keys in cohort_keys.items():
        print(f"  {role:>20s}: {len(keys):>5d} platekey(s) loaded")

    # ---- 2. LabKey lookup for the small cohorts -----------------------------
    all_lookup_keys = sorted(
        {k for keys in cohort_keys.values() for k in keys}
    )
    print(f"\nLooking up {len(all_lookup_keys)} platekeys in LabKey "
          f"transcriptome_file_paths_and_types ...")
    pk_to_pid = lookup_participants(all_lookup_keys)
    n_hit = sum(1 for k in all_lookup_keys if k in pk_to_pid)
    print(f"  Participant_id resolved for {n_hit}/{len(all_lookup_keys)} "
          f"platekeys")

    missing = [k for k in all_lookup_keys if k not in pk_to_pid]
    if missing:
        print(f"  WARNING: {len(missing)} platekey(s) had no participant_id "
              f"in transcriptome_file_paths_and_types:")
        for k in missing[:10]:
            print(f"    {k}")
        if len(missing) > 10:
            print(f"    ... ({len(missing) - 10} more)")

    # ---- 3. 6k whole-blood cohort from manifest -----------------------------
    wb_rows: list[dict] = []
    if not args.include_6k:
        print("\nSkipping 6k WB control cohort (pass --include-6k to add it).")
    elif WB_6K_MANIFEST.exists():
        wb = pd.read_csv(WB_6K_MANIFEST, dtype = str)
        # Expected columns: at minimum rnaseq_platekey, participant_id
        cols = set(wb.columns)
        pk_col  = "rnaseq_platekey" if "rnaseq_platekey" in cols else (
                  "platekey"        if "platekey"        in cols else None)
        pid_col = "participant_id"  if "participant_id"  in cols else None
        if pk_col and pid_col:
            wb_clean = wb[[pk_col, pid_col]].dropna()
            wb_clean = wb_clean.drop_duplicates(subset = [pk_col])
            print(f"\n6k WB control manifest: {len(wb_clean)} rows from "
                  f"{WB_6K_MANIFEST}")
            wb_rows = [
                dict(cohort_role     = "wb_6k_control",
                     group           = "control",
                     rnaseq_platekey = row[pk_col],
                     participant_id  = row[pid_col])
                for _, row in wb_clean.iterrows()
            ]
        else:
            print(f"WARNING: {WB_6K_MANIFEST} missing platekey or "
                  f"participant_id columns; skipping 6k cohort")
    else:
        print(f"\nWARNING: {WB_6K_MANIFEST} not found; skipping 6k cohort")

    # ---- 4. Assemble output -------------------------------------------------
    role_to_group = {
        "dev_case":         "case",
        "dev_control":      "control",
        "biallelic_RNU4_2": "carrier",
        "hd_renu_case":     "case",
    }
    rows: list[dict] = []
    for role, keys in cohort_keys.items():
        for pk in keys:
            rows.append(dict(
                cohort_role     = role,
                group           = role_to_group.get(role, ""),
                rnaseq_platekey = pk,
                participant_id  = pk_to_pid.get(pk, ""),
            ))
    rows.extend(wb_rows)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents = True, exist_ok = True)
    with out_path.open("w", newline = "") as fh:
        w = csv.DictWriter(fh, fieldnames = [
            "cohort_role", "group", "rnaseq_platekey", "participant_id"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {out_path}")

    # ---- 5. Per-cohort summary printed to stdout ----------------------------
    df_out = pd.DataFrame(rows)
    print("\nPer-cohort participant counts:")
    counts = (df_out
              .assign(has_pid = df_out["participant_id"].fillna("") != "")
              .groupby("cohort_role")
              .agg(n_rows         = ("rnaseq_platekey", "count"),
                   n_unique_pids  = ("participant_id", "nunique"),
                   n_missing_pid  = ("has_pid", lambda s: int((~s).sum())))
              .reset_index()
              .sort_values("cohort_role"))
    print(counts.to_string(index = False))

    # ---- 6. De-duplicated participant_ids txt (airlock payload) -------------
    if args.participant_txt.upper() != "NONE":
        txt_path = Path(args.participant_txt)
        txt_path.parent.mkdir(parents = True, exist_ok = True)
        pids = (
            df_out["participant_id"]
            .dropna()
            .map(str)
            .map(str.strip)
        )
        pids = sorted({p for p in pids if p and p.lower() != "nan"})
        txt_path.write_text("\n".join(pids) + "\n")
        print(f"\nWrote {len(pids)} unique participant_ids to {txt_path}")

    print("\nNOTE: participant_ids are GEL-protected. Outputs are GEL-side "
          "only; airlock if they need to leave the RE.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
