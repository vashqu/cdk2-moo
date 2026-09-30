#!/usr/bin/env python3
"""
Stage 1: download raw CDK2 bioactivity data from ChEMBL and cache it.

Writes ONE file: data/raw/chembl_<TARGET>_activities.csv

This is the raw layer. No filtering, no deduplication, no unit conversion.
Censored values ('IC50 > 10000 nM'), non-nM units and curator-flagged records
are all kept, because every one of those is a curation decision and curation
happens in stage 2 against this cached file. Re-running curation must never
require another API crawl.

Run:
    python scripts/01_fetch_chembl.py
    python scripts/01_fetch_chembl.py --force   # re-download over the cache
"""

import argparse

import pandas as pd

from cdk2moo import config
from cdk2moo.data import fetch_activities

# Downloaded regardless of what we end up using, so that revisiting the
# IC50-only decision costs nothing. See notebooks/explore_target.ipynb.
DOWNLOAD_TYPES = ["IC50", "Ki", "Kd"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="re-download even if the cache file exists")
    args = ap.parse_args()

    if config.TARGET_CHEMBL_ID is None:
        raise SystemExit(
            "config.TARGET_CHEMBL_ID is None.\n"
            "Run notebooks/explore_target.ipynb and set the verified ID first."
        )

    out_path = config.RAW_DIR / f"chembl_{config.TARGET_CHEMBL_ID}_activities.csv"

    if out_path.exists() and not args.force:
        df = pd.read_csv(out_path)
        print(f"Cache exists: {out_path}")
        print(f"  {len(df)} records, {df['molecule_chembl_id'].nunique()} unique compounds")
        print("  Use --force to re-download.")
        return

    print(f"Target: {config.TARGET_CHEMBL_ID} (UniProt {config.TARGET_UNIPROT})")
    print(f"Types:  {', '.join(DOWNLOAD_TYPES)}\n")

    df = fetch_activities(config.TARGET_CHEMBL_ID, DOWNLOAD_TYPES)

    df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")

    # --- What did we actually get? ------------------------------------
    print(f"\n{len(df)} records, "
          f"{df['molecule_chembl_id'].nunique()} unique compounds, "
          f"{df['document_chembl_id'].nunique()} source documents\n")

    print("Records per measurement type:")
    print(df["standard_type"].value_counts().to_string())

    print("\nRelation (censoring):")
    print(df["standard_relation"].value_counts(dropna=False).to_string())

    print("\nUnits:")
    print(df["standard_units"].value_counts(dropna=False).head(10).to_string())

    print("\nAssay type (B = binding, F = functional):")
    print(df["assay_type"].value_counts(dropna=False).to_string())

    flagged = df["data_validity_comment"].notna().sum()
    print(f"\nCurator-flagged records: {flagged}")
    if flagged:
        print(df["data_validity_comment"].value_counts().to_string())

    missing_smiles = df["canonical_smiles"].isna().sum()
    print(f"\nRecords with no SMILES: {missing_smiles}")

    print("\nNext: stage 2 curation.")


if __name__ == "__main__":
    main()

