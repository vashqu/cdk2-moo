"""
Curating one measurement set end to end (the IC50 set, or the Ki/Kd set).

This wires together the steps defined in curate.py and standardize.py and
prints the funnel as it goes, so the script that calls it stays a short list of
"what to curate". Writes nothing to disk.
"""

import pandas as pd

from cdk2moo import config
from cdk2moo.curate import (
    filter_records,
    add_pactivity,
    check_against_pchembl,
    aggregate_by_structure,
)
from cdk2moo.standardize import standardize_frame


def print_funnel(funnel):
    """Print each filtering step with how many records it removed."""
    prev = None
    for step, n in funnel:
        delta = "" if prev is None else f"   (-{prev - n})"
        print(f"  {n:>6}  {step}{delta}")
        prev = n


def curate_set(raw, keep_types, assay_types, tautomer=True, label=""):
    """Filter -> pActivity -> standardize -> aggregate. Returns (agg, report)."""
    print(f"\n--- {label} ---")

    df, funnel = filter_records(
        raw,
        keep_types=keep_types,
        keep_relation=config.KEEP_RELATION,
        keep_units=config.KEEP_UNITS,
        drop_curator_flagged=config.DROP_CURATOR_FLAGGED,
        mutant_regex=config.MUTANT_REGEX,
        keep_assay_types=assay_types,
    )
    print_funnel(funnel)

    if not len(df):
        return pd.DataFrame(), {"funnel": funnel}

    df = add_pactivity(df)

    n_cmp, n_bad, max_diff = check_against_pchembl(df)
    print(f"\n  pActivity vs ChEMBL's pchembl_value: {n_cmp} comparable, "
          f"{n_bad} disagree by >0.01 (max diff {max_diff:.4f})")
    if n_bad:
        print("  WARNING: disagreement suggests a unit-handling bug. Investigate "
              "before trusting these labels.")

    df = standardize_frame(df, canonical_tautomer=tautomer)

    n_before = len(df)
    fails = df[~df["std_ok"]]
    if len(fails):
        print("\n  standardization failures:")
        print(fails["std_fail_reason"].value_counts().to_string())
    df = df[df["std_ok"]].copy()
    funnel.append(("standardized ok", len(df)))
    print(f"\n  {len(df)} records survived standardization (-{n_before - len(df)})")

    agg = aggregate_by_structure(df)
    funnel.append(("unique structures", len(agg)))

    n_ids = df["molecule_chembl_id"].nunique()
    print(f"\n  {n_ids} unique ChEMBL IDs -> {len(agg)} unique structures")
    if n_ids > len(agg):
        print(f"  {n_ids - len(agg)} ID(s) collapsed into another structure "
              "(salts, tautomers, charge states). These would have leaked "
              "across a train/test split if we had deduplicated by ID.")

    return agg, {"funnel": funnel}
