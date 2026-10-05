#!/usr/bin/env python3
"""
Stage 3b: which training molecules may start a GA run (matched initialization for both preparation policies).

Reads : processed/cdk2_ic50_curated.csv
Writes: processed/eligible_starts.csv   one row per curated molecule: eligibility for the primary (20-39 heavy atoms) and the
                                        relaxed (15-50) window, with the reason each excluded molecule was excluded

The rule is defined once under the stricter corrected policy (cdk2moo/starts.py): inside the size window, unchanged by standardization,
closed-shell, and unique after standardization. Legacy and corrected GA runs, and every constrained run and its unconstrained twin, draw
their starting molecules from this list, so paired comparisons start from identical populations.

Run (inside a campaign scope):
    CDK2_CAMPAIGN=<id> CDK2_SCOPE=shared python scripts/03b_eligible_starts.py
"""

import argparse

import pandas as pd

from cdk2moo import config
from cdk2moo.starts import eligible_rows


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    out = config.PROCESSED_DIR / "eligible_starts.csv"
    if out.exists() and not args.overwrite:
        raise SystemExit(f"{out} exists; pass --overwrite to replace it")

    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    cache = {}
    smiles, atoms = df["std_smiles"].tolist(), df["n_heavy_atoms"].tolist()
    primary = eligible_rows(smiles, atoms, (config.GA_MIN_HEAVY_ATOMS, config.GA_MAX_HEAVY_ATOMS), cache)
    relaxed = eligible_rows(smiles, atoms, (config.GA_RELAXED_MIN_HEAVY_ATOMS, config.GA_RELAXED_MAX_HEAVY_ATOMS), cache)
    table = pd.DataFrame({"row": range(len(df)), "inchikey": df["inchikey"],
                          "eligible_primary": [r["eligible"] for r in primary], "reason_primary": [r["reason"] for r in primary],
                          "eligible_relaxed": [r["eligible"] for r in relaxed], "reason_relaxed": [r["reason"] for r in relaxed],
                          "tautomer_status": [r["tautomer_status"] for r in primary]})
    table.to_csv(out, index=False)
    print(f"{len(df)} curated molecules")
    for label, column, reason in [("primary window 20-39", "eligible_primary", "reason_primary"),
                                  ("relaxed window 15-50", "eligible_relaxed", "reason_relaxed")]:
        print(f"  {label}: {int(table[column].sum())} eligible; excluded by reason: "
              f"{table.loc[~table[column], reason].value_counts().to_dict()}")
    print(f"  tautomer enumeration status among molecules reaching standardization: {table['tautomer_status'].value_counts().to_dict()}")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
