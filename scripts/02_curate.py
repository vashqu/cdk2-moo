#!/usr/bin/env python3
"""
Stage 2: curate raw ChEMBL records into one labelled row per molecule.

Reads : data/raw/chembl_<TARGET>_activities.csv
Writes: data/processed/cdk2_ic50_curated.csv     <- the modelling dataset
        data/processed/cdk2_kikd_external.csv    <- Ki/Kd-only compounds
        results/02_curation_funnel.json          <- what was dropped, and why

The second file is compounds measured against the same protein on a different
readout, with no IC50 at all. They never enter training and serve as an extra
external check later.

Run:
    python scripts/02_curate.py
    python scripts/02_curate.py --no-tautomer   # faster, less thorough
"""

import argparse
import json

import pandas as pd

from cdk2moo import config
from cdk2moo.curate import replicate_spread
from cdk2moo.curate_set import curate_set


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-tautomer", action="store_true",
                    help="skip canonical tautomer generation (faster, less thorough)")
    args = ap.parse_args()
    tautomer = not args.no_tautomer

    raw_path = config.RAW_DIR / f"chembl_{config.TARGET_CHEMBL_ID}_activities.csv"
    if not raw_path.exists():
        raise SystemExit(f"Missing {raw_path}. Run scripts/01_fetch_chembl.py first.")

    raw = pd.read_csv(raw_path)
    print(f"Loaded {len(raw)} raw records from {raw_path.name}")

    # ---- the modelling dataset ---------------------------------------
    agg, report = curate_set(raw, config.KEEP_TYPES, config.KEEP_ASSAY_TYPES,
                             tautomer, "IC50 (modelling set)")

    spread = replicate_spread(agg)
    print("\n  Disagreement between repeated records (diagnostic only):")
    for k, v in spread.items():
        print(f"    {k}: {v:.3f}" if isinstance(v, float) else f"    {k}: {v}")
    print("\n  This is NOT a noise floor: most replicates are the same "
          "measurement curated twice. Cite the literature figure (~0.5 log "
          "units) as external. See decisions.md D-09.")

    print("\n  pActivity distribution:")
    print(agg["pactivity"].describe().to_string())

    out_path = config.PROCESSED_DIR / "cdk2_ic50_curated.csv"
    agg.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}  ({len(agg)} molecules)")

    # ---- external set: Ki/Kd compounds with no IC50 ------------------
    # No assay-type filter here: ChEMBL labels the 395-compound kinome Kd panel
    # (CHEMBL1201862) as F although Kd is a binding readout. See D-11.
    ext, _ = curate_set(raw, ["Ki", "Kd"], None, tautomer, "Ki/Kd (external check)")
    if len(ext):
        ext = ext[~ext["inchikey"].isin(set(agg["inchikey"]))].copy()
        ext_path = config.PROCESSED_DIR / "cdk2_kikd_external.csv"
        ext.to_csv(ext_path, index=False)
        print(f"\nWrote {ext_path}  ({len(ext)} molecules with no IC50)")

    report["replicate_spread"] = spread
    report["n_molecules"] = int(len(agg))
    report["n_external"] = int(len(ext))
    with open(config.RESULTS_DIR / "02_curation_funnel.json", "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    print("\nNext: stage 3, featurization and splits.")


if __name__ == "__main__":
    main()