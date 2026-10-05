#!/usr/bin/env python3
"""
Stage 12: read-only audits of train/test leakage and of the saved docking and PoseBusters artifacts.

Reads : processed/cdk2_ic50_curated.csv, cdk2_splits.csv, ecfp4.npy, results/05e_holdout_sets.csv of the selected campaign scope, plus its docking export
        results/08d_export/{docking_scores.csv, poses_seed*.sdf} and results/posebusters_audit/posebusters_checks.csv
Writes (inside --out-dir, default results/audit/; no split, score, pose or result is changed):
        leakage_by_seed.csv, leakage_summary.csv, control_overlap.csv, score_pose_coverage.csv,
        posebusters_coverage.json

What each overlap means is defined in cdk2moo/audit.py. Splits are used exactly as stored; the script verifies that the stored
split rows line up with the curated molecules before using them, and stops with a diagnostic if not.

`primary_document` is the first non-null document met while aggregating a structure's records; it is not a verified earliest
publication, and the paper split is not a temporal split.

Run:
    python scripts/12_audit_leakage_and_artifacts.py
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.audit import (check_alignment, key_overlap, fingerprint_keys, connectivity_keys, document_overlap)
from cdk2moo.posebusters_audit import REQUIRED_CHECKS, read_sdf_records

SPLIT_COLUMNS = {"random": "random_split", "scaffold": "scaffold_split", "paper": "paper_split"}


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(config.RESULTS_DIR / "audit"))
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    out = Path(args.out_dir)
    if out.exists() and any(out.iterdir()) and not args.overwrite:
        sys.exit(f"{out} is not empty; choose another --out-dir or pass --overwrite.")
    out.mkdir(parents=True, exist_ok=True)

    curated = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    splits = pd.read_csv(config.PROCESSED_DIR / "cdk2_splits.csv")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    keys = {"identity_full": curated["inchikey"].to_numpy(dtype=object),
            "identity_connectivity": connectivity_keys(curated["inchikey"]),
            "fingerprint": fingerprint_keys(fps)}

    # ---- A. leakage across every stored split ---------------------------------------------------
    rows = []
    for seed in sorted(splits["seed"].unique()):
        one = check_alignment(splits[splits["seed"] == seed], curated)
        scaffolds = one["scaffold"].to_numpy(dtype=object)
        for name, column in SPLIT_COLUMNS.items():
            train = np.where(one[column] == "train")[0]
            test = np.where(one[column] == "test")[0]
            row = {"seed": seed, "split": name, "n_train": len(train), "n_test": len(test)}
            for label, k in keys.items():
                row[f"{label}_test_molecules"], row[f"{label}_distinct_keys"] = key_overlap(train, test, k)
            row["scaffold_test_molecules"], row["scaffold_distinct"] = key_overlap(train, test, scaffolds)
            docs = document_overlap(train, test, curated["all_documents"].fillna("").to_numpy(),
                                    curated["primary_document"].fillna("").to_numpy())
            row["any_document_test_molecules"], row["primary_document_test_molecules"] = docs["any_document"], docs["primary_document"]
            rows.append(row)
    by_seed = pd.DataFrame(rows)
    by_seed.to_csv(out / "leakage_by_seed.csv", index=False)
    cols = [c for c in by_seed.columns if c.endswith("test_molecules")]
    summary = by_seed.groupby("split")[cols].agg(["median", "min", "max"])
    summary.to_csv(out / "leakage_summary.csv")
    print(f"{len(curated)} curated molecules; {splits['seed'].nunique()} split seeds; every stored split verified to align with the curated rows.")
    print("\nTest molecules (of ~403 per split) overlapping the training set, median [min, max] over seeds:")
    print(f"  {'overlap type':<34}" + "".join(f"{s:>22}" for s in SPLIT_COLUMNS))
    for c in cols:
        cells = [f"{summary.loc[s, (c, 'median')]:.0f} [{summary.loc[s, (c, 'min')]:.0f}, {summary.loc[s, (c, 'max')]:.0f}]" for s in SPLIT_COLUMNS]
        print(f"  {c.replace('_test_molecules', '').replace('_', ' '):<34}" + "".join(f"{x:>22}" for x in cells))
    print(f"  (identity_full is 0 by construction of the curation: duplicates were merged by InChIKey. identical fingerprints "
          f"between distinct molecules arise from stereoisomers and ECFP4 collisions.)")

    # ---- B. positive controls (held-out vs remaining training) -------------------------------------
    roles = pd.read_csv(config.RESULTS_DIR / "05e_holdout_sets.csv")
    one = roles.assign(inchikey=curated["inchikey"]).sort_values("row")
    scaffolds = roles.sort_values("row")["scaffold"].to_numpy(dtype=object)
    control_rows = []
    for control in ("molecule", "scaffold"):
        held_all = np.where(roles.sort_values("row")[f"{control}_control"] == "held")[0]
        train = np.where(roles.sort_values("row")[f"{control}_control"] == "train")[0]
        actives = held_all[roles.sort_values("row")["is_top_decile"].to_numpy()[held_all]]
        for group, idx in (("held-out actives", actives), ("all held-out molecules", held_all)):
            row = {"control": control, "group": group, "n": len(idx)}
            for label, k in keys.items():
                row[f"{label}_overlap_with_remaining_training"] = key_overlap(train, idx, k)[0]
            row["scaffold_overlap_with_remaining_training"] = key_overlap(train, idx, scaffolds)[0]
            row["any_document_overlap"] = document_overlap(train, idx, curated["all_documents"].fillna("").to_numpy())["any_document"]
            control_rows.append(row)
    controls = pd.DataFrame(control_rows)
    controls.to_csv(out / "control_overlap.csv", index=False)
    print("\nPositive controls: held-out molecules that overlap the remaining training set")
    print(controls.to_string(index=False))

    # ---- C. score-to-pose coverage by Vina seed (historical docking artifacts) --------------------------
    scores = pd.read_csv(config.RESULTS_DIR / "08d_export" / "docking_scores.csv")
    coverage = []
    for seed in sorted(scores["vina_seed"].unique()):
        pose_file = config.RESULTS_DIR / "08d_export" / f"poses_seed{seed}.sdf"
        records = read_sdf_records(pose_file) if pose_file.exists() else []
        pose_smiles = [r["smiles"] for r in records if r["smiles"] is not None]
        s = scores[scores["vina_seed"] == seed]
        ok = set(s[s["status"] == "ok"]["smiles"])
        coverage.append({"vina_seed": seed, "score_rows": len(s), "ok_scores": len(ok), "pose_records": len(records),
                         "records_without_identity": len(records) - len(pose_smiles), "distinct_pose_identities": len(set(pose_smiles)),
                         "duplicate_pose_identities": len(pose_smiles) - len(set(pose_smiles)),
                         "duplicate_score_rows": int(s.duplicated("smiles").sum()),
                         "ok_scores_without_pose": len(ok - set(pose_smiles)), "poses_without_ok_score": len(set(pose_smiles) - ok),
                         "non_ok_status_counts": json.dumps(s[s["status"] != "ok"]["status"].value_counts().to_dict())})
    coverage = pd.DataFrame(coverage)
    coverage.to_csv(out / "score_pose_coverage.csv", index=False)
    print("\nScore-to-pose coverage by Vina seed (this scope's docking export)")
    print(coverage.drop(columns="non_ok_status_counts").to_string(index=False))
    print("  non-ok statuses:", {int(r.vina_seed): r.non_ok_status_counts for r in coverage.itertuples()})

    # ---- D. PoseBusters check coverage ---------------------------------------------------------------------
    checks = pd.read_csv(config.RESULTS_DIR / "posebusters_audit" / "posebusters_checks.csv")
    report = {"required_checks": len(REQUIRED_CHECKS), "audit_poses": int(len(checks)),
              "required_checks_present_as_columns": int(sum(c in checks.columns for c in REQUIRED_CHECKS)),
              "status_counts": checks["status"].value_counts().to_dict(),
              "missing_by_check": {c: int(checks[c].isna().sum()) for c in REQUIRED_CHECKS if checks[c].isna().any()},
              "failed_by_check": {c: int((checks[c] == False).sum()) for c in REQUIRED_CHECKS if (checks[c] == False).any()}}   # noqa: E712
    with open(out / "posebusters_coverage.json", "w") as fh:
        json.dump(report, fh, indent=2)
    print("\nPoseBusters required-check coverage:")
    for k, v in report.items():
        print(f"  {k}: {v}")

    # ---- E. terminology ---------------------------------------------------------------------------------------
    multi = int((curated["n_documents"] > 1).sum())
    print(f"\n`primary_document` = first non-null document met when aggregating a structure's records (order of the downloaded file); "
          f"{multi} of {len(curated)} molecules were reported in more than one document. Not a verified earliest publication; no temporal split exists.")
    print(f"\nWrote audit files to {out}")


if __name__ == "__main__":
    main()
