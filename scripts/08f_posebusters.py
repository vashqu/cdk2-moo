#!/usr/bin/env python3
"""
Stage 8f: are the docked poses physically sensible? (PoseBusters, corrected audit)

Reads : results/08d_export/poses_seed42.sdf (best pose of every docked molecule of this scope), structures/receptor_H.pdb (the receptor the poses
        were docked into), results/08c_docking_sets.csv, results/08d_export/docking_scores.csv (to know how many poses to expect)
Writes (all inside --out-dir, default results/posebusters_audit/ of the selected scope; an existing audit is never replaced):
        posebusters_checks.csv    one row per SDF record: state of every required check, status, counts
        posebusters_summary.json  expected / parsed / audited / passed / failed / unevaluable
        posebusters_by_set.csv    per-set counts with denominators
        posebusters_audit.png     per-set pass / fail / unevaluable

What changed from the first version of this script, and why: see cdk2moo/posebusters_audit.py. In short,
the first version picked check columns by dtype, which silently dropped any check that had a missing
value (it lost `internal_energy`), and it could not tell a missing result from a pass. The historical
outputs (results/08f_posebusters_by_molecule.csv, figures/08_posebusters.png, figures/final/fig5) are kept
as they were and should be read with that caveat.

A pose "passes" only if all 22 required PoseBusters checks are explicitly True. A docking engine
places ligands without overlaps by construction, so passing is weak evidence of a correct pose.

Run:
    python scripts/08f_posebusters.py
    python scripts/08f_posebusters.py --limit 20 --out-dir /tmp/pb_smoke     # smoke test
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from posebusters import PoseBusters

from cdk2moo import config
from cdk2moo.ga_plots import INK, GRID
from cdk2moo.posebusters_audit import (REQUIRED_CHECKS, read_sdf_records, run_posebusters,
                                       aggregate_status, summarize, attach_to_molecules)

def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--pose-file", default=str(config.RESULTS_DIR / "08d_export" / f"poses_seed{config.RANDOM_SEED}.sdf"))
    ap.add_argument("--protein", default=str(config.STRUCTURES_DIR / "receptor_H.pdb"))
    ap.add_argument("--out-dir", default=str(config.RESULTS_DIR / "posebusters_audit"))
    ap.add_argument("--limit", type=int, default=None, help="audit only the first N records (smoke test)")
    ap.add_argument("--overwrite", action="store_true", help="allow replacing files in --out-dir")
    args = ap.parse_args()

    out = Path(args.out_dir)
    targets = [out / n for n in ("posebusters_checks.csv", "posebusters_summary.json", "posebusters_by_set.csv", "posebusters_audit.png")]
    if any(t.exists() for t in targets) and not args.overwrite:
        sys.exit(f"{out} already holds audit files; use another --out-dir or --overwrite.")
    out.mkdir(parents=True, exist_ok=True)

    records = read_sdf_records(args.pose_file)
    if args.limit:
        records = records[:args.limit]
    scores = pd.read_csv(config.RESULTS_DIR / "08d_export" / "docking_scores.csv")
    expected = int(((scores["vina_seed"] == config.RANDOM_SEED) & (scores["status"] == "ok")).sum())
    parsed = sum(r["mol"] is not None for r in records)
    print(f"{len(records)} SDF records read ({parsed} parsed), {expected} poses expected from the scores file"
          + (" (smoke test: limited)" if args.limit else ""))

    report = run_posebusters(records, args.protein, PoseBusters(config="dock"))
    table = aggregate_status(report, records)
    summary = summarize(table, expected if not args.limit else len(records))
    table.to_csv(out / "posebusters_checks.csv", index=False)
    with open(out / "posebusters_summary.json", "w") as fh:
        json.dump({**summary, "required_checks": list(REQUIRED_CHECKS)}, fh, indent=2)
    print("\nSummary (denominator for every rate below = records audited):")
    for key, value in summary.items():
        print(f"  {key:<22}{value}")
    if not args.limit and summary["records_read"] != expected:
        print(f"  WARNING: {summary['records_read']} records read but {expected} expected")

    # ---- per set, with denominators --------------------------------------------------------
    sets = pd.read_csv(config.RESULTS_DIR / "08c_docking_sets.csv")[["smiles", "set"]]
    SETS = list(dict.fromkeys(sets["set"]))
    merged = attach_to_molecules(sets, table)
    rows = []
    print(f"\n  {'set':<17}{'n':>5}{'pass':>6}{'fail':>6}{'unevaluable':>13}{'no audit':>10}{'pass rate':>11}   (rate = pass / n in set)")
    for name in SETS:
        d = merged[merged["set"] == name]
        counts = d["status"].value_counts()
        row = {"set": name, "n": len(d), "pass": int(counts.get("pass", 0)), "fail": int(counts.get("fail", 0)),
               "unevaluable": int(counts.get("unevaluable", 0)), "no_audit_result": int(counts.get("no_audit_result", 0))}
        row["pass_rate"] = row["pass"] / row["n"]
        rows.append(row)
        print(f"  {name:<17}{row['n']:>5}{row['pass']:>6}{row['fail']:>6}{row['unevaluable']:>13}"
              f"{row['no_audit_result']:>10}{100 * row['pass_rate']:>10.0f}%")
    by_set = pd.DataFrame(rows)
    by_set.to_csv(out / "posebusters_by_set.csv", index=False)

    failing = table["failed_checks"].str.split("|").explode()
    missing = table["missing_checks"].str.split("|").explode()
    print("\nChecks that failed (poses):", failing[failing != ""].value_counts().to_dict() or "none")
    print("Checks that were missing (poses):", missing[missing != ""].value_counts().to_dict() or "none")

    vina = scores[(scores["vina_seed"] == config.RANDOM_SEED) & (scores["status"] == "ok")][["smiles", "vina_score"]]
    joined = vina.merge(table[["smiles", "status"]], on="smiles", how="left", validate="many_to_one")
    joined["status"] = joined["status"].fillna("no_audit_result")
    print("\nMedian Vina score by status (kcal/mol):",
          {k: round(v, 2) for k, v in joined.groupby("status")["vina_score"].median().items()})

    fig, ax = plt.subplots(figsize=(9, 4.8), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    left = pd.Series(0.0, index=by_set.index)
    for column, color in [("pass", "#2a78d6"), ("fail", "#e34948"), ("unevaluable", "#eda100"), ("no_audit_result", "#9a9a94")]:
        share = 100 * by_set[column] / by_set["n"]
        ax.barh(range(len(by_set)), share, left=left, color=color, height=0.6, label=column)
        left += share
    ax.set_yticks(range(len(by_set)))
    ax.set_yticklabels(by_set["set"], color=INK, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of molecules in the set (denominator = n in set)", color=INK)
    ax.grid(True, axis="x", color=GRID, linewidth=0.8)
    ax.legend(frameon=False, fontsize=8, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "posebusters_audit.png", dpi=150)
    print(f"\nWrote {out}/posebusters_checks.csv, posebusters_summary.json, posebusters_by_set.csv, posebusters_audit.png")


if __name__ == "__main__":
    main()
