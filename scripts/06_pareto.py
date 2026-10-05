#!/usr/bin/env python3
"""
Stage 6: Pareto analysis of the generated molecules over three objectives:
predicted pActivity (real surrogate, maximise), QED (maximise), SA (minimise).

Reads : results/05_ga_populations_scored.csv, results/04_test_predictions.csv,
        data/processed/cdk2_ic50_curated.csv
Writes: results/06_pareto_front.csv            <- the pooled front, with similarity columns
        figures/06_pareto_fronts.png           <- three 2-D views with overlays
        figures/06_pareto_structures.png       <- four front molecules drawn out

Generated molecules: unique SMILES among population members with
birth_generation > 0, pooled over all generations and seeds of an arm.

Overlays, and what they are NOT:
  known actives      top decile of the training set by MEASURED pActivity. They
                     are in-sample (the surrogate was trained on them), not the
                     held-out positive control, which is a separate experiment
                     not yet built. Their activity axis is measured; the
                     generated molecules' is predicted.
  random reference   300 random training molecules (measured). This is a stand-in
                     for the planned random-ChEMBL floor; no external ChEMBL
                     sample has been downloaded.

Each front molecule also carries its raw ECFP4 max-Tanimoto and size percentile,
so the front can be read against Stage 4's error-vs-similarity curve.

Run:
    python scripts/06_pareto.py
"""

import argparse

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.objectives import qed_and_sa, ARMS
from cdk2moo.pareto import non_dominated, dominated_by_any
from cdk2moo.pareto_plots import plot_fronts, draw_structures


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=config.RANDOM_SEED)
    args = ap.parse_args()
    print(f"seed = {args.seed} (random reference sample only)")

    pop = pd.read_csv(config.RESULTS_DIR / "05_ga_populations_scored.csv")
    gen = pop[pop["birth_generation"] > 0]
    train = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")

    # ---- unique generated molecules per arm, and each arm's own front ------
    by_arm = {}
    print("Generated molecules (birth_generation > 0), unique per arm, pooled over seeds and generations")
    for arm in ARMS:
        d = (gen[gen["arm"] == arm].drop_duplicates("smiles")
             .rename(columns={"pred_real": "pred"}).reset_index(drop=True))
        objectives = np.column_stack([d["pred"], d["qed"], -d["sa"]])
        d["on_arm_front"] = non_dominated(objectives)
        by_arm[arm] = d
        print(f"  {arm:<16} {len(d):>6} molecules, own Pareto front {int(d['on_arm_front'].sum()):>4}")

    # ---- the pooled front: all arms together, one row per distinct molecule --
    pooled = (gen.rename(columns={"pred_real": "pred"})
              .groupby("smiles")
              .agg(pred=("pred", "first"), qed=("qed", "first"), sa=("sa", "first"),
                   max_tanimoto=("max_tanimoto", "first"),
                   size_pctile=("size_pctile", "first"),
                   n_heavy_atoms=("n_heavy_atoms", "first"),
                   arms=("arm", lambda a: "|".join(sorted(set(a)))))
              .reset_index())
    obj = np.column_stack([pooled["pred"], pooled["qed"], -pooled["sa"]])
    front = pooled[non_dominated(obj)].sort_values("pred", ascending=False).reset_index(drop=True)
    print(f"\nPooled: {len(pooled)} distinct generated molecules; Pareto front has {len(front)}")

    print("  which arm(s) produced the front molecules:")
    for arm in ARMS:
        n = front["arms"].str.contains(arm).sum()
        print(f"    {arm:<16} {n:>4} of {len(front)}")
    print(f"  front ranges: predicted pActivity {front['pred'].min():.2f}-{front['pred'].max():.2f}, "
          f"QED {front['qed'].min():.2f}-{front['qed'].max():.2f}, "
          f"SA {front['sa'].min():.2f}-{front['sa'].max():.2f}, "
          f"heavy atoms {front['n_heavy_atoms'].min()}-{front['n_heavy_atoms'].max()}")

    # ---- reference sets ----------------------------------------------------
    qed, sa = qed_and_sa(list(train["std_smiles"]))
    train["qed"], train["sa"], train["pred"] = qed, sa, train["pactivity"]
    cutoff = train["pactivity"].quantile(0.9)
    actives = train[train["pactivity"] >= cutoff]
    reference = train.sample(300, random_state=args.seed)
    print(f"\nKnown actives: {len(actives)} training molecules with measured pActivity >= {cutoff:.2f}"
          f" (in-sample); random reference: {len(reference)} training molecules")

    front_pts = np.column_stack([front["pred"], front["qed"], -front["sa"]])
    act_pts = np.column_stack([actives["pred"], actives["qed"], -actives["sa"]])
    ref_pts = np.column_stack([reference["pred"], reference["qed"], -reference["sa"]])
    print("  Trade-off comparison (generated = PREDICTED activity, reference = MEASURED; indicative only):")
    print(f"    known actives dominated by at least one front molecule: "
          f"{int(dominated_by_any(act_pts, front_pts).sum())}/{len(actives)}")
    print(f"    random training molecules dominated by a front molecule: "
          f"{int(dominated_by_any(ref_pts, front_pts).sum())}/{len(reference)}")
    print(f"    front molecules dominated by a known active: "
          f"{int(dominated_by_any(front_pts, act_pts).sum())}/{len(front)}")
    print(f"    known actives: median QED {actives['qed'].median():.2f}, median SA {actives['sa'].median():.2f}; "
          f"front: median QED {front['qed'].median():.2f}, median SA {front['sa'].median():.2f}")

    # ---- where is the front, relative to what Stage 4 says we can trust? ---
    print("\nFront molecules against Stage 4's error-vs-similarity curve (ECFP4 Tanimoto, 2048-bit):")
    p4 = pd.read_csv(config.RESULTS_DIR / "04_test_predictions.csv")
    p4["abs_err"] = (p4["pred_rf"] - p4["y_true"]).abs()
    edges = [0, 0.4, 0.6, 0.8, 1.01]
    names = ["<0.4", "0.4-0.6", "0.6-0.8", ">=0.8"]
    p4["bin"] = pd.cut(p4["max_tanimoto_to_train"], edges, right=False, labels=names)
    front["bin"] = pd.cut(front["max_tanimoto"], edges, right=False, labels=names)
    mae = p4.groupby("bin", observed=True)["abs_err"].mean()
    share = front["bin"].value_counts(normalize=True)
    print(f"  {'raw max Tanimoto bin':<22}{'Stage 4 MAE':>12}{'share of front':>16}")
    for nm in names:
        print(f"  {nm:<22}{mae.get(nm, float('nan')):>12.2f}{100 * share.get(nm, 0):>15.0f}%")
    print(f"  front median raw max Tanimoto {front['max_tanimoto'].median():.2f}, "
          f"median size percentile {front['size_pctile'].median():.1f}")
    print(f"  front molecules predicted >= 8 : {int((front['pred'] >= 8).sum())}; "
          f"of those with raw similarity < 0.6: {int(((front['pred'] >= 8) & (front['max_tanimoto'] < 0.6)).sum())}")
    front.drop(columns="bin").to_csv(config.RESULTS_DIR / "06_pareto_front.csv", index=False,
                                     float_format="%.4f")

    # ---- four front molecules to draw -------------------------------------
    # Rule (fixed before looking): the front's best on each objective, plus the
    # best on the GA's own multi-objective fitness (geometric mean of the scaled terms).
    span = config.ACTIVITY_HIGH - config.ACTIVITY_LOW
    fit = (np.clip((front["pred"] - config.ACTIVITY_LOW) / span, 0, 1)
           * front["qed"] * (10 - front["sa"]) / 9) ** (1 / 3)
    picks = [("highest predicted pActivity", front["pred"].idxmax()),
             ("highest QED", front["qed"].idxmax()),
             ("lowest SA score", front["sa"].idxmin()),
             ("best balanced (GA fitness)", fit.idxmax())]
    chosen, legends = [], []
    for label, i in picks:
        r = front.loc[i]
        chosen.append(r["smiles"])
        legends.append(f"{label}: pred {r['pred']:.2f}, QED {r['qed']:.2f}, SA {r['sa']:.2f}, "
                       f"sim {r['max_tanimoto']:.2f} (pctile {r['size_pctile']:.0f}), {r['arms']}")
    print("\nFront molecules drawn:")
    for s, lg in zip(chosen, legends):
        print(f"  {lg}\n    {s}")

    plot_fronts(by_arm, front, actives, reference,
                config.FIGURES_DIR / "06_pareto_fronts.png",
                "Pareto view of generated molecules: predicted pActivity (real surrogate), QED, SA score")
    draw_structures(chosen, legends, config.FIGURES_DIR / "06_pareto_structures.png")
    print("\nWrote results/06_pareto_front.csv, figures/06_pareto_fronts.png, figures/06_pareto_structures.png")


if __name__ == "__main__":
    main()
