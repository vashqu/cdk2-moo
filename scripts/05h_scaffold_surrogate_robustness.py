#!/usr/bin/env python3
"""
Stage 5h: is the Stage 5 picture stable if the surrogate is trained on a scaffold
split's training set (1,613 molecules) instead of all 2,016?

Reads : results/05_ga_populations.csv                    (primary: surrogate on all 2,016)
        results/05f_scaffold_split[_sNN]_populations.csv (secondary: five scaffold splits,
                                                           split seeds 42-46, 5 GA seeds each)
        results/04_surrogate_metrics.csv
Writes: results/05h_robustness_summary.csv

Same GA, arms, GA seeds and size window everywhere. In each experiment "similarity" is
ECFP4 (Morgan r=2, 2048-bit) max-Tanimoto to THAT experiment's training molecules, and real
predictions come from THAT experiment's surrogate, so absolute values are not comparable
across experiments; the qualitative claims below are.

Six claims were fixed before the extra splits (43-46) were run. Rule, also fixed in
advance: a claim is STABLE if it holds in the primary run and in at least 4 of the 5
scaffold splits; otherwise it is FRAGILE and is reported as such. Only generated
molecules (birth_generation > 0) of the final generation are used; differences are
paired by GA seed (same seed = same start molecules).

Run:
    python scripts/05h_scaffold_surrogate_robustness.py
"""

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.features import murcko_scaffold
from cdk2moo.objectives import ARMS

SPLIT_SEEDS = config.SPLIT_SEEDS[:5]
EXPERIMENTS = {"primary": "05_ga_populations.csv"}
for split_seed in SPLIT_SEEDS:
    suffix = "" if split_seed == SPLIT_SEEDS[0] else f"_s{split_seed}"
    EXPERIMENTS[f"split {split_seed}"] = f"05f_scaffold_split{suffix}_populations.csv"
SPLITS = [e for e in EXPERIMENTS if e != "primary"]
COLUMNS = ["pred_real", "max_tanimoto", "qed", "sa", "atoms", "n_scaffolds", "drift_lt_0.4"]


def main():
    config.require_campaign()
    # ---- one number per (experiment, arm, GA seed) -----------------------------------
    rows = []
    for exp, filename in EXPERIMENTS.items():
        pop = pd.read_csv(config.RESULTS_DIR / filename)
        final = pop[(pop["generation"] == pop["generation"].max()) & (pop["birth_generation"] > 0)]
        scaffold_of = {s: murcko_scaffold(s) for s in final["smiles"].unique()}
        for (arm, seed), g in final.groupby(["arm", "seed"]):
            rows.append({"experiment": exp, "arm": arm, "seed": seed,
                         "pred_real": g["pred_real"].median(), "max_tanimoto": g["max_tanimoto"].median(),
                         "qed": g["qed"].median(), "sa": g["sa"].median(), "atoms": g["n_heavy_atoms"].median(),
                         "n_scaffolds": g["smiles"].map(scaffold_of).nunique(),
                         "drift_lt_0.4": (g["max_tanimoto"] < 0.4).mean()})
    stats = pd.DataFrame(rows)
    stats.to_csv(config.RESULTS_DIR / "05h_robustness_summary.csv", index=False)

    values = {}                      # (experiment, arm, column) -> array over GA seeds
    for (exp, arm), d in stats.groupby(["experiment", "arm"]):
        d = d.sort_values("seed")
        for col in COLUMNS:
            values[(exp, arm, col)] = d[col].to_numpy(dtype=float)

    # ---- key statistics: primary vs the scaffold-trained splits --------------------------
    print("Generation 50, generated molecules.")
    print("  primary: median [min, max] over 5 GA seeds.")
    print(f"  scaffold-trained: median [min, max] over the {len(SPLITS)} splits of each split's median.")
    for col, label in [("pred_real", "median real prediction"), ("max_tanimoto", "median raw max Tanimoto"),
                       ("qed", "median QED"), ("sa", "median SA"), ("atoms", "median heavy atoms"),
                       ("n_scaffolds", "distinct scaffolds per 100"), ("drift_lt_0.4", "share with raw similarity < 0.4")]:
        print(f"\n  {label}")
        print(f"    {'arm':<16}{'primary':>26}{'scaffold-trained':>30}")
        for arm in ARMS:
            v = values[("primary", arm, col)]
            per_split = np.array([np.median(values[(sp, arm, col)]) for sp in SPLITS])
            print(f"    {arm:<16}{np.median(v):>8.2f} [{v.min():.2f}, {v.max():.2f}]"
                  f"{np.median(per_split):>12.2f} [{per_split.min():.2f}, {per_split.max():.2f}]")

    # ---- the six claims, per experiment -------------------------------------------------
    labels = ["multi_real raises real prediction above druglike_only in every GA seed",
              "multi_real raw similarity >= druglike_only (median paired difference)",
              "activity_only ends closer to training than druglike_only in every GA seed",
              "multi_scrambled beats druglike_only on REAL prediction in at most 2 GA seeds",
              "multi_real and activity_only keep fewer distinct scaffolds than druglike_only",
              "druglike_only drifts at least as far as activity_only (share < 0.4)"]
    holds = {}                       # (experiment, claim index) -> bool
    d_multi_sim_all = []             # paired multi_real - druglike similarity, every (split, GA seed)
    d_act_scaf_all = []
    for exp in EXPERIMENTS:
        base_pred = values[(exp, "druglike_only", "pred_real")]
        base_sim = values[(exp, "druglike_only", "max_tanimoto")]
        base_scaf = np.median(values[(exp, "druglike_only", "n_scaffolds")])
        base_drift = np.median(values[(exp, "druglike_only", "drift_lt_0.4")])
        d_multi_pred = values[(exp, "multi_real", "pred_real")] - base_pred
        d_multi_sim = values[(exp, "multi_real", "max_tanimoto")] - base_sim
        d_act_sim = values[(exp, "activity_only", "max_tanimoto")] - base_sim
        d_scr_pred = values[(exp, "multi_scrambled", "pred_real")] - base_pred
        n_multi = np.median(values[(exp, "multi_real", "n_scaffolds")])
        n_act = np.median(values[(exp, "activity_only", "n_scaffolds")])
        holds[(exp, 0)] = int((d_multi_pred > 0).sum()) == len(d_multi_pred)
        holds[(exp, 1)] = np.median(d_multi_sim) >= 0
        holds[(exp, 2)] = int((d_act_sim > 0).sum()) == len(d_act_sim)
        holds[(exp, 3)] = int((d_scr_pred > 0).sum()) <= 2
        holds[(exp, 4)] = n_multi < base_scaf and n_act < base_scaf
        holds[(exp, 5)] = base_drift >= np.median(values[(exp, "activity_only", "drift_lt_0.4")])
        if exp != "primary":
            d_multi_sim_all += list(d_multi_sim)
            d_act_scaf_all += list(values[(exp, "activity_only", "n_scaffolds")] - values[(exp, "druglike_only", "n_scaffolds")])

    print("\nDo the Stage 5 claims hold?  (paired by GA seed, arm minus druglike_only)")
    print(f"  {'claim':<78}{'primary':>9}{'splits holding':>16}{'verdict':>10}")
    stable_all = True
    for k, label in enumerate(labels):
        n_hold = sum(holds[(sp, k)] for sp in SPLITS)
        stable = holds[("primary", k)] and n_hold >= 4
        stable_all = stable_all and stable
        print(f"  {label:<78}{'holds' if holds[('primary', k)] else 'FAILS':>9}{n_hold:>13}/{len(SPLITS)}"
              f"{'STABLE' if stable else 'FRAGILE':>10}")

    print("\nThe two claims that were marginal with one split, per split:")
    print(f"  {'':<10}{'multi_real minus druglike, raw similarity':>44}{'scaffolds per 100: multi_real / activity_only / druglike':>58}")
    for exp in EXPERIMENTS:
        d = np.median(values[(exp, "multi_real", "max_tanimoto")] - values[(exp, "druglike_only", "max_tanimoto")])
        print(f"  {exp:<10}{d:>+44.3f}"
              f"{np.median(values[(exp, 'multi_real', 'n_scaffolds')]):>26.0f} /"
              f"{np.median(values[(exp, 'activity_only', 'n_scaffolds')]):>5.0f} /"
              f"{np.median(values[(exp, 'druglike_only', 'n_scaffolds')]):>5.0f}")
    q = np.percentile(d_multi_sim_all, [25, 50, 75])
    print(f"\n  Pooled over all {len(d_multi_sim_all)} (split, GA seed) pairs, multi_real minus druglike_only raw similarity: "
          f"median {q[1]:+.3f}, IQR [{q[0]:+.3f}, {q[2]:+.3f}]; positive in {int((np.array(d_multi_sim_all) > 0).sum())}/{len(d_multi_sim_all)}")
    q = np.percentile(d_act_scaf_all, [25, 50, 75])
    print(f"  Pooled activity_only minus druglike_only distinct scaffolds per 100: median {q[1]:+.0f}, IQR [{q[0]:+.0f}, {q[2]:+.0f}]; "
          f"negative (activity_only less diverse) in {int((np.array(d_act_scaf_all) < 0).sum())}/{len(d_act_scaf_all)}")

    print("\n  Verdict: " + ("every claim is STABLE under the pre-declared rule."
                            if stable_all else "at least one claim is FRAGILE; state it as fragile (see table)."))

    m = pd.read_csv(config.RESULTS_DIR / "04_surrogate_metrics.csv")
    acc = m[(m["seed"].isin(SPLIT_SEEDS)) & (m["split"] == "scaffold") & (m["model"] == "rf")]
    print(f"\n  Context: scaffold-trained surrogates' accuracy on their held-out molecules (Stage 4, seeds {SPLIT_SEEDS[0]}-{SPLIT_SEEDS[-1]}): "
          f"R2 median {acc['r2'].median():.2f} [{acc['r2'].min():.2f}, {acc['r2'].max():.2f}], "
          f"Spearman median {acc['spearman'].median():.2f}")
    print("\nWrote results/05h_robustness_summary.csv")


if __name__ == "__main__":
    main()
