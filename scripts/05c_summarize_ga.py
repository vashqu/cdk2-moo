#!/usr/bin/env python3
"""
Stage 5c: summarise the GA runs: trajectories, final-generation spread, and
paired contrasts against the druglike_only baseline.

Reads : results/05_ga_populations_scored.csv   (stage 5b: populations + size percentile)
Writes: results/05_ga_per_run.csv       <- per arm x seed x generation statistics
        results/05_ga_final_summary.csv <- final generation, per arm x seed
        figures/05_ga_trajectories.png
        figures/05_ga_final_generation.png

The contrast that matters is NOT "did similarity fall" (it falls in every arm,
including druglike_only, which never sees the surrogate) but whether an
activity-driven arm moves further than druglike_only, arm by arm and seed by
seed. All four arms of one seed share the same starting molecules, so that
comparison is paired.

Run:
    python scripts/05c_summarize_ga.py

Headline statistics describe GENERATED chemistry only: a molecule counts if its
birth_generation > 0. Starting molecules are training molecules (similarity 1.0,
in-sample predictions), so they are shown only as generation 0, the explicit
starting reference. From generation 1 on, each statistic is taken over the
population members born after generation 0.

Similarity is always ECFP4 (Morgan radius 2, 2048-bit, binary) Tanimoto to the
2,016 training molecules; size_pctile is that value as a percentile among
size-matched training molecules (leave-one-out reference; cdk2moo/similarity.py).
(The pilot run predates the birth_generation column, so this script does not
process it; its figures were made from all population members.)
"""

import argparse

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.ga_plots import plot_trajectories, plot_final_generation

STATS = ["median_pred_real", "best_pred_real", "median_pred_scrambled",
         "median_max_tanimoto", "median_size_pctile", "median_qed", "median_sa",
         "median_heavy_atoms"]


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", default="05_ga_populations_scored.csv")
    ap.add_argument("--tag", default="", help="suffix for output names")
    args = ap.parse_args()
    tag = f"_{args.tag}" if args.tag else ""

    pop = pd.read_csv(config.RESULTS_DIR / args.pop)
    last = int(pop["generation"].max())
    print(f"{args.pop}: {pop['arm'].nunique()} arms, {pop['seed'].nunique()} seeds, "
          f"{last} generations, {len(pop)} rows")

    n_start_left = (pop["birth_generation"] == 0) & (pop["generation"] > 0)
    print(f"population members that are still starting molecules after generation 0: "
          f"{int(n_start_left.sum())} of {int((pop['generation'] > 0).sum())} rows "
          f"(excluded from every statistic below)")
    pop = pop[(pop["generation"] == 0) | (pop["birth_generation"] > 0)]

    per_run = pop.groupby(["arm", "seed", "generation"]).agg(
        n_generated=("smiles", "size"),
        median_size_pctile=("size_pctile", "median"),
        best_pred_real=("pred_real", "max"),
        median_pred_real=("pred_real", "median"),
        median_pred_scrambled=("pred_scrambled", "median"),
        median_max_tanimoto=("max_tanimoto", "median"),
        median_qed=("qed", "median"),
        median_sa=("sa", "median"),
        median_heavy_atoms=("n_heavy_atoms", "median"),
        median_fitness=("fitness", "median"),
    ).reset_index()
    # The prediction of the surrogate each arm optimises: real for multi_real and
    # activity_only, scrambled for multi_scrambled, none for druglike_only.
    per_run["median_followed_pred"] = np.select(
        [per_run["arm"].isin(["multi_real", "activity_only"]), per_run["arm"] == "multi_scrambled"],
        [per_run["median_pred_real"], per_run["median_pred_scrambled"]], default=np.nan)
    per_run.to_csv(config.RESULTS_DIR / f"05_ga_per_run{tag}.csv", index=False)

    # Sanity: arms of one seed must start from identical populations.
    start = per_run[per_run["generation"] == 0].groupby("seed")["median_pred_real"].nunique()
    assert (start == 1).all(), "arms did not share a starting population"

    first = per_run[per_run["generation"] == 0].set_index(["arm", "seed"])
    final = per_run[per_run["generation"] == last].set_index(["arm", "seed"])
    final.reset_index().to_csv(config.RESULTS_DIR / f"05_ga_final_summary{tag}.csv", index=False)

    arms = [a for a in ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]
            if a in pop["arm"].unique()]

    # ---- how many generated molecules is each statistic based on? ---------
    g = pop[pop["generation"] == last]
    print(f"\nGeneration {last}: generated molecules (birth_generation > 0) per run, of {config.GA_POP_SIZE}")
    for arm in arms:
        a = g[g["arm"] == arm]
        n = a.groupby("seed").size()
        uniq = a["smiles"].nunique()
        print(f"  {arm:<16} per run min {n.min()}, median {int(n.median())}; "
              f"unique SMILES across 5 runs {uniq} of {len(a)}")

    # ---- final generation, per arm: median [min, max] over seeds ---------
    print(f"\nGeneration {last}: per-seed population statistic, median [min, max] over seeds")
    print(f"  {'':<16}" + "".join(f"{s.replace('median_', 'med ').replace('best_', 'best '):>24}" for s in STATS))
    for arm in arms:
        cells = []
        for s in STATS:
            v = final.loc[arm, s].to_numpy()
            cells.append(f"{np.median(v):.3f} [{v.min():.2f},{v.max():.2f}]")
        print(f"  {arm:<16}" + "".join(f"{c:>24}" for c in cells))

    # ---- change from start, per arm --------------------------------------
    print(f"\nChange from generation 0 to {last}, median [min, max] over seeds")
    changes = STATS
    print(f"  {'':<16}" + "".join(f"{s.replace('median_', 'med ').replace('best_', 'best '):>24}" for s in changes))
    for arm in arms:
        cells = []
        for s in changes:
            d = (final.loc[arm, s] - first.loc[arm, s]).to_numpy()
            cells.append(f"{np.median(d):+.3f} [{d.min():+.2f},{d.max():+.2f}]")
        print(f"  {arm:<16}" + "".join(f"{c:>24}" for c in cells))

    # ---- paired contrast against druglike_only ---------------------------
    # Same seed, same start molecules, same operators: the only difference is
    # the fitness function. Positive = the arm moved that statistic MORE than
    # druglike_only did.
    if "druglike_only" in arms:
        print(f"\nPaired contrast at generation {last}: arm minus druglike_only, same seed")
        print("  (how much of the movement is attributable to the activity objective)")
        cols = ["median_pred_real", "median_max_tanimoto", "median_size_pctile",
                "median_qed", "median_sa", "median_heavy_atoms"]
        print(f"  {'':<16}" + "".join(f"{c.replace('median_', 'med '):>24}" for c in cols)
              + f"{'seeds arm>base on pred_real':>30}")
        for arm in arms:
            if arm == "druglike_only":
                continue
            cells = []
            for c in cols:
                d = (final.loc[arm, c].sort_index() - final.loc["druglike_only", c].sort_index()).to_numpy()
                cells.append(f"{np.median(d):+.3f} [{d.min():+.2f},{d.max():+.2f}]")
            d = (final.loc[arm, "median_pred_real"].sort_index()
                 - final.loc["druglike_only", "median_pred_real"].sort_index()).to_numpy()
            print(f"  {arm:<16}" + "".join(f"{c:>24}" for c in cells)
                  + f"{int((d > 0).sum()):>27}/{len(d)}")

    # ---- figures -----------------------------------------------------------
    note = "20-39 heavy atoms"
    plot_trajectories(per_run, config.FIGURES_DIR / f"05_ga_trajectories{tag}.png",
                      f"GA trajectories, 5 seeds, {note}; generation 0 = starting molecules, later generations = molecules born after gen 0; median over seeds, band = min-max")
    plot_final_generation(final.reset_index(), config.FIGURES_DIR / f"05_ga_final_generation{tag}.png",
                          f"Generation {last}, one dot per GA seed, bar = median ({note})")
    print(f"\nWrote figures/05_ga_trajectories{tag}.png and figures/05_ga_final_generation{tag}.png")


if __name__ == "__main__":
    main()
