#!/usr/bin/env python3
"""
Stage 5: run the genetic algorithm, four arms x five seeds.

Reads : data/processed/cdk2_ic50_curated.csv, ecfp4.npy
Writes: results/05_ga_populations.csv  for the PRIMARY protocol only (legacy policy, 20-39 heavy atoms, default
        arms, seeds and generations), plus <output>.manifest.json (settings, package versions, platform, git commit).

Where output goes (see cdk2moo/run_outputs.py):
  * The primary file name is reserved. --relaxed, --policy corrected, or a non-default --arms/--seeds/--generations
    derive a different default name (e.g. 05_ga_relaxed_populations.csv, 05_ga_seeds-42_gen-3_populations.csv).
  * --out chooses a destination explicitly; an existing file is never replaced unless --overwrite.

Surrogates: ONE real random forest and ONE label-scrambled random forest, each trained once on all 2,016 curated
molecules with fixed seeds, reused by every GA seed and arm. GA seeds therefore vary only the search.

Arms (see objectives.fitness):
  multi_real       activity x QED x SA, real surrogate       <- the main arm
  multi_scrambled  activity x QED x SA, scrambled surrogate  <- H4 control
  activity_only    activity alone                            <- ablation
  druglike_only    QED x SA, no surrogate                    <- ablation

For one seed, all four arms start from the SAME 100 training molecules (drawn at random from those inside the size
window) and use identical operators and settings, so arms differ only in what they optimise.

Size window: primary runs use 20-39 heavy atoms (the central 5th-95th percentile of the training molecules' heavy-atom
counts; the training set itself spans 5-78). --relaxed uses 15-50 and is a stress test, never a primary result.

Preparation policy (see cdk2moo/candidates.py): "legacy" (default) scores edited molecules as they are, as in every
historical run. "corrected" standardizes each candidate like the training molecules, accepts only closed-shell
structures, and scores that same representation. It changes the search: it is a new experiment variant and its
outputs are never comparable with legacy outputs as if they were one experiment.

This stage only produces raw trajectories. Summary and figures are scripts/05b_size_conditioned_similarity.py then
scripts/05c_summarize_ga.py (they read the primary file); the H1-H4 analysis is a later stage.

Run:
    python scripts/05_run_ga.py                                  # primary protocol (refuses if the file exists)
    python scripts/05_run_ga.py --relaxed                        # writes 05_ga_relaxed_populations.csv
    python scripts/05_run_ga.py --policy corrected               # writes 05_ga_corrected_populations.csv
    python scripts/05_run_ga.py --seeds 42 --generations 3       # smoke test, derived non-primary name
    python scripts/05_run_ga.py --out my_run.csv --overwrite     # explicit destination
"""

import argparse

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.candidates import POLICIES
from cdk2moo.ga import run_ga
from cdk2moo.objectives import ARMS
from cdk2moo.run_outputs import choose_output, commit_csv, run_tokens, write_manifest
from cdk2moo.starts import load_eligible, select_start_population
from cdk2moo.surrogate import fit_forest, scramble_labels

PRIMARY_NAME = "05_ga_populations.csv"


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser(description="GA runs; see the file docstring for output naming and policies.")
    ap.add_argument("--seeds", type=int, nargs="+", default=config.GA_SEEDS)
    ap.add_argument("--generations", type=int, default=config.GA_N_GENERATIONS)
    ap.add_argument("--arms", nargs="+", default=ARMS, choices=ARMS)
    ap.add_argument("--relaxed", action="store_true", help="15-50 heavy-atom window: a stress test, not primary")
    ap.add_argument("--policy", choices=POLICIES, default=config.default_policy(),
                    help="candidate preparation policy (default: the scope's policy, else legacy)")
    ap.add_argument("--out", default=None, help="output CSV (default derived from the settings; see docstring)")
    ap.add_argument("--overwrite", action="store_true", help="allow replacing an existing output file")
    args = ap.parse_args()
    config.check_policy(args.policy)
    policy_token = "legacy" if args.policy == config.default_policy() else args.policy   # a scope's own policy needs no name token

    size_range = ((config.GA_RELAXED_MIN_HEAVY_ATOMS, config.GA_RELAXED_MAX_HEAVY_ATOMS) if args.relaxed
                  else (config.GA_MIN_HEAVY_ATOMS, config.GA_MAX_HEAVY_ATOMS))
    tokens = run_tokens(args.relaxed, policy_token, args.arms, ARMS, args.seeds, config.GA_SEEDS,
                        args.generations, config.GA_N_GENERATIONS)
    out_path = choose_output(config.RESULTS_DIR, PRIMARY_NAME, args.out, tokens, args.overwrite)   # fail before any work
    print(f"seeds = {args.seeds}, generations = {args.generations}, size window = {size_range}, "
          f"policy = {args.policy}, arms = {args.arms}")
    print(f"output: {out_path}")

    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    y = df["pactivity"].to_numpy()

    print(f"Training ONE real and ONE scrambled surrogate on all {len(df)} molecules "
          f"(seeds {config.SURROGATE_SEED} and {config.SCRAMBLE_SEED})")
    real = fit_forest(fps, y, config.SURROGATE_SEED)
    scrambled = fit_forest(fps, scramble_labels(y, config.SCRAMBLE_SEED), config.SCRAMBLE_SEED)

    eligible = load_eligible(df, args.relaxed)
    print(f"{len(eligible)} of {len(df)} training molecules are eligible as starting molecules (size window, unchanged by standardization, "
          "closed-shell, unique; the same list for both preparation policies)")

    runs, rejection_counts = [], {}
    for seed in args.seeds:
        _, start = select_start_population(eligible, np.arange(len(df)), df["std_smiles"].tolist(), config.GA_POP_SIZE, seed)
        for arm in args.arms:
            hist = run_ga(start, arm, real, scrambled, fps, seed, size_range, args.generations, policy=args.policy)
            runs.append(hist)
            rejection_counts[f"{arm}/seed{seed}"] = hist.attrs.get("rejects", {})
            last = hist[hist["generation"] == args.generations]
            first = hist[hist["generation"] == 0]
            print(f"  seed {seed}  {arm:<16} pred_real {first['pred_real'].mean():.2f} -> "
                  f"{last['pred_real'].mean():.2f}   max Tanimoto {first['max_tanimoto'].mean():.2f} -> "
                  f"{last['max_tanimoto'].mean():.2f}", flush=True)

    out = pd.concat(runs, ignore_index=True)
    commit_csv(out, out_path, index=False, float_format="%.5f")
    inputs = {"curated": config.PROCESSED_DIR / "cdk2_ic50_curated.csv", "fingerprints": config.PROCESSED_DIR / "ecfp4.npy",
              "eligible_starts": config.PROCESSED_DIR / "eligible_starts.csv"}
    manifest = write_manifest(out_path, {
        "seeds": args.seeds, "generations": args.generations, "arms": args.arms, "relaxed": args.relaxed,
        "size_window": list(size_range), "policy": args.policy, "population_size": config.GA_POP_SIZE,
        "children_per_generation": config.GA_N_CHILDREN, "mutation_probability": config.GA_MUTATION_PROB,
        "surrogate_seed": config.SURROGATE_SEED, "scramble_seed": config.SCRAMBLE_SEED,
        "rf_trees": config.RF_N_TREES, "rf_max_features": config.RF_MAX_FEATURES,
        "n_training_molecules": int(len(df)), "training_subset": "all curated molecules", "rejection_counts": rejection_counts},
        config.PROJECT_ROOT, inputs=inputs)
    print(f"\nWrote {out_path} ({len(out)} rows) and {manifest.name}")


if __name__ == "__main__":
    main()
