#!/usr/bin/env python3
"""
Stage 10a: the H5 mitigation arm. Rerun the GA with a hard similarity floor to the training set.

Reads : data/processed/cdk2_ic50_curated.csv, ecfp4.npy, results/05e_holdout_sets.csv
Writes: results/10_constrained_<experiment>_populations.csv   (experiment = primary, molecule, scaffold)

For each experiment, arm in {multi_real, activity_only}, floor tau in {0.4, 0.5, 0.6, 0.7} and
GA seed 42-46: the same GA as Stage 5 (same surrogates, start molecules, operators and size
window), except that any molecule whose raw ECFP4 max-Tanimoto to the experiment's training set
is below tau has fitness 0 and cannot survive. The unconstrained twin of every run already exists
(Stage 5 and 5f), so comparisons are paired by seed.

Experiments:
  primary   surrogate trained on all 2,016 molecules (D-18): measures how much predicted activity
            the floor costs.
  molecule, scaffold   surrogate trained without the held-out actives (D-26): measures whether
            the floor helps recover them.
Only the columns needed downstream are saved.

Run:
    python scripts/10_run_constrained_ga.py
    python scripts/10_run_constrained_ga.py --experiments scaffold --taus 0.6 --seeds 42 --generations 3   # smoke test (derived name)
    python scripts/10_run_constrained_ga.py --policy corrected   # corrected variant, separate output names
"""

import argparse

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.candidates import POLICIES
from cdk2moo.ga import run_ga
from cdk2moo.run_outputs import choose_output, commit_csv, run_tokens, write_manifest
from cdk2moo.starts import load_eligible, select_start_population
from cdk2moo.surrogate import fit_forest, scramble_labels

DEFAULT_TAUS = [0.4, 0.5, 0.6, 0.7]
DEFAULT_ARMS = ["multi_real", "activity_only"]

KEEP = ["smiles", "pred_real", "qed", "sa", "max_tanimoto", "n_heavy_atoms", "fitness",
        "birth_generation", "generation", "arm", "seed", "min_similarity"]


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--experiments", nargs="+", default=["primary", "molecule", "scaffold"])
    ap.add_argument("--taus", type=float, nargs="+", default=DEFAULT_TAUS)
    ap.add_argument("--arms", nargs="+", default=DEFAULT_ARMS)
    ap.add_argument("--seeds", type=int, nargs="+", default=config.GA_SEEDS)
    ap.add_argument("--generations", type=int, default=config.GA_N_GENERATIONS)
    ap.add_argument("--tag", default="", help="extra label added to the derived output name")
    ap.add_argument("--policy", choices=POLICIES, default=config.default_policy(),
                    help="candidate preparation policy (default: the scope's policy, else legacy)")
    ap.add_argument("--overwrite", action="store_true", help="allow replacing an existing output file")
    args = ap.parse_args()
    config.check_policy(args.policy)
    policy_token = "legacy" if args.policy == config.default_policy() else args.policy   # a scope's own policy needs no name token
    extra = ([] if args.taus == DEFAULT_TAUS else ["floors-" + "+".join(str(t) for t in args.taus)]) + ([args.tag] if args.tag else [])
    tokens = run_tokens(False, policy_token, args.arms, DEFAULT_ARMS, args.seeds, config.GA_SEEDS, args.generations,
                        config.GA_N_GENERATIONS, extra=extra)
    out_paths = {e: choose_output(config.RESULTS_DIR, f"10_constrained_{e}_populations.csv", None, tokens, args.overwrite)
                 for e in args.experiments}                       # fail before any work if a destination is taken

    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    y = df["pactivity"].to_numpy()
    roles = pd.read_csv(config.RESULTS_DIR / "05e_holdout_sets.csv")
    size_range = (config.GA_MIN_HEAVY_ATOMS, config.GA_MAX_HEAVY_ATOMS)
    sizes = df["n_heavy_atoms"].to_numpy()

    for experiment in args.experiments:
        if experiment == "primary":
            train = np.arange(len(df))
        else:
            train = np.where(roles[f"{experiment}_control"] == "train")[0]
        print(f"\n{experiment}: surrogate trained on {len(train)} molecules; floors {args.taus}", flush=True)
        real = fit_forest(fps[train], y[train], config.SURROGATE_SEED)
        scrambled = fit_forest(fps[train], scramble_labels(y[train], config.SCRAMBLE_SEED), config.SCRAMBLE_SEED)
        eligible = load_eligible(df)

        runs, rejection_counts = [], {}
        for seed in args.seeds:
            _, start = select_start_population(eligible, train, df["std_smiles"].tolist(), config.GA_POP_SIZE, seed)
            for arm in args.arms:
                for tau in args.taus:
                    hist = run_ga(start, arm, real, scrambled, fps[train], seed, size_range,
                                  args.generations, min_similarity=tau, policy=args.policy)
                    rejection_counts[f"{arm}/tau{tau}/seed{seed}"] = hist.attrs.get("rejects", {})
                    runs.append(hist[KEEP])
                    last = hist[(hist["generation"] == args.generations) & (hist["birth_generation"] > 0)]
                    print(f"  seed {seed} {arm:<14} tau {tau}: {len(last):>3} generated; pred_real {last['pred_real'].median():.2f}, "
                          f"similarity {last['max_tanimoto'].median():.2f}", flush=True)
        out = pd.concat(runs, ignore_index=True)
        commit_csv(out, out_paths[experiment], index=False, float_format="%.5f")
        inputs = {"curated": config.PROCESSED_DIR / "cdk2_ic50_curated.csv", "fingerprints": config.PROCESSED_DIR / "ecfp4.npy",
                  "eligible_starts": config.PROCESSED_DIR / "eligible_starts.csv", "holdout_assignments": config.RESULTS_DIR / "05e_holdout_sets.csv"}
        write_manifest(out_paths[experiment], {"experiment": experiment, "seeds": args.seeds, "generations": args.generations,
                                               "arms": args.arms, "floors": args.taus, "policy": args.policy,
                                               "size_window": list(size_range), "n_training_molecules": int(len(train)),
                                               "training_subset": ("all curated molecules" if experiment == "primary" else f"rows with role 'train' in column {experiment}_control"),
                                               "rejection_counts": rejection_counts}, config.PROJECT_ROOT, inputs=inputs)
        print(f"Wrote {out_paths[experiment]} ({len(out)} rows) and its manifest")


if __name__ == "__main__":
    main()
