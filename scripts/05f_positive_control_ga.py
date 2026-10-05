#!/usr/bin/env python3
"""
Stage 5f: run the Stage 5 GA unchanged, but with a surrogate trained WITHOUT the
held-out molecules of each positive control.

Reads : data/processed/cdk2_ic50_curated.csv, ecfp4.npy, results/05e_holdout_sets.csv
Writes: results/05f_<control>_populations.csv   (control = molecule, scaffold: the positive
        controls; scaffold_split, scaffold_split_s43..s46: the robustness runs on the
        Stage 3 scaffold splits' training sets, evaluated by 05h)

Same GA, arms, seeds, size window (20-39 heavy atoms), generations and
operators as scripts/05_run_ga.py. What changes: the real and scrambled forests
are trained on the control's remaining molecules only, starting molecules are drawn
from those, and "similarity to the training set" in the output means similarity
to those remaining molecules. Held-out molecules are used only afterwards, in 05g.

Run:
    python scripts/05f_positive_control_ga.py
    python scripts/05f_positive_control_ga.py --controls scaffold --seeds 42 --generations 3   # smoke test (derived name)
    python scripts/05f_positive_control_ga.py --policy corrected   # corrected variant, separate output names
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


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls", nargs="+", default=["molecule", "scaffold"])
    ap.add_argument("--seeds", type=int, nargs="+", default=config.GA_SEEDS)
    ap.add_argument("--generations", type=int, default=config.GA_N_GENERATIONS)
    ap.add_argument("--tag", default="", help="extra label added to the derived output name")
    ap.add_argument("--policy", choices=POLICIES, default=config.default_policy(),
                    help="candidate preparation policy (default: the scope's policy, else legacy)")
    ap.add_argument("--overwrite", action="store_true", help="allow replacing an existing output file")
    args = ap.parse_args()
    config.check_policy(args.policy)
    policy_token = "legacy" if args.policy == config.default_policy() else args.policy   # a scope's own policy needs no name token

    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    y = df["pactivity"].to_numpy()
    roles = pd.read_csv(config.RESULTS_DIR / "05e_holdout_sets.csv")
    size_range = (config.GA_MIN_HEAVY_ATOMS, config.GA_MAX_HEAVY_ATOMS)

    for control in args.controls:
        tokens = run_tokens(False, policy_token, None, None, args.seeds, config.GA_SEEDS, args.generations,
                            config.GA_N_GENERATIONS, extra=[args.tag] if args.tag else [])
        out_path = choose_output(config.RESULTS_DIR, f"05f_{control}_populations.csv", None, tokens, args.overwrite)
        print(f"output: {out_path}")
        train = np.where(roles[f"{control}_control"] == "train")[0]
        print(f"\n{control} control: surrogate trained on {len(train)} molecules "
              f"(seeds {config.SURROGATE_SEED} and {config.SCRAMBLE_SEED}), size window {size_range}")
        real = fit_forest(fps[train], y[train], config.SURROGATE_SEED)
        scrambled = fit_forest(fps[train], scramble_labels(y[train], config.SCRAMBLE_SEED),
                               config.SCRAMBLE_SEED)
        eligible = load_eligible(df)

        runs, rejection_counts = [], {}
        for seed in args.seeds:
            _, start = select_start_population(eligible, train, df["std_smiles"].tolist(), config.GA_POP_SIZE, seed)
            for arm in ARMS:
                hist = run_ga(start, arm, real, scrambled, fps[train], seed, size_range,
                              args.generations, policy=args.policy)
                rejection_counts[f"{arm}/seed{seed}"] = hist.attrs.get("rejects", {})
                runs.append(hist)
                last = hist[hist["generation"] == args.generations]
                print(f"  seed {seed}  {arm:<16} pred_real -> {last['pred_real'].mean():.2f}, "
                      f"similarity to remaining training -> {last['max_tanimoto'].mean():.2f}", flush=True)
        out = pd.concat(runs, ignore_index=True)
        commit_csv(out, out_path, index=False, float_format="%.5f")
        inputs = {"curated": config.PROCESSED_DIR / "cdk2_ic50_curated.csv", "fingerprints": config.PROCESSED_DIR / "ecfp4.npy",
                  "eligible_starts": config.PROCESSED_DIR / "eligible_starts.csv", "holdout_assignments": config.RESULTS_DIR / "05e_holdout_sets.csv"}
        write_manifest(out_path, {"control": control, "seeds": args.seeds, "generations": args.generations,
                                  "policy": args.policy, "size_window": list(size_range), "n_training_molecules": int(len(train)),
                                  "surrogate_seed": config.SURROGATE_SEED, "scramble_seed": config.SCRAMBLE_SEED,
                                  "training_subset": f"rows with role 'train' in column {control}_control of 05e_holdout_sets.csv",
                                  "rejection_counts": rejection_counts}, config.PROJECT_ROOT, inputs=inputs)
        print(f"Wrote {out_path} ({len(out)} rows) and its manifest")


if __name__ == "__main__":
    main()
