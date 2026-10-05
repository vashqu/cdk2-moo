#!/usr/bin/env python3
"""
Stage 3: fingerprints, scaffolds, and the random, scaffold and paper splits.

Reads : data/processed/cdk2_ic50_curated.csv
Writes: data/processed/ecfp4.npy          <- fingerprints, row i = curated row i
        data/processed/cdk2_splits.csv    <- LONG format: one row per molecule per
                                             seed, with which side of each of the
                                             three splits it is on
        results/03_splits.json            <- overlap statistics, every seed

Every later stage must load the split from cdk2_splits.csv and never re-split,
so that all surrogates and all GA arms are judged on identical test sets.

Run:
    python scripts/03_features_and_splits.py
    python scripts/03_features_and_splits.py --seeds 42 43
"""

import argparse
import json

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.features import ecfp4, murcko_scaffold
from cdk2moo.splits import random_split, group_split, overlap_report


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=config.SPLIT_SEEDS,
                    help="split seeds (default: the predetermined list in config)")
    args = ap.parse_args()
    print(f"seeds = {args.seeds}")

    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    print(f"Loaded {len(df)} molecules")

    # ---- features ----------------------------------------------------
    fps = ecfp4(df["std_smiles"], config.FP_RADIUS, config.FP_BITS)
    np.save(config.PROCESSED_DIR / "ecfp4.npy", fps)
    bits_on = fps.sum(axis=1)
    print(f"ECFP4: {fps.shape}, on-bits per molecule: median {np.median(bits_on):.0f}, "
          f"min {bits_on.min()}, max {bits_on.max()}")

    scaffolds = [murcko_scaffold(s) for s in df["std_smiles"]]
    df["scaffold"] = scaffolds

    # ---- scaffold statistics ----------------------------------------
    sizes = df["scaffold"].value_counts()
    n_acyclic = int((df["scaffold"] == "").sum())
    print(f"\nScaffolds: {len(sizes)} unique for {len(df)} molecules")
    print(f"  singletons (scaffold used by one molecule): {int((sizes == 1).sum())}")
    print(f"  acyclic molecules (empty scaffold, lumped into one group): {n_acyclic}")
    print("  largest scaffold groups:")
    for scaffold, n in sizes.head(5).items():
        print(f"    {n:>4}  {scaffold if scaffold else '(empty)'}")

    primary_docs = df["primary_document"].to_numpy()
    all_docs = df["all_documents"].to_numpy()
    pactivity = df["pactivity"].to_numpy()

    report = {
        "seeds": args.seeds,
        "n_molecules": int(len(df)),
        "n_scaffolds": int(len(sizes)),
        "n_singleton_scaffolds": int((sizes == 1).sum()),
        "n_acyclic": n_acyclic,
        "largest_scaffold_size": int(sizes.iloc[0]),
        "by_seed": {},
    }
    long_frames = []

    # ---- one set of three splits per seed ---------------------------
    for seed in args.seeds:
        splits = {
            "random": random_split(len(df), config.TEST_FRAC, seed),
            "scaffold": group_split(scaffolds, config.TEST_FRAC, seed),
            "paper": group_split(primary_docs, config.TEST_FRAC, seed),
        }
        out = pd.DataFrame({"row": np.arange(len(df)), "inchikey": df["inchikey"],
                            "scaffold": scaffolds, "seed": seed})
        report["by_seed"][seed] = {}

        for name, (train, test) in splits.items():
            # A molecule must not be on both sides.
            assert not set(train) & set(test)
            assert len(train) + len(test) == len(df)
            out[f"{name}_split"] = "train"
            out.loc[test, f"{name}_split"] = "test"
            report["by_seed"][seed][name] = overlap_report(
                train, test, scaffolds, fps, primary_docs, all_docs, pactivity)

        # The defining property of a scaffold split: no scaffold on both sides.
        train, test = splits["scaffold"]
        shared = {scaffolds[i] for i in train} & {scaffolds[i] for i in test}
        assert not shared, f"scaffold split leaks (seed {seed}): {shared}"

        train, test = splits["paper"]
        shared = {primary_docs[i] for i in train} & {primary_docs[i] for i in test}
        assert not shared, f"paper split leaks (seed {seed}): {shared}"

        long_frames.append(out)

    # ---- how different are the splits, across seeds? ----------------
    # Median [min, max] over seeds. The spread matters: the paper split's label
    # shift (test pActivity mean minus train) varies with the draw.
    print("\nTest set vs train set, median [min, max] over seeds:")
    print(f"  {'':<44}{'random':>22}{'scaffold':>22}{'paper':>22}")
    stats = {
        "n_test": lambda r: r["n_test"],
        "test molecules with scaffold in train": lambda r: r["test_molecules_with_scaffold_in_train"],
        "test molecules with ANY paper in train": lambda r: r["test_molecules_with_any_paper_in_train"],
        "median max Tanimoto to train": lambda r: r["median_max_tanimoto_to_train"],
        "frac test with max Tanimoto > 0.7": lambda r: r["frac_test_with_nn_above_0.7"],
        "pActivity shift (test mean - train mean)": lambda r: r["pactivity_mean_test"] - r["pactivity_mean_train"],
    }
    for label, get in stats.items():
        cells = []
        for name in ("random", "scaffold", "paper"):
            v = np.array([get(report["by_seed"][sd][name]) for sd in args.seeds])
            cells.append(f"{np.median(v):.3g} [{v.min():.3g}, {v.max():.3g}]")
        print(f"  {label:<44}" + "".join(f"{c:>22}" for c in cells))

    # ---- save --------------------------------------------------------
    pd.concat(long_frames).to_csv(config.PROCESSED_DIR / "cdk2_splits.csv", index=False)
    with open(config.RESULTS_DIR / "03_splits.json", "w") as fh:
        json.dump(report, fh, indent=2)
    print("\nWrote cdk2_splits.csv, ecfp4.npy, results/03_splits.json")
    print("Next: stage 4, surrogate.")


if __name__ == "__main__":
    main()
