#!/usr/bin/env python3
"""
Stage 4: train the random-forest surrogate on each split, with controls.

Reads : data/processed/cdk2_ic50_curated.csv, ecfp4.npy, cdk2_splits.csv
Writes: results/04_surrogate_metrics.csv     <- one row per seed x split x model
        results/04_test_predictions.csv      <- one row per test molecule per
                                                seed x split, with its
                                                max-Tanimoto to that train set

Four models are scored on every test set, so no number appears without a control:

  rf          the surrogate: random forest on ECFP4
  scrambled   same forest, trained on shuffled training labels. Should score ~0.
  mean        predicts the training-set mean for everything. The R2 = 0 line.
  size_only   forest on heavy-atom count alone. If it beats "scrambled" by much,
              part of what the surrogate 'learns' is just molecule size
              (bigger molecules tend to bind tighter) and not chemistry.

The per-molecule predictions file is what the applicability-domain audit
(stage 7) will analyse: error as a function of similarity to the training set.

Run:
    python scripts/04_train_surrogate.py
    python scripts/04_train_surrogate.py --seeds 42 43
"""

import argparse

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.features import max_tanimoto
from cdk2moo.surrogate import fit_forest, scramble_labels, regression_metrics

SPLITS = ["random", "scaffold", "paper"]


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, nargs="+", default=config.SPLIT_SEEDS)
    args = ap.parse_args()
    print(f"seeds = {args.seeds}")

    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    splits = pd.read_csv(config.PROCESSED_DIR / "cdk2_splits.csv")
    y = df["pactivity"].to_numpy()
    size = df[["n_heavy_atoms"]].to_numpy(dtype=float)

    metric_rows = []
    pred_rows = []

    for seed in args.seeds:
        this = splits[splits["seed"] == seed].sort_values("row")
        assert len(this) == len(df), f"no splits for seed {seed}; rerun stage 3"
        assert (this["inchikey"].to_numpy() == df["inchikey"].to_numpy()).all()

        for split in SPLITS:
            is_train = (this[f"{split}_split"] == "train").to_numpy()
            train, test = np.where(is_train)[0], np.where(~is_train)[0]

            rf = fit_forest(fps[train], y[train], seed)
            scrambled = fit_forest(fps[train], scramble_labels(y[train], seed), seed)
            # One input feature, so max_features must be 1.0 (all of it).
            size_only = fit_forest(size[train], y[train], seed, max_features=1.0)

            preds = {
                "rf": rf.predict(fps[test]),
                "scrambled": scrambled.predict(fps[test]),
                "mean": np.full(len(test), y[train].mean()),
                "size_only": size_only.predict(size[test]),
            }
            for model, p in preds.items():
                m = regression_metrics(y[test], p)
                metric_rows.append({"seed": seed, "split": split, "model": model,
                                    "n_test": len(test), **m})

            pred_rows.append(pd.DataFrame({
                "seed": seed, "split": split, "row": test,
                "y_true": y[test],
                "pred_rf": preds["rf"],
                "pred_scrambled": preds["scrambled"],
                "pred_size_only": preds["size_only"],
                "max_tanimoto_to_train": max_tanimoto(fps[test], fps[train]),
            }))
        print(f"  seed {seed} done")

    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(config.RESULTS_DIR / "04_surrogate_metrics.csv", index=False)
    pd.concat(pred_rows).to_csv(config.RESULTS_DIR / "04_test_predictions.csv",
                                index=False)

    # ---- summary: median [IQR] over seeds ---------------------------
    print(f"\nTest-set metrics, median [25th, 75th percentile] over {len(args.seeds)} seeds")
    for split in SPLITS:
        print(f"\n  {split} split")
        print(f"    {'model':<11}{'R2':>24}{'RMSE':>24}{'Spearman':>24}")
        for model in ["rf", "scrambled", "mean", "size_only"]:
            sub = metrics[(metrics["split"] == split) & (metrics["model"] == model)]
            cells = []
            for col in ["r2", "rmse", "spearman"]:
                v = sub[col].dropna()
                if len(v):
                    q1, med, q3 = np.percentile(v, [25, 50, 75])
                    cells.append(f"{med:.3f} [{q1:.3f}, {q3:.3f}]")
                else:
                    cells.append("n/a")
            print(f"    {model:<11}" + "".join(f"{c:>24}" for c in cells))

    print("\nWrote results/04_surrogate_metrics.csv and results/04_test_predictions.csv")
    print("Next: inspect. Stage 5 (GA) needs a decision on which split trains the GA surrogate.")


if __name__ == "__main__":
    main()
