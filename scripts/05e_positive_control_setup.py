#!/usr/bin/env python3
"""
Stage 5e: define the two held-out positive controls and test the control
surrogates on the molecules they never saw, BEFORE any GA is run.

Reads : data/processed/cdk2_ic50_curated.csv, ecfp4.npy
Writes: results/05e_holdout_sets.csv        <- each molecule's role in each control (and in
                                               the scaffold-split robustness runs)
        results/05e_heldout_predictions.csv <- control surrogate on its held-out actives
        results/05e_setup_report.json

Selection rules are fixed in config and in cdk2moo/holdout.py (see there for the
reasoning); nothing here is chosen after seeing a GA result.

  molecule control: the 204 top-decile actives (measured pActivity >= the 90th
                    percentile) are removed from training; siblings stay.
  scaffold control: every scaffold family with >= 5 top-decile members is removed
                    entirely (all members, potent or not); other potent molecules stay.

Run:
    python scripts/05e_positive_control_setup.py
"""

import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import KFold

from cdk2moo import config
from cdk2moo.features import max_tanimoto, murcko_scaffold
from cdk2moo.holdout import molecule_holdout, scaffold_holdout, top_decile
from cdk2moo.similarity import loo_max_tanimoto, size_conditioned_percentile
from cdk2moo.surrogate import fit_forest


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=config.SURROGATE_SEED)
    args = ap.parse_args()
    print(f"seed = {args.seed} (surrogate and cross-validation folds)")

    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    y = df["pactivity"].to_numpy()
    sizes = df["n_heavy_atoms"].to_numpy()
    scaffolds = np.array([murcko_scaffold(s) for s in df["std_smiles"]])
    is_top, cutoff = top_decile(y, config.HOLDOUT_QUANTILE)
    print(f"{len(df)} molecules; top decile = measured pActivity >= {cutoff:.3f} ({int(is_top.sum())} molecules, "
          f"{len(set(scaffolds[is_top]))} scaffolds)")

    controls = {}
    tr, held, _ = molecule_holdout(y, config.HOLDOUT_QUANTILE)
    controls["molecule"] = (tr, held, None)
    tr, held, _, families = scaffold_holdout(y, scaffolds, config.HOLDOUT_QUANTILE,
                                             config.HOLDOUT_MIN_TOP_PER_FAMILY)
    controls["scaffold"] = (tr, held, families)

    roles = pd.DataFrame({"row": np.arange(len(df)), "inchikey": df["inchikey"],
                          "scaffold": scaffolds, "pactivity": y, "is_top_decile": is_top})
    # Not a positive control: the training set of the Stage 3 scaffold split (default
    # seed). Recorded here so the secondary robustness GA (05f --controls scaffold_split)
    # trains its surrogate on exactly that split's training molecules.
    split = pd.read_csv(config.PROCESSED_DIR / "cdk2_splits.csv")
    split = split[split["seed"] == config.RANDOM_SEED].sort_values("row")
    roles["scaffold_split_control"] = np.where(split["scaffold_split"] == "test", "held", "train")
    # The same for the next four predetermined split seeds (43-46), so the robustness
    # run can use five splits, one per GA seed, as the primary run does.
    all_splits = pd.read_csv(config.PROCESSED_DIR / "cdk2_splits.csv")
    for extra_seed in config.SPLIT_SEEDS[1:5]:
        one = all_splits[all_splits["seed"] == extra_seed].sort_values("row")
        roles[f"scaffold_split_s{extra_seed}_control"] = np.where(one["scaffold_split"] == "test", "held", "train")
    report = {"top_decile_cutoff": cutoff, "n_top_decile": int(is_top.sum())}
    prediction_rows = []

    for name, (train, held, families) in controls.items():
        held_top = held[is_top[held]]                       # the held-out ACTIVES
        roles[f"{name}_control"] = "train"
        roles.loc[held, f"{name}_control"] = "held"
        print(f"\n===== {name.upper()}-LEVEL CONTROL =====")
        print(f"  training {len(train)}; held out {len(held)} molecules, of which {len(held_top)} are top-decile actives")
        print(f"  held-out actives: {len(set(scaffolds[held_top]))} distinct scaffolds")
        if families is not None:
            print(f"  held-out scaffold FAMILIES (>= {config.HOLDOUT_MIN_TOP_PER_FAMILY} top-decile members each): {len(families)}; "
                  f"{len(held) - len(held_top)} held-out members are not top decile")
            for scaf, (n, nt) in sorted(families.items(), key=lambda kv: -kv[1][1]):
                print(f"    {nt:>3} actives / {n:>3} members  {scaf[:66]}")
        q = np.percentile(y[held_top], [0, 25, 50, 75, 100])
        print(f"  held-out actives' measured pActivity: min {q[0]:.2f}, quartiles {q[1]:.2f}/{q[2]:.2f}/{q[3]:.2f}, max {q[4]:.2f}")
        print(f"  remaining training labels: max {y[train].max():.2f}, 90th percentile {np.percentile(y[train], 90):.2f}, "
              f"top-decile actives still in training: {int(is_top[train].sum())}")

        # ---- how supported are the held-out actives, after the holdout? -------
        sim = max_tanimoto(fps[held_top], fps[train])
        loo = loo_max_tanimoto(fps[train])
        pct, _, _ = size_conditioned_percentile(sim, sizes[held_top], sizes[train], loo)
        print(f"  nearest-training ECFP4 Tanimoto of held-out actives: median {np.median(sim):.2f}, "
              f"quartiles {np.percentile(sim, 25):.2f}/{np.percentile(sim, 75):.2f}; "
              f"{100 * (sim >= 0.6).mean():.0f}% >= 0.6, {100 * (sim >= 0.8).mean():.0f}% >= 0.8; "
              f"median size percentile {np.median(pct):.0f}")

        # ---- the control surrogate on the molecules it never saw --------------
        forest = fit_forest(fps[train], y[train], args.seed)
        pred = forest.predict(fps[held_top])
        err = pred - y[held_top]
        rho = spearmanr(y[held_top], pred)[0]
        # Out-of-fold predictions for the remaining molecules: a fair reference for "a typical molecule".
        oof = np.zeros(len(train))
        for fit_rows, test_rows in KFold(5, shuffle=True, random_state=args.seed).split(train):
            f = fit_forest(fps[train][fit_rows], y[train][fit_rows], args.seed)
            oof[test_rows] = f.predict(fps[train][test_rows])
        labels = np.r_[np.ones(len(pred)), np.zeros(len(oof))]
        auc = roc_auc_score(labels, np.r_[pred, oof])
        print(f"  control surrogate on held-out actives (out-of-sample): MAE {np.abs(err).mean():.2f}, "
              f"RMSE {np.sqrt((err ** 2).mean()):.2f}, bias {err.mean():+.2f} (negative = under-predicts), "
              f"Spearman within the narrow held-out range {rho:.2f}")
        print(f"    predicted pActivity of held-out actives: median {np.median(pred):.2f} (measured median {np.median(y[held_top]):.2f}); "
              f"{100 * (pred >= 8).mean():.0f}% predicted >= 8, {100 * (pred >= 7).mean():.0f}% >= 7")
        print(f"    discrimination: AUC {auc:.2f} for ranking held-out actives above 5-fold out-of-fold predictions of the "
              f"remaining molecules; {100 * (pred > np.percentile(oof, 75)).mean():.0f}% of held-out actives exceed the remainder's 75th percentile")

        prediction_rows.append(pd.DataFrame({"control": name, "row": held_top, "pactivity": y[held_top],
                                             "pred": pred, "max_tanimoto_to_train": sim, "size_pctile": pct}))
        report[name] = {
            "n_train": int(len(train)), "n_held": int(len(held)), "n_held_actives": int(len(held_top)),
            "n_held_scaffolds": int(len(set(scaffolds[held_top]))),
            "n_families": None if families is None else len(families),
            "held_pactivity_quartiles": q.tolist(), "train_label_max": float(y[train].max()),
            "nn_train_similarity_median": float(np.median(sim)),
            "frac_nn_ge_0.6": float((sim >= 0.6).mean()), "frac_nn_ge_0.8": float((sim >= 0.8).mean()),
            "surrogate_mae": float(np.abs(err).mean()), "surrogate_bias": float(err.mean()),
            "surrogate_auc_vs_remainder": float(auc)}

    roles.to_csv(config.RESULTS_DIR / "05e_holdout_sets.csv", index=False)
    pd.concat(prediction_rows).to_csv(config.RESULTS_DIR / "05e_heldout_predictions.csv", index=False, float_format="%.4f")
    with open(config.RESULTS_DIR / "05e_setup_report.json", "w") as fh:
        json.dump(report, fh, indent=2)
    print("\nWrote results/05e_holdout_sets.csv, 05e_heldout_predictions.csv, 05e_setup_report.json")
    print("Next: scripts/05f_positive_control_ga.py")


if __name__ == "__main__":
    main()
