#!/usr/bin/env python3
"""
Stage 7a: how does the surrogate's error depend on similarity to its training set?

Reads : results/04_test_predictions.csv (out-of-sample predictions, 3 splits x 10 seeds),
        data/processed/cdk2_splits.csv, cdk2_ic50_curated.csv, ecfp4.npy
Writes: results/07_stage4_ad.csv         <- Stage 4 predictions + size percentile
        results/07_reliability_bins.csv  <- error per region, both measures
        figures/07_reliability_curves.png

Both similarity measures are computed against the training set of the split the
molecule was tested on: raw ECFP4 max-Tanimoto (already in the Stage 4 file) and
the size-conditioned percentile, whose reference is that split's OWN training
molecules, leave-one-out.

Error is reported only where at least config.AD_MIN_UNIQUE distinct molecules
back it. The same molecule is in the test set of several seeds, so rows are not
independent evidence; the distinct-molecule count is what limits the claim.

Run:
    python scripts/07_ad_reliability.py
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cdk2moo import config
from cdk2moo.applicability import add_region, error_table, large_error_auc
from cdk2moo.ad_plots import plot_reliability
from cdk2moo.similarity import loo_max_tanimoto, size_conditioned_percentile


def main():
    config.require_campaign()
    p4 = pd.read_csv(config.RESULTS_DIR / "04_test_predictions.csv")
    splits = pd.read_csv(config.PROCESSED_DIR / "cdk2_splits.csv")
    mols = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    sizes = mols["n_heavy_atoms"].to_numpy()

    # ---- size percentile for every Stage 4 test molecule ---------------------
    parts = []
    for (seed, split), g in p4.groupby(["seed", "split"]):
        this = splits[splits["seed"] == seed].sort_values("row")
        train = np.where((this[f"{split}_split"] == "train").to_numpy())[0]
        loo = loo_max_tanimoto(fps[train])                        # train vs train, no self
        g = g.copy()
        g["n_heavy_atoms"] = sizes[g["row"].to_numpy()]
        g["size_pctile"], g["stratum_halfwidth"], _ = size_conditioned_percentile(
            g["max_tanimoto_to_train"], g["n_heavy_atoms"], sizes[train], loo)
        parts.append(g)
    p4 = pd.concat(parts, ignore_index=True)
    p4["mol"] = p4["row"]                          # a molecule's identity = its curated row
    p4["err"] = p4["pred_rf"] - p4["y_true"]
    p4["abs_err"] = p4["err"].abs()
    p4 = add_region(p4, "max_tanimoto_to_train", config.AD_RAW_EDGES, config.AD_RAW_NAMES, "raw_region")
    p4 = add_region(p4, "size_pctile", config.AD_PCT_EDGES, config.AD_PCT_NAMES, "pct_region")
    p4.to_csv(config.RESULTS_DIR / "07_stage4_ad.csv", index=False, float_format="%.5f")
    print(f"{len(p4)} Stage 4 test predictions, {p4['mol'].nunique()} distinct molecules; "
          f"stratum widened beyond +/-{config.SIZE_STRATUM_HALFWIDTH} for "
          f"{int((p4['stratum_halfwidth'] > config.SIZE_STRATUM_HALFWIDTH).sum())} rows")

    # ---- does either measure track error? ------------------------------------
    print("\nRank correlation of |error| with similarity (negative = error rises as similarity falls),")
    print(f"and AUC for flagging |error| >= {config.AD_LARGE_ERROR} (0.5 = useless, 1.0 = perfect):")
    print(f"  {'':<10}{'rho raw':>10}{'rho pctile':>12}{'AUC raw':>10}{'AUC pctile':>12}{'large-error share':>20}")
    for name in ["random", "scaffold", "paper", "pooled"]:
        d = p4 if name == "pooled" else p4[p4["split"] == name]
        rho_raw = spearmanr(d["abs_err"], d["max_tanimoto_to_train"])[0]
        rho_pct = spearmanr(d["abs_err"], d["size_pctile"])[0]
        share = (d["abs_err"] >= config.AD_LARGE_ERROR).mean()
        print(f"  {name:<10}{rho_raw:>10.2f}{rho_pct:>12.2f}"
              f"{large_error_auc(d, 'max_tanimoto_to_train'):>10.3f}"
              f"{large_error_auc(d, 'size_pctile'):>12.3f}{100 * share:>19.0f}%")

    # ---- error by region, each measure ----------------------------------------
    outputs = []
    for label, col in [("raw similarity region", "raw_region"), ("size-percentile region", "pct_region")]:
        t = error_table(p4, [col])
        t.insert(0, "measure", label)
        outputs.append(t.rename(columns={col: "region"}))
        print(f"\nStage 4 error by {label} (pooled over splits and seeds):")
        print(f"  {'region':<26}{'rows':>6}{'distinct mols':>15}{'MAE':>7}{'RMSE':>7}{'bias':>7}  supported")
        for _, r in t.iterrows():
            print(f"  {str(r[col]):<26}{r['n_rows']:>6}{r['n_unique']:>15}{r['mae']:>7.2f}"
                  f"{r['rmse']:>7.2f}{r['bias']:>+7.2f}  {'yes' if r['supported'] else 'NO (< ' + str(config.AD_MIN_UNIQUE) + ' molecules)'}")
    for split in ["random", "scaffold", "paper"]:
        t = error_table(p4[p4["split"] == split], ["raw_region"])
        t.insert(0, "measure", f"raw similarity region, {split} split")
        outputs.append(t.rename(columns={"raw_region": "region"}))

    # ---- does the size percentile add anything beyond raw similarity? -----------
    # Inside one raw-similarity region, does error still depend on the percentile?
    print("\nInside each raw region: MAE by size-percentile region (distinct molecules in brackets).")
    print("If the percentile carried extra information, MAE would change across a row.")
    both = error_table(p4, ["raw_region", "pct_region"])
    both.insert(0, "measure", "raw region x percentile region")
    outputs.append(both.rename(columns={"raw_region": "region", "pct_region": "region2"}))
    print(f"  {'raw region':<12}" + "".join(f"{n[:22]:>26}" for n in config.AD_PCT_NAMES))
    for rn in config.AD_RAW_NAMES:
        cells = []
        for pn in config.AD_PCT_NAMES:
            r = both[(both["raw_region"] == rn) & (both["pct_region"] == pn)]
            if r.empty:
                cells.append("-")
            elif not r.iloc[0]["supported"]:
                cells.append(f"n/a ({int(r.iloc[0]['n_unique'])})")
            else:
                cells.append(f"{r.iloc[0]['mae']:.2f} ({int(r.iloc[0]['n_unique'])})")
        print(f"  {rn:<12}" + "".join(f"{c:>26}" for c in cells))

    pd.concat(outputs, ignore_index=True).to_csv(config.RESULTS_DIR / "07_reliability_bins.csv",
                                                  index=False, float_format="%.4f")
    plot_reliability(p4, config.FIGURES_DIR / "07_reliability_curves.png",
                     "Surrogate error against similarity to its training set (Stage 4, out-of-sample; 8 equal-count bins per curve)")
    print("\nWrote results/07_stage4_ad.csv, results/07_reliability_bins.csv, figures/07_reliability_curves.png")


if __name__ == "__main__":
    main()
