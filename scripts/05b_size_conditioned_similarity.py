#!/usr/bin/env python3
"""
Stage 5b: add a size-conditioned similarity percentile to the GA populations.

Reads : data/processed/cdk2_ic50_curated.csv, ecfp4.npy
        results/05_ga_populations.csv       (never modified)
Writes: results/05_training_loo_similarity.csv   <- reference: each training
                                                    molecule's nearest-neighbour
                                                    Tanimoto to the OTHER 2,015
        results/05_ga_populations_scored.csv     <- populations + size_pctile,
                                                    stratum_halfwidth, stratum_n

Both measures are kept: max_tanimoto (raw ECFP4 / Morgan radius 2, 2048-bit,
binary Tanimoto to all 2,016 training molecules) and size_pctile (that value as
a percentile among size-matched training molecules; see cdk2moo/similarity.py).

Run:
    python scripts/05b_size_conditioned_similarity.py
"""

import argparse

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.similarity import loo_max_tanimoto, size_conditioned_percentile


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", default="05_ga_populations.csv")
    ap.add_argument("--out", default="05_ga_populations_scored.csv")
    args = ap.parse_args()

    train = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")

    # ---- the reference distribution: leave-one-out, per training molecule ----
    loo = loo_max_tanimoto(fps)
    ref = pd.DataFrame({"inchikey": train["inchikey"],
                        "n_heavy_atoms": train["n_heavy_atoms"],
                        "loo_max_tanimoto": loo})
    ref.to_csv(config.RESULTS_DIR / "05_training_loo_similarity.csv", index=False)

    print("Training set, leave-one-out nearest-neighbour ECFP4 Tanimoto (2048-bit, radius 2)")
    print(f"  all molecules: median {np.median(loo):.2f}, "
          f"10th/90th percentile {np.percentile(loo, 10):.2f}/{np.percentile(loo, 90):.2f}")
    print(f"  identical fingerprint to another training molecule: {(loo >= 0.999).sum()}")
    print("  by size (this is the size dependence the percentile removes):")
    for lo, hi in [(5, 19), (20, 24), (25, 29), (30, 34), (35, 39), (40, 80)]:
        sel = ref["n_heavy_atoms"].between(lo, hi)
        print(f"    {lo:>2}-{hi:<2} heavy atoms  n={int(sel.sum()):>4}  median LOO similarity "
              f"{ref.loc[sel, 'loo_max_tanimoto'].median():.2f}")

    # ---- annotate the generated molecules ----------------------------------
    pop = pd.read_csv(config.RESULTS_DIR / args.pop)
    pct, width, n_ref = size_conditioned_percentile(
        pop["max_tanimoto"], pop["n_heavy_atoms"],
        train["n_heavy_atoms"].to_numpy(), loo)
    pop["size_pctile"] = pct
    pop["stratum_halfwidth"] = width
    pop["stratum_n"] = n_ref
    pop.to_csv(config.RESULTS_DIR / args.out, index=False, float_format="%.5f")

    print(f"\nStrata: half-width used {config.SIZE_STRATUM_HALFWIDTH} unless fewer than "
          f"{config.SIZE_STRATUM_MIN_N} reference molecules")
    print("  half-width actually used (rows):", pop["stratum_halfwidth"].value_counts().sort_index().to_dict())
    print(f"  reference molecules per stratum: min {pop['stratum_n'].min()}, median "
          f"{int(pop['stratum_n'].median())}, max {pop['stratum_n'].max()}")
    print(f"\nWrote results/{args.out} ({len(pop)} rows)")


if __name__ == "__main__":
    main()
