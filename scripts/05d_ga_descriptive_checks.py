#!/usr/bin/env python3
"""
Stage 5d: descriptive checks on two things Stage 5 surfaced, with no reruns.

  1. Diversity collapse: multi_real and activity_only end with ~21 distinct
     scaffolds out of 100. When does it happen, is it the same few scaffolds in
     every seed, and are those scaffolds the ones known potent training
     molecules have?
  2. Boundary piling: many generated molecules sit exactly at the lower size
     limit (20 heavy atoms). Is that the window capping an objective's own
     preference for small molecules?

Reads : results/05_ga_populations_scored.csv, data/processed/cdk2_ic50_curated.csv, ecfp4.npy
Only generated molecules (birth_generation > 0) are counted, except generation 0.
Prints everything; run with `> results/05_descriptive_checks.txt` to keep it.

Run:
    python scripts/05d_ga_descriptive_checks.py
"""

import itertools

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED

from cdk2moo import config
from cdk2moo.features import ecfp4, murcko_scaffold
from cdk2moo.objectives import sascorer

ARMS = ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]


def main():
    config.require_campaign()
    pop = pd.read_csv(config.RESULTS_DIR / "05_ga_populations_scored.csv")
    pop = pop[(pop["generation"] == 0) | (pop["birth_generation"] > 0)].copy()
    train = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    train_fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    last = int(pop["generation"].max())

    # Scaffolds are slow-ish to compute, so do each unique SMILES once.
    unique_smiles = pop["smiles"].unique()
    scaffold_of = {s: murcko_scaffold(s) for s in unique_smiles}
    pop["scaffold"] = pop["smiles"].map(scaffold_of)
    train["scaffold"] = train["std_smiles"].map(murcko_scaffold)

    # ================= 1. diversity collapse =================
    print("1. DIVERSITY: distinct scaffolds per population of 100 (median over seeds)")
    marks = [0, 1, 3, 5, 10, 20, 30, 50]
    print(f"   {'gen':>4}" + "".join(f"{a:>18}" for a in ARMS))
    counts = pop.groupby(["arm", "seed", "generation"])["scaffold"].nunique().reset_index()
    for g in marks:
        row = [counts[(counts["arm"] == a) & (counts["generation"] == g)]["scaffold"].median()
               for a in ARMS]
        print(f"   {g:>4}" + "".join(f"{v:>18.0f}" for v in row))
    print("   (generation 0 = 100 random training molecules; later rows count only molecules born after gen 0)")

    print("\n   Do different seeds collapse onto the SAME scaffolds? Mean Jaccard overlap of the")
    print("   final-generation scaffold sets across the 10 pairs of seeds (0 = disjoint, 1 = identical):")
    final = pop[pop["generation"] == last]
    for arm in ARMS:
        sets = {s: set(g["scaffold"]) for s, g in final[final["arm"] == arm].groupby("seed")}
        jac = [len(sets[a] & sets[b]) / len(sets[a] | sets[b])
               for a, b in itertools.combinations(sorted(sets), 2)]
        print(f"     {arm:<16} {np.mean(jac):.2f}")

    print("\n   Top scaffolds in the final populations (5 seeds pooled, 500 molecules per arm),")
    print("   with how many training molecules share that scaffold and their mean measured pActivity:")
    for arm in ["multi_real", "activity_only"]:
        top = final[final["arm"] == arm]["scaffold"].value_counts().head(4)
        print(f"     {arm}")
        for scaf, n in top.items():
            t = train[train["scaffold"] == scaf]
            tinfo = (f"{len(t)} training molecules, mean pActivity {t['pactivity'].mean():.2f}"
                     if len(t) else "no training molecule has this scaffold")
            print(f"       {n:>3}/{len(final[final['arm'] == arm])}  {scaf[:58]:<58}  {tinfo}")

    print("\n   Nearest training molecule of each final generated molecule (by ECFP4 Tanimoto):")
    top_decile = train["pactivity"].quantile(0.9)
    print(f"   (training top decile = measured pActivity >= {top_decile:.2f}; training median {train['pactivity'].median():.2f})")
    print(f"   {'arm':<16}{'median pAct of NN':>20}{'NN is top-decile':>20}{'distinct NN':>14}{'of':>5}")
    for arm in ARMS:
        a = final[final["arm"] == arm]
        fps = ecfp4(list(a["smiles"]), config.FP_RADIUS, config.FP_BITS).astype(np.float32)
        tf = train_fps.astype(np.float32)
        shared = fps @ tf.T
        sim = shared / (fps.sum(1)[:, None] + tf.sum(1)[None, :] - shared)
        nn = sim.argmax(axis=1)
        nn_act = train["pactivity"].to_numpy()[nn]
        print(f"   {arm:<16}{np.median(nn_act):>20.2f}{100 * (nn_act >= top_decile).mean():>19.0f}%"
              f"{len(set(nn)):>14}{len(nn):>5}")

    # ================= 2. boundary piling =================
    print("\n2. BOUNDARY PILING: share of generated molecules at the window edges (20 or 39 heavy atoms)")
    print(f"   {'gen':>4}" + "".join(f"{a:>18}" for a in ARMS))
    for g in [1, 5, 10, 30, 50]:
        row = []
        for a in ARMS:
            s = pop[(pop["arm"] == a) & (pop["generation"] == g)]
            at20 = 100 * (s["n_heavy_atoms"] == 20).mean()
            at39 = 100 * (s["n_heavy_atoms"] == 39).mean()
            row.append(f"{at20:.0f}% @20, {at39:.0f}% @39")
        print(f"   {g:>4}" + "".join(f"{c:>18}" for c in row))

    print("\n   Does the drug-likeness objective itself prefer small molecules? Training molecules,")
    print("   median QED and SA score by heavy-atom count (no GA involved):")
    mols = [Chem.MolFromSmiles(s) for s in train["std_smiles"]]
    train["qed"] = [QED.qed(m) for m in mols]
    train["sa"] = [sascorer.calculateScore(m) for m in mols]
    for lo, hi in [(20, 20), (21, 22), (23, 25), (26, 30), (31, 35), (36, 39)]:
        t = train[train["n_heavy_atoms"].between(lo, hi)]
        print(f"     {lo}-{hi} atoms  n={len(t):>3}  median QED {t['qed'].median():.2f}  median SA {t['sa'].median():.2f}  "
              f"median measured pActivity {t['pactivity'].median():.2f}")

    print("\n   Within the final multi_real generation, size vs the two things it trades off:")
    m = final[final["arm"] == "multi_real"]
    print(f"     Spearman(heavy atoms, real prediction) = {m['n_heavy_atoms'].corr(m['pred_real'], method='spearman'):+.2f}, "
          f"Spearman(heavy atoms, QED) = {m['n_heavy_atoms'].corr(m['qed'], method='spearman'):+.2f}")
    at_edge = m["n_heavy_atoms"] <= 21
    print(f"     molecules with <=21 atoms: {int(at_edge.sum())}/{len(m)}; their median real prediction "
          f"{m.loc[at_edge, 'pred_real'].median():.2f} vs {m.loc[~at_edge, 'pred_real'].median():.2f} for larger ones")


if __name__ == "__main__":
    main()
