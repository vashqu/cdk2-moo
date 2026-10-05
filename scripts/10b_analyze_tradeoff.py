#!/usr/bin/env python3
"""
Stage 10b: H5, what does a similarity floor cost, and what does it buy?

Reads : results/10_constrained_<experiment>_populations.csv  (floors 0.4-0.7)
        results/05_ga_populations.csv, 05f_molecule_populations.csv, 05f_scaffold_populations.csv
        (the unconstrained twins, same seeds and start molecules), results/05e_holdout_sets.csv, ecfp4.npy
Writes: results/10_tradeoff_v2.csv, figures/10_tradeoff_v2.png   (suffix: --out-suffix; the historical
        results/10_tradeoff.csv and figures/10_tradeoff.png are kept and use the older labels "recovery" and "NN share")

Cost (H1 side): the real surrogate's predicted pActivity of generated molecules, final generation.
Benefit side, measured on the held-out-actives controls (D-26) by FINGERPRINT PROXIMITY, not by the surrogate:
  proximity share  share of generated final-generation molecules with ECFP4 (Morgan r=2, 2048-bit) Tanimoto >= 0.6 to
                   some held-out active. Closeness in fingerprint space; it is not recovery of the molecule and not
                   evidence of activity (the earlier outputs called this "recovery").
  nearest is held  share whose nearest neighbour among all 2,016 molecules is a held-out active (ties: lowest index).
  proximity recall share of held-out actives with at least one distinct generated molecule at proximity >= 0.6, any generation.
Exchange rate: percentage points of proximity share gained per 1.0 pActivity of predicted activity given up, versus the
unconstrained twin of the same seed. Undefined when the predicted activity does not fall.

Only generated molecules (birth_generation > 0) of the final generation are used, except recall.

Run:
    python scripts/10b_analyze_tradeoff.py
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.features import ecfp4
from cdk2moo.ga_plots import ARM_STYLE, INK, MUTED, GRID

TAUS = [None, 0.4, 0.5, 0.6, 0.7]
ARMS = ["multi_real", "activity_only"]
BASELINE_FILE = {"primary": "05_ga_populations.csv", "molecule": "05f_molecule_populations.csv",
                 "scaffold": "05f_scaffold_populations.csv"}


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-suffix", default="v2")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    csv_path = config.RESULTS_DIR / f"10_tradeoff_{args.out_suffix}.csv"
    fig_path = config.FIGURES_DIR / f"10_tradeoff_{args.out_suffix}.png"
    for path in (csv_path, fig_path):
        if path.exists() and not args.overwrite:
            raise SystemExit(f"{path} exists; choose another --out-suffix or pass --overwrite.")
    fps = np.load(config.PROCESSED_DIR / "ecfp4.npy").astype(np.float32)
    roles = pd.read_csv(config.RESULTS_DIR / "05e_holdout_sets.csv")
    rows = []

    for experiment in ["primary", "molecule", "scaffold"]:
        base = pd.read_csv(config.RESULTS_DIR / BASELINE_FILE[experiment])
        base = base[base["arm"].isin(ARMS)].assign(min_similarity=np.nan)
        cons = pd.read_csv(config.RESULTS_DIR / f"10_constrained_{experiment}_populations.csv")
        allpop = pd.concat([base[cons.columns.intersection(base.columns)], cons], ignore_index=True)
        if experiment != "primary":
            held = np.where((roles[f"{experiment}_control"] == "held") & roles["is_top_decile"])[0]
        for arm in ARMS:
            for tau in TAUS:
                sub = allpop[(allpop["arm"] == arm) & (allpop["min_similarity"].isna() if tau is None
                                                       else np.isclose(allpop["min_similarity"].fillna(-1), tau))]
                gen = sub[sub["birth_generation"] > 0]
                final = gen[gen["generation"] == gen["generation"].max()]
                for seed, g in final.groupby("seed"):
                    row = {"experiment": experiment, "arm": arm, "tau": np.nan if tau is None else tau, "seed": seed,
                           "n_generated": len(g), "pred_real": g["pred_real"].median(), "sim": g["max_tanimoto"].median(),
                           "qed": g["qed"].median(), "sa": g["sa"].median(), "atoms": g["n_heavy_atoms"].median()}
                    if experiment != "primary" and len(g):
                        f = ecfp4(list(g["smiles"]), config.FP_RADIUS, config.FP_BITS).astype(np.float32)
                        sh = f @ fps.T
                        sim = sh / (f.sum(1)[:, None] + fps.sum(1)[None, :] - sh)
                        row["proximity_share"] = float((sim[:, held].max(axis=1) >= config.REDISCOVERY_SIM).mean())
                        row["nearest_is_heldout"] = float(np.isin(sim.argmax(axis=1), held).mean())
                    rows.append(row)
                if experiment != "primary":                     # recall: pooled seeds, any generation
                    smiles = list(gen["smiles"].unique())
                    f = ecfp4(smiles, config.FP_RADIUS, config.FP_BITS).astype(np.float32)
                    sh = f @ fps[held].T
                    sim = sh / (f.sum(1)[:, None] + fps[held].sum(1)[None, :] - sh)
                    rows.append({"experiment": experiment, "arm": arm, "tau": np.nan if tau is None else tau,
                                 "seed": -1, "recall": float((sim >= config.REDISCOVERY_SIM).any(axis=0).mean())})
    res = pd.DataFrame(rows)
    res.to_csv(csv_path, index=False)
    per_seed = res[res["seed"] >= 0]
    recall = {(r["experiment"], r["arm"], None if np.isnan(r["tau"]) else round(float(r["tau"]), 1)): r["recall"]
              for _, r in res[res["seed"] == -1].iterrows()}

    # ---- A. what the floor costs (primary surrogate) ---------------------------------
    print("A. COST: primary surrogate (all 2,016 molecules), generation 50, generated molecules")
    print("   median over 5 seeds [min, max]; change is paired with the same seed's unconstrained run")
    for arm in ARMS:
        print(f"\n   {arm}")
        print(f"   {'floor':<7}{'predicted pActivity':>24}{'change vs none':>24}{'raw similarity':>18}{'QED':>7}{'SA':>6}{'atoms':>7}{'n/seed':>8}")
        d = per_seed[(per_seed["experiment"] == "primary") & (per_seed["arm"] == arm)]
        none = d[d["tau"].isna()].set_index("seed")["pred_real"]
        for tau in sorted(d["tau"].unique(), key=lambda t: -1 if np.isnan(t) else t):
            g = d[(d["tau"].isna() if np.isnan(tau) else d["tau"] == tau)].set_index("seed")
            delta = (g["pred_real"] - none).dropna()
            tau_text = "none" if np.isnan(tau) else f"{tau:.1f}"
            print(f"   {tau_text:<7}{g['pred_real'].median():>9.2f} [{g['pred_real'].min():.2f}, {g['pred_real'].max():.2f}]"
                  f"{delta.median():>12.2f} [{delta.min():+.2f}, {delta.max():+.2f}]{g['sim'].median():>18.2f}"
                  f"{g['qed'].median():>7.2f}{g['sa'].median():>6.2f}{g['atoms'].median():>7.0f}{g['n_generated'].median():>8.0f}")

    # ---- B-C. what it buys on the held-out controls, and the exchange rate -----------------
    exchange = []
    for experiment in ["molecule", "scaffold"]:
        print(f"\nB. TRANSFER: {experiment}-level held-out-actives control (surrogate never saw the actives)")
        for arm in ARMS:
            print(f"\n   {arm}")
            print(f"   {'floor':<7}{'predicted':>11}{'proximity share >= 0.6':>28}{'nearest is held-out':>22}{'proximity recall':>18}{'sim':>7}")
            d = per_seed[(per_seed["experiment"] == experiment) & (per_seed["arm"] == arm)]
            none = d[d["tau"].isna()].set_index("seed")
            for tau in sorted(d["tau"].unique(), key=lambda t: -1 if np.isnan(t) else t):
                g = d[(d["tau"].isna() if np.isnan(tau) else d["tau"] == tau)].set_index("seed")
                tau_text = "none" if np.isnan(tau) else f"{tau:.1f}"
                print(f"   {tau_text:<7}{g['pred_real'].median():>11.2f}"
                      f"{100 * g['proximity_share'].median():>14.0f}% [{100 * g['proximity_share'].min():.0f}, {100 * g['proximity_share'].max():.0f}]"
                      f"{100 * g['nearest_is_heldout'].median():>20.0f}%{100 * recall[(experiment, arm, None if np.isnan(tau) else round(float(tau), 1))]:>17.0f}%{g['sim'].median():>7.2f}")
                if not np.isnan(tau):
                    d_prox = 100 * (g["proximity_share"] - none["proximity_share"])
                    d_pred = g["pred_real"] - none["pred_real"]
                    exchange.append({"experiment": experiment, "arm": arm, "tau": tau,
                                     "d_proximity_pp": d_prox.median(), "d_pred": d_pred.median()})
    print("\nC. EXCHANGE RATE: percentage points of proximity share gained per 1.0 pActivity of predicted activity given up")
    print("   (median paired change over seeds; 'n/a' where predicted activity does not fall by more than 0.05)")
    print(f"   {'experiment':<10}{'arm':<15}{'floor':>6}{'proximity change':>17}{'predicted change':>18}{'pp per pActivity lost':>24}")
    for e in exchange:
        rate = f"{e['d_proximity_pp'] / -e['d_pred']:+.1f}" if e["d_pred"] < -0.05 else "n/a"
        print(f"   {e['experiment']:<10}{e['arm']:<15}{e['tau']:>6.1f}{e['d_proximity_pp']:>+14.1f} pp{e['d_pred']:>+18.2f}{rate:>24}")

    # ---- figure ------------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.4), facecolor="#fcfcfb")
    for ax in axes:
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    for ax, experiment in zip(axes[:2], ["molecule", "scaffold"]):
        for arm in ARMS:
            d = per_seed[(per_seed["experiment"] == experiment) & (per_seed["arm"] == arm)]
            xs, ys, lo, hi, taus = [], [], [], [], []
            for tau in sorted(d["tau"].unique(), key=lambda t: -1 if np.isnan(t) else t):
                g = d[d["tau"].isna() if np.isnan(tau) else d["tau"] == tau]
                xs.append(g["pred_real"].median()); ys.append(100 * g["proximity_share"].median())
                lo.append(100 * g["proximity_share"].min()); hi.append(100 * g["proximity_share"].max()); taus.append(tau)
            st = ARM_STYLE[arm]
            ax.errorbar(xs, ys, yerr=[np.array(ys) - lo, np.array(hi) - ys], color=st["color"], ls=st["ls"],
                        marker="o", linewidth=2, capsize=3, label=arm)
            for x, y, t in zip(xs, ys, taus):
                t_text = "none" if np.isnan(t) else f"{t:.1f}"
                ax.annotate(t_text, (x, y), textcoords="offset points", xytext=(5, 5), fontsize=8, color=MUTED)
        ax.set_xlabel("predicted pActivity of generated molecules (control surrogate)", color=INK)
        ax.set_ylabel("generated molecules with ECFP4 proximity >= 0.6 to a held-out active (%)", color=INK)
        ax.set_title(f"{experiment}-level control; labels = similarity floor", color=INK, fontsize=10)
    for arm in ARMS:
        d = per_seed[(per_seed["experiment"] == "primary") & (per_seed["arm"] == arm)]
        order = sorted(d["tau"].unique(), key=lambda t: -1 if np.isnan(t) else t)
        axes[2].plot(range(len(order)), [d[d["tau"].isna() if np.isnan(t) else d["tau"] == t]["pred_real"].median() for t in order],
                     color=ARM_STYLE[arm]["color"], ls=ARM_STYLE[arm]["ls"], marker="o", linewidth=2, label=arm)
    axes[2].set_xticks(range(len(TAUS)))
    axes[2].set_xticklabels(["none"] + [f"{t:.1f}" for t in TAUS[1:]])
    axes[2].set_xlabel("similarity floor", color=INK)
    axes[2].set_ylabel("predicted pActivity (primary surrogate)", color=INK)
    axes[2].set_title("what the floor costs in predicted activity", color=INK, fontsize=10)
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("H5: similarity floor against predicted activity and proximity to held-out actives", color=INK, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(fig_path, dpi=150)
    print(f"\nWrote {csv_path.name} and {fig_path.name}")


if __name__ == "__main__":
    main()
