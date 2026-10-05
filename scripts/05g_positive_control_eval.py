#!/usr/bin/env python3
"""
Stage 5g: how close did the GA get to the held-out actives? (corrected evaluation, v2 outputs)

Reads : results/05f_<control>_populations.csv, results/05e_holdout_sets.csv,
        data/processed/cdk2_ic50_curated.csv
Writes: results/05g_positive_control_summary_v2.csv, figures/05g_positive_control_v2.png  (suffix: --out-suffix)

The historical outputs (results/05g_positive_control_summary.csv, figures/05g_positive_control.png) are kept
unchanged. They used labels that overstated what was measured; this version defines every quantity separately
(cdk2moo/recovery.py has the full definitions):
  proximity            ECFP4 Tanimoto to the nearest held-out active, shares at 0.5/0.6/0.7/0.8. Closeness in
                       fingerprint space only.
  fingerprint identity ECFP4 bit vector equal to a held-out active's. ECFP4 ignores chirality: NOT molecular identity.
  exact identity       standardized structure equal to a held-out active's: std_isomeric (the project's standardization,
                       which strips stereo at tautomerizable centres), std_flat (no stereo), strict_stereo (stereo kept,
                       tautomers not merged).
  nearest is held-out  the nearest reference molecule among all 2,016 is a held-out active. Ties: lowest index (tie share and
                       "any tied neighbour" reported). The "reference prevalence" column is the share of the 2,016 that
                       are held-out actives; it is NOT a chance rate or a null model.
  neighbour annotation among queries with a neighbour at >= 0.6, the share whose nearest reference molecule has MEASURED
                       pActivity >= the top-decile cutoff: a label copied from the neighbour, nothing measured for the query.
Self-matches: the starting molecules and the "remaining training" baseline are reference molecules, so their own identity is
excluded from their nearest-neighbour search; generated molecules are not excluded (reproducing a reference is a rediscovery).
Counts are population-weighted (every occurrence, seeds pooled) and distinct (each standardized isomeric structure once).

Groups: each arm's final generation (generated molecules only), the starting molecules, and the remaining training molecules in
the size window (what siblings alone achieve). Recall of held-out actives pools every distinct generated molecule of all
generations and seeds.

Run:
    python scripts/05g_positive_control_eval.py
    python scripts/05g_positive_control_eval.py --out-suffix v2b      # another output name
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.features import murcko_scaffold
from cdk2moo.ga_plots import ARM_STYLE, INK, MUTED, GRID
from cdk2moo.objectives import ARMS
from cdk2moo.recovery import build_reference, evaluate_group, proximity_to_held, recall_of_held_actives

BASELINE_STYLE = {"start molecules": ("#9a9a94", ":"), "remaining training": ("#52514e", "--")}


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-suffix", default="v2")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    csv_path = config.RESULTS_DIR / f"05g_positive_control_summary_{args.out_suffix}.csv"
    fig_path = config.FIGURES_DIR / f"05g_positive_control_{args.out_suffix}.png"
    for path in (csv_path, fig_path):
        if path.exists() and not args.overwrite:
            raise SystemExit(f"{path} exists; choose another --out-suffix or pass --overwrite.")

    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    sizes = df["n_heavy_atoms"].to_numpy()
    roles = pd.read_csv(config.RESULTS_DIR / "05e_holdout_sets.csv")
    cutoff = roles.loc[roles["is_top_decile"], "pactivity"].min()
    cache = {}
    reference = build_reference(df["std_smiles"], df["pactivity"], cache)
    thresholds = config.REDISCOVERY_SIMS
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), facecolor="#fcfcfb")

    for ax, control in zip(axes, ["molecule", "scaffold"]):
        held_role = roles[f"{control}_control"] == "held"
        held_actives = np.where(held_role & roles["is_top_decile"])[0]
        train = np.where(~held_role)[0]
        pop = pd.read_csv(config.RESULTS_DIR / f"05f_{control}_populations.csv")
        print(f"\n===== {control.upper()}-LEVEL CONTROL: {len(held_actives)} held-out actives (measured pActivity >= {cutoff:.2f}) =====")

        final = pop[(pop["generation"] == pop["generation"].max()) & (pop["birth_generation"] > 0)]
        groups = {arm: (list(final[final["arm"] == arm]["smiles"]), False) for arm in ARMS}
        groups["start molecules"] = (list(pop[(pop["generation"] == 0) & (pop["arm"] == "multi_real")]["smiles"]), True)
        in_window = train[(sizes[train] >= config.GA_MIN_HEAVY_ATOMS) & (sizes[train] <= config.GA_MAX_HEAVY_ATOMS)]
        groups["remaining training"] = (list(df["std_smiles"].iloc[in_window]), True)

        print(f"  {'group':<20}{'n (distinct)':>14}{'median':>8}" + "".join(f"{'>=' + str(t):>7}" for t in thresholds)
              + f"{'fp-ident':>10}{'exact iso/strict':>18}{'NN held':>9}{'(prevalence)':>14}{'annot.':>8}")
        for name, (smiles, exclude_self) in groups.items():
            m = evaluate_group(smiles, reference, held_actives, thresholds, cutoff, config.REDISCOVERY_SIM,
                               exclude_self=exclude_self, cache=cache)
            print(f"  {name:<20}{m['n_weighted']:>7} ({m['n_distinct']:>3}){m['median_proximity_to_held_actives_weighted']:>8.2f}"
                  + "".join(f"{100 * m[f'proximity_ge_{t}_weighted']:>6.0f}%" for t in thresholds)
                  + f"{100 * m['fingerprint_identical_weighted']:>9.1f}%"
                  + f"{100 * m['exact_std_isomeric_weighted']:>10.1f}%/{100 * m['exact_strict_stereo_weighted']:.1f}%"
                  + f"{100 * m['nearest_is_heldout_weighted']:>8.0f}%{100 * m['reference_prevalence_of_heldout']:>13.0f}%"
                  + f"{100 * m['neighbour_annotation_potent_weighted']:>7.0f}%")
            rows.append({"control": control, "group": name, "self_excluded": exclude_self, **m})
            xs = np.sort(proximity_to_held(smiles, reference, held_actives))
            style = ARM_STYLE.get(name)
            color, ls, label = (style["color"], style["ls"], style["label"]) if style else (*BASELINE_STYLE[name], name)
            ax.plot(xs, np.arange(1, len(xs) + 1) / len(xs), color=color, ls=ls, linewidth=2, label=label)

        print(f"\n  Recall: share of held-out actives with >= 1 distinct generated molecule at proximity >= {config.REDISCOVERY_SIM} "
              "(all generations, seeds pooled; proximity, not exact recovery)")
        held_scaf = {murcko_scaffold(s) for s in df["std_smiles"].iloc[held_actives]}
        for arm in ARMS:
            generated = list(pop[(pop["arm"] == arm) & (pop["birth_generation"] > 0)]["smiles"])
            share, hits, n_distinct = recall_of_held_actives(generated, reference, held_actives, config.REDISCOVERY_SIM)
            scaffolds = {murcko_scaffold(s) for s in dict.fromkeys(generated)}
            exact_scaf = held_scaf & scaffolds
            print(f"  {arm:<17} {hits:>3} of {len(held_actives)} actives ({100 * share:.0f}%); distinct generated molecules {n_distinct} (distinct standardized isomeric identity); "
                  f"scaffold equality (not molecule identity; small scaffolds match by chance) for {len(exact_scaf)} of {len(held_scaf)} held-out scaffolds")
            rows.append({"control": control, "group": arm + " (recall)", "recall_share": share, "actives_hit": hits,
                         "n_actives": len(held_actives), "n_distinct_generated": n_distinct, "scaffold_equality_hits": len(exact_scaf), "n_held_scaffolds": len(held_scaf)})
        ceiling = proximity_to_held(df["std_smiles"].iloc[train], reference, held_actives)
        sibling = reference["fps"][held_actives] @ reference["fps"][train].T
        sibling = sibling / (reference["fps"][held_actives].sum(1)[:, None] + reference["fps"][train].sum(1)[None, :] - sibling)
        print(f"  (siblings alone: {100 * (sibling.max(axis=1) >= config.REDISCOVERY_SIM).mean():.0f}% of held-out actives already have a "
              f"remaining-training molecule within {config.REDISCOVERY_SIM})")

        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color=GRID, linewidth=0.8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.axvline(config.REDISCOVERY_SIM, color=MUTED, linewidth=0.8, linestyle=":")
        ax.set_xlabel("ECFP4 Tanimoto proximity to the nearest held-out active", color=INK)
        ax.set_ylabel("cumulative share of molecules", color=INK)
        ax.set_title(f"{control}-level holdout ({len(held_actives)} held-out actives)", color=INK, fontsize=11)
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle(f"Proximity of generated molecules to held-out actives (fingerprint closeness; not identity, not measured activity). "
                 f"Campaign {config.CAMPAIGN}, {config.SCOPE} policy; curves weight every occurrence of a molecule, 5 GA seeds pooled.",
                 color=INK, fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(fig_path, dpi=150)
    pd.DataFrame(rows).to_csv(csv_path, index=False, float_format="%.4f")
    print(f"\nWrote {csv_path.name} and {fig_path.name}")


if __name__ == "__main__":
    main()
