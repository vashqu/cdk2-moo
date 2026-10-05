#!/usr/bin/env python3
"""
Stage 7c: descriptive structural-alert audit (PAINS and Brenk) of the generated molecules.

Reads : results/07_generated_regions.csv, results/06_pareto_front.csv,
        results/07_stage4_ad.csv, data/processed/cdk2_ic50_curated.csv
Writes: results/07_generated_alerts.csv   <- alert counts and names, every distinct generated molecule
        results/07_front_alerts.csv       <- same for the Pareto front, with its AD class
        figures/07_generated_on_domain.png

Alerts are RECORDED, never used to remove molecules. Three concepts stay separate:
  graph validity     every generated molecule sanitizes by construction (100% valid)
  structural alerts  PAINS / Brenk substructure matches (cdk2moo/alerts.py)
  domain support     similarity-based class A / B / C from stage 7b

Baselines: the training set and its top decile show what alert rates look like for
real CDK2 inhibitors. Many real kinase inhibitors carry Brenk alerts, so the useful
question is not "zero alerts?" but "more than known chemistry, and more where support is weak?"

Run:
    python scripts/07c_structural_alerts.py
"""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from cdk2moo import config
from cdk2moo.alerts import alert_hits, alert_row
from cdk2moo.ad_plots import plot_generated_on_domain


def main():
    config.require_campaign()
    gen = pd.read_csv(config.RESULTS_DIR / "07_generated_regions.csv")
    front = pd.read_csv(config.RESULTS_DIR / "06_pareto_front.csv")
    train = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    p4 = pd.read_csv(config.RESULTS_DIR / "07_stage4_ad.csv")

    # ---- alerts for every distinct generated molecule and for the baselines ---
    pains_names, brenk_names = [], []
    for smiles in gen["smiles"]:
        p, b = alert_hits(smiles)
        pains_names.append("|".join(p))
        brenk_names.append("|".join(b))
    gen["pains"], gen["brenk"] = pains_names, brenk_names
    gen["n_pains"] = gen["pains"].map(lambda s: len(s.split("|")) if s else 0)
    gen["n_brenk"] = gen["brenk"].map(lambda s: len(s.split("|")) if s else 0)
    gen["any_alert"] = (gen["n_pains"] + gen["n_brenk"]) > 0
    gen.to_csv(config.RESULTS_DIR / "07_generated_alerts.csv", index=False, float_format="%.4f")

    baseline = {}
    cutoff = train["pactivity"].quantile(0.9)
    for name, sub in [("training set (2,016)", train),
                      ("known actives, top decile", train[train["pactivity"] >= cutoff])]:
        hits = [alert_hits(s) for s in sub["std_smiles"]]
        baseline[name] = (np.mean([len(p) > 0 for p, b in hits]), np.mean([len(b) > 0 for p, b in hits]),
                          np.mean([len(p) + len(b) > 0 for p, b in hits]))

    print("Graph validity: 100% of generated molecules sanitize (required to be generated at all).")
    print("That says nothing about whether they are sensible; see the alerts and the examples below.\n")

    print(f"  {'set':<46}{'n':>7}{'PAINS':>10}{'Brenk':>10}{'either':>10}{'Brenk/mol':>10}")
    for name, (pa, br, ei) in baseline.items():
        print(f"  {name:<46}{'':>7}{100 * pa:>9.0f}%{100 * br:>9.0f}%{100 * ei:>9.0f}%")
    for arm in ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]:
        print(alert_row(f"generated, {arm} (distinct)", gen[gen["arms"].str.contains(arm)]))
    print(alert_row("generated, all distinct", gen))

    # ---- enrichment: are high-prediction, weakly supported molecules special? --
    print("\nAlert rate by applicability-domain class (all distinct generated molecules):")
    for cls in sorted(gen["candidate_class"].unique()):
        print(alert_row(cls, gen[gen["candidate_class"] == cls]))
    a, b = gen[gen["candidate_class"].str.startswith("A")], gen[gen["candidate_class"].str.startswith("B")]
    table = [[int(b["any_alert"].sum()), int((~b["any_alert"]).sum())],
             [int(a["any_alert"].sum()), int((~a["any_alert"]).sum())]]
    odds, pval = fisher_exact(table)
    print(f"  B (weak support) vs A (supported), any alert: odds ratio {odds:.2f}, Fisher p = {pval:.3g}"
          f" (the same molecule recurs in related forms, so p is optimistic)")
    print("\nAlert rate by raw-similarity region, generated molecules predicted >= 8 only:")
    hi = gen[gen["pred"] >= 8]
    for rn in config.AD_RAW_NAMES:
        d = hi[hi["raw_region"] == rn]
        if len(d):
            print(alert_row(f"raw similarity {rn}", d))

    # ---- the Pareto front ------------------------------------------------------
    front = front.merge(gen[["smiles", "pains", "brenk", "n_pains", "n_brenk", "any_alert",
                             "candidate_class", "stage4_can_estimate_error"]], on="smiles")
    front.to_csv(config.RESULTS_DIR / "07_front_alerts.csv", index=False, float_format="%.4f")
    print(f"\nPareto front ({len(front)} molecules):")
    print(f"  {'group':<46}{'n':>7}{'PAINS':>10}{'Brenk':>10}{'either':>10}{'Brenk/mol':>10}")
    print(alert_row("front, all", front))
    for cls in sorted(front["candidate_class"].unique()):
        print(alert_row("front, " + cls, front[front["candidate_class"] == cls]))
    print("  most frequent alerts on the front:")
    names = pd.Series([n for s in pd.concat([front["pains"].fillna(""), front["brenk"].fillna("")])
                       for n in s.split("|") if n])
    for n, c in names.value_counts().head(6).items():
        print(f"    {c:>3}  {n}")

    print("\nClass B front molecules (high prediction, weak support): alerts found")
    for _, r in front[front["candidate_class"].str.startswith("B")].sort_values("pred", ascending=False).iterrows():
        found = "|".join(x for x in [r["pains"] if isinstance(r["pains"], str) else "",
                                      r["brenk"] if isinstance(r["brenk"], str) else ""] if x)
        print(f"  pred {r['pred']:.2f} raw {r['max_tanimoto']:.2f}  alerts: {found if found else 'none'}")

    print("\nDo the alerts catch the implausible molecule drawn in Stage 6?")
    top = front.sort_values("pred", ascending=False).iloc[0]
    print(f"  highest-predicted front molecule: {top['smiles']}")
    print(f"    PAINS: {top['pains'] if isinstance(top['pains'], str) and top['pains'] else 'none'}")
    print(f"    Brenk: {top['brenk'] if isinstance(top['brenk'], str) and top['brenk'] else 'none'}")

    plot_generated_on_domain(p4, front, config.FIGURES_DIR / "07_generated_on_domain.png",
                             "Pareto-front molecules on the surrogate's validated domain")
    print("\nWrote results/07_generated_alerts.csv, results/07_front_alerts.csv, figures/07_generated_on_domain.png")


if __name__ == "__main__":
    main()
