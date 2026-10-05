#!/usr/bin/env python3
"""
Stage 7b: where do the generated molecules sit relative to what Stage 4 can vouch for?

Reads : results/05_ga_populations_scored.csv, results/06_pareto_front.csv,
        results/07_stage4_ad.csv
Writes: results/07_generated_regions.csv  <- every distinct generated molecule with
                                             its similarity region, prediction band,
                                             and whether Stage 4 can estimate error there

Regions are defined on the raw ECFP4 max-Tanimoto (config.AD_RAW_*), with the size
percentile reported alongside. Prediction bands: <6, 6-7, 7-8, >=8 pActivity.

What this stage does NOT claim: that the surrogate is reliable anywhere Stage 4
has fewer than config.AD_MIN_UNIQUE distinct validation molecules, nor that low
similarity means the surrogate was exploited. Similarity falls in druglike_only as
much as anywhere (its GA never sees the surrogate), so a drift away from the
training set is not, by itself, evidence of exploitation.

Run:
    python scripts/07b_ad_generated_overlay.py
"""

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.applicability import add_region, error_table

ARMS = ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]
PRED_EDGES = [0.0, 6.0, 7.0, 8.0, 10.0]
PRED_NAMES = ["<6", "6-7", "7-8", ">=8"]


def main():
    config.require_campaign()
    pop = pd.read_csv(config.RESULTS_DIR / "05_ga_populations_scored.csv")
    gen = pop[pop["birth_generation"] > 0]          # generated molecules only (D-21)
    front = pd.read_csv(config.RESULTS_DIR / "06_pareto_front.csv")
    p4 = pd.read_csv(config.RESULTS_DIR / "07_stage4_ad.csv")

    # ---- what Stage 4 can vouch for, cell by cell -----------------------------
    p4 = add_region(p4, "max_tanimoto_to_train", config.AD_RAW_EDGES, config.AD_RAW_NAMES, "raw_region")
    p4 = add_region(p4, "pred_rf", PRED_EDGES, PRED_NAMES, "pred_band")
    cells = error_table(p4, ["raw_region", "pred_band"])
    print("Stage 4 validation coverage: distinct test molecules (MAE, bias) by similarity region x prediction band")
    print(f"'n/a' = fewer than {config.AD_MIN_UNIQUE} distinct molecules: error cannot be estimated there")
    print(f"  {'similarity':<11}" + "".join(f"{b:>24}" for b in PRED_NAMES))
    for rn in config.AD_RAW_NAMES:
        row = []
        for pn in PRED_NAMES:
            c = cells[(cells["raw_region"] == rn) & (cells["pred_band"] == pn)]
            if c.empty:
                row.append("none")
            elif not c.iloc[0]["supported"]:
                row.append(f"n/a  (n={int(c.iloc[0]['n_unique'])})")
            else:
                c = c.iloc[0]
                row.append(f"{c['mae']:.2f} {c['bias']:+.2f} (n={int(c['n_unique'])})")
        print(f"  {rn:<11}" + "".join(f"{c:>24}" for c in row))
    print("  (cell entries: MAE, bias = prediction minus measured)")

    # ---- every distinct generated molecule, labelled --------------------------
    pooled = (gen.groupby("smiles")
              .agg(pred=("pred_real", "first"), qed=("qed", "first"), sa=("sa", "first"),
                   max_tanimoto=("max_tanimoto", "first"), size_pctile=("size_pctile", "first"),
                   n_heavy_atoms=("n_heavy_atoms", "first"),
                   arms=("arm", lambda a: "|".join(sorted(set(a)))))
              .reset_index())
    pooled = add_region(pooled, "max_tanimoto", config.AD_RAW_EDGES, config.AD_RAW_NAMES, "raw_region")
    pooled = add_region(pooled, "size_pctile", config.AD_PCT_EDGES, config.AD_PCT_NAMES, "pct_region")
    pooled = add_region(pooled, "pred", PRED_EDGES, PRED_NAMES, "pred_band")
    pooled["on_front"] = pooled["smiles"].isin(set(front["smiles"]))
    key = cells.set_index(["raw_region", "pred_band"])["supported"]
    pooled["stage4_can_estimate_error"] = [
        bool(key.get((r, b), False)) for r, b in zip(pooled["raw_region"], pooled["pred_band"])]
    # Three clearly defined classes of candidate, on the raw measure:
    pooled["candidate_class"] = np.select(
        [(pooled["pred"] >= 8) & (pooled["max_tanimoto"] >= 0.6),
         (pooled["pred"] >= 8) & (pooled["max_tanimoto"] < 0.6)],
        ["A: high prediction, supported (raw >= 0.6)",
         "B: high prediction, weak support (raw < 0.6)"],
        default="C: prediction < 8")
    pooled.to_csv(config.RESULTS_DIR / "07_generated_regions.csv", index=False, float_format="%.4f")

    # ---- the final generation of each arm, by region ---------------------------
    final = gen[gen["generation"] == gen["generation"].max()].copy()
    final = add_region(final, "max_tanimoto", config.AD_RAW_EDGES, config.AD_RAW_NAMES, "raw_region")
    final = add_region(final, "size_pctile", config.AD_PCT_EDGES, config.AD_PCT_NAMES, "pct_region")
    print("\nGeneration 50, generated molecules, 5 seeds pooled: share in each raw-similarity region")
    print(f"  {'arm':<16}" + "".join(f"{n:>10}" for n in config.AD_RAW_NAMES) + f"{'n':>6}")
    for arm in ARMS:
        a = final[final["arm"] == arm]
        shares = a["raw_region"].value_counts(normalize=True)
        print(f"  {arm:<16}" + "".join(f"{100 * shares.get(n, 0):>9.0f}%" for n in config.AD_RAW_NAMES) + f"{len(a):>6}")
    print("  ...and in each size-percentile region")
    print(f"  {'arm':<16}" + "".join(f"{n[:16]:>18}" for n in config.AD_PCT_NAMES))
    for arm in ARMS:
        a = final[final["arm"] == arm]
        shares = a["pct_region"].value_counts(normalize=True)
        print(f"  {arm:<16}" + "".join(f"{100 * shares.get(n, 0):>17.0f}%" for n in config.AD_PCT_NAMES))
    weak = {arm: (final[final["arm"] == arm]["max_tanimoto"] < 0.4).mean() for arm in ARMS}
    print(f"\n  DRIFT IS NOT EXPLOITATION: share with raw similarity < 0.4: "
          + ", ".join(f"{a} {100 * v:.0f}%" for a, v in weak.items()))
    print("  druglike_only never sees the surrogate yet drifts as far; low similarity alone proves nothing about exploitation.")

    # ---- candidate classes --------------------------------------------------
    print("\nCandidate classes among distinct generated molecules (all generations, 5 seeds pooled):")
    print(f"  {'arm':<16}{'distinct':>9}{'A supported':>13}{'B weak support':>16}{'C pred<8':>10}")
    for arm in ARMS:
        a = pooled[pooled["arms"].str.contains(arm)]
        cls = a["candidate_class"].str[0].value_counts()
        print(f"  {arm:<16}{len(a):>9}{cls.get('A', 0):>13}{cls.get('B', 0):>16}{cls.get('C', 0):>10}")
    f = pooled[pooled["on_front"]]
    cls = f["candidate_class"].str[0].value_counts()
    print(f"  {'Pareto front':<16}{len(f):>9}{cls.get('A', 0):>13}{cls.get('B', 0):>16}{cls.get('C', 0):>10}")

    print("\nPareto-front molecules by Stage 4 support:")
    print(f"  front molecules in a (similarity x prediction) cell where Stage 4 can estimate error: "
          f"{int(f['stage4_can_estimate_error'].sum())}/{len(f)}")
    hi = f[f["pred"] >= 8]
    print(f"  of the {len(hi)} front molecules predicted >= 8: "
          f"{int((hi['stage4_can_estimate_error']).sum())} sit in a supported cell, "
          f"{int((~hi['stage4_can_estimate_error']).sum())} do not")
    print("  class B = high prediction and weak support (no error claim possible):")
    b = f[f["candidate_class"].str.startswith("B")].sort_values("pred", ascending=False)
    print(f"    {'pred':>5} {'QED':>5} {'SA':>5} {'rawsim':>7} {'pctile':>7} {'atoms':>6}  arms / SMILES")
    for _, r in b.iterrows():
        print(f"    {r['pred']:>5.2f} {r['qed']:>5.2f} {r['sa']:>5.2f} {r['max_tanimoto']:>7.2f} "
              f"{r['size_pctile']:>7.1f} {r['n_heavy_atoms']:>6}  {r['arms']}  {r['smiles']}")

    print("\nAgreement of the two measures for generated molecules (distinct, pooled): rows = raw region, cols = percentile region")
    ct = pd.crosstab(pooled["raw_region"], pooled["pct_region"])
    print(ct.to_string())
    print("\nWrote results/07_generated_regions.csv")


if __name__ == "__main__":
    main()
