#!/usr/bin/env python3
"""
Stage 8e: what do the Vina scores say, and is Vina usable as an orthogonal signal? (campaign version; per experiment group)

Reads : results/08c_docking_sets.csv, results/08d_export/docking_scores.csv, results/05*_populations.csv (predictions of generated molecules),
        results/posebusters_audit/posebusters_checks.csv (for the sensitivity analysis)
Writes: results/08e_docking_by_molecule.csv   every (molecule, set) row; jobs that failed or are missing keep their status
        results/08e_denominators.csv          per set: rows, completed, failed by status, with no job
        results/08e_set_summary.csv           per set: heavy atoms, Vina median and IQR (completed jobs)
        results/08e_validity_gate.csv         per group: can Vina discriminate the comparator actives from matched decoys and random molecules?
        results/08e_h3_endpoints.csv          per group and arm: P_sup against the comparator actives, raw and size-stratified, cluster-bootstrap intervals
        results/08e_pb_sensitivity.csv        the same endpoints excluding poses that did not pass all PoseBusters checks
        results/08e_replicate_noise.csv       engine noise over the three Vina seeds
        figures/08_docking_<group>.png

Definitions are in the campaign plan (section 3) and cdk2moo/docking_stats.py. "Better" = more negative Vina score. P_sup(a over b) is the chance that a random
molecule of a scores better than a random molecule of b (0.5 = no difference). Groups are never pooled: P = primary surrogate (comparator actives were IN its training
data), M = molecule-level holdout, S = scaffold-level holdout (comparators held out). H3 is assessed only in M and S, and only if the validity gate passes there.

Run (inside a policy scope, after 08d and 08f):
    CDK2_CAMPAIGN=<id> CDK2_SCOPE=legacy python scripts/08e_docking_analysis.py
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cdk2moo import config, docking_stats as ds
from cdk2moo.figures_validation import fig_docking_group

GROUPS = {"P": ("training_top_decile", "random_training", "in-sample comparator (the surrogate saw these actives)"),
          "M": ("heldout_actives", "random_remaining", "held-out comparator (molecule-level holdout)"),
          "S": ("heldout_actives", "random_remaining", "held-out comparator (scaffold-level holdout)")}
ARMS = ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]
POPULATIONS = {"P": "05_ga_populations.csv", "M": "05f_molecule_populations.csv", "S": "05f_scaffold_populations.csv"}


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=ds.N_BOOT)
    args = ap.parse_args()
    R, seed0 = config.RESULTS_DIR, config.RANDOM_SEED
    outputs = [R / n for n in ("08e_docking_by_molecule.csv", "08e_denominators.csv", "08e_set_summary.csv", "08e_validity_gate.csv",
                               "08e_h3_endpoints.csv", "08e_pb_sensitivity.csv", "08e_replicate_noise.csv")]
    for path in outputs:
        if path.exists():
            raise SystemExit(f"{path} exists; outputs are never replaced in place")

    sets = pd.read_csv(R / "08c_docking_sets.csv")
    scores = pd.read_csv(R / "08d_export" / "docking_scores.csv")
    first = scores[scores["vina_seed"] == seed0][["smiles", "status", "vina_score", "n_heavy_atoms", "n_rotatable"]]
    df = sets.merge(first, on="smiles", how="left", validate="many_to_one")
    df["status"] = df["status"].fillna("no_job")
    for group, name in POPULATIONS.items():
        pred = pd.read_csv(R / name, usecols=["smiles", "pred_real"]).drop_duplicates("smiles").set_index("smiles")["pred_real"]
        df.loc[df["group"] == group, "pred_real"] = df.loc[df["group"] == group, "smiles"].map(pred)
    pb_path = R / "posebusters_audit" / "posebusters_checks.csv"
    pb = pd.read_csv(pb_path)[["smiles", "status"]].rename(columns={"status": "pb_status"}) if pb_path.exists() else None
    df = df.merge(pb, on="smiles", how="left", validate="many_to_one") if pb is not None else df.assign(pb_status=np.nan)
    df.to_csv(outputs[0], index=False)
    ok = df[df["status"] == "ok"]

    # ---- denominators: failures and missing jobs stay visible ----------------------------------------------
    denominators = df.groupby("set").agg(rows=("smiles", "size"), completed=("status", lambda s: int((s == "ok").sum())),
                                         prep_failed=("status", lambda s: int((s == "prep_failed").sum())),
                                         dock_failed=("status", lambda s: int((s == "dock_failed").sum())),
                                         pose_export_failed=("status", lambda s: int((s == "pose_export_failed").sum())),
                                         no_job=("status", lambda s: int((s == "no_job").sum()))).reset_index()
    denominators.to_csv(outputs[1], index=False)
    print("Denominators (rows in set; completed jobs; failures by kind; jobs missing):")
    print(denominators.to_string(index=False))

    summary = ok.groupby(["group", "set", "role"]).agg(n_completed=("vina_score", "size"), heavy_atoms_median=("n_heavy_atoms", "median"),
                                                       vina_median=("vina_score", "median"), vina_q25=("vina_score", lambda s: s.quantile(0.25)),
                                                       vina_q75=("vina_score", lambda s: s.quantile(0.75))).reset_index()
    summary.to_csv(outputs[2], index=False)

    gate_rows, h3_rows, sens_rows = [], [], []
    for group, (actives_name, random_name, role_text) in GROUPS.items():
        g = ok[ok["group"] == group]
        A, D, Rr = (g[g["set"] == f"{group}:{n}"] for n in (actives_name, "decoys", random_name))
        print(f"\n===== GROUP {group}: {role_text}; {len(A)} comparator actives, {len(D)} decoys, {len(Rr)} random molecules completed =====")
        slope = np.polyfit(g["n_heavy_atoms"], g["vina_score"], 1)[0]
        print(f"  size dependence over all {len(g)} docked molecules of the group: {slope:+.3f} kcal/mol per heavy atom, "
              f"Spearman {spearmanr(g['n_heavy_atoms'], g['vina_score'])[0]:+.2f}")

        # validity gate: can Vina separate the comparator actives from matched decoys / random molecules?
        row = {"group": group, "comparator_role": role_text, "n_actives": len(A), "n_decoys": len(D), "n_random": len(Rr)}
        for label, other in (("decoys", D), ("random", Rr)):
            cmp = ds.compare_samples(A["vina_score"], A["n_heavy_atoms"], A["cluster"], other["vina_score"], other["n_heavy_atoms"], other["cluster"], args.n_boot, seed0)
            row.update({f"auc_vs_{label}": cmp["p_sup"], f"auc_vs_{label}_low": cmp["p_sup_low"], f"auc_vs_{label}_high": cmp["p_sup_high"],
                        f"auc_vs_{label}_stratified": cmp["p_sup_stratified"], f"auc_vs_{label}_strat_low": cmp["p_sup_strat_low"],
                        f"auc_vs_{label}_strat_high": cmp["p_sup_strat_high"]})
        ref = pd.concat([A, Rr])
        ref = ref[ref["pactivity"].notna()]
        rho = lambda ia, ib: spearmanr(-ref["vina_score"].to_numpy()[ia], ref["pactivity"].to_numpy()[ia])[0] if len(ia) > 3 else np.nan   # noqa: E731
        low, high, _ = ds.cluster_bootstrap(rho, ref["vina_score"], ref["cluster"], [0], ["x"], args.n_boot, seed0)
        row.update({"spearman_negvina_vs_measured": float(rho(np.arange(len(ref)), None)), "spearman_low": low, "spearman_high": high, "n_with_measured": len(ref)})
        passed, label = ds.validity_gate(row["auc_vs_decoys_low"], row["auc_vs_decoys_strat_low"], row["auc_vs_decoys"], row["auc_vs_decoys_stratified"])
        row.update({"gate_passed": passed, "gate_label": label})
        gate_rows.append(row)
        print(f"  gate: AUC actives vs decoys {row['auc_vs_decoys']:.2f} [{row['auc_vs_decoys_low']:.2f}, {row['auc_vs_decoys_high']:.2f}], stratified "
              f"{row['auc_vs_decoys_stratified']:.2f} [{row['auc_vs_decoys_strat_low']:.2f}, {row['auc_vs_decoys_strat_high']:.2f}] -> {label.upper()}")
        print(f"        AUC actives vs random {row['auc_vs_random']:.2f} [{row['auc_vs_random_low']:.2f}, {row['auc_vs_random_high']:.2f}]; "
              f"Spearman(-Vina, measured pActivity) {row['spearman_negvina_vs_measured']:+.2f} [{low:+.2f}, {high:+.2f}] (n={len(ref)})")

        # H3 endpoint: does each optimized set score better than the comparator actives?
        for sample_name in ARMS + ([] if group != "P" else ["pareto_front"]) + ["decoys", random_name]:
            S = g[g["set"] == f"{group}:{sample_name}"]
            cmp = ds.compare_samples(S["vina_score"], S["n_heavy_atoms"], S["cluster"], A["vina_score"], A["n_heavy_atoms"], A["cluster"], args.n_boot, seed0)
            if group == "P":
                decision = "in-sample comparison (not H3)"
            elif not passed:
                decision = "not assessable (validity gate failed or undefined)"
            else:
                decision = "exceeded" if cmp["p_sup_low"] > 0.5 else "not exceeded"
            h3_rows.append({"group": group, "sample": sample_name, "is_generated_arm": sample_name in ARMS or sample_name == "pareto_front", "decision": decision, **cmp})
            if pb is not None:
                keep = S[S["pb_status"] == "pass"]
                keepA = A[A["pb_status"] == "pass"]
                sens = ds.compare_samples(keep["vina_score"], keep["n_heavy_atoms"], keep["cluster"], keepA["vina_score"], keepA["n_heavy_atoms"], keepA["cluster"], args.n_boot, seed0)
                sens_rows.append({"group": group, "sample": sample_name, "n_a_all": len(S), "n_a_pass": len(keep), "n_b_all": len(A), "n_b_pass": len(keepA), **sens})
        print("  P_sup of each sample over the comparator actives (0.5 = same; 95% cluster-bootstrap interval; stratified = within heavy-atom strata):")
        for r in [x for x in h3_rows if x["group"] == group]:
            print(f"    {r['sample']:<22} n={r['n_a']:>3}  {r['p_sup']:.2f} [{r['p_sup_low']:.2f}, {r['p_sup_high']:.2f}]   stratified "
                  f"{r['p_sup_stratified']:.2f} [{r['p_sup_strat_low']:.2f}, {r['p_sup_strat_high']:.2f}] (pairs {r['pairs_used']}/{r['pairs_total']})   -> {r['decision']}")

        fig_docking_group(g, group, role_text, config.CAMPAIGN, config.SCOPE, config.FIGURES_DIR / f"08_docking_{group}.png")

    pd.DataFrame(gate_rows).to_csv(outputs[3], index=False)
    pd.DataFrame(h3_rows).to_csv(outputs[4], index=False)
    pd.DataFrame(sens_rows).to_csv(outputs[5], index=False)

    # ---- engine noise over the three Vina seeds ------------------------------------------------------------------
    rep = scores[scores["status"] == "ok"].groupby("smiles").filter(lambda t: t["vina_seed"].nunique() == 3)
    if len(rep):
        wide = rep.pivot(index="smiles", columns="vina_seed", values="vina_score")
        noise = pd.DataFrame({"n_molecules": [len(wide)], "median_sd": [float(wide.std(axis=1).median())],
                              "max_range": [float((wide.max(axis=1) - wide.min(axis=1)).max())],
                              "spearman_seed_pairs": ["; ".join(f"{a}-{b}: {spearmanr(wide[a], wide[b])[0]:.2f}"
                                                                for i, a in enumerate(wide.columns) for b in list(wide.columns)[i + 1:])]})
        noise.to_csv(outputs[6], index=False)
        print(f"\nEngine noise: {len(wide)} molecules docked with 3 seeds: median SD {noise['median_sd'][0]:.3f} kcal/mol, max range {noise['max_range'][0]:.2f}")
    else:
        pd.DataFrame().to_csv(outputs[6], index=False)
    print("\nWrote results/08e_* and figures/08_docking_<group>.png")


if __name__ == "__main__":
    main()
