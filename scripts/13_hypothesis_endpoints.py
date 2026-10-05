#!/usr/bin/env python3
"""
Stage 13: the pre-specified endpoints for H1, H2, H4 and H5 (attenuation and proximity), for the selected policy scope.

Reads : results/05_ga_populations_scored.csv (primary GA with size percentile), results/05f_molecule_populations.csv, 05f_scaffold_populations.csv,
        results/10_constrained_{primary,molecule,scaffold}_populations.csv, results/05e_holdout_sets.csv, processed/ (curated data, fingerprints)
Writes: results/13_h1_runs.csv, 13_h1_summary.csv, 13_h2_runs.csv, 13_h2_summary.csv, 13_h4_runs.csv, 13_h4_summary.csv, 13_h5_runs.csv, 13_h5_summary.csv,
        13_final_generation.csv

Definitions are in the campaign plan, section 3 (src/cdk2moo/endpoints.py has the trajectory statistics). A RUN is one arm x one GA seed; n = 5 runs per arm
and condition; comparisons are paired by seed. Statistics about generated chemistry use molecules born after generation 0; a generation needs at least 10 such
molecules to enter step and trend statistics. Nothing here pools policies, experiments or training subsets.

H5 transfer in the plan's AMENDED, fingerprint sense is the share of the run's final generated molecules with ECFP4 Tanimoto >= 0.6 to a held-out active. It is
a closeness measure, not activity. H5 transfer in H3 terms needs docking and is handled by the docking stages (and only if the validity gate passes).

Run (inside a policy scope, after the GA stages):
    CDK2_CAMPAIGN=<id> CDK2_SCOPE=legacy python scripts/13_hypothesis_endpoints.py
"""

import numpy as np
import pandas as pd

from cdk2moo import config, endpoints as ep
from cdk2moo.features import ecfp4
from cdk2moo.recovery import tanimoto_matrix
from cdk2moo.surrogate import fit_forest, scramble_labels

ARMS = ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]
FOLLOWED = {"multi_real": "pred_real", "activity_only": "pred_real", "multi_scrambled": "pred_scrambled", "druglike_only": None}


def main():
    config.require_campaign()
    R = config.RESULTS_DIR
    names = ["13_h1_runs", "13_h1_summary", "13_h2_runs", "13_h2_summary", "13_h4_runs", "13_h4_summary", "13_h5_runs", "13_h5_summary", "13_final_generation"]
    for n in names:
        if (R / f"{n}.csv").exists():
            raise SystemExit(f"{R / (n + '.csv')} exists; outputs are never replaced in place")
    print(f"policy scope: {config.SCOPE}; campaign {config.CAMPAIGN}")

    pop = pd.read_csv(R / "05_ga_populations_scored.csv")
    seeds = sorted(pop["seed"].unique())
    last = int(pop["generation"].max())
    tables = {}
    for arm in ARMS:
        for seed in seeds:
            run = pop[(pop["arm"] == arm) & (pop["seed"] == seed)]
            tables[(arm, seed)] = {c: ep.generation_table(run, c) for c in ("pred_real", "pred_scrambled", "max_tanimoto", "size_pctile", "qed", "sa", "n_heavy_atoms")}

    # ---- H1: does predicted activity rise (and monotonically)? -----------------------------------------------------
    rows = [{"arm": a, "seed": s, **ep.trend_statistics(tables[(a, s)]["pred_real"])} for a in ARMS for s in seeds]
    h1 = pd.DataFrame(rows)
    h1.to_csv(R / "13_h1_runs.csv", index=False)
    summary = []
    for arm in ARMS:
        d = h1[h1["arm"] == arm]
        summary.append({"arm": arm, "runs": len(d), "net_rise_median": d["net_rise"].median(), "net_rise_min": d["net_rise"].min(), "net_rise_max": d["net_rise"].max(),
                        "runs_with_rise": int((d["net_rise"] > 0).sum()), "runs_spearman_ge_0.9": int((d["spearman_generation"] >= 0.9).sum()),
                        "runs_literally_monotone": int((d["fraction_nondecreasing"] == 1).sum()), "fraction_nondecreasing_median": d["fraction_nondecreasing"].median(),
                        "rises_by_plan": bool((d["net_rise"] > 0).all() and (d["spearman_generation"] >= 0.9).sum() >= 4),
                        "monotone_by_plan": bool((d["fraction_nondecreasing"] == 1).all())})
    pd.DataFrame(summary).to_csv(R / "13_h1_summary.csv", index=False)
    print("\nH1  predicted pActivity (real surrogate) of generated molecules over generations; 5 runs per arm")
    for r in summary:
        print(f"  {r['arm']:<16} net rise median {r['net_rise_median']:+.2f} [{r['net_rise_min']:+.2f}, {r['net_rise_max']:+.2f}]; rises in {r['runs_with_rise']}/{r['runs']} runs, "
              f"Spearman >= 0.9 in {r['runs_spearman_ge_0.9']}, literally monotone in {r['runs_literally_monotone']}  -> rises: {r['rises_by_plan']}, monotone: {r['monotone_by_plan']}")

    # ---- H2: does similarity to the training set fall as predicted activity climbs? ---------------------------------------
    rows = []
    for arm in ARMS:
        for seed in seeds:
            t = tables[(arm, seed)]
            rho, n = ep.association_across_generations(t["pred_real"], t["max_tanimoto"])
            base = tables[("druglike_only", seed)]
            final = lambda table, col: [r for r in table[col] if r["generation"] == last][0]["median"]    # noqa: E731
            rows.append({"arm": arm, "seed": seed, "spearman_activity_vs_similarity": rho, "n_generations": n,
                         "final_similarity": final(t, "max_tanimoto"), "baseline_final_similarity": final(base, "max_tanimoto"),
                         "delta_similarity_vs_druglike": final(t, "max_tanimoto") - final(base, "max_tanimoto"),
                         "final_size_percentile": final(t, "size_pctile"), "delta_size_percentile_vs_druglike": final(t, "size_pctile") - final(base, "size_pctile"),
                         "final_heavy_atoms": final(t, "n_heavy_atoms"), "delta_heavy_atoms_vs_druglike": final(t, "n_heavy_atoms") - final(base, "n_heavy_atoms")})
    h2 = pd.DataFrame(rows)
    h2.to_csv(R / "13_h2_runs.csv", index=False)
    summary = []
    for arm in ARMS:
        d = h2[h2["arm"] == arm]
        summary.append({"arm": arm, "runs_negative_literal": int((d["spearman_activity_vs_similarity"] < 0).sum()),
                        "spearman_median": d["spearman_activity_vs_similarity"].median(),
                        "delta_similarity_median": d["delta_similarity_vs_druglike"].median(), "delta_similarity_min": d["delta_similarity_vs_druglike"].min(),
                        "delta_similarity_max": d["delta_similarity_vs_druglike"].max(), "runs_farther_than_baseline": int((d["delta_similarity_vs_druglike"] < 0).sum()),
                        "delta_size_percentile_median": d["delta_size_percentile_vs_druglike"].median(), "delta_heavy_atoms_median": d["delta_heavy_atoms_vs_druglike"].median(),
                        "literal_by_plan": bool((d["spearman_activity_vs_similarity"] < 0).sum() >= 4),
                        "attributable_by_plan": bool((d["delta_similarity_vs_druglike"] < 0).sum() >= 4)})
    pd.DataFrame(summary).to_csv(R / "13_h2_summary.csv", index=False)
    print("\nH2  (literal) Spearman between predicted activity and similarity across generations 1-50; (attributable) final similarity minus the druglike_only baseline of the same seed")
    for r in summary:
        print(f"  {r['arm']:<16} Spearman median {r['spearman_median']:+.2f}, negative in {r['runs_negative_literal']}/5 -> literal: {r['literal_by_plan']};  "
              f"delta similarity median {r['delta_similarity_median']:+.3f} [{r['delta_similarity_min']:+.3f}, {r['delta_similarity_max']:+.3f}], farther in "
              f"{r['runs_farther_than_baseline']}/5 -> attributable: {r['attributable_by_plan']}  (delta heavy atoms {r['delta_heavy_atoms_median']:+.1f})")

    # ---- H4: do real and scrambled surrogates give similar-looking curves, and does the real prediction move? ------------------
    curated = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    train_fps = np.load(config.PROCESSED_DIR / "ecfp4.npy")
    y = curated["pactivity"].to_numpy()
    real = fit_forest(train_fps, y, config.SURROGATE_SEED)
    scrambled = fit_forest(train_fps, scramble_labels(y, config.SCRAMBLE_SEED), config.SCRAMBLE_SEED)
    sd = {"pred_real": float(np.std(real.predict(train_fps))), "pred_scrambled": float(np.std(scrambled.predict(train_fps)))}
    rows = []
    for arm in ARMS:
        column = FOLLOWED[arm]
        for seed in seeds:
            followed = ep.trend_statistics(tables[(arm, seed)][column or "pred_real"])
            base = [r for r in tables[("druglike_only", seed)]["pred_real"] if r["generation"] == last][0]["median"]
            real_final = [r for r in tables[(arm, seed)]["pred_real"] if r["generation"] == last][0]["median"]
            rows.append({"arm": arm, "seed": seed, "followed": column or "none (real shown)", "net_rise_followed": followed["net_rise"],
                         "spearman_followed": followed["spearman_generation"], "surrogate_sd_on_training": sd[column or "pred_real"],
                         "standardized_rise": followed["net_rise"] / sd[column or "pred_real"], "real_prediction_minus_druglike": real_final - base})
    h4 = pd.DataFrame(rows)
    h4.to_csv(R / "13_h4_runs.csv", index=False)
    summary = []
    for arm in ARMS:
        d = h4[h4["arm"] == arm]
        summary.append({"arm": arm, "followed": d["followed"].iloc[0], "net_rise_median": d["net_rise_followed"].median(), "standardized_rise_median": d["standardized_rise"].median(),
                        "runs_spearman_ge_0.9": int((d["spearman_followed"] >= 0.9).sum()), "real_minus_druglike_median": d["real_prediction_minus_druglike"].median(),
                        "real_minus_druglike_min": d["real_prediction_minus_druglike"].min(), "real_minus_druglike_max": d["real_prediction_minus_druglike"].max(),
                        "runs_real_above_baseline": int((d["real_prediction_minus_druglike"] > 0).sum())})
    summary = pd.DataFrame(summary)
    summary.to_csv(R / "13_h4_summary.csv", index=False)
    mr, ms = summary.set_index("arm").loc["multi_real"], summary.set_index("arm").loc["multi_scrambled"]
    shape = bool(mr["runs_spearman_ge_0.9"] >= 4 and ms["runs_spearman_ge_0.9"] >= 4)
    print(f"\nH4  followed-prediction curves (surrogate SD over its own training molecules: real {sd['pred_real']:.3f}, scrambled {sd['pred_scrambled']:.3f})")
    for _, r in summary.iterrows():
        print(f"  {r['arm']:<16} follows {r['followed']:<18} net rise {r['net_rise_median']:+.2f} (standardized {r['standardized_rise_median']:+.2f}), Spearman >= 0.9 in {r['runs_spearman_ge_0.9']}/5; "
              f"real prediction minus druglike_only: {r['real_minus_druglike_median']:+.2f} [{r['real_minus_druglike_min']:+.2f}, {r['real_minus_druglike_max']:+.2f}], above in {r['runs_real_above_baseline']}/5")
    print(f"  similar-looking curve by the plan's shape criterion (Spearman >= 0.9 in >= 4/5 runs for both arms): {shape}; "
          f"standardized-rise ratio scrambled / real = {ms['standardized_rise_median'] / mr['standardized_rise_median']:.2f} (no threshold)")

    # ---- H5: attenuation and (amended, fingerprint) proximity transfer -------------------------------------------------------
    roles = pd.read_csv(R / "05e_holdout_sets.csv")
    unconstrained = {"primary": R / "05_ga_populations.csv", "molecule": R / "05f_molecule_populations.csv", "scaffold": R / "05f_scaffold_populations.csv"}
    rows = []
    for experiment, path in unconstrained.items():
        base = pd.read_csv(path)
        cons = pd.read_csv(R / f"10_constrained_{experiment}_populations.csv")
        held = np.where((roles[f"{experiment}_control"] == "held") & roles["is_top_decile"])[0] if experiment != "primary" else None
        held_fps = train_fps[held] if held is not None else None
        for arm in ("multi_real", "activity_only"):
            for floor in sorted(cons["min_similarity"].dropna().unique()):
                for seed in seeds:
                    twin = base[(base["arm"] == arm) & (base["seed"] == seed) & (base["generation"] == last) & (base["birth_generation"] > 0)]
                    run = cons[(cons["arm"] == arm) & (cons["seed"] == seed) & np.isclose(cons["min_similarity"], floor) & (cons["generation"] == last) & (cons["birth_generation"] > 0)]
                    row = {"experiment": experiment, "arm": arm, "floor": floor, "seed": seed, "n_generated_constrained": len(run), "n_generated_twin": len(twin),
                           "pred_constrained": run["pred_real"].median(), "pred_twin": twin["pred_real"].median(),
                           "delta_pred": run["pred_real"].median() - twin["pred_real"].median(),
                           "delta_qed": run["qed"].median() - twin["qed"].median(), "delta_heavy_atoms": run["n_heavy_atoms"].median() - twin["n_heavy_atoms"].median()}
                    if held_fps is not None:
                        share = lambda frame: float((tanimoto_matrix(ecfp4(frame["smiles"].tolist()).astype(np.float32), held_fps.astype(np.float32)).max(axis=1) >= config.REDISCOVERY_SIM).mean())   # noqa: E731
                        row.update({"proximity_constrained": share(run), "proximity_twin": share(twin), "delta_proximity_pp": 100 * (share(run) - share(twin))})
                    rows.append(row)
    h5 = pd.DataFrame(rows)
    h5.to_csv(R / "13_h5_runs.csv", index=False)
    summary = []
    for (experiment, arm, floor), d in h5.groupby(["experiment", "arm", "floor"]):
        row = {"experiment": experiment, "arm": arm, "floor": floor, "seeds": len(d), "delta_pred_median": d["delta_pred"].median(), "delta_pred_min": d["delta_pred"].min(),
               "delta_pred_max": d["delta_pred"].max(), "runs_attenuated": int((d["delta_pred"] < 0).sum()), "attenuation_by_plan": bool((d["delta_pred"] < 0).sum() >= 4),
               "delta_qed_median": d["delta_qed"].median(), "delta_heavy_atoms_median": d["delta_heavy_atoms"].median()}
        if "delta_proximity_pp" in d:
            point, low, high, defined = ep.bootstrap_ratio(d["delta_proximity_pp"].fillna(0), d["delta_pred"])
            row.update({"delta_proximity_pp_median": d["delta_proximity_pp"].median(), "delta_proximity_pp_min": d["delta_proximity_pp"].min(),
                        "delta_proximity_pp_max": d["delta_proximity_pp"].max(), "exchange_pp_per_pActivity": point, "exchange_low": low, "exchange_high": high,
                        "exchange_defined_share_of_resamples": defined})
        summary.append(row)
    summary = pd.DataFrame(summary)
    summary.to_csv(R / "13_h5_summary.csv", index=False)
    print("\nH5  constrained minus unconstrained twin (paired by seed); attenuation = predicted activity lower in >= 4/5 seeds")
    for r in summary.itertuples():
        prox = (f"; proximity {r.delta_proximity_pp_median:+.1f} pp [{r.delta_proximity_pp_min:+.1f}, {r.delta_proximity_pp_max:+.1f}], exchange "
                f"{'undefined' if not np.isfinite(r.exchange_pp_per_pActivity) else format(r.exchange_pp_per_pActivity, '+.1f') + ' pp per pActivity'}") if hasattr(r, "delta_proximity_pp_median") and np.isfinite(r.delta_proximity_pp_median) else ""
        print(f"  {r.experiment:<9} {r.arm:<14} floor {r.floor:.1f}: delta predicted {r.delta_pred_median:+.2f} [{r.delta_pred_min:+.2f}, {r.delta_pred_max:+.2f}], attenuated in {r.runs_attenuated}/{r.seeds}; "
              f"delta QED {r.delta_qed_median:+.2f}, atoms {r.delta_heavy_atoms_median:+.0f}{prox}")

    # ---- final-generation description, per arm (secondary) -----------------------------------------------------------------------
    final_rows = []
    for arm in ARMS:
        for seed in seeds:
            final = pop[(pop["arm"] == arm) & (pop["seed"] == seed) & (pop["generation"] == last) & (pop["birth_generation"] > 0)]
            final_rows.append({"arm": arm, "seed": seed, "n_generated": len(final), **{c: final[c].median() for c in ("pred_real", "pred_scrambled", "max_tanimoto", "size_pctile", "qed", "sa", "n_heavy_atoms")},
                               "tautomer_incomplete_share": float((final["tautomer_status"].fillna("") .isin(["MaxTransformsReached", "MaxTautomersReached", "Canceled"])).mean())})
    pd.DataFrame(final_rows).to_csv(R / "13_final_generation.csv", index=False)
    print("\nWrote results/13_*.csv")


if __name__ == "__main__":
    main()
