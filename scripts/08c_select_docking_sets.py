#!/usr/bin/env python3
"""
Stage 8c: choose the molecules to dock, by experiment group. Rules are fixed by seed; nothing is chosen after seeing a docking score.

Reads : processed/cdk2_ic50_curated.csv, results/05e_holdout_sets.csv, results/05_ga_populations.csv (primary GA),
        results/05f_molecule_populations.csv, results/05f_scaffold_populations.csv, results/06_pareto_front.csv  (the selected scope's own)
Writes: results/08c_docking_sets.csv        one row per (molecule, set), with group, role, policy, cluster and replicate flag
        results/08c_selection_report.json   counts, pool sizes, decoy matching quality

Groups P, M and S, and why they exist, are explained in cdk2moo/docksets.py: only M and S have a comparator that the surrogate that drove
the GA never saw (held-out actives). P's actives are the top decile of the training set, which that surrogate DID see ("in_sample_top_decile").

Sets of 100 molecules each (the Pareto front, and group S's held-out actives, are their own sizes):
  P:training_top_decile  P:random_training  P:decoys  P:multi_real P:multi_scrambled P:activity_only P:druglike_only  P:pareto_front
  M:heldout_actives      M:random_remaining M:decoys  M:<four arms>
  S:heldout_actives (all 93)  S:random_remaining  S:decoys  S:<four arms>
Generated sets are random draws of distinct molecules from the final generation (birth_generation > 0) of each arm, five seeds pooled.
Decoys are property-matched to the group's sampled actives (heavy atoms, MW, cLogP, donors, acceptors, rotatable bonds) from molecules of the group's
druglike_only and multi_scrambled arms (any generation), each with ECFP4 Tanimoto < 0.5 to EVERY active of the group; nothing guarantees a decoy is inactive.
Replicates: 6 molecules from each of the actives set, the random set, multi_real, activity_only and decoys of each group are flagged for two extra Vina seeds.
`cluster` is the unit of replication used later for uncertainty: the GA seed for generated molecules and decoys, the Murcko scaffold for reference molecules.

Run (inside a policy scope):
    CDK2_CAMPAIGN=<id> CDK2_SCOPE=legacy python scripts/08c_select_docking_sets.py
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.docksets import ARMS, GROUP_ROLE, keyed_rng, match_decoys
from cdk2moo.features import ecfp4, murcko_scaffold

TRAINING_DESCRIPTION = {"P": "all curated molecules (the comparator actives were in the training data)",
                        "M": "curated molecules without the 204 top-decile actives",
                        "S": "curated molecules without 9 scaffold families (260 molecules, 93 top-decile actives)"}


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=config.RANDOM_SEED)
    ap.add_argument("--out", default=None, help="sets CSV (default results/08c_docking_sets.csv); an existing file is never replaced")
    args = ap.parse_args()
    out = Path(args.out) if args.out else config.RESULTS_DIR / "08c_docking_sets.csv"
    if out.exists():
        raise SystemExit(f"{out} exists; docking sets are never replaced in place (use --out for another file)")
    policy = config.SCOPE_POLICY or "legacy"

    R = config.RESULTS_DIR
    df = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    roles = pd.read_csv(R / "05e_holdout_sets.csv")
    if len(roles) != len(df) or not (roles["inchikey"].to_numpy() == df["inchikey"].to_numpy()).all():
        raise SystemExit("holdout assignments do not line up with the curated molecules")
    pactivity = df["pactivity"].to_numpy()
    smiles_of = df["std_smiles"].tolist()
    scaffold_of = [murcko_scaffold(s) for s in smiles_of]
    cutoff = np.quantile(pactivity, config.HOLDOUT_QUANTILE)
    top = np.where(pactivity >= cutoff)[0]
    held_m = np.where((roles["molecule_control"] == "held") & roles["is_top_decile"])[0]
    held_s = np.where((roles["scaffold_control"] == "held") & roles["is_top_decile"])[0]
    train_m = np.where(roles["molecule_control"] == "train")[0]
    train_s = np.where(roles["scaffold_control"] == "train")[0]

    populations = {"P": pd.read_csv(R / "05_ga_populations.csv"), "M": pd.read_csv(R / "05f_molecule_populations.csv"),
                   "S": pd.read_csv(R / "05f_scaffold_populations.csv")}
    front = pd.read_csv(R / "06_pareto_front.csv")
    rows, report = [], {"policy": policy, "seed": args.seed, "top_decile_cutoff": float(cutoff)}

    # ---- reference sets: identical in every scope (they depend on the data, the holdout and the seed only) ----
    actives = {"P": np.sort(keyed_rng(args.seed, "P:training_top_decile").choice(top, 100, replace=False)),
               "M": np.sort(keyed_rng(args.seed, "M:heldout_actives").choice(held_m, 100, replace=False)),
               "S": np.sort(held_s)}
    others_p = np.setdiff1d(np.arange(len(df)), actives["P"])
    random_ids = {"P": np.sort(keyed_rng(args.seed, "P:random_training").choice(others_p, 100, replace=False)),
                  "M": np.sort(keyed_rng(args.seed, "M:random_remaining").choice(train_m, 100, replace=False)),
                  "S": np.sort(keyed_rng(args.seed, "S:random_remaining").choice(train_s, 100, replace=False))}
    reference_sets = [("P", "training_top_decile", actives["P"], GROUP_ROLE["P"]), ("M", "heldout_actives", actives["M"], GROUP_ROLE["M"]),
                      ("S", "heldout_actives", actives["S"], GROUP_ROLE["S"]), ("P", "random_training", random_ids["P"], "random_reference"),
                      ("M", "random_remaining", random_ids["M"], "random_reference"), ("S", "random_remaining", random_ids["S"], "random_reference")]
    for group, name, ids, role in reference_sets:
        for i in ids:
            rows.append({"smiles": smiles_of[i], "group": group, "set": f"{group}:{name}", "role": role, "policy": "shared",
                         "pactivity": pactivity[i], "source": f"curated row {i}", "cluster": scaffold_of[i], "replicate": False})

    # ---- generated sets and decoys, from this scope's GA populations --------------------------------------
    all_active_ids = {"P": top, "M": np.where((roles["molecule_control"] == "held") & roles["is_top_decile"])[0], "S": held_s}
    for group, pop in populations.items():
        gen = pop[pop["birth_generation"] > 0]
        final = gen[gen["generation"] == gen["generation"].max()]
        taken = set()
        for arm in ARMS:
            pool = final[final["arm"] == arm].drop_duplicates("smiles").reset_index(drop=True)
            pick = pool.iloc[np.sort(keyed_rng(args.seed, group, arm).choice(len(pool), 100, replace=False))]
            taken |= set(pick["smiles"])
            for _, r in pick.iterrows():
                rows.append({"smiles": r["smiles"], "group": group, "set": f"{group}:{arm}", "role": "generated", "policy": policy,
                             "pactivity": np.nan, "source": f"arm {arm}, seed {r['seed']}, final generation", "cluster": f"seed{r['seed']}",
                             "replicate": False})
        if group == "P":
            for _, r in front.iterrows():
                rows.append({"smiles": r["smiles"], "group": "P", "set": "P:pareto_front", "role": "generated", "policy": policy,
                             "pactivity": np.nan, "source": f"pooled front ({r['arms']})", "cluster": "front", "replicate": False})
        pool = (gen[gen["arm"].isin(["druglike_only", "multi_scrambled"])].drop_duplicates("smiles")
                .pipe(lambda d: d[~d["smiles"].isin(taken)]).reset_index(drop=True))
        chosen, smd = match_decoys([smiles_of[i] for i in actives[group]], pool["smiles"].tolist(),
                                   ecfp4([smiles_of[i] for i in all_active_ids[group]]), keyed_rng(args.seed, group, "decoys"))
        for j in chosen:
            rows.append({"smiles": pool["smiles"][j], "group": group, "set": f"{group}:decoys", "role": "decoy", "policy": policy,
                         "pactivity": np.nan, "source": f"matched decoy, arm {pool['arm'][j]}, seed {pool['seed'][j]}",
                         "cluster": f"seed{pool['seed'][j]}", "replicate": False})
        report[f"{group}_decoy_standardized_mean_differences"] = {k: float(v) for k, v in smd.items()}
        report[f"{group}_decoy_pool_after_exclusions"] = int(len(pool))

    sets = pd.DataFrame(rows)
    sets["surrogate_training"] = sets["group"].map(TRAINING_DESCRIPTION)
    for group in "PMS":
        for name in [f"{group}:{n}" for n in (["training_top_decile", "random_training"] if group == "P" else
                     ["heldout_actives", "random_remaining"]) + ["multi_real", "activity_only", "decoys"]]:
            idx = sets.index[sets["set"] == name]
            sets.loc[keyed_rng(args.seed, name, "replicate").choice(idx, 6, replace=False), "replicate"] = True
    sets.to_csv(out, index=False)
    report["set_sizes"] = sets.groupby("set").size().to_dict()
    report["distinct_molecules"] = int(sets["smiles"].nunique())
    report["replicate_molecules"] = int(sets.loc[sets["replicate"], "smiles"].nunique())
    (out.with_name(out.stem.replace("docking_sets", "selection_report") + ".json")).write_text(json.dumps(report, indent=2))
    print(f"policy scope: {policy}; groups P (all-data surrogate), M (molecule-level holdout), S (scaffold-level holdout)")
    print(sets.groupby(["group", "set", "role"]).size().to_string())
    print(f"\n{len(sets)} rows, {report['distinct_molecules']} distinct molecules; {report['replicate_molecules']} replicate molecules")
    print("decoy matching (standardized mean differences, decoys minus actives):")
    for group in "PMS":
        print(f"  {group}: " + ", ".join(f"{k} {v:+.2f}" for k, v in report[f"{group}_decoy_standardized_mean_differences"].items()))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
