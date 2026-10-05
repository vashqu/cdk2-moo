#!/usr/bin/env python3
"""
Stage 11b: Table 1, a per-set summary of everything measured on the molecules that were docked, for one policy scope.

Reads : results/08c_docking_sets.csv, 08e_docking_by_molecule.csv, 05e_holdout_sets.csv, data/processed/cdk2_ic50_curated.csv, ecfp4.npy
Writes: results/table1.csv, figures/final/table1.md

One row per docked set (named GROUP:set; groups P, M, S as in Stage 8c). Columns, medians over the set unless noted:
  n            molecules in the set; every other denominator below is stated, never silently the completed jobs only
  n_vina       molecules whose docking job succeeded (Vina is the median over these)
  pActivity    MEASURED for sets drawn from the curated data, PREDICTED by the group's own surrogate for generated molecules and decoys; the
               "basis" column says which, and the two are not comparable
  Vina         best Vina score, kcal/mol (more negative = better); a docking score is not a binding free energy
  QED, SA      drug-likeness (0-1, higher better) and synthetic accessibility (1-10, lower easier)
  MW, cLogP    molecular weight and Crippen logP
  diversity    mean over pairs of (1 - ECFP4 Tanimoto) within the set
  %PAINS       share of the set with a PAINS alert (PAINS = substructures that tend to give false positives in biochemical assays)
  %PB pass     share of the set's molecules whose docked pose passes all 22 PoseBusters checks; molecules without an audit result count as
               not passing (pass / n)
  max Tanimoto raw ECFP4 max-Tanimoto to the group's training molecules (P: all 2,016; M, S: the training part of that holdout). For sets drawn
               from the training data the nearest OTHER training molecule is used (the molecule itself is excluded).

Run (inside a policy scope, after the docking analysis):
    CDK2_CAMPAIGN=<id> CDK2_SCOPE=legacy python scripts/11b_table1.py
"""

import argparse

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

from cdk2moo import config
from cdk2moo.alerts import alert_hits
from cdk2moo.features import ecfp4
from cdk2moo.objectives import qed_and_sa


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    R = config.RESULTS_DIR
    csv_path, md_path = R / "table1.csv", config.FIGURES_DIR / "final" / "table1.md"
    for path in (csv_path, md_path):
        if path.exists() and not args.overwrite:
            raise SystemExit(f"{path} exists and is not replaced; pass --overwrite")
    sets = pd.read_csv(R / "08c_docking_sets.csv")
    docked = pd.read_csv(R / "08e_docking_by_molecule.csv")
    roles = pd.read_csv(R / "05e_holdout_sets.csv")
    curated = pd.read_csv(config.PROCESSED_DIR / "cdk2_ic50_curated.csv")
    train_fps = np.load(config.PROCESSED_DIR / "ecfp4.npy").astype(np.float32)
    index_of = {smiles: i for i, smiles in enumerate(curated["std_smiles"])}
    in_training = {"P": np.ones(len(curated), dtype=bool), "M": (roles["molecule_control"] == "train").to_numpy(),
                   "S": (roles["scaffold_control"] == "train").to_numpy()}

    rows = []
    for name in dict.fromkeys(sets["set"]):
        s = sets[sets["set"] == name].reset_index(drop=True)
        d = docked[docked["set"] == name]
        smiles = list(s["smiles"])
        qed, sa = qed_and_sa(smiles)
        mols = [Chem.MolFromSmiles(x) for x in smiles]
        fps = ecfp4(smiles, config.FP_RADIUS, config.FP_BITS).astype(np.float32)
        shared = fps @ fps.T
        tanimoto = shared / (fps.sum(1)[:, None] + fps.sum(1)[None, :] - shared)
        diversity = 1 - (tanimoto.sum() - len(smiles)) / (len(smiles) * (len(smiles) - 1))
        pains = np.mean([len(alert_hits(x)[0]) > 0 for x in smiles])

        # similarity to the group's training molecules, leaving a molecule out of its own comparison
        reference = train_fps[in_training[name[0]]]
        reference_rows = np.where(in_training[name[0]])[0]
        shared = fps @ reference.T
        sim = shared / (fps.sum(1)[:, None] + reference.sum(1)[None, :] - shared)
        for k, x in enumerate(smiles):
            if s["role"][k] in ("in_sample_top_decile", "random_reference") and x in index_of:
                sim[k, reference_rows == index_of[x]] = -1.0
        measured = s["pactivity"].notna().all()
        rows.append({
            "set": name, "n": len(s), "n_vina": int((d["status"] == "ok").sum()),
            "pActivity": s["pactivity"].median() if measured else d["pred_real"].median(), "basis": "measured" if measured else "predicted",
            "Vina": d["vina_score"].median(), "QED": np.median(qed), "SA": np.median(sa),
            "MW": np.median([Descriptors.MolWt(m) for m in mols]), "cLogP": np.median([Crippen.MolLogP(m) for m in mols]),
            "diversity": diversity, "pct_PAINS": 100 * pains,
            "pct_PB_pass": 100 * (d["pb_status"] == "pass").sum() / len(s), "max_Tanimoto": np.median(sim.max(axis=1))})
    table = pd.DataFrame(rows)
    table.to_csv(csv_path, index=False, float_format="%.3f")

    fmt = table.copy()
    for col, digits in [("pActivity", 2), ("Vina", 2), ("QED", 2), ("SA", 2), ("MW", 0), ("cLogP", 2),
                        ("diversity", 2), ("pct_PAINS", 0), ("pct_PB_pass", 0), ("max_Tanimoto", 2)]:
        fmt[col] = fmt[col].map(lambda v: f"{v:.{digits}f}")
    fmt.columns = ["set", "n", "n with Vina", "pActivity", "basis", "Vina", "QED", "SA", "MW", "cLogP", "diversity", "% PAINS",
                   "% PB pass", "max Tanimoto"]
    print(fmt.to_string(index=False))
    lines = ["| " + " | ".join(fmt.columns) + " |", "|" + "---|" * len(fmt.columns)]
    lines += ["| " + " | ".join(row) + " |" for row in fmt.astype(str).to_numpy()]
    note = (f"\n\nTable 1. Campaign {config.CAMPAIGN}, {config.SCOPE} policy scope. Per-set summary of docked molecules (medians; diversity = mean "
            "pairwise 1 - ECFP4 Tanimoto). Groups: P = surrogate trained on all 2,016 molecules (its comparator actives were in training); "
            "M = molecule-level holdout; S = scaffold-level holdout. pActivity is measured for sets drawn from the curated data and predicted "
            "by the group's surrogate otherwise; not comparable. Vina in kcal/mol (more negative = better), median over the n with Vina "
            "molecules. % PB pass = passes all 22 PoseBusters checks / n in set. Max Tanimoto is raw ECFP4 (Morgan r=2, 2048-bit) to the group's "
            "training molecules, excluding the molecule itself for sets drawn from training. random_training / random_remaining are "
            "stand-ins for a random ChEMBL sample, not a ChEMBL sample.\n")
    md_path.write_text("\n".join(lines) + note)
    print(f"\nWrote {csv_path.name} and {md_path.name}")


if __name__ == "__main__":
    main()
