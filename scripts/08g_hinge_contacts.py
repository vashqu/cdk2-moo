#!/usr/bin/env python3
"""
Stage 8g: how close do docked ligand N and O atoms come to the kinase hinge backbone? (distance only)

Reads : results/08d_export/poses_seed42.sdf, structures/receptor_H.pdb, structures/4kd1_1QK_crystal.sdf, results/08c_docking_sets.csv,
        results/08d_export/docking_scores.csv
Writes: results/08g_hinge_contacts.csv   one row per docked molecule with a pose: nearest ligand N/O distance to the hinge atoms
        results/08g_hinge_by_set.csv     per set: molecules in the set, with a pose, with a contact, and the two shares

A kinase's ATP pocket is anchored by the "hinge", the strand joining the kinase's two lobes. ATP-competitive inhibitors, including dinaciclib, usually make
hydrogen bonds to the hinge backbone. This script measures ONLY DISTANCE: the smallest distance between a ligand nitrogen or oxygen atom and the backbone
atoms Glu81 O, Leu83 N and Leu83 O (the pairs dinaciclib uses). A distance <= 3.5 A is called "hinge proximity". It does not check donor/acceptor roles,
angles or hydrogen positions, so it is not a hydrogen-bond detector, and a contact with a non-polar-looking pose is not evidence of binding.
Two shares are given, with their denominators: contacts / molecules in the set (jobs without a pose count as no contact), and contacts / molecules with a pose.

Run (inside a policy scope, after 08d):
    CDK2_CAMPAIGN=<id> CDK2_SCOPE=legacy python scripts/08g_hinge_contacts.py
"""

import numpy as np
import pandas as pd
from rdkit import Chem

from cdk2moo import config
from cdk2moo.hinge import hinge_atoms, nearest_polar_distance

CUTOFF = 3.5


def main():
    config.require_campaign()
    R = config.RESULTS_DIR
    for path in (R / "08g_hinge_contacts.csv", R / "08g_hinge_by_set.csv"):
        if path.exists():
            raise SystemExit(f"{path} exists; outputs are never replaced in place")
    protein = Chem.MolFromPDBFile(str(config.STRUCTURES_DIR / "receptor_H.pdb"), removeHs=False)
    hinge_xyz = hinge_atoms(protein)

    crystal = Chem.MolFromMolFile(str(config.STRUCTURES_DIR / "4kd1_1QK_crystal.sdf"))
    print(f"crystal dinaciclib: nearest ligand N/O to the 3 hinge atoms {nearest_polar_distance(crystal, hinge_xyz):.2f} A (reference)")

    rows = []
    for mol in Chem.SDMolSupplier(str(R / "08d_export" / f"poses_seed{config.RANDOM_SEED}.sdf"), removeHs=False):
        if mol is not None:
            rows.append({"smiles": mol.GetProp("smiles"), "hinge_distance": nearest_polar_distance(mol, hinge_xyz)})
    hinge = pd.DataFrame(rows)
    if hinge["smiles"].duplicated().any():
        raise SystemExit("duplicate poses in the export")
    hinge["hinge_contact"] = hinge["hinge_distance"] <= CUTOFF
    hinge.to_csv(R / "08g_hinge_contacts.csv", index=False)

    sets = pd.read_csv(R / "08c_docking_sets.csv")
    df = sets.merge(hinge, on="smiles", how="left", validate="many_to_one")
    df["has_pose"] = df["hinge_distance"].notna()
    df["hinge_contact"] = df["hinge_contact"].fillna(False).astype(bool)
    table = df.groupby("set", sort=False).agg(n_in_set=("smiles", "size"), n_with_pose=("has_pose", "sum"), n_contact=("hinge_contact", "sum"),
                                              median_distance=("hinge_distance", "median")).reset_index()
    table["share_of_set"] = table["n_contact"] / table["n_in_set"]
    table["share_of_posed"] = table["n_contact"] / table["n_with_pose"]
    table.to_csv(R / "08g_hinge_by_set.csv", index=False)
    print(f"\nLigand N/O within {CUTOFF} A of Glu81 O / Leu83 N / Leu83 O (distance only; not a verified hydrogen bond)")
    print(table.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print(f"\nWrote {R / '08g_hinge_contacts.csv'} and 08g_hinge_by_set.csv")


if __name__ == "__main__":
    main()
