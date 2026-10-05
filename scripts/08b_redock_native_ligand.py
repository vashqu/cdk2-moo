#!/usr/bin/env python3
"""
Stage 8b: redock dinaciclib into 4KD1 and compare with the crystal pose.

Reads : data/structures/receptor.pdbqt, receptor_box.txt, 4kd1_1QK_crystal.sdf
Writes: results/08b_redocking.json
        data/structures/redock_poses_seed<N>.sdf   <- every docked pose, for inspection

The test. Take the ligand that the crystal shows in the pocket, forget where it was,
let Vina search for the best placement, and measure how far the top-ranked pose
lies from the crystal pose: the root-mean-square deviation (RMSD) over heavy atoms,
in Angstrom, accounting for symmetry (a benzene ring flipped 180 degrees is the same
molecule). Convention: <= 2 A means the docking setup can reproduce a known binding
mode. If it cannot, no docking score downstream means anything, and that must be
reported as such (CLAUDE.md control 7).

Also reported, because a failure has two possible causes:
  - Vina's score of the CRYSTAL pose itself. If the crystal pose scores WORSE than
    the best docked pose, the scoring function prefers a wrong pose (a scoring
    problem). If it scores better but the search missed it, that is a search problem.
  - Three different Vina random seeds (engine noise, control 9).
Protonation: the ligand is docked as deposited, a 1-hydroxypyridinium cation (N-OH).

Run:
    python scripts/08b_redock_native_ligand.py
"""

import argparse
import json
import time

import numpy as np
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from vina import Vina

from cdk2moo import config


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=config.RANDOM_SEED)
    args = ap.parse_args()
    seeds = [args.seed, args.seed + 1, args.seed + 2]
    print(f"Vina seeds = {seeds}, exhaustiveness {config.DOCK_EXHAUSTIVENESS}, box {config.DOCK_BOX_SIZE} A")

    box = {k: v for k, *v in (line.split() for line in open(config.STRUCTURES_DIR / "receptor_box.txt"))}
    centre = [float(x) for x in box["center"]]
    size = [float(x) for x in box["size"]]

    # ---- ligand: crystal coordinates, hydrogens added, converted to PDBQT --------
    crystal = Chem.MolFromMolFile(str(config.STRUCTURES_DIR / "4kd1_1QK_crystal.sdf"), removeHs=False)
    reference = Chem.RemoveHs(crystal)
    crystal_h = Chem.AddHs(crystal, addCoords=True)
    setup = MoleculePreparation().prepare(crystal_h)[0]
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise SystemExit(f"Meeko could not prepare the ligand: {error}")
    print(f"Ligand: {reference.GetNumHeavyAtoms()} heavy atoms, {setup.flexibility_model['rigid_body_connectivity'] and len(setup.flexibility_model['rigid_body_connectivity'])} rigid groups")

    receptor = str(config.STRUCTURES_DIR / "receptor.pdbqt")

    # ---- the crystal pose, scored as it is and after local relaxation ------------
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(receptor)
    v.set_ligand_from_string(pdbqt)
    v.compute_vina_maps(center=centre, box_size=size)
    crystal_score = float(v.score()[0])
    crystal_relaxed = float(v.optimize()[0])
    print(f"Vina score of the crystal pose: {crystal_score:.2f} kcal/mol as deposited, "
          f"{crystal_relaxed:.2f} after local relaxation")

    # ---- redocking, three seeds --------------------------------------------------
    report = {"crystal_score": crystal_score, "crystal_score_relaxed": crystal_relaxed, "seeds": {}}
    top_rmsd = []
    for seed in seeds:
        v = Vina(sf_name="vina", seed=seed, verbosity=0)
        v.set_receptor(receptor)
        v.set_ligand_from_string(pdbqt)
        v.compute_vina_maps(center=centre, box_size=size)
        start = time.time()
        v.dock(exhaustiveness=config.DOCK_EXHAUSTIVENESS, n_poses=config.DOCK_N_POSES)
        seconds = time.time() - start
        energies = v.energies(n_poses=config.DOCK_N_POSES)[:, 0]
        poses = PDBQTMolecule(v.poses(n_poses=config.DOCK_N_POSES), is_dlg=False, skip_typing=True)
        mol = RDKitMolCreate.from_pdbqt_mol(poses)[0]
        mol_heavy = Chem.RemoveHs(mol)

        rmsds = [rdMolAlign.CalcRMS(mol_heavy, reference, prbId=i) for i in range(mol_heavy.GetNumConformers())]
        writer = Chem.SDWriter(str(config.STRUCTURES_DIR / f"redock_poses_seed{seed}.sdf"))
        for i in range(mol.GetNumConformers()):
            mol.SetProp("vina_score", f"{energies[i]:.3f}")
            mol.SetProp("rmsd_to_crystal", f"{rmsds[i]:.3f}")
            writer.write(mol, confId=i)
        writer.close()

        top_rmsd.append(rmsds[0])
        print(f"  seed {seed}: {seconds:.0f} s; top pose {energies[0]:.2f} kcal/mol, RMSD {rmsds[0]:.2f} A; "
              f"best RMSD among {len(rmsds)} poses {min(rmsds):.2f} A (pose {int(np.argmin(rmsds)) + 1}, "
              f"score {energies[int(np.argmin(rmsds))]:.2f})")
        report["seeds"][seed] = {"seconds": seconds, "top_score": float(energies[0]), "top_rmsd": float(rmsds[0]),
                                 "best_rmsd": float(min(rmsds)), "scores": energies.tolist(), "rmsds": rmsds}

    median = float(np.median(top_rmsd))
    print(f"\nTop-ranked pose RMSD to the crystal pose: median {median:.2f} A over {len(seeds)} seeds "
          f"(range {min(top_rmsd):.2f}-{max(top_rmsd):.2f})")
    if median <= config.REDOCK_MAX_RMSD:
        print(f"  <= {config.REDOCK_MAX_RMSD} A: the setup reproduces the crystal binding mode.")
    else:
        print(f"  > {config.REDOCK_MAX_RMSD} A: REDOCKING FAILED. No downstream docking result may be trusted. "
              f"Diagnose, or fall back to {config.PDB_ID_FALLBACK}.")
    report["top_rmsd_median"] = median
    report["passed"] = bool(median <= config.REDOCK_MAX_RMSD)
    with open(config.RESULTS_DIR / "08b_redocking.json", "w") as fh:
        json.dump(report, fh, indent=2)
    print("Wrote results/08b_redocking.json")


if __name__ == "__main__":
    main()
