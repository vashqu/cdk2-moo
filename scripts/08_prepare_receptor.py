#!/usr/bin/env python3
"""
Stage 8a: prepare the CDK2 receptor (PDB 4KD1) for Vina docking.

Writes: data/structures/4KD1.pdb                 <- downloaded (raw, never edited)
        data/structures/4kd1_1QK_crystal.sdf     <- dinaciclib's crystal pose
        data/structures/receptor_H.pdb           <- protein only, hydrogens added
        data/structures/receptor.pdbqt           <- what Vina reads
        data/structures/receptor_box.txt         <- search box (centre, size)

Steps and what each means:
  1. Keep the protein (chain A). Drop the ligand 1QK, the cryoprotectant EDO and all
     waters. Waters are dropped because docking programs have no good model of which
     crystal waters stay when a different ligand binds (a limitation, not a fix).
  2. Add hydrogens at pH 7.4 (pdbfixer/OpenMM). X-ray structures do not show
     hydrogens, but hydrogen bonds are the physics docking scores. Histidine
     protonation (which ring nitrogen carries the H) is chosen by OpenMM's
     heuristic; no manual flipping.
  3. Convert to PDBQT with Meeko (the arm64-compatible tool; ADFR/MGLTools are not
     used). PDBQT is Vina's input format: atoms plus a partial charge and an
     "atom type" that says what interaction the atom can make (H-bond donor,
     acceptor, aromatic carbon...).
  4. The search box is a cube centred on the crystal ligand's atoms.

4KD1 has 298 modelled residues and no gaps, so no loop rebuilding is done.

Run:
    python scripts/08_prepare_receptor.py
"""

import shutil
import subprocess
import urllib.request

import numpy as np
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem

from cdk2moo import config


def main():
    config.require_campaign()
    pdb_path = config.STRUCTURES_DIR / f"{config.PDB_ID}.pdb"
    sdf_path = config.STRUCTURES_DIR / f"{config.PDB_ID.lower()}_{config.LIGAND_CODE}_crystal.sdf"
    for path, url in [(pdb_path, f"https://files.rcsb.org/download/{config.PDB_ID}.pdb"), (sdf_path, config.LIGAND_SDF_URL)]:
        if path.exists():
            continue
        frozen = config.INPUT_STRUCTURES_DIR / path.name
        if frozen.exists() and frozen != path:
            shutil.copy2(frozen, path)                       # the campaign's frozen structural input
        elif config.CAMPAIGN == "historical":
            urllib.request.urlretrieve(url, path)
        else:
            raise SystemExit(f"{path.name} not found in the campaign inputs ({config.INPUT_STRUCTURES_DIR}); a campaign never "
                             "downloads structures in place of its frozen inputs.")
    print(f"Structure {pdb_path.name}, crystal ligand {sdf_path.name}")

    # ---- 1-2. protein only, with hydrogens -----------------------------------
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.removeChains([i for i, c in enumerate(fixer.topology.chains()) if c.id != "A"])
    fixer.removeHeterogens(keepWater=False)          # drops 1QK, EDO and waters
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    print(f"  missing residues: {len(fixer.missingResidues)}, missing atoms: "
          f"{sum(len(v) for v in fixer.missingAtoms.values())}, non-standard residues: {len(fixer.nonstandardResidues)}")
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(config.DOCK_PH)
    receptor_h = config.STRUCTURES_DIR / "receptor_H.pdb"
    with open(receptor_h, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    n_atoms = fixer.topology.getNumAtoms()
    print(f"  protein with hydrogens: {n_atoms} atoms")

    # ---- 4. box centred on the crystal ligand ---------------------------------
    crystal = Chem.MolFromMolFile(str(sdf_path), removeHs=False)
    xyz = crystal.GetConformer().GetPositions()
    centre = xyz.mean(axis=0)
    extent = xyz.max(axis=0) - xyz.min(axis=0)
    print(f"  crystal ligand: {crystal.GetNumHeavyAtoms()} heavy atoms, centre "
          f"({centre[0]:.1f}, {centre[1]:.1f}, {centre[2]:.1f}), extent "
          f"{extent[0]:.1f} x {extent[1]:.1f} x {extent[2]:.1f} A; box {config.DOCK_BOX_SIZE} A cube")

    # ---- 3. PDBQT with Meeko --------------------------------------------------
    base = config.STRUCTURES_DIR / "receptor"
    cmd = ["mk_prepare_receptor.py", "--read_pdb", str(receptor_h), "-o", str(base), "-p", "-v",
           "--box_center", *[f"{c:.3f}" for c in centre], "--box_size", *[str(config.DOCK_BOX_SIZE)] * 3,
           "--default_altloc", "A"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout[-600:] if result.stdout else "", result.stderr[-800:] if result.returncode else "")
    if result.returncode != 0:
        raise SystemExit("Meeko receptor preparation failed; see output above.")
    with open(config.STRUCTURES_DIR / "receptor_box.txt", "w") as fh:
        fh.write(f"center {centre[0]:.3f} {centre[1]:.3f} {centre[2]:.3f}\nsize {config.DOCK_BOX_SIZE} "
                 f"{config.DOCK_BOX_SIZE} {config.DOCK_BOX_SIZE}\n")
    print("\nWrote receptor_H.pdb, receptor.pdbqt, receptor_box.txt in data/structures/")
    print("Next: scripts/08b_redock_native_ligand.py")


if __name__ == "__main__":
    main()
