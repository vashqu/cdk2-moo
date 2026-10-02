"""
Docking one molecule with Vina: build a 3-D ligand, dock it, return the best pose.

A SMILES string has no 3-D shape, so each ligand is first given one: RDKit embeds a
plausible 3-D conformer (ETKDG, a distance-geometry method that respects known
bond-length and torsion preferences) and relaxes it with a force field (MMFF).
Vina then searches over position, orientation and the ligand's rotatable bonds, so
the starting conformer matters little. Unspecified stereocentres are assigned
whatever the embedding picks (seeded); stereo and protonation are NOT enumerated,
which is a stated limitation.

Each molecule gets its own freshly seeded Vina object and pre-computed grid maps, so
a molecule's score does not depend on which molecules were docked before it.
"""

import numpy as np
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from vina import Vina

from cdk2moo import config

MAP_PREFIX = str(config.STRUCTURES_DIR / "receptor_maps" / "maps")


def write_grid_maps():
    """Compute Vina's affinity grids once for the receptor and save them for reuse."""
    box = {k: v for k, *v in (line.split() for line in open(config.STRUCTURES_DIR / "receptor_box.txt"))}
    (config.STRUCTURES_DIR / "receptor_maps").mkdir(exist_ok=True)
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(str(config.STRUCTURES_DIR / "receptor.pdbqt"))
    v.compute_vina_maps(center=[float(x) for x in box["center"]], box_size=[float(x) for x in box["size"]],
                        force_even_voxels=True)
    v.write_maps(MAP_PREFIX)


def prepare_ligand(smiles, seed):
    """
    SMILES -> (PDBQT string, RDKit molecule with hydrogens, rotatable bond count),
    or None if embedding or Meeko preparation fails.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    try:
        setup = MoleculePreparation().prepare(mol)[0]
        pdbqt, ok, _ = PDBQTWriterLegacy.write_string(setup)
    except Exception:
        return None
    if not ok:
        return None
    return pdbqt, mol, rdMolDescriptors.CalcNumRotatableBonds(mol)


def dock(pdbqt, seed):
    """
    Dock a prepared ligand. Returns (best score in kcal/mol, RDKit molecule holding
    the best pose). Score is Vina's estimate of binding energy: more negative is better.
    """
    v = Vina(sf_name="vina", seed=seed, verbosity=0)
    v.load_maps(MAP_PREFIX)
    v.set_ligand_from_string(pdbqt)
    v.dock(exhaustiveness=config.DOCK_EXHAUSTIVENESS, n_poses=config.DOCK_N_POSES)
    best = float(v.energies(n_poses=1)[0][0])
    poses = PDBQTMolecule(v.poses(n_poses=1), is_dlg=False, skip_typing=True)
    return best, RDKitMolCreate.from_pdbqt_mol(poses)[0]
