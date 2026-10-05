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

Grid maps are only reused after validation against a manifest of file hashes (see
write_grid_maps / validate_grid_maps). The directory data/structures/receptor_maps, made by the
first docking run, has no manifest: it is legacy and is not adopted for new runs.
"""

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem, __version__ as rdkit_version
from rdkit.Chem import AllChem, rdMolDescriptors
import vina as vina_package
from vina import Vina

from cdk2moo import config
from cdk2moo.dock_jobs import ResumeRefused, sha256_of

LEGACY_MAP_PREFIX = str(config.STRUCTURES_DIR / "receptor_maps" / "maps")
MAPS_MANIFEST = "maps_manifest.json"


def read_box(path=None):
    """The docking box written by scripts/08_prepare_receptor.py in the selected campaign scope: (centre, size) as float lists."""
    with open(path or config.STRUCTURES_DIR / "receptor_box.txt") as fh:
        box = {k: v for k, *v in (line.split() for line in fh)}
    return [float(x) for x in box["center"]], [float(x) for x in box["size"]]


def write_grid_maps(directory, receptor=None, box=None):
    """
    Compute Vina's affinity grids for the receptor and save them with a manifest.

    The manifest records the receptor hash, box and the size and SHA-256 of every map file, so a
    later run can tell a complete, matching set of maps from a partial or stale one.
    """
    receptor = Path(receptor or config.STRUCTURES_DIR / "receptor.pdbqt")
    centre, size = box or read_box()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(str(receptor))
    v.compute_vina_maps(center=centre, box_size=size, force_even_voxels=True)
    v.write_maps(str(directory / "maps"))
    files = {p.name: {"bytes": p.stat().st_size, "sha256": sha256_of(p)}
             for p in sorted(directory.glob("maps.*"))}
    manifest = {"receptor_sha256": sha256_of(receptor), "box_center": centre, "box_size": size,
                "force_even_voxels": True, "files": files}
    (directory / MAPS_MANIFEST).write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def validate_grid_maps(directory, receptor=None, box=None):
    """
    Accept a grid-map directory only if its manifest exists and every listed file is present with
    the recorded size and hash, and the manifest's receptor and box match the current ones.
    Raises ResumeRefused with the reason otherwise.
    """
    receptor = Path(receptor or config.STRUCTURES_DIR / "receptor.pdbqt")
    directory = Path(directory)
    manifest_path = directory / MAPS_MANIFEST
    if not manifest_path.exists():
        raise ResumeRefused(f"{directory} has no {MAPS_MANIFEST}: legacy or unverified grid maps, not reused")
    manifest = json.loads(manifest_path.read_text())
    if not manifest.get("files"):
        raise ResumeRefused(f"{directory}: manifest lists no map files")
    for name, info in manifest["files"].items():
        path = directory / name
        if not path.exists():
            raise ResumeRefused(f"{directory}: grid map {name} is missing (incomplete directory)")
        if path.stat().st_size != info["bytes"] or sha256_of(path) != info["sha256"]:
            raise ResumeRefused(f"{directory}: grid map {name} does not match its recorded hash")
    centre, size = box or read_box()
    if manifest["receptor_sha256"] != sha256_of(receptor):
        raise ResumeRefused(f"{directory}: maps were made from a different receptor file")
    if manifest["box_center"] != centre or manifest["box_size"] != size:
        raise ResumeRefused(f"{directory}: maps were made for a different search box")
    return manifest


def provenance(prep_seed, receptor=None, box=None):
    """Everything that could change a docking result, for the run manifest."""
    receptor = Path(receptor or config.STRUCTURES_DIR / "receptor.pdbqt")
    centre, size = box or read_box()
    import meeko
    return {"receptor_sha256": sha256_of(receptor), "box_center": centre, "box_size": size,
            "scoring_function": "vina", "exhaustiveness": config.DOCK_EXHAUSTIVENESS,
            "n_poses": config.DOCK_N_POSES, "vina_cpu": 0, "force_even_voxels": True,
            "ligand_prep_seed": int(prep_seed), "ligand_prep": "AddHs + ETKDGv3 + MMFF(500) + Meeko",
            "vina_version": getattr(vina_package, "__version__", "unknown"),
            "meeko_version": getattr(meeko, "__version__", "unknown"), "rdkit_version": rdkit_version,
            "platform": platform.platform(), "python": platform.python_version()}


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


def dock(pdbqt, seed, maps_prefix):
    """
    Dock a prepared ligand using the maps at `maps_prefix`. Returns (best score in kcal/mol, RDKit
    molecule holding the best pose, or None if the pose could not be converted). Score is Vina's
    estimate of binding energy: more negative is better.
    """
    v = Vina(sf_name="vina", seed=seed, verbosity=0)
    v.load_maps(maps_prefix)
    v.set_ligand_from_string(pdbqt)
    v.dock(exhaustiveness=config.DOCK_EXHAUSTIVENESS, n_poses=config.DOCK_N_POSES)
    best = float(v.energies(n_poses=1)[0][0])
    try:
        poses = PDBQTMolecule(v.poses(n_poses=1), is_dlg=False, skip_typing=True)
        pose = RDKitMolCreate.from_pdbqt_mol(poses)[0]
    except Exception:
        pose = None
    return best, pose
