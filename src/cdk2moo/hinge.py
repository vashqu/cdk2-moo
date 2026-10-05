"""
Distance-only hinge proximity (see scripts/08g_hinge_contacts.py for the biology and the limits of this measure).

`hinge_atoms` picks the three backbone atoms dinaciclib hydrogen-bonds to (Glu81 O, Leu83 N, Leu83 O); `nearest_polar_distance` returns the smallest distance
from any ligand nitrogen or oxygen atom to them. Nothing here looks at donor/acceptor roles, angles or hydrogens.
"""

import numpy as np

HINGE = {("GLU", 81): ["O"], ("LEU", 83): ["N", "O"]}


def hinge_atoms(protein):
    """Coordinates (3 x 3 array) of the hinge backbone atoms of an RDKit protein molecule with PDB residue info; raises if not exactly three are found."""
    xyz = protein.GetConformer().GetPositions()
    found = []
    for atom in protein.GetAtoms():
        info = atom.GetPDBResidueInfo()
        residue = (info.GetResidueName().strip(), info.GetResidueNumber())
        if residue in HINGE and info.GetName().strip() in HINGE[residue]:
            found.append(xyz[atom.GetIdx()])
    if len(found) != 3:
        raise ValueError(f"expected 3 hinge backbone atoms, found {len(found)}")
    return np.array(found)


def nearest_polar_distance(mol, hinge_xyz):
    """Smallest distance (A) from a ligand N or O atom to the hinge atoms; infinity if the ligand has no N or O."""
    polar = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ("N", "O")]
    if not polar:
        return float("inf")
    pos = mol.GetConformer().GetPositions()
    return float(np.linalg.norm(pos[polar][:, None, :] - hinge_xyz[None, :, :], axis=2).min())
