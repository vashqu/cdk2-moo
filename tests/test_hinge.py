"""Tests for the distance-only hinge measure. Run: python -m unittest discover tests"""

import unittest

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from cdk2moo import hinge


def ligand_with_oxygen_at(x, y, z):
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    conf = mol.GetConformer()
    oxygen = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "O"][0]
    shift = np.array([x, y, z]) - np.array(conf.GetAtomPosition(oxygen))
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, np.array(conf.GetAtomPosition(i)) + shift)
    return mol


class HingeTests(unittest.TestCase):
    HINGE = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0]])

    def test_distance_is_to_the_nearest_hinge_atom(self):
        mol = ligand_with_oxygen_at(2.0, 0.0, 0.0)
        self.assertAlmostEqual(hinge.nearest_polar_distance(mol, self.HINGE), 2.0, places=3)

    def test_a_ligand_without_n_or_o_is_infinitely_far(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCC"))
        AllChem.EmbedMolecule(mol, randomSeed=1)
        self.assertEqual(hinge.nearest_polar_distance(mol, self.HINGE), float("inf"))

    def test_hinge_atoms_requires_exactly_three(self):
        mol = Chem.MolFromSmiles("CC")
        AllChem.EmbedMolecule(mol, randomSeed=1)
        with self.assertRaises(Exception):
            hinge.hinge_atoms(mol)


if __name__ == "__main__":
    unittest.main()
