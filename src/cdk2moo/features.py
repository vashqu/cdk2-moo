"""
Molecular features: fingerprints, scaffolds, and fingerprint similarity.

Fingerprint (ECFP4)
-------------------
A machine-learning model cannot read a molecule; it needs a fixed-length list
of numbers. ECFP4 ("extended-connectivity fingerprint, diameter 4", also called
a Morgan fingerprint of radius 2) builds one like this:

  1. Give every atom a label from its own properties (element, charge, number
     of bonds...).
  2. Repeat twice: update each atom's label by combining it with its
     neighbours' labels. After 2 rounds each label describes the atom plus
     everything within 2 bonds of it - a small substructure.
  3. Hash every such substructure to a number between 0 and 2047 and switch that
     bit on.

The result is a 2048-long 0/1 vector. Two molecules sharing many substructures
share many on-bits. Different substructures can hash to the same bit
("collisions"); with 2048 bits and drug-sized molecules this is mild but real.
Chirality is not encoded, so two stereoisomers get identical fingerprints.

Bemis-Murcko scaffold
---------------------
The scaffold is a molecule's skeleton: take the molecule, delete every side
chain, and keep only the ring systems plus the chains of atoms that link rings
together. (Atoms double-bonded to a ring or linker atom, like the C=O of a
ring ketone, are kept.) Example: the drug-like molecule
"CCOc1ccc(cc1)C(=O)Nc1ccccc1" becomes "O=C(c1ccccc1)Nc1ccccc1" - two benzene
rings joined by an amide linker; the ethoxy tail disappears.

Why it matters here: medicinal chemists make series of analogues that share one
scaffold and differ in their side chains, and papers publish them together.
Two molecules with the same scaffold are therefore usually close relatives. A
random split puts relatives on both sides, so a model can score well by
recognising the series. Grouping by scaffold keeps whole series on one side.
A molecule with no ring at all has an empty scaffold ("").
"""

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


def ecfp4(smiles_list, radius=2, n_bits=2048):
    """
    Return an (n_molecules, n_bits) uint8 array of ECFP4 fingerprints.

    Input must be the standardized SMILES from stage 2, so that a molecule has
    one fingerprint regardless of how it was originally drawn.
    """
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius,
                                                          fpSize=n_bits)
    fps = np.zeros((len(smiles_list), n_bits), dtype=np.uint8)
    for i, smiles in enumerate(smiles_list):
        fps[i] = generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
    return fps


def murcko_scaffold(smiles):
    """
    Return the Bemis-Murcko scaffold of one molecule as a SMILES string.

    Stereochemistry is removed so that two stereoisomers share a scaffold.
    Returns "" for molecules with no ring.
    """
    return MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(smiles),
                                               includeChirality=False)


def max_tanimoto(query_fps, reference_fps):
    """
    For each query molecule, the highest Tanimoto similarity to any reference.

    Tanimoto = (bits both molecules have on) / (bits either has on). It is 1 for
    identical fingerprints and 0 for none shared. The "max to reference" is a
    nearest-neighbour distance: how close is the closest molecule the model has
    seen? This is the quantity we later use to ask whether a generated molecule
    is inside or outside the training distribution.
    """
    q = query_fps.astype(np.float32)
    r = reference_fps.astype(np.float32)
    shared = q @ r.T
    either = q.sum(axis=1)[:, None] + r.sum(axis=1)[None, :] - shared
    return (shared / either).max(axis=1)
