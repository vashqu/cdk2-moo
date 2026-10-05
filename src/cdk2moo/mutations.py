"""
The genetic algorithm's two ways of inventing new molecules: mutation and
crossover. Both edit the molecular GRAPH (atoms are nodes, bonds are edges);
they know no chemistry beyond "RDKit refuses molecules with impossible
valences". A proposed molecule that RDKit cannot sanitize is simply discarded,
so every edit here is a gamble, and most gambles fail. That is expected.

Chemically these edits are naive: swapping a carbon for a nitrogen in a ring,
or bolting an oxygen onto a random atom, often gives something nobody could
make or that would be unstable. The GA does not care, and neither does the
surrogate. That is part of what the project examines.

Atoms are edited in Kekulé form (aromatic rings drawn with explicit alternating
single/double bonds) so that RDKit can re-detect aromaticity afterwards.
"""

import numpy as np
from rdkit import Chem, RDLogger

from cdk2moo import config

# Failed edits are expected and handled by returning None; RDKit's valence
# messages for them are noise.
RDLogger.DisableLog("rdApp.*")


def count_radicals(mol):
    """Number of atoms carrying unpaired (radical) electrons. 0 means closed-shell as RDKit perceives it."""
    return sum(1 for atom in mol.GetAtoms() if atom.GetNumRadicalElectrons() > 0)


def rejection_reason(mol, size_range, closed_shell=False):
    """
    Why a proposed molecule is not accepted, or None if it is.

    Checks in order: it sanitizes (legal valences), it is one connected piece, its heavy-atom count is in
    the size window, and, only when closed_shell=True, no atom carries radical electrons.
    closed_shell=False is the LEGACY acceptance used by every historical run (radicals were admitted; saved
    activity-only populations contain some). closed_shell=True is the corrected variant's policy. Neither
    establishes that a molecule is stable, synthesizable or sensible; they only exclude structures RDKit
    sees as broken, fragmented, out of range or open-shell.
    """
    if mol is None:
        return "no_molecule"
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return "sanitize"
    if len(Chem.GetMolFrags(mol)) != 1:
        return "disconnected"
    low, high = size_range
    if not low <= mol.GetNumHeavyAtoms() <= high:
        return "size"
    if closed_shell and count_radicals(mol) > 0:
        return "radical"
    return None


def is_valid(mol, size_range, closed_shell=False):
    """True if rejection_reason finds nothing wrong (legacy acceptance unless closed_shell=True)."""
    return rejection_reason(mol, size_range, closed_shell) is None


def _record(rejects, reason):
    if rejects is not None:
        rejects[reason] = rejects.get(reason, 0) + 1


def _largest_fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    return max(frags, key=lambda m: m.GetNumAtoms())


def _append_atom(mol, rng):
    """Attach a new atom by a single bond to a random atom."""
    i = int(rng.integers(mol.GetNumAtoms()))
    symbol = config.GA_ALLOWED_ELEMENTS[int(rng.integers(len(config.GA_ALLOWED_ELEMENTS)))]
    new = mol.AddAtom(Chem.Atom(symbol))
    mol.AddBond(i, new, Chem.BondType.SINGLE)
    return mol


def _delete_atom(mol, rng):
    """Remove a random atom; if that splits the molecule keep the biggest piece."""
    mol.RemoveAtom(int(rng.integers(mol.GetNumAtoms())))
    return _largest_fragment(mol)


def _change_atom(mol, rng):
    """Swap one atom's element for C, N, O or S (e.g. a ring CH becomes N)."""
    i = int(rng.integers(mol.GetNumAtoms()))
    atom = mol.GetAtomWithIdx(i)
    atom.SetAtomicNum([6, 7, 8, 16][int(rng.integers(4))])
    atom.SetNoImplicit(False)
    return mol


def _insert_atom(mol, rng, aromatic_bonds):
    """Break a bond and put a new C, N, O or S atom in the gap."""
    options = [b for b in mol.GetBonds() if b.GetIdx() not in aromatic_bonds]
    if not options:
        return None
    b = options[int(rng.integers(len(options)))]
    i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    mol.RemoveBond(i, j)
    new = mol.AddAtom(Chem.Atom([6, 7, 8, 16][int(rng.integers(4))]))
    mol.AddBond(i, new, Chem.BondType.SINGLE)
    mol.AddBond(new, j, Chem.BondType.SINGLE)
    return mol


def _flip_bond_order(mol, rng, aromatic_bonds):
    """Turn a non-aromatic single bond into a double bond, or the reverse."""
    options = [b for b in mol.GetBonds() if b.GetIdx() not in aromatic_bonds
               and b.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE)]
    if not options:
        return None
    b = options[int(rng.integers(len(options)))]
    if b.GetBondType() == Chem.BondType.SINGLE:
        b.SetBondType(Chem.BondType.DOUBLE)
    else:
        b.SetBondType(Chem.BondType.SINGLE)
    return mol


def _close_ring(mol, rng):
    """Bond two atoms that are 4 or 5 bonds apart, making a 5- or 6-membered ring."""
    dist = Chem.GetDistanceMatrix(mol)
    pairs = np.argwhere(np.triu((dist == 4) | (dist == 5)))
    if not len(pairs):
        return None
    i, j = pairs[int(rng.integers(len(pairs)))]
    mol.AddBond(int(i), int(j), Chem.BondType.SINGLE)
    return mol


def _open_ring(mol, rng, aromatic_bonds):
    """Break a non-aromatic bond that sits in a ring."""
    options = [b for b in mol.GetBonds()
               if b.IsInRing() and b.GetIdx() not in aromatic_bonds]
    if not options:
        return None
    b = options[int(rng.integers(len(options)))]
    mol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    return mol


def mutate(mol, rng, size_range, max_tries=20, closed_shell=False, rejects=None):
    """
    Return a valid mutated copy of `mol`, or None if every try failed.

    Picks one of seven edits at random and retries with another draw if the
    result is rejected (see rejection_reason). `rejects`, if given, is a dict that accumulates
    rejection reasons by count (no per-molecule logging).
    """
    for _ in range(max_tries):
        aromatic_bonds = {b.GetIdx() for b in mol.GetBonds() if b.GetIsAromatic()}
        new = Chem.RWMol(mol)
        Chem.Kekulize(new, clearAromaticFlags=True)
        edit = int(rng.integers(7))
        if edit == 0:
            new = _append_atom(new, rng)
        elif edit == 1:
            new = _delete_atom(new, rng)
        elif edit == 2:
            new = _change_atom(new, rng)
        elif edit == 3:
            new = _insert_atom(new, rng, aromatic_bonds)
        elif edit == 4:
            new = _flip_bond_order(new, rng, aromatic_bonds)
        elif edit == 5:
            new = _close_ring(new, rng)
        else:
            new = _open_ring(new, rng, aromatic_bonds)
        reason = "edit_not_applicable" if new is None else rejection_reason(new, size_range, closed_shell)
        if reason is None:
            return new
        _record(rejects, reason)
    return None


def _cut_halves(mol, rng):
    """
    Cut a random acyclic single bond between two non-terminal atoms and return
    the two halves, each carrying a dummy atom (atomic number 0) where the cut
    was. Returns None if the molecule has no such bond.
    """
    options = [b.GetIdx() for b in mol.GetBonds()
               if not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE
               and b.GetBeginAtom().GetDegree() > 1 and b.GetEndAtom().GetDegree() > 1]
    if not options:
        return None
    bond = options[int(rng.integers(len(options)))]
    cut = Chem.FragmentOnBonds(mol, [bond], addDummies=True)
    return Chem.GetMolFrags(cut, asMols=True, sanitizeFrags=False)


def crossover(mol_a, mol_b, rng, size_range, max_tries=10, closed_shell=False, rejects=None):
    """
    Return a valid child built from one half of each parent, or None.

    Each parent is cut at a random acyclic single bond; one half from each is
    joined across their cut points. The child inherits a chunk of each parent,
    which is how the GA recombines building blocks that scored well.
    """
    for _ in range(max_tries):
        halves_a = _cut_halves(mol_a, rng)
        halves_b = _cut_halves(mol_b, rng)
        if halves_a is None or halves_b is None:
            return None
        part_a = halves_a[int(rng.integers(2))]
        part_b = halves_b[int(rng.integers(2))]
        combo = Chem.RWMol(Chem.CombineMols(part_a, part_b))
        dummies = [a.GetIdx() for a in combo.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummies) != 2:
            continue
        anchor_a = combo.GetAtomWithIdx(dummies[0]).GetNeighbors()[0].GetIdx()
        anchor_b = combo.GetAtomWithIdx(dummies[1]).GetNeighbors()[0].GetIdx()
        combo.AddBond(anchor_a, anchor_b, Chem.BondType.SINGLE)
        for d in sorted(dummies, reverse=True):   # remove high index first
            combo.RemoveAtom(d)
        reason = rejection_reason(combo, size_range, closed_shell)
        if reason is None:
            return combo
        _record(rejects, reason)
    return None
