"""
Choosing the molecules to dock, organised by EXPERIMENT GROUP so that every comparison states what the surrogate did and did not see.

Groups (each is docked and analysed separately; groups are never pooled):

  P  primary experiment. The surrogate was trained on ALL curated molecules. Its "known actives" comparator is the top decile of the
     TRAINING set, so those actives were IN the surrogate's training data: an in-sample comparator, labelled `in_sample_top_decile`,
     which is NOT held out. It can show whether docking tracks activity at all, not whether optimized molecules beat unseen actives.
  M  molecule-level control. The surrogate was trained WITHOUT the 204 top-decile actives (siblings stay). Comparator: those held-out actives.
  S  scaffold-level control. The surrogate was trained without 9 whole scaffold families (260 molecules; 93 of them top-decile actives).
     Comparator: the 93 held-out actives.

H3 ("Vina scores of optimized molecules do not exceed those of held-out known actives") needs a comparator the surrogate never saw; only groups M
and S provide one. Group P is reported as an in-sample comparison.

Sets inside a group: <G>:<comparator actives>, <G>:random_*, <G>:decoys, <G>:<arm> for the four GA arms, and for P the pooled Pareto front.
Reference sets (actives, random molecules) depend only on the curated data, the holdout assignment and the seed, so they are IDENTICAL in
the legacy and corrected scopes and are docked once. Generated sets and decoys come from that scope's GA populations.
Every draw uses a generator keyed by (seed, set), so a set does not depend on the order in which sets are drawn.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

from cdk2moo.features import ecfp4, max_tanimoto

ARMS = ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]
GROUP_ROLE = {"P": "in_sample_top_decile", "M": "heldout_active", "S": "heldout_active"}
SET_KEY = {"P:training_top_decile": 1, "P:random_training": 2, "M:heldout_actives": 3, "M:random_remaining": 4, "S:random_remaining": 6}
PROPERTY_NAMES = ["heavy atoms", "MW", "cLogP", "HBD", "HBA", "rotatable bonds"]


def keyed_rng(seed, *key):
    """A generator that depends only on the seed and the (string or integer) key, never on how many draws were made before."""
    numbers = [seed] + [sum(ord(c) * 31 ** i for i, c in enumerate(str(k))) % (2 ** 31) for k in key]
    return np.random.default_rng(numbers)


def properties(smiles_list):
    """The properties that dominate docking scores, for decoy matching: heavy atoms, MW, cLogP, H-bond donors and acceptors, rotatable bonds."""
    out = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        out.append([m.GetNumHeavyAtoms(), Descriptors.MolWt(m), Descriptors.MolLogP(m), Lipinski.NumHDonors(m),
                    Lipinski.NumHAcceptors(m), rdMolDescriptors.CalcNumRotatableBonds(m)])
    return np.array(out, dtype=float)


def match_decoys(active_smiles, pool_smiles, all_active_fps, rng, max_similarity=0.5):
    """
    One decoy per active: the nearest unused pool molecule in z-scored property space, restricted to pool molecules whose ECFP4 Tanimoto to EVERY
    active in `all_active_fps` is below `max_similarity`. Actives are served in a random order so none gets first pick systematically.
    Returns (chosen pool indices, standardized mean differences per property). Raises if the pool is too small.
    """
    pool_fps = ecfp4(list(pool_smiles))
    far = np.where(max_tanimoto(pool_fps, all_active_fps) < max_similarity)[0]
    if len(far) < len(active_smiles):
        raise ValueError(f"only {len(far)} decoy candidates for {len(active_smiles)} actives")
    p_pool = properties([pool_smiles[i] for i in far])
    p_act = properties(list(active_smiles))
    scale = p_pool.std(axis=0)
    scale[scale == 0] = 1.0
    used = np.zeros(len(far), dtype=bool)
    chosen = []
    for k in rng.permutation(len(active_smiles)):
        dist = np.sqrt((((p_pool - p_act[k]) / scale) ** 2).sum(axis=1))
        dist[used] = np.inf
        j = int(np.argmin(dist))
        used[j] = True
        chosen.append(int(far[j]))
    smd = (p_pool[[int(np.where(far == c)[0][0]) for c in chosen]].mean(axis=0) - p_act.mean(axis=0)) / scale
    return chosen, dict(zip(PROPERTY_NAMES, smd))
