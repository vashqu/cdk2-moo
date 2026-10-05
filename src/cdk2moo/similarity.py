"""
Size-conditioned similarity: is a molecule as close to the training set as a
training molecule of the same size typically is?

The problem. The raw measure used everywhere else is the highest Tanimoto
similarity between a molecule's ECFP4 fingerprint (Morgan radius 2, 2048 bits,
binary, no chirality) and the fingerprints of the 2,016 training molecules.
That number depends on molecule size: a small molecule switches on few bits,
so it has little to share with anything, and scores low even if it is an
obvious fragment of a known compound.

The remedy, with no invented formula. Ask a relative question. For every
TRAINING molecule, compute its nearest-neighbour Tanimoto to the OTHER 2,015
training molecules (leave-one-out; comparing a molecule with itself would give
1.0 and mean nothing). That gives an empirical distribution of "how close is a
typical training molecule of size h to the rest of the training set". A
generated molecule's raw max-Tanimoto is then reported as a PERCENTILE within
the distribution for training molecules of similar size:

    percentile = 100 x (fraction of size-matched training molecules whose
                        leave-one-out nearest-neighbour Tanimoto <= this value)

50 = as close to the training set as a typical training molecule of that size;
10 = closer than only 10% of them, i.e. unusually novel for its size;
90 = unusually well covered.

Caveat: the reference is the training set's own local density. Training
molecules come in analogue series, so their leave-one-out similarities are
high, and a "typical" percentile is therefore a demanding bar.
"""

import numpy as np

from cdk2moo import config


def loo_max_tanimoto(fps):
    """
    For each molecule, its highest Tanimoto to any OTHER molecule in `fps`.

    Two different molecules can share an identical fingerprint (for example
    stereoisomers, which ECFP4 cannot tell apart); they score 1.0 against each
    other, which is genuine and is counted.
    """
    f = fps.astype(np.float32)
    shared = f @ f.T
    either = f.sum(axis=1)[:, None] + f.sum(axis=1)[None, :] - shared
    sim = shared / either
    np.fill_diagonal(sim, -1.0)   # never compare a molecule with itself
    return sim.max(axis=1)


def size_conditioned_percentile(sim, sizes, ref_sizes, ref_sim):
    """
    Percentile of each raw similarity within its size-matched reference.

    sim, sizes         raw max-Tanimoto and heavy-atom count of the molecules to score
    ref_sizes, ref_sim heavy-atom count and leave-one-out similarity of the training set

    Stratum rule (config): heavy atoms within +/- SIZE_STRATUM_HALFWIDTH, widened
    one atom at a time until at least SIZE_STRATUM_MIN_N reference molecules fall
    inside. Returns (percentile, half-width used, reference count).
    """
    sim, sizes = np.asarray(sim), np.asarray(sizes)
    pct = np.zeros(len(sim))
    used = np.zeros(len(sim), dtype=int)
    n_ref = np.zeros(len(sim), dtype=int)

    for h in np.unique(sizes):
        width = config.SIZE_STRATUM_HALFWIDTH
        while True:
            in_stratum = np.abs(ref_sizes - h) <= width
            if in_stratum.sum() >= config.SIZE_STRATUM_MIN_N or width > 100:
                break
            width += 1
        reference = np.sort(ref_sim[in_stratum])
        rows = sizes == h
        pct[rows] = 100.0 * np.searchsorted(reference, sim[rows], side="right") / len(reference)
        used[rows] = width
        n_ref[rows] = len(reference)
    return pct, used, n_ref
