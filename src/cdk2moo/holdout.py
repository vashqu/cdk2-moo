"""
Held-out positive controls: can the pipeline find potency it was never shown?

The main GA experiment asks what an optimizer does to a surrogate trained on all
the data. A positive control turns the question round: take away some known
potent molecules, train the surrogate on the rest, run the same GA, and ask
whether the search finds chemistry near the molecules that were taken away.
Measured activities of the held-out molecules are the ground truth, so success
is judged against experiment and not against the surrogate's own prediction.

Two holdouts, answering different questions. Both use rules fixed here, from the
data alone, before any GA output exists; nothing is chosen after seeing results.

  molecule level   Remove the top-decile actives themselves. Close analogues and
                   scaffold-mates stay in training, so the surrogate has seen
                   their siblings. A GA success here shows interpolation: it can
                   get to a potent molecule by staying near potent siblings.

  scaffold level   Remove every molecule belonging to the scaffold families
                   (see features.py for scaffolds) that hold several top-decile
                   actives, potent or not. The surrogate never sees those ring
                   systems at all. A success here would be evidence of genuine
                   chemical generalization; a failure is informative, not a bug.

Note what a held-out surrogate cannot do: a random forest predicts averages of
training labels, so it can never predict above the highest label it was trained
on. Success must therefore be measured by similarity to held-out actives, never
by the predicted activity alone.
"""

import numpy as np


def top_decile(pactivity, quantile):
    """
    Boolean mask of the most potent molecules and the cutoff used.

    A molecule is "top decile" if its measured pActivity is >= the `quantile`
    quantile of the dataset (>= so that ties at the cutoff are all included).
    """
    cutoff = float(np.quantile(pactivity, quantile))
    return np.asarray(pactivity) >= cutoff, cutoff


def molecule_holdout(pactivity, quantile):
    """
    Hold out exactly the top-decile molecules.

    Returns (train_idx, held_idx, cutoff). Deterministic: no sampling.
    """
    is_top, cutoff = top_decile(pactivity, quantile)
    return np.where(~is_top)[0], np.where(is_top)[0], cutoff


def scaffold_holdout(pactivity, scaffolds, quantile, min_top_per_family):
    """
    Hold out whole scaffold families that contain several top-decile actives.

    A family is held out if at least `min_top_per_family` of its members are top
    decile; then EVERY member is held out, potent or not. Families with fewer top
    actives stay in training, including their potent members. The threshold picks
    out potent SERIES rather than one-off potent molecules.

    Returns (train_idx, held_idx, cutoff, families) where `families` maps each held
    scaffold to (n_members, n_top_decile). Deterministic: no sampling.
    """
    scaffolds = np.asarray(scaffolds)
    is_top, cutoff = top_decile(pactivity, quantile)

    families = {}
    for scaffold in np.unique(scaffolds):
        members = scaffolds == scaffold
        n_top = int((members & is_top).sum())
        if n_top >= min_top_per_family:
            families[scaffold] = (int(members.sum()), n_top)

    held = np.isin(scaffolds, list(families))
    return np.where(~held)[0], np.where(held)[0], cutoff, families
