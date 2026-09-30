"""
Train/test splits, and the statistics that say how different they really are.

Three splits over the same molecules, the same test fraction and the same seed:

  random split    every molecule is assigned independently. Test molecules
                  typically have close relatives (same scaffold, same paper) in
                  the training set. It measures interpolation.

  scaffold split  all molecules sharing a Bemis-Murcko scaffold (see
                  features.py) go to the same side. Test molecules belong to
                  ring systems the model has never seen. It is a much harder
                  test and a better proxy for what a generative optimizer
                  asks of a surrogate: predictions on new chemistry.

  paper split     same idea, but grouped by the paper a molecule first appeared
                  in. Analogue series are published together, and a new scaffold
                  is often one ring swapped in a known series, so the scaffold
                  split still leaves near-relatives across the divide. The
                  paper split is the strictest of the three.

The gaps between them are the point, so all are kept and reported together.
"""

import numpy as np

from cdk2moo.features import max_tanimoto


def random_split(n, test_frac, seed):
    """Return (train_idx, test_idx), sorted arrays of row indices."""
    order = np.random.default_rng(seed).permutation(n)
    n_test = round(n * test_frac)
    return np.sort(order[n_test:]), np.sort(order[:n_test])


def group_split(group_labels, test_frac, seed):
    """
    Return (train_idx, test_idx), keeping every group on one side.

    `group_labels` has one label per molecule (its scaffold, or its paper).

    This is the "balanced" scaffold split: a plain split that puts the largest
    groups in train and the singletons in test makes the test set artificially
    all-singletons, and is deterministic so the seed would do nothing. Instead:

      1. Groups bigger than half the test set go to train first (one giant
         series must not decide the test set's size).
      2. The remaining groups are shuffled with the seed and added to train
         until it holds (1 - test_frac) of the molecules; the rest is test.

    Group sizes are lumpy, so the test fraction is approximate, and it can only
    be >= the requested fraction.
    """
    groups = {}
    for i, label in enumerate(group_labels):
        groups.setdefault(label, []).append(i)

    n = len(group_labels)
    n_train_target = n - round(n * test_frac)
    half_test = round(n * test_frac) / 2

    big = [g for g in groups.values() if len(g) > half_test]
    small = [g for g in groups.values() if len(g) <= half_test]

    train, test = [], []
    for g in big:
        train += g
    for k in np.random.default_rng(seed).permutation(len(small)):
        g = small[k]
        if len(train) + len(g) <= n_train_target:
            train += g
        else:
            test += g
    return np.sort(train), np.sort(test)


def overlap_report(train, test, scaffolds, fps, primary_docs, all_docs, pactivity):
    """
    Measure how much of the test set has a close relative in the training set.

    Four different yardsticks, because no single one settles it:
      scaffolds  - does the test molecule's ring skeleton appear in train?
      papers     - did a training molecule come from the same paper? Counted
                   twice: by the test molecule's first paper only, and by ANY
                   paper that measured it (strict; catches reference compounds
                   that appear in many papers)
      max Tanimoto - similarity to the nearest training molecule
      pActivity  - are the two sets comparable in label distribution?
    """
    train_scaffolds = {scaffolds[i] for i in train}
    test_scaffolds = {scaffolds[i] for i in test}
    train_primary = {primary_docs[i] for i in train}
    train_any = set()
    for i in train:
        train_any |= set(all_docs[i].split("|"))

    in_train = [scaffolds[i] in train_scaffolds for i in test]
    primary_in_train = [primary_docs[i] in train_primary for i in test]
    any_in_train = [bool(set(all_docs[i].split("|")) & train_any) for i in test]
    nn = max_tanimoto(fps[test], fps[train])

    return {
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "test_fraction": float(len(test) / (len(train) + len(test))),
        "n_test_scaffolds": len(test_scaffolds),
        "test_scaffolds_seen_in_train": len(test_scaffolds & train_scaffolds),
        "test_molecules_with_scaffold_in_train": int(sum(in_train)),
        "test_molecules_with_primary_paper_in_train": int(sum(primary_in_train)),
        "test_molecules_with_any_paper_in_train": int(sum(any_in_train)),
        "median_max_tanimoto_to_train": float(np.median(nn)),
        "frac_test_with_nn_above_0.7": float(np.mean(nn > 0.7)),
        "pactivity_mean_train": float(np.mean(pactivity[train])),
        "pactivity_mean_test": float(np.mean(pactivity[test])),
        "pactivity_sd_train": float(np.std(pactivity[train])),
        "pactivity_sd_test": float(np.std(pactivity[test])),
    }
