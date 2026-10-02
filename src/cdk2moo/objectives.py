"""
Scoring generated molecules, and turning scores into GA fitness.

Every generated molecule gets the same panel of numbers, whichever arm made it:

  pred_real       predicted pActivity from the real surrogate
  pred_scrambled  predicted pActivity from the label-scrambled surrogate
  qed             drug-likeness, 0-1 (see below)
  sa              synthetic accessibility, 1 (easy) to 10 (hard)
  max_tanimoto    similarity to the nearest TRAINING molecule (all 2,016)
  n_heavy_atoms   size

Recording all of them for every arm is deliberate: the arm that followed the
scrambled surrogate still gets a pred_real, which is how we ask whether it
"improved" by the real surrogate's standards.

QED (quantitative estimate of drug-likeness) combines eight properties (weight,
logP, hydrogen-bond donors/acceptors, polar surface area, rotatable bonds,
aromatic rings, structural alerts) into one number by matching their
distributions in known oral drugs. It therefore rewards resembling existing
drugs, which is partly circular.

SA score (Ertl & Schuffenhauer) is a heuristic: molecules built from fragments
that are common in known compounds score easy, unusual fragments, stereocentres
and macrocycles score hard. It is not a synthetic route.
"""

import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDConfig
from rdkit.Chem import QED

from cdk2moo import config
from cdk2moo.features import ecfp4, max_tanimoto

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402  (ships with RDKit's Contrib directory)

ARMS = ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]


def qed_and_sa(smiles_list):
    """Return (QED list, SA-score list) for a list of SMILES."""
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    return [QED.qed(m) for m in mols], [sascorer.calculateScore(m) for m in mols]


def score_molecules(smiles_list, real_forest, scrambled_forest, train_fps):
    """Return a DataFrame with one row per SMILES and the columns listed above."""
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    qed, sa = qed_and_sa(smiles_list)
    fps = ecfp4(smiles_list, config.FP_RADIUS, config.FP_BITS)
    return pd.DataFrame({
        "smiles": smiles_list,
        "pred_real": real_forest.predict(fps),
        "pred_scrambled": scrambled_forest.predict(fps),
        "qed": qed,
        "sa": sa,
        "max_tanimoto": max_tanimoto(fps, train_fps),
        "n_heavy_atoms": [m.GetNumHeavyAtoms() for m in mols],
    })


def fitness(scores, arm):
    """
    The number the GA maximises, per arm. Each term is scaled to 0-1 and terms
    are combined by GEOMETRIC mean, so a molecule must be decent on every term:
    a zero on any one sinks the whole score (the standard choice in molecular
    optimisation benchmarks such as GuacaMol).

    Arms differ in exactly one respect each:
      multi_real       activity (real surrogate) x QED x SA     <- the main arm
      multi_scrambled  activity (scrambled surrogate) x QED x SA  <- control (H4)
      activity_only    activity (real surrogate) alone          <- ablation
      druglike_only    QED x SA, no surrogate at all            <- ablation
    """
    span = config.ACTIVITY_HIGH - config.ACTIVITY_LOW
    act_real = np.clip((scores["pred_real"] - config.ACTIVITY_LOW) / span, 0, 1)
    act_scrambled = np.clip((scores["pred_scrambled"] - config.ACTIVITY_LOW) / span, 0, 1)
    qed = scores["qed"]
    sa = (10 - scores["sa"]) / 9

    if arm == "multi_real":
        return (act_real * qed * sa) ** (1 / 3)
    if arm == "multi_scrambled":
        return (act_scrambled * qed * sa) ** (1 / 3)
    if arm == "activity_only":
        return act_real
    if arm == "druglike_only":
        return (qed * sa) ** 0.5
    raise ValueError(f"unknown arm: {arm}")
