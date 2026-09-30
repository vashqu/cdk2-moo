"""Applicability domain: how far a molecule sits from what the model saw.

Will live here:

- maximum Tanimoto similarity of a molecule to the surrogate's training set
- nearest-neighbour summaries over a whole generation or molecule set
- empirical surrogate error as a function of that similarity, measured on
  held-out REAL molecules, which is what calibrates the reward-hacking claim

This module carries H2: if predicted activity climbs while max Tanimoto to
training falls, the optimizer is walking out of the domain where the
surrogate has any evidence.
"""
