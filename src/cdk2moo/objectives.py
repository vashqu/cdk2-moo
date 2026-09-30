"""Optimization objectives behind one interface.

Will live here:

- predicted pActivity, from a loaded ``surrogate`` model
- QED, a 0-1 heuristic score of how 'drug-like' a molecule's properties are
- synthetic accessibility (SA), a 1-10 heuristic of how hard a molecule looks
  to make, where lower is easier
- a common wrapper so the GA scores single-objective and multi-objective arms
  through the same call, enabling the ablation arms

Sign conventions (which objectives are maximised) are fixed here once so no
caller has to remember them.
"""
