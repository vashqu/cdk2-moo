"""Thin wrapper around the vendored graph-based genetic algorithm.

Will live here:

- configuration of population size, generations, mutation and crossover rates
- the per-generation loop, scoring the population through ``objectives``
- optional similarity constraint to the training distribution, which is the
  mitigation arm for H5
- per-generation logging of every molecule with its scores, written to disk

The GA implementation itself is not written here: it is vendored verbatim
under ``vendor`` so that third-party code stays unedited and attributable.
"""
