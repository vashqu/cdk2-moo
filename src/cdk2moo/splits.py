"""Train/test splitting.

Will live here:

- a random split (the optimistic baseline)
- a Bemis-Murcko scaffold split, which keeps molecules sharing a core
  skeleton on the same side and so measures generalisation to new chemotypes
- selection of the held-out known actives (top-decile pActivity, excluded
  from training) used as a reference set downstream

Both splits are reported side by side in every surrogate table. The scaffold
split is expected to look worse; that is the point, not a defect.
"""
