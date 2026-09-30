"""The activity surrogate: ECFP4 -> RandomForest.

Will live here:

- fit / predict / save / load for a RandomForest regressor on pActivity
- the label-scrambled control model (labels permuted before fitting), which
  is the baseline for H4
- honest evaluation under both splits from ``splits``

The surrogate is an instrument, not a result. It is deliberately not tuned
for a flattering R^2; a mediocre surrogate is the expected and acceptable
case, and its error structure is what stages 4 and 5 interrogate.
"""
