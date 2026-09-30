"""Non-dominated sorting over the objective triple.

Will live here:

- Pareto dominance test and non-dominated front extraction over
  (pActivity, QED, SA)
- front ranking, for overlaying successive fronts on one figure
- helpers to place reference sets (known actives, random ChEMBL) in the same
  objective space as the GA output

Plain non-dominated sorting only. NSGA-II selection is named as future work.
"""
