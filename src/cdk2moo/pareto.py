"""
Pareto dominance.

With several objectives there is usually no single best molecule. Molecule A
DOMINATES molecule B if A is at least as good on every objective and strictly
better on at least one; B is then never a sensible choice. The molecules that no
other molecule dominates form the Pareto front: each is the best available
answer for some way of weighing the objectives, and moving along the front
trades one objective against another.
"""

import numpy as np


def non_dominated(points, chunk=400):
    """
    Boolean mask of the rows of `points` that nothing else dominates.

    `points` is (n, k) with every column to be MAXIMISED (negate a column that
    should be minimised). Identical rows do not dominate each other.
    """
    points = np.asarray(points, dtype=float)
    keep = np.ones(len(points), dtype=bool)
    for start in range(0, len(points), chunk):
        block = points[start:start + chunk]                       # candidates
        at_least_as_good = (points[None, :, :] >= block[:, None, :]).all(axis=2)
        strictly_better = (points[None, :, :] > block[:, None, :]).any(axis=2)
        keep[start:start + chunk] = ~(at_least_as_good & strictly_better).any(axis=1)
    return keep


def dominated_by_any(points, front):
    """For each row of `points`, is it dominated by at least one row of `front`? (all maximised)"""
    points, front = np.asarray(points, float), np.asarray(front, float)
    at_least_as_good = (front[None, :, :] >= points[:, None, :]).all(axis=2)
    strictly_better = (front[None, :, :] > points[:, None, :]).any(axis=2)
    return (at_least_as_good & strictly_better).any(axis=1)
