"""Sanity checks on molecules and poses.

Will live here:

- PoseBusters pass rate, which checks whether a docked pose is physically
  sensible (bond lengths, clashes, geometry) rather than merely well scored
- PAINS substructure flags, for fragments that show up as false positives
  across unrelated assays
- the size confound check: docking scores correlate with heavy-atom count, so
  a heavy-atom-count-only 'scoring function' is included as a baseline
- internal diversity of a molecule set

Pose-level checks read files produced by ``docking`` but do not import it, so
this module survives the deletion of stage 6.
"""
