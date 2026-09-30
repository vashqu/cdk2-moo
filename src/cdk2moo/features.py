"""Molecular featurisation.

Will live here:

- Morgan/ECFP4 bit fingerprints, using ``config.MORGAN_RADIUS`` and
  ``config.MORGAN_N_BITS``
- RDKit physicochemical descriptors (MW, cLogP, TPSA, rotatable bonds,
  heavy-atom count) used both as surrogate inputs and as confound checks
- Bemis-Murcko scaffold extraction, consumed by ``splits``

One SMILES in, one fixed-width vector out. No fitting, no state.
"""
