"""Stage 1b: turning raw ChEMBL rows into a clean molecule table.

Will live here:

- salt stripping (keeping the largest organic fragment of a salt form)
- neutralization of charged groups to a consistent protonation state
- canonical tautomer selection
- canonical SMILES generation and de-duplication on that canonical form
- aggregation of replicate assay values, and conversion to pActivity

Ordering matters for leakage: de-duplication and replicate aggregation happen
HERE, before any split is drawn in ``splits``. Aggregating after splitting
would put the same molecule on both sides.
"""
