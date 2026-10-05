"""
Matched, deterministic initialization of the genetic algorithm.

Why this exists. A GA run starts from a population of training molecules. The two candidate-preparation policies (candidates.py)
treat starting molecules differently: the corrected policy standardizes and validates them, and refuses to start if that would
merge molecules or reject any. To compare the policies, and to pair a constrained run with its unconstrained twin, every run of
a given seed must start from the SAME molecules. So the eligibility rule is defined once, under the stricter (corrected) policy,
and used by both policies:

  eligible  heavy-atom count inside the size window, AND the stored structure is unchanged by standardization (idempotent), AND it
            is closed-shell and passes the corrected policy, AND its standardized string is unique among eligible molecules.

`eligible_rows` is computed once per campaign (scripts/03b_eligible_starts.py) and saved with the reason every excluded molecule was
excluded; `select_start_population` then draws the starting molecules from the eligible rows allowed in an experiment (for example,
only the remaining training molecules of a holdout control) with a seeded generator, so the draw depends only on the seed, the
eligible set and the allowed set.
"""

import numpy as np
from rdkit import Chem

from cdk2moo.candidates import prepare_candidate


def eligible_rows(smiles_list, n_heavy_atoms, size_range, cache=None):
    """
    Decide eligibility for each molecule. Returns a list of dicts (row, eligible, reason, tautomer_status), one per input, in
    input order. Reasons: size, standardization_failed, not_idempotent, rejected_<reason>, duplicate_after_standardization.
    """
    cache = {} if cache is None else cache
    low, high = size_range
    rows, seen = [], {}
    for row, (smiles, atoms) in enumerate(zip(smiles_list, n_heavy_atoms)):
        record = {"row": row, "eligible": False, "reason": "", "tautomer_status": ""}
        if not low <= atoms <= high:
            record["reason"] = "size"
        else:
            std, reason, status = prepare_candidate(Chem.MolFromSmiles(smiles), "corrected", size_range, cache=cache)
            record["tautomer_status"] = status or ""
            if std is None:
                record["reason"] = f"rejected_{reason}"
            elif std != Chem.MolToSmiles(Chem.MolFromSmiles(smiles)):
                record["reason"] = "not_idempotent"
            elif std in seen:
                record["reason"] = "duplicate_after_standardization"
            else:
                seen[std] = row
                record["eligible"] = True
        rows.append(record)
    return rows


def select_start_population(eligible_row_ids, allowed_row_ids, smiles_list, n, seed):
    """
    Draw n starting molecules: the eligible rows that are also allowed in this experiment, sorted by row, drawn with a generator seeded by
    `seed`. Raises if fewer than n are available (the population is never silently smaller). Returns (row ids, SMILES list).
    """
    pool = np.array(sorted(set(int(r) for r in eligible_row_ids) & set(int(r) for r in allowed_row_ids)))
    if len(pool) < n:
        raise ValueError(f"only {len(pool)} eligible starting molecules for a population of {n}")
    chosen = np.sort(np.random.default_rng(seed).choice(pool, n, replace=False))
    return chosen, [smiles_list[i] for i in chosen]


def load_eligible(curated, relaxed=False):
    """
    Read processed/eligible_starts.csv of the selected campaign scope and return the eligible row ids for the primary (or relaxed) size
    window. Verifies that the file lines up with the curated molecules (same length, same InChIKey at every row) before using it.
    """
    import pandas as pd
    from cdk2moo import config
    table = pd.read_csv(config.PROCESSED_DIR / "eligible_starts.csv")
    if len(table) != len(curated) or not (table["inchikey"].to_numpy() == curated["inchikey"].to_numpy()).all():
        raise ValueError("eligible_starts.csv does not line up with the curated molecules; rerun scripts/03b_eligible_starts.py")
    column = "eligible_relaxed" if relaxed else "eligible_primary"
    return table.loc[table[column], "row"].to_numpy()
