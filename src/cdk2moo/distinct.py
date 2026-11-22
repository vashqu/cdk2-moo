"""
The counting unit of the analysis: DISTINCT STANDARDIZED MOLECULES.

The campaign plan names distinct standardized molecules as the primary unit. Under the `corrected` preparation policy this holds by
construction: every generated candidate is standardized like the training data (candidates.py) and its canonical SMILES is what is stored, and
the GA never admits a SMILES already present in the population (ga.py). So within one population (one arm, seed, floor and generation) each row
is a different standardized molecule, and "rows" and "distinct standardized molecules" are the same thing.

The `legacy` policy does not give that guarantee (it stores candidates as edited, so two rows can be one compound), which is why population
tables from it must never be summarized as if they were distinct molecules.

require_distinct_standardized() turns this from an assumption into a check at the point where an analysis loads populations: it refuses a table
that is not entirely `corrected`, or that holds the same SMILES twice in one population.
"""

import pandas as pd

POPULATION_KEYS = ["arm", "seed", "generation"]


def require_distinct_standardized(population, name="population table"):
    """
    Raise ValueError unless every row of `population` was prepared under the corrected policy and no SMILES occurs twice within one population.
    Returns a small dict of counts for logging. Tables without a `min_similarity` column (unconstrained runs) and with one (constrained runs, NaN or
    a floor value) are both handled; the floor is part of the population's identity.
    """
    if "policy" not in population.columns:
        raise ValueError(f"{name}: no `policy` column, so it cannot be shown that the molecules are standardized")
    other = sorted(set(population["policy"].dropna().unique()) - {"corrected"})
    if other or population["policy"].isna().any():
        raise ValueError(f"{name}: rows prepared under policy {other or 'missing'}; only corrected populations are distinct standardized molecules")
    keys = [k for k in POPULATION_KEYS if k in population.columns]
    if "min_similarity" in population.columns:
        keys.append("min_similarity")
    frame = population.assign(**{k: population[k].astype(object).where(population[k].notna(), "none") for k in keys})
    duplicated = int(frame.duplicated(keys + ["smiles"]).sum())
    if duplicated:
        raise ValueError(f"{name}: {duplicated} rows repeat a SMILES within one population; rows would not be distinct molecules")
    return {"rows": int(len(population)), "populations": int(frame.groupby(keys).ngroups)}
