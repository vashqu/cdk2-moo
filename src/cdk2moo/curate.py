"""
Curation: turning raw ChEMBL activity records into one labelled row per
molecule.

Every filter here is explicit and counted, so the funnel printed by
scripts/02_curate.py is an auditable account of what was discarded and why.
Nothing in this module touches the network.
"""

import numpy as np
import pandas as pd


def filter_records(df, keep_types, keep_relation="=", keep_units="nM",
                   drop_curator_flagged=True, mutant_regex=None,
                   keep_assay_types=None):
    """
    Apply the record-level filters and return (filtered_df, funnel).

    `funnel` is a list of (step, n_remaining) so the caller can print exactly
    where records were lost.
    """
    funnel = [("raw records", len(df))]

    df = df[df["standard_type"].isin(keep_types)]
    funnel.append((f"type in {keep_types}", len(df)))

    # Censored values ('IC50 > 10000 nM') carry real information - the compound
    # is weak - but there is no number to regress on. Dropping them biases what
    # remains toward actives; that bias is a limitation, not a bug, and it is
    # why the surrogate has little experience of weak compounds.
    df = df[df["standard_relation"] == keep_relation]
    funnel.append((f"relation == '{keep_relation}'", len(df)))

    df = df[df["standard_units"] == keep_units]
    funnel.append((f"units == '{keep_units}'", len(df)))

    if drop_curator_flagged:
        df = df[df["data_validity_comment"].isna()]
        funnel.append(("no curator flag", len(df)))

    # A mutant CDK2 is a different protein: the hinge residues (F82, L83, H84)
    # that an ATP-competitive inhibitor hydrogen-bonds to are changed, so the
    # measurement does not describe wild-type binding. Up to 3.8 log units
    # apart within one paper for the same compound.
    if mutant_regex is not None:
        is_mutant = df["assay_description"].fillna("").str.contains(
            mutant_regex, case=False, regex=True)
        df = df[~is_mutant]
        funnel.append(("not a mutant-CDK2 assay", len(df)))

    # Dropped for cleanliness, not effect: too few records to matter.
    if keep_assay_types is not None:
        df = df[df["assay_type"].isin(keep_assay_types)]
        funnel.append((f"assay_type in {keep_assay_types}", len(df)))

    df = df[df["standard_value"].notna() & (df["standard_value"] > 0)]
    funnel.append(("positive value present", len(df)))

    df = df[df["canonical_smiles"].notna()]
    funnel.append(("SMILES present", len(df)))

    return df.copy(), funnel


def add_pactivity(df, value_col="standard_value"):
    """
    Convert a concentration in nM to pActivity = -log10(concentration in M).

        1 nM   -> 9
        1 uM   -> 6
        1 mM   -> 3

    Higher is better. The log scale is the natural one: binding free energy is
    proportional to the log of the dissociation constant, so a fixed step in
    pActivity is a fixed step in energy. It also makes the error structure
    roughly symmetric, which matters for how we aggregate replicates below.
    """
    out = df.copy()
    out["pactivity"] = 9.0 - np.log10(out[value_col].astype(float))
    return out


def check_against_pchembl(df, tol=0.01):
    """
    Cross-check our pActivity against ChEMBL's own `pchembl_value`.

    ChEMBL computes the same quantity independently. If our numbers disagree on
    rows where theirs exists, our unit handling is wrong - a bug that would
    otherwise be invisible, since a systematically shifted label still trains a
    model that looks fine.

    Returns (n_compared, n_disagreeing, max_abs_diff).
    """
    sub = df[df["pchembl_value"].notna()].copy()
    if not len(sub):
        return 0, 0, float("nan")
    diff = (sub["pactivity"] - sub["pchembl_value"].astype(float)).abs()
    return len(sub), int((diff > tol).sum()), float(diff.max())


def aggregate_by_structure(df, key="inchikey"):
    """
    Collapse replicate measurements to one row per unique structure.

    Aggregation is by STANDARDIZED STRUCTURE, not by molecule_chembl_id.
    ChEMBL gives a free base and its hydrochloride salt different IDs; after
    standardization they are the same molecule. Deduplicating by ID would leave
    both in the dataset, and a random split would put one in train and one in
    test - the model would then "predict" a compound it had memorised.

    The label is the MEDIAN pActivity across replicates:

      - median, not mean, because a single measurement made at a different ATP
        concentration can sit a log unit away from the others and would drag a
        mean with it;
      - of pActivity, not of the raw nM value, because averaging in linear
        concentration space lets the weakest measurement dominate.

    The spread across replicates is retained for inspection, but it is NOT a
    usable estimate of label noise: most "replicates" here are the same
    measurement curated twice, so the spread is ~0 (see decisions.md D-09).
    """
    g = df.groupby(key)

    out = pd.DataFrame({
        "pactivity": g["pactivity"].median(),
        "n_measurements": g["pactivity"].size(),
        "pactivity_std": g["pactivity"].std(),
        "pactivity_range": g["pactivity"].max() - g["pactivity"].min(),
        "std_smiles": g["std_smiles"].first(),
        "mw": g["mw"].first(),
        "n_heavy_atoms": g["n_heavy_atoms"].first(),
        "n_documents": g["document_chembl_id"].nunique(),
        # Kept so a compound can be traced back to ChEMBL by hand.
        "molecule_chembl_ids": g["molecule_chembl_id"].apply(
            lambda s: "|".join(sorted(set(s.dropna())))
        ),
        # The paper a compound first came from. Analogue series are published
        # together, so this is how we detect them when splitting.
        "primary_document": g["document_chembl_id"].apply(
            lambda s: s.dropna().iloc[0] if s.notna().any() else None
        ),
        # Every paper that measured the compound. Reference inhibitors are
        # re-measured in many papers, so a split by primary_document alone can
        # still leave a test molecule's other papers on the train side.
        "all_documents": g["document_chembl_id"].apply(
            lambda s: "|".join(sorted(set(s.dropna())))
        ),
        "assay_types": g["assay_type"].apply(
            lambda s: "|".join(sorted(set(s.dropna())))
        ),
    }).reset_index()

    return out


def replicate_spread(agg):
    """
    Summarise how much repeated records of the same compound disagree.

    This is a diagnostic, not a noise floor. In principle the within-compound SD
    would bound any model's RMSE, but here most replicates are duplicate
    curation of one measurement, so the median SD is ~0 and says nothing about
    real inter-lab noise (decisions.md D-09). The mean SD and max range are
    dominated by the few compounds that do disagree, e.g. across mutant assays.
    """
    rep = agg[agg["n_measurements"] > 1]
    if not len(rep):
        return {"n_replicated": 0}
    return {
        "n_replicated": int(len(rep)),
        "frac_replicated": float(len(rep) / len(agg)),
        "median_within_compound_sd": float(rep["pactivity_std"].median()),
        "mean_within_compound_sd": float(rep["pactivity_std"].mean()),
        "median_range": float(rep["pactivity_range"].median()),
        "max_range": float(rep["pactivity_range"].max()),
    }
