"""
Evaluating generated molecules against held-out reference molecules, with every term defined.

Separate, plainly named quantities (the first evaluation used labels that claimed more than the numbers show):

  proximity            ECFP4 (Morgan r=2, 2048-bit, binary, no chirality) Tanimoto to the nearest held-out active, and the share at or
                       above each threshold. Closeness in fingerprint space; not recovery of the molecule, not evidence of activity.
  fingerprint identity the ECFP4 bit vector equals a held-out active's. ECFP4 ignores stereochemistry, so stereoisomers (and some
                       different molecules) collide. NOT molecular identity.
  exact identity       SMILES strings of STANDARDIZED structures are equal, under three keys:
                         std_isomeric   the project's standardization (cleanup, largest fragment, neutralise, canonical tautomer),
                                        isomeric SMILES. RDKit's tautomer step removes sp3 and double-bond stereo where tautomerism could
                                        interconvert it (e.g. the alpha carbon of an amino acid), so some stereoisomers merge here.
                         std_flat       the same without stereochemistry (connectivity-level identity).
                         strict_stereo  cleanup, largest fragment and neutralise only: stereo kept everywhere, tautomers NOT merged.
  nearest reference    the reference molecule with the highest Tanimoto to the query. Ties go to the lowest reference index; the share of
                       queries with a tie, and the share for which ANY tied neighbour is a held-out active, are reported.
  reference prevalence the share of reference molecules that are held-out actives: a property of the reference set, NOT a chance or
                       null rate (nearest-neighbour assignments are not uniform).
  neighbour annotation among queries with a reference neighbour at >= the proximity threshold, the share whose nearest reference
                       molecule has MEASURED pActivity >= a cutoff. A label of the neighbour copied onto the query; nothing was measured
                       for the query.

TWO COUNTING DEFINITIONS, with different identity and representation (never mixed):
  weighted   every query occurrence counts once, and fingerprint-derived values use the RAW SMILES as supplied. This is the historical
             definition (the same molecule found by several seeds counts once per occurrence).
  distinct   each standardized isomeric structure counts once. The IDENTITY used to deduplicate is the std_isomeric key, and the
             REPRESENTATION used for every fingerprint-derived value is that same standardized SMILES, so a distinct metric is a function
             of the set of structures only: it cannot depend on query order or on which raw form of a structure came first. A query whose
             standardization fails is its own distinct identity (represented by its raw string) and is counted in the failure columns.

Self-matches. A reference molecule identical to the query (same std_isomeric key) is excluded from the nearest-neighbour search when
exclude_self=True (starting molecules and the no-search "remaining training" baseline are reference molecules). Generated molecules are
not excluded: reproducing a reference molecule is reported as a rediscovery. A query left with no reference after exclusion is counted in
`no_reference_neighbour` and left out of the neighbour-based shares, whose denominator `n_with_reference` is reported.
"""

import numpy as np
from rdkit import Chem

from cdk2moo.features import ecfp4
from cdk2moo.standardize import standardize_with_status, stereo_preserving_form

NAN = float("nan")


def standardized_keys(smiles_list, cache=None):
    """
    Per input SMILES: (std_isomeric, std_flat, strict_stereo, tautomer_status). A part that cannot be computed is None
    (status None when standardization failed before enumeration). `cache` maps SMILES to the tuple.
    """
    cache = {} if cache is None else cache
    keys = []
    for smiles in smiles_list:
        if smiles not in cache:
            raw = Chem.MolFromSmiles(smiles)
            mol, _, status = standardize_with_status(Chem.Mol(raw) if raw is not None else None)
            strict = stereo_preserving_form(Chem.Mol(raw) if raw is not None else None)
            cache[smiles] = (None if mol is None else Chem.MolToSmiles(mol),
                             None if mol is None else Chem.MolToSmiles(mol, isomericSmiles=False),
                             None if strict is None else Chem.MolToSmiles(strict), status)
        keys.append(cache[smiles])
    return keys


def build_reference(smiles_list, pactivity, cache=None):
    """The reference set: fingerprints, the identity keys and measured pActivity, in a fixed order."""
    keys = standardized_keys(smiles_list, cache)
    return {"smiles": list(smiles_list), "fps": ecfp4(list(smiles_list)).astype(np.float32),
            "iso": np.array([k[0] for k in keys], dtype=object), "flat": np.array([k[1] for k in keys], dtype=object),
            "strict": np.array([k[2] for k in keys], dtype=object), "pactivity": np.asarray(pactivity, dtype=float)}


def tanimoto_matrix(query_fps, reference_fps):
    """Tanimoto of every query row against every reference row (binary fingerprints)."""
    shared = query_fps @ reference_fps.T
    either = query_fps.sum(axis=1)[:, None] + reference_fps.sum(axis=1)[None, :] - shared
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(either > 0, shared / either, 0.0)


def nearest_reference(similarity, exclude=None):
    """
    Nearest reference per query. `exclude` is a list (one entry per query) of reference indices that may not be chosen.
    Returns (index, similarity, tied, available): index is the lowest index among ties, `tied` marks every reference sharing the
    maximum, and `available` is False for a query with no allowed reference left (its other outputs are then meaningless).
    """
    similarity = similarity.copy()
    if exclude is not None:
        for row, columns in enumerate(exclude):
            similarity[row, columns] = -1.0
    best = similarity.max(axis=1) if similarity.shape[1] else np.full(len(similarity), -1.0)
    available = best >= 0
    index = similarity.argmax(axis=1) if similarity.shape[1] else np.zeros(len(similarity), dtype=int)
    tied = similarity >= best[:, None] - 1e-9
    return index, best, tied, available


def _share(mask, rows):
    return float(np.mean(mask[rows])) if len(rows) else NAN


def _metrics(rows_smiles, rows_iso, rows_flat, rows_strict, reference, held, thresholds, annotation_cutoff, proximity_threshold, exclude_self, suffix):
    """All fingerprint-based and identity-based metrics for one list of query strings; keys end with `suffix`."""
    sim = tanimoto_matrix(ecfp4(rows_smiles).astype(np.float32), reference["fps"])
    sim_held = sim[:, held].max(axis=1)
    exclude = None
    if exclude_self:
        exclude = [np.where(reference["iso"] == k)[0] if k is not None else np.array([], dtype=int) for k in rows_iso]
    nn, nn_sim, tied, available = nearest_reference(sim, exclude)
    is_held = np.zeros(len(reference["smiles"]), dtype=bool)
    is_held[held] = True
    held_sets = {n: set(reference[n][held]) - {None} for n in ("iso", "flat", "strict")}
    everyone = np.arange(len(rows_smiles))
    with_ref = np.where(available)[0]
    near = np.where(available & (nn_sim >= proximity_threshold))[0]
    out = {f"median_proximity_to_held_actives_{suffix}": float(np.median(sim_held)),
           f"fingerprint_identical_{suffix}": _share(np.isclose(sim_held, 1.0), everyone),
           f"exact_std_isomeric_{suffix}": _share(np.array([k is not None and k in held_sets["iso"] for k in rows_iso]), everyone),
           f"exact_std_flat_{suffix}": _share(np.array([k is not None and k in held_sets["flat"] for k in rows_flat]), everyone),
           f"exact_strict_stereo_{suffix}": _share(np.array([k is not None and k in held_sets["strict"] for k in rows_strict]), everyone),
           f"nearest_is_heldout_{suffix}": _share(is_held[nn], with_ref),
           f"nearest_tie_share_{suffix}": _share(tied.sum(axis=1) > 1, with_ref),
           f"any_tied_nearest_is_heldout_{suffix}": _share((tied & is_held[None, :]).any(axis=1), with_ref),
           f"neighbour_annotation_potent_{suffix}": _share(reference["pactivity"][nn] >= annotation_cutoff, near),
           f"n_with_reference_{suffix}": int(len(with_ref)), f"no_reference_neighbour_{suffix}": int(len(everyone) - len(with_ref)),
           f"n_with_reference_neighbour_above_threshold_{suffix}": int(len(near))}
    for t in thresholds:
        out[f"proximity_ge_{t}_{suffix}"] = _share(sim_held >= t, everyone)
    return out


def evaluate_group(queries, reference, held_idx, thresholds, annotation_cutoff, proximity_threshold=0.6, exclude_self=False, cache=None):
    """
    Metrics for one group of query SMILES against held-out reference molecules `held_idx` (indices into `reference`).

    Returns a dict with `*_weighted` values (every occurrence, raw SMILES as supplied: the historical definition) and `*_distinct`
    values (each standardized isomeric structure once, represented by its standardized SMILES: order-invariant); see the module text.
    Also returns counts, standardization failures, incomplete tautomer enumerations and the reference prevalence. An empty group returns
    status "empty" with undefined (NaN) values instead of failing.
    """
    queries = list(queries)
    held = np.asarray(held_idx, dtype=int)
    out = {"n_weighted": len(queries), "reference_prevalence_of_heldout": float(np.isin(np.arange(len(reference["smiles"])), held).mean())}
    if not queries:
        return {**out, "status": "empty", "n_distinct": 0, "n_standardization_failed": 0, "n_distinct_standardization_failed": 0,
                "n_tautomer_incomplete": 0}
    if len(held) == 0:
        raise ValueError("no held-out actives to compare with")
    cache = {} if cache is None else cache
    keys = standardized_keys(queries, cache)
    iso = [k[0] for k in keys]
    failed = [k[0] is None for k in keys]
    incomplete = [k[3] is not None and k[3] != "Completed" for k in keys]

    # distinct identities, in sorted order so nothing depends on how the queries were ordered
    identity = [k[0] if k[0] is not None else "__failed__:" + q for k, q in zip(keys, queries)]
    distinct_ids = sorted(set(identity))
    representative = {}
    for ident, q, k in zip(identity, queries, keys):
        representative.setdefault(ident, q if k[0] is None else k[0])     # standardized SMILES; raw string only for failures
    rep_smiles = [representative[i] for i in distinct_ids]
    rep_keys = standardized_keys(rep_smiles, cache)

    weighted = _metrics(queries, iso, [k[1] for k in keys], [k[2] for k in keys], reference, held, thresholds, annotation_cutoff,
                        proximity_threshold, exclude_self, "weighted")
    distinct = _metrics(rep_smiles, [k[0] for k in rep_keys], [k[1] for k in rep_keys], [k[2] for k in rep_keys], reference, held,
                        thresholds, annotation_cutoff, proximity_threshold, exclude_self, "distinct")
    return {**out, "status": "ok", "n_distinct": len(distinct_ids), "n_standardization_failed": int(sum(failed)),
            "n_distinct_standardization_failed": int(sum(i.startswith("__failed__:") for i in distinct_ids)),
            "n_tautomer_incomplete": int(sum(incomplete)), **weighted, **distinct}


def proximity_to_held(queries, reference, held_idx):
    """ECFP4 Tanimoto of each query (raw SMILES as supplied) to its nearest held-out active, as an array (for plots and tables)."""
    queries = list(queries)
    if not queries:
        return np.array([])
    sim = tanimoto_matrix(ecfp4(queries).astype(np.float32), reference["fps"])
    return sim[:, np.asarray(held_idx)].max(axis=1)


def recall_of_held_actives(generated, reference, held_idx, threshold, cache=None):
    """
    Share of held-out actives with at least one DISTINCT generated structure at proximity >= threshold, and how many.

    Distinct means distinct std_isomeric identity (a structure whose standardization fails is its own identity), represented by its
    standardized SMILES, so the result does not depend on the order of `generated`. Counts proximity, not exact recovery.
    Returns (share, hits, n_distinct_generated).
    """
    generated = list(generated)
    if not generated:
        return NAN, 0, 0
    keys = standardized_keys(generated, {} if cache is None else cache)
    reps = {}
    for q, k in zip(generated, keys):
        reps.setdefault(k[0] if k[0] is not None else "__failed__:" + q, q if k[0] is None else k[0])
    smiles = [reps[i] for i in sorted(reps)]
    sim = tanimoto_matrix(ecfp4(smiles).astype(np.float32), reference["fps"][np.asarray(held_idx)])
    hit = (sim >= threshold).any(axis=0)
    return float(hit.mean()), int(hit.sum()), len(smiles)
