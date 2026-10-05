"""
Read-only leakage and artifact audits over the stored splits, fingerprints and docking files.

Nothing here changes a split or a result. Each function answers one narrow question, so a reader can see
what "overlap" means in each case:

  identity overlap          a test molecule shares a standardized-structure key with a training molecule. Two keys:
                            the full InChIKey (stereochemistry included, as standardized) and its first block, the
                            connectivity layer (stereo-agnostic).
  fingerprint overlap       a test molecule's ECFP4 bit vector (Morgan r=2, 2048-bit, no chirality) is bit-for-bit equal to
                            some training molecule's. A different molecule can collide (stereoisomers do); this is not identity.
  scaffold overlap          a test molecule's Bemis-Murcko scaffold string appears among training scaffolds.
  document overlap          a test molecule was reported in at least one document that also reported a training molecule,
                            using ALL documents of each molecule (the `all_documents` column), not only the primary one.

`primary_document` (stage 2) is the first non-null document met when the records of a structure are aggregated, in the
order of the downloaded file. It was described elsewhere as the "first paper"; that is not verified: no publication
date was checked, so it is NOT the earliest publication and the paper split is NOT a temporal split.
"""

import numpy as np


def check_alignment(split_rows, curated, id_column="inchikey"):
    """
    Return the split table for ONE seed ordered by row, after verifying it lines up with the curated table.

    Requires a `row` column covering 0..n-1 exactly once and an identity key equal to the curated table's at the same row.
    Raises ValueError with the reason if keys are missing or misaligned (for example if inputs were reordered separately).
    """
    for column in ("row", id_column):
        if column not in split_rows.columns:
            raise ValueError(f"split table lacks the identity column {column!r}")
        if column not in curated.columns and column != "row":
            raise ValueError(f"curated table lacks the identity column {column!r}")
    if split_rows[id_column].isna().any() or curated[id_column].isna().any():
        raise ValueError(f"missing {id_column} values; cannot align splits to molecules")
    ordered = split_rows.sort_values("row").reset_index(drop=True)
    if len(ordered) != len(curated) or not (ordered["row"].to_numpy() == np.arange(len(curated))).all():
        raise ValueError("split rows do not cover 0..n-1 exactly once for the curated molecules")
    if not (ordered[id_column].to_numpy() == curated[id_column].to_numpy()).all():
        raise ValueError(f"split table and curated table disagree on {id_column} at the same row (misaligned inputs)")
    return ordered


def key_overlap(train_idx, test_idx, keys):
    """
    (test molecules whose key occurs in training, distinct shared keys) for any per-molecule key array.
    """
    keys = np.asarray(keys, dtype=object)
    train_keys = set(keys[np.asarray(train_idx)])
    test_keys = keys[np.asarray(test_idx)]
    hits = [k in train_keys for k in test_keys]
    return int(sum(hits)), len({k for k, h in zip(test_keys, hits) if h})


def fingerprint_keys(fps):
    """One hashable key per fingerprint row (the packed bits), so equal bit vectors have equal keys."""
    return np.array([row.tobytes() for row in np.packbits(np.asarray(fps, dtype=np.uint8), axis=1)], dtype=object)


def connectivity_keys(inchikeys):
    """First block of each InChIKey: the connectivity layer, which ignores stereochemistry."""
    return np.array([str(k).split("-")[0] for k in inchikeys], dtype=object)


def document_overlap(train_idx, test_idx, documents, primary=None):
    """
    Test molecules sharing at least one document with training: using every document ("pipe-joined" string per molecule),
    and, if `primary` is given, using only each molecule's primary document, for comparison.
    """
    lists = [set(str(d).split("|")) - {""} for d in documents]
    train_docs = set().union(*[lists[i] for i in train_idx]) if len(train_idx) else set()
    any_hits = sum(bool(lists[i] & train_docs) for i in test_idx)
    out = {"any_document": int(any_hits)}
    if primary is not None:
        train_primary = {primary[i] for i in train_idx}
        out["primary_document"] = int(sum(primary[i] in train_primary for i in test_idx))
    return out
