"""
PoseBusters audit of docked poses, with three explicit states per molecule.

PoseBusters runs a battery of physical-plausibility checks on a docked pose. Each check is
True (passes), False (fails) or missing (NaN: the check could not be evaluated).

Where missing values come from. In PoseBusters 0.6.5 a module that cannot run returns NaN
for its outputs (for example `energy_ratio.py` returns NaN when the molecule has no
conformer, does not sanitize, lacks force-field parameters, or cannot be converted to
InChI). A check that is merely *inapplicable* is NOT reported as NaN: with no cofactors or
waters in the receptor, the cofactor and water checks return True. So in this audit a
missing value always means "not evaluated", never "not applicable".

Why this module exists. The first audit (scripts/08f, results/08f_posebusters_by_molecule.csv)
chose check columns with `dtype == bool`. A check column holding True/False and any NaN is
object-typed, so it silently disappeared from the audit; `internal_energy` was dropped
that way (the saved report has 21 of the 22 checks), and a pose with a missing result could not
be told from a clean pass. The functions below name the required checks explicitly and keep
pass, fail and unevaluable apart.

Molecule-level status:
  fail         at least one required check is explicitly False
  unevaluable  no explicit failure, but at least one required check is missing
  pass         every required check is explicitly True
`pb_pass` is True only for status "pass".
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

# The 22 outputs of PoseBusters(config="dock") in version 0.6.5, by their reported names.
# Frozen here on purpose; tests/test_posebusters_audit.py compares it with the installed API.
REQUIRED_CHECKS = (
    "mol_pred_loaded", "mol_cond_loaded", "sanitization", "inchi_convertible",
    "all_atoms_connected", "no_radicals", "bond_lengths", "bond_angles",
    "internal_steric_clash", "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy", "protein-ligand_maximum_distance",
    "minimum_distance_to_protein", "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors", "minimum_distance_to_waters",
    "volume_overlap_with_protein", "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors", "volume_overlap_with_waters",
)

_SMILES_PROPERTY = re.compile(r">\s*<smiles>[^\n]*\n([^\n]*)\n")   # RDKit writes ">  <smiles>  (1) "


def read_sdf_records(path):
    """
    Read every record of an SDF file, keeping its position even if it cannot be parsed.

    Returns a list of dicts: position (0-based record number), pose_id (stable, derived from the
    position), smiles (the `smiles` property read from the raw text, so a record that RDKit
    cannot parse keeps its identity), and mol (RDKit molecule with hydrogens kept, or None).
    """
    text = Path(path).read_text()
    # Split on the record delimiter line, taking its newline with it, so every block starts at the
    # molecule-name line (which may legitimately be blank).
    blocks = [b for b in re.split(r"\$\$\$\$[ \t]*\r?\n?", text) if b.strip()]
    records = []
    for position, block in enumerate(blocks):
        found = _SMILES_PROPERTY.search(block)
        mol = Chem.MolFromMolBlock(block, removeHs=False)
        pose_id = f"pose{position:06d}"
        if mol is not None:
            mol.SetProp("_Name", pose_id)
        records.append({"position": position, "pose_id": pose_id,
                        "smiles": found.group(1).strip() if found else None, "mol": mol})
    return records


def run_posebusters(records, protein_path, buster):
    """
    Run PoseBusters on the parsed records and return its report indexed by pose_id.

    The pose id travels through PoseBusters in the molecule name, so the mapping back to the
    record never depends on row order. Raises if an id is repeated or any parsed pose is missing
    from the report.
    """
    parsed = [r for r in records if r["mol"] is not None]
    if not parsed:
        return pd.DataFrame()
    report = buster.bust([r["mol"] for r in parsed], None, str(protein_path), full_report=False)
    report = report.reset_index()
    report = report.rename(columns={"molecule": "pose_id"}).drop(columns=["file", "position"], errors="ignore")
    if report["pose_id"].duplicated().any():
        raise ValueError("PoseBusters returned a repeated pose id: "
                         f"{report.loc[report['pose_id'].duplicated(), 'pose_id'].tolist()[:5]}")
    missing = {r["pose_id"] for r in parsed} - set(report["pose_id"])
    if missing:
        raise ValueError(f"PoseBusters report lacks {len(missing)} parsed poses, e.g. {sorted(missing)[:5]}")
    return report.set_index("pose_id")


def aggregate_status(report, records, required=REQUIRED_CHECKS):
    """
    One row per SDF record with the state of every required check and an overall status.

    `report` is the PoseBusters report indexed by pose_id (parsed records only); `records` is the
    full record list, so an unparsed record still gets a row, marked as having every check missing.
    Raises ValueError if a required check is absent from the report as a column.
    """
    ids = [r["pose_id"] for r in records]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate pose ids among SDF records")
    absent = [c for c in required if len(report) and c not in report.columns]
    if absent:
        raise ValueError(f"required PoseBusters checks absent from the report: {absent}. "
                         "Not treating them as passes. Check the PoseBusters version and config.")
    table = report.reindex(ids)                       # unparsed records become all-missing rows
    checks = pd.DataFrame(index=ids)
    for name in required:
        column = table[name] if name in table.columns else pd.Series(np.nan, index=ids)
        checks[name] = pd.array(column.where(column.notna(), None), dtype="boolean")
    failed = (checks == False).fillna(False).astype(bool)      # noqa: E712  (nullable booleans)
    missing = checks.isna()
    n_failed = failed.sum(axis=1)
    n_missing = missing.sum(axis=1)
    status = np.where(n_failed > 0, "fail", np.where(n_missing > 0, "unevaluable", "pass"))

    out = pd.DataFrame({
        "pose_id": ids,
        "record_position": [r["position"] for r in records],
        "smiles": [r["smiles"] for r in records],
        "parsed": [r["mol"] is not None for r in records],
        "status": status,
        "pb_pass": status == "pass",
        "n_failed": n_failed.to_numpy(),
        "n_missing": n_missing.to_numpy(),
        "failed_checks": ["|".join(c for c in required if failed.loc[i, c]) for i in ids],
        "missing_checks": ["|".join(c for c in required if missing.loc[i, c]) for i in ids],
    })
    return pd.concat([out, checks.reset_index(drop=True)], axis=1)


def summarize(status_table, expected):
    """Counts with explicit denominators: expected, parsed, audited, passed, failed, unevaluable."""
    counts = status_table["status"].value_counts().to_dict()
    return {
        "expected_poses": int(expected),
        "records_read": int(len(status_table)),
        "parsed": int(status_table["parsed"].sum()),
        "audited": int((status_table["parsed"]).sum()),
        "passed": int(counts.get("pass", 0)),
        "failed": int(counts.get("fail", 0)),
        "unevaluable": int(counts.get("unevaluable", 0)),
        "denominator_for_rates": int(len(status_table)),
    }


def attach_to_molecules(molecules, status_table, key="smiles"):
    """
    Add audit status to a molecule table without dropping any molecule.

    `molecules` may list a molecule more than once (it can belong to several sets); the audit
    must hold exactly one row per key. Molecules with no audit row get status "no_audit_result"
    (counted, never treated as a pass). Raises if the audit keys are not unique.
    """
    if status_table[key].isna().any():
        raise ValueError("audit rows without an identity key; cannot join to molecules")
    if status_table[key].duplicated().any():
        dup = status_table.loc[status_table[key].duplicated(), key].tolist()[:3]
        raise ValueError(f"audit keys are not unique, e.g. {dup}")
    keep = ["status", "pb_pass", "n_failed", "n_missing", "failed_checks", "missing_checks"]
    merged = molecules.merge(status_table[[key] + keep], on=key, how="left", validate="many_to_one")
    merged["status"] = merged["status"].fillna("no_audit_result")
    merged["pb_pass"] = merged["pb_pass"].fillna(False).astype(bool)
    return merged
