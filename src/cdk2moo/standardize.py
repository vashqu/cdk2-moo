"""
Structure standardization.

Two molecules that are chemically the same must end up as the same string, or
deduplication silently fails and the same compound lands on both sides of a
train/test split. ChEMBL stores what the depositing paper reported, so the same
compound appears as a free base in one entry and a hydrochloride salt in
another, drawn with different tautomers, sometimes with a counter-ion or a
solvent of crystallisation attached.

The pipeline below reduces each of those to one canonical form.
"""

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

# RDKit is loud about every valence problem it meets; we handle failures by
# returning None and counting them, so the log noise is not useful.
RDLogger.DisableLog("rdApp.*")

# Elements that occur in ordinary organic medicinal chemistry. Anything else
# (a transition metal, an unusual main-group element) is a organometallic or a
# curation artifact, and our fingerprints and docking treatment are not valid
# for it.
ORGANIC_ELEMENTS = {
    "H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Se", "Br", "I",
}

_uncharger = rdMolStandardize.Uncharger()
_tautomer_enumerator = rdMolStandardize.TautomerEnumerator()


def standardize_mol(mol, canonical_tautomer=True):
    """
    Reduce a molecule to a canonical form. Returns (mol, None) or (None, reason).

    Steps, in order:

    1. Cleanup            - sanitize, fix common drawing conventions (e.g. a
                            nitro group drawn with a pentavalent nitrogen is
                            rewritten to the charge-separated form).

    2. FragmentParent     - keep only the largest organic fragment. This is
                            "salt stripping": a hydrochloride salt is stored as
                            two disconnected fragments, the drug and a chloride
                            ion. The chloride is not part of the molecule that
                            binds the protein.

    3. Uncharge           - neutralise where chemistry allows. A carboxylate
                            drawn as COO- in one paper and COOH in another is
                            the same compound at the same pH; we pick one.

    4. Canonical tautomer - tautomers are isomers that interconvert by moving a
                            proton and a double bond (keto/enol, amide/imidic
                            acid). They are genuinely the same substance in
                            solution, but different SMILES and different
                            fingerprints. RDKit scores all tautomers by a fixed
                            rule set and returns the top one, so that any
                            tautomer of a molecule maps to the same output.
                            This is the slowest step, roughly 10-50 ms/molecule.
    """
    if mol is None:
        return None, "unparseable"

    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = _uncharger.uncharge(mol)
        if canonical_tautomer:
            mol = _tautomer_enumerator.Canonicalize(mol)
    except Exception as exc:
        return None, f"standardization_error: {type(exc).__name__}"

    if mol is None or mol.GetNumAtoms() == 0:
        return None, "empty_after_standardization"

    symbols = {a.GetSymbol() for a in mol.GetAtoms()}
    if not symbols <= ORGANIC_ELEMENTS:
        bad = ",".join(sorted(symbols - ORGANIC_ELEMENTS))
        return None, f"inorganic: {bad}"

    return mol, None


def standardize_with_status(mol):
    """
    The same steps as standardize_mol (cleanup, salt stripping, neutralising, canonical tautomer, element
    check), but also reports whether tautomer enumeration ran to completion.

    Returns (mol, None, status) or (None, reason, status); `status` is the name of RDKit's
    TautomerEnumeratorStatus ("Completed", "MaxTransformsReached", "MaxTautomersReached", "Canceled") or
    None if enumeration was not reached. Only "Completed" shows the search for the canonical tautomer
    was exhaustive; otherwise the result is a canonical choice among the tautomers found so far, and two
    tautomers of one compound are NOT guaranteed to map to the same structure. In the 2,016 curated
    training molecules enumeration was incomplete for 23 (1.1%), all "MaxTransformsReached".
    """
    if mol is None:
        return None, "unparseable", None
    status = None
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = _uncharger.uncharge(mol)
        enumeration = _tautomer_enumerator.Enumerate(mol)
        status = str(enumeration.status).split(".")[-1]
        mol = _tautomer_enumerator.PickCanonical(enumeration)
    except Exception as exc:
        return None, f"standardization_error: {type(exc).__name__}", status
    if mol is None or mol.GetNumAtoms() == 0:
        return None, "empty_after_standardization", status
    symbols = {a.GetSymbol() for a in mol.GetAtoms()}
    if not symbols <= ORGANIC_ELEMENTS:
        return None, f"inorganic: {','.join(sorted(symbols - ORGANIC_ELEMENTS))}", status
    return mol, None, status


def stereo_preserving_form(mol):
    """
    Cleanup, largest fragment and neutralising ONLY: no tautomer canonicalization, so stereochemistry is kept as given.

    Why it exists: RDKit's tautomer canonicalization removes sp3 and double-bond stereochemistry at positions where
    tautomerism could interconvert it (by default, e.g. the alpha carbon of an amino acid), so the standard
    pipeline above merges those stereoisomers. This form is the strict alternative for identity checks; it does not
    merge tautomers. Returns a molecule or None.
    """
    if mol is None:
        return None
    try:
        mol = _uncharger.uncharge(rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(mol)))
    except Exception:
        return None
    return mol if mol is not None and mol.GetNumAtoms() > 0 else None


def standardize_smiles(smiles, canonical_tautomer=True):
    """
    Standardize one SMILES string.

    Returns a dict of the standardized structure and a few descriptors, or a
    dict with `ok=False` and a reason. Never raises: bad input is data, not an
    exception, and we want to count it.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return {"ok": False, "reason": "missing_smiles"}

    mol, reason = standardize_mol(Chem.MolFromSmiles(smiles),
                                  canonical_tautomer=canonical_tautomer)
    if mol is None:
        return {"ok": False, "reason": reason}

    return {
        "ok": True,
        "reason": None,
        # The canonical SMILES is our structure identifier for joins.
        "std_smiles": Chem.MolToSmiles(mol),
        # InChIKey is a hash of the structure: a fixed-length string that is
        # the same for the same molecule regardless of how it was drawn. It is
        # the safer key to deduplicate on.
        "inchikey": Chem.MolToInchiKey(mol),
        "mw": Descriptors.MolWt(mol),
        "n_heavy_atoms": mol.GetNumHeavyAtoms(),
    }


def standardize_frame(df, smiles_col="canonical_smiles",
                      canonical_tautomer=True, verbose=True):
    """
    Standardize a DataFrame's SMILES column, adding the result as new columns.

    Rows that fail keep a reason in `std_fail_reason` rather than being dropped
    here: dropping is the caller's decision and should be visible in the funnel.
    """
    import pandas as pd

    # Standardize each UNIQUE input string once, not each row. The tautomer
    # step is the expensive part and the same SMILES recurs across records.
    uniq = df[smiles_col].dropna().unique()
    if verbose:
        print(f"  standardizing {len(uniq)} unique input SMILES "
              f"(tautomer canonicalization: {canonical_tautomer}) ...")

    cache = {}
    for i, smi in enumerate(uniq, 1):
        cache[smi] = standardize_smiles(smi, canonical_tautomer=canonical_tautomer)
        if verbose and i % 500 == 0:
            print(f"    {i}/{len(uniq)}")

    empty = {"ok": False, "reason": "missing_smiles"}
    results = [cache.get(s, empty) if isinstance(s, str) else empty
               for s in df[smiles_col]]

    res = pd.DataFrame(results, index=df.index)
    out = df.copy()
    out["std_ok"] = res["ok"]
    out["std_fail_reason"] = res["reason"]
    for col in ("std_smiles", "inchikey", "mw", "n_heavy_atoms"):
        out[col] = res[col] if col in res.columns else None

    return out