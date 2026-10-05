"""
How a proposed molecule becomes the molecule that is scored: the two preparation policies.

legacy     (every historical experiment) The edited graph is accepted if it sanitizes, is one piece and is
           in the size window; its canonical SMILES is scored as it is. No salt stripping, neutralising or
           tautomer canonicalization is applied, although the TRAINING molecules went through all of them
           (src/cdk2moo/standardize.py). A generated molecule can therefore be scored in a different
           representation than the training set it is compared with.

corrected  The edited graph is standardized exactly as the training molecules were (cleanup, largest
           fragment, neutralise, canonical tautomer, organic-element check), then validated again
           (sanitizes, one piece, size window, NO radical electrons), and its canonical SMILES is what is
           deduplicated, scored (activity, similarity, QED, SA) and stored. Everything downstream sees one
           representation. This changes the search and is a new experiment variant, never a recalculation of a
           legacy run.

Limits worth stating. Passing these checks does not make a molecule stable or synthesizable. Tautomer
canonicalization is only complete when RDKit's enumeration finishes; when it stops at its transform or
tautomer limit the outcome is flagged ("tautomer_incomplete") and two tautomers of one compound are not
guaranteed to map to the same structure (23 of 2,016 training molecules hit this limit).
"""

from rdkit import Chem

from cdk2moo.mutations import rejection_reason
from cdk2moo.standardize import standardize_with_status

POLICIES = ("legacy", "corrected")


def prepare_candidate(mol, policy, size_range, cache=None, rejects=None):
    """
    Turn an accepted-by-the-edit molecule into the molecule to score.

    Returns (smiles, reason, tautomer_status). `smiles` is the representation to score and store, or None
    if the policy rejects the molecule, in which case `reason` says why. `cache` (a dict) avoids
    re-standardizing the same proposal within a run; `rejects` (a dict) accumulates reason counts.
    """
    if policy not in POLICIES:
        raise ValueError(f"unknown preparation policy {policy!r}; choose from {POLICIES}")
    proposal = Chem.MolToSmiles(mol)
    if policy == "legacy":
        return proposal, None, None
    if cache is not None and proposal in cache:
        result = cache[proposal]
    else:
        result = _standardize_and_validate(proposal, size_range)
        if cache is not None:
            cache[proposal] = result
    smiles, reason, status = result
    if rejects is not None:
        if reason is not None:
            rejects[reason] = rejects.get(reason, 0) + 1
        elif status != "Completed":
            rejects["tautomer_incomplete"] = rejects.get("tautomer_incomplete", 0) + 1
    return result


def _standardize_and_validate(proposal, size_range):
    """The corrected policy for one proposal SMILES: standardize, validate, canonical SMILES."""
    std, reason, status = standardize_with_status(Chem.MolFromSmiles(proposal))
    if std is None:
        return None, f"standardize_failed: {reason}", status
    reason = rejection_reason(std, size_range, closed_shell=True)
    if reason is not None:
        return None, f"after_standardization_{reason}", status
    smiles = Chem.MolToSmiles(std)
    again = Chem.MolFromSmiles(smiles)                       # the stored string must rebuild the same molecule
    if again is None or Chem.MolToSmiles(again) != smiles:
        return None, "unstable_representation", status
    return smiles, None, status
