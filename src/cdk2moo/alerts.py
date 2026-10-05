"""
Medicinal-chemistry structural alerts, using RDKit's published filter sets.

An alert is a substructure that, in past experience, tends to cause trouble:
PAINS (pan-assay interference compounds) are motifs that show up as false
positives in many biochemical assays (by reacting, aggregating, fluorescing or
chelating); Brenk alerts are groups that are chemically reactive, unstable or
toxic-prone. An alert is a flag for a human, NOT proof of a bad molecule (many
approved kinase inhibitors carry some), and its absence is not proof of a good
one. Nothing here removes molecules; it only records what matched.

Three concepts are kept separate throughout this project:
  graph validity     RDKit can build the molecule with legal valences
  structural alerts  the substructure filters in this module
  domain support     similarity to the training set (applicability.py)
A molecule can be valid, alert-free and far from training data, or valid,
well-supported and full of alerts.
"""

from rdkit import Chem
from rdkit.Chem import FilterCatalog


def _catalog(which):
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(which)
    return FilterCatalog.FilterCatalog(params)


_pains = _catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_brenk = _catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)


def alert_hits(smiles):
    """Return (PAINS alert names, Brenk alert names) that match one molecule."""
    mol = Chem.MolFromSmiles(smiles)
    pains = sorted({m.GetDescription() for m in _pains.GetMatches(mol)})
    brenk = sorted({m.GetDescription() for m in _brenk.GetMatches(mol)})
    return pains, brenk


def alert_row(label, d):
    """One printed line of alert rates (share of molecules with a PAINS hit, a Brenk hit, either; mean Brenk hits) for the molecules in d."""
    return (f"  {label:<46}{len(d):>7}{100 * (d['n_pains'] > 0).mean():>9.0f}%"
            f"{100 * (d['n_brenk'] > 0).mean():>9.0f}%{100 * d['any_alert'].mean():>9.0f}%"
            f"{d['n_brenk'].mean():>10.2f}")
