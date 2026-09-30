import pandas as pd
"""
ChEMBL access: finding the right target, and counting what data exists for it.

Nothing here writes to disk or makes curation decisions. Fetching and curation
live in later stages.
"""

from chembl_webresource_client.new_client import new_client


def search_targets(query):
    """
    Search ChEMBL for targets matching a name string.

    Returns a list of dicts. Note that a name search is fuzzy: querying 'CDK2'
    will also return protein complexes (CDK2/cyclin A), non-human orthologues,
    and cell lines. Filtering is the caller's job.
    """
    rows = []
    for hit in new_client.target.search(query):
        rows.append({
            "chembl_id": hit.get("target_chembl_id"),
            "pref_name": hit.get("pref_name"),
            "organism": hit.get("organism"),
            "target_type": hit.get("target_type"),
            "accessions": _accessions(hit),
        })
    return rows


def _accessions(hit):
    """
    Pull the UniProt accession(s) out of a target record.

    A ChEMBL target can be made of several protein components: a SINGLE
    PROTEIN has one, a PROTEIN COMPLEX has several. The accession is the
    authoritative identity of each component.
    """
    out = []
    for comp in hit.get("target_components") or []:
        acc = comp.get("accession")
        if acc:
            out.append(acc)
    return out


def match_by_accession(rows, uniprot, organism=None):
    """
    Select targets whose protein component is exactly the given UniProt
    accession, and which are a single protein rather than a complex.

    This is the identity check that matters. Matching on pref_name would be
    fragile; matching on accession is not.
    """
    out = []
    for r in rows:
        if r["target_type"] != "SINGLE PROTEIN":
            continue
        if uniprot not in r["accessions"]:
            continue
        if organism and (r["organism"] or "").lower() != organism.lower():
            continue
        out.append(r)
    return out


def count_activities(target_id, measurement_types, units="nM", relation="="):
    """
    Count bioactivity records for a target, per measurement type.

    Reports both the raw count and the count that survives two basic filters:

      relation == '='  drops censored values such as 'IC50 > 10000 nM', which
                       carry real information but cannot be regressed on. Note
                       this biases the surviving set toward active compounds.

      units == 'nM'    drops records ChEMBL could not normalise, which are
                       usually a different kind of measurement or a curation
                       oddity. Cheap insurance against unit-conversion bugs.

    Each count is a separate API call, so this takes a minute or two.
    """
    api = new_client.activity
    total = len(api.filter(target_chembl_id=target_id))

    per_type = {}
    for t in measurement_types:
        n_all = len(api.filter(target_chembl_id=target_id, standard_type=t))
        n_clean = len(api.filter(
            target_chembl_id=target_id,
            standard_type=t,
            standard_units=units,
            standard_relation=relation,
        ))
        per_type[t] = {"all": n_all, "clean": n_clean}

    return {"total": total, "per_type": per_type}

# --------------------------------------------------------------------------
# Raw download
# --------------------------------------------------------------------------

ACTIVITY_FIELDS = [
    "activity_id",
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_type",
    "standard_value",
    "standard_units",
    "standard_relation",
    "pchembl_value",
    "assay_chembl_id",
    "assay_type",
    "assay_description",
    "document_chembl_id",
    "target_chembl_id",
    "data_validity_comment",
    "potential_duplicate",
    "standard_flag",
]


def fetch_activities(target_id, standard_types, fields=None, verbose=True):
    """
    Download every activity record for a target and set of measurement types.

    NO filtering beyond the measurement type. Censored values, odd units and
    curator-flagged records all come through, because this is the raw layer:
    decisions belong downstream where they can be revisited without another
    API crawl.
    """
    fields = fields or ACTIVITY_FIELDS
    api = new_client.activity

    frames = []
    for t in standard_types:
        query = api.filter(target_chembl_id=target_id, standard_type=t).only(fields)
        if verbose:
            print(f"  fetching {t} ...", end="", flush=True)
        df = pd.DataFrame.from_records(list(query))
        if verbose:
            print(f" {len(df)} records")
        if len(df):
            frames.append(df)

    if not frames:
        return pd.DataFrame(columns=fields)

    out = pd.concat(frames, ignore_index=True)
    for f in fields:
        if f not in out.columns:
            out[f] = None
    return out[fields]
