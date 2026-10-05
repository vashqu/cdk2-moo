"""
Validation of a docking job store: is every committed job internally consistent, and is each pose really that job's pose?

A readable SDF is not evidence that it belongs to the job whose record points at it, so reconciliation checks, for every record:

  schema          the record has exactly the expected fields, with the right types, an allowed status, and finite numbers where
                  a number is required (an "ok" job needs a finite score)
  identity chain  the job id recomputed from (molecule, Vina seed, preparation seed) equals the file name, the record's job_id and
                  the job_id stored inside the pose
  pose artifact   exactly one pose record, whose SHA-256 equals the hash committed in the record
  pose metadata   the pose's `smiles` and `vina_score` properties equal the record's
  pose chemistry  the pose is the same molecule as the job's input (see pose_identity)

Findings are sorted into distinct classes, never collapsed:
  complete, failed (a valid record of a recorded failure), incomplete_records (malformed JSON or schema violations),
  id_mismatch, missing_poses, unreadable_poses, pose_mismatch (hash, count, metadata or chemistry), orphan_poses,
  duplicate_keys, and `unattributable` (a file whose job identity cannot be recovered).
"""

import hashlib
import json
import math
from pathlib import Path

from rdkit import Chem

OK, PREP_FAILED, DOCK_FAILED, POSE_FAILED = "ok", "prep_failed", "dock_failed", "pose_export_failed"
STATUSES = (OK, PREP_FAILED, DOCK_FAILED, POSE_FAILED)
# field -> accepted types; None is allowed where noted in validate_record
SCHEMA = {"job_id": str, "smiles": str, "vina_seed": int, "prep_seed": int, "status": str, "vina_score": (int, float, type(None)),
          "n_heavy_atoms": (int, type(None)), "n_rotatable": (int, type(None)), "seconds": (int, float), "error": str,
          "pose_sha256": (str, type(None))}
ANOMALY_CLASSES = ["orphan_poses", "missing_poses", "unreadable_poses", "pose_mismatch", "incomplete_records", "id_mismatch",
                   "duplicate_keys", "unattributable"]


def job_id(smiles, vina_seed, prep_seed):
    """Stable identifier of one job: hash of molecule, Vina seed and ligand-preparation seed."""
    text = json.dumps([smiles, int(vina_seed), int(prep_seed)])
    return hashlib.sha1(text.encode()).hexdigest()[:16]


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_record(record):
    """List of schema problems (empty if the record is well formed). Booleans are not accepted as numbers."""
    if not isinstance(record, dict):
        return ["record is not an object"]
    problems = [f"missing field {k}" for k in SCHEMA if k not in record]
    problems += [f"unexpected field {k}" for k in record if k not in SCHEMA]
    for key, types in SCHEMA.items():
        if key in record and (isinstance(record[key], bool) or not isinstance(record[key], types)):
            problems.append(f"field {key} has type {type(record[key]).__name__}")
    if problems:
        return problems
    if record["status"] not in STATUSES:
        problems.append(f"status {record['status']!r} is not one of {STATUSES}")
    for key in ("vina_score", "seconds"):
        value = record[key]
        if value is not None and not math.isfinite(value):
            problems.append(f"{key} is not finite")
    if record["status"] == OK:
        if record["vina_score"] is None:
            problems.append("ok job without a score")
        if not record["pose_sha256"]:
            problems.append("ok job without a pose hash")
    elif record["pose_sha256"] is not None:
        problems.append("failed job carries a pose hash")
    return problems


def read_pose(path):
    """(molecules, ok): every record of an SDF file, with hydrogens kept. A record RDKit cannot parse is None."""
    try:
        return list(Chem.SDMolSupplier(str(path), removeHs=False)), True
    except Exception:
        return [], False


def pose_identity(input_smiles, pose):
    """
    Is the pose the same molecule as the job's input? Returns (ok, reason).

    Definition. Hydrogens are ignored on both sides (RemoveHs; the docked pose has only polar hydrogens). The two structures must
    have the same InChIKey connectivity block (formula, connections, hydrogen layer), which is independent of 3-D coordinates. When
    the input SMILES specifies EVERY stereocentre, the full InChIKey including the stereo layer must also match. When the input leaves
    a stereocentre unspecified, the 3-D embedding chooses it, so stereochemistry is not compared (connectivity only). Preparation may
    change the representation (explicit hydrogens, 3-D coordinates, atom order); none of that changes this identity. Fingerprint
    equality is never used as identity.
    """
    inp = Chem.MolFromSmiles(input_smiles)
    if inp is None or pose is None:
        return False, "structure could not be built"
    key_in, key_pose = Chem.MolToInchiKey(inp), Chem.MolToInchiKey(Chem.RemoveHs(pose))
    if not key_in or not key_pose:
        return False, "InChIKey unavailable"
    if key_in.split("-")[0] != key_pose.split("-")[0]:
        return False, "pose has different connectivity from the job's molecule"
    centres = Chem.FindMolChiralCenters(inp, includeUnassigned=True, useLegacyImplementation=False)
    if centres and all(label != "?" for _, label in centres) and key_in != key_pose:
        return False, "stereochemistry specified in the input differs from the pose"
    return True, ""


def check_pose(job_dir, jid, record):
    """The problems with one ok job's pose artifact, as (class, reason) or None if it is consistent with its record."""
    pose_path = Path(job_dir) / f"{jid}.sdf"
    if not pose_path.exists():
        return "missing_poses", "no pose file"
    mols, readable = read_pose(pose_path)
    if not readable or not mols:
        return "unreadable_poses", "file holds no readable record"
    if len(mols) != 1:
        return "pose_mismatch", f"{len(mols)} pose records, expected exactly 1"
    mol = mols[0]
    if mol is None:
        return "unreadable_poses", "pose record cannot be parsed"
    if sha256_of(pose_path) != record["pose_sha256"]:
        return "pose_mismatch", "pose file hash differs from the committed hash"
    # Read properties as raw strings. GetPropsAsDict would turn a job id made only of digits (or digits with one "e") into a number:
    # "4445476191389595" -> 4445476191389595.0 and "84868095015190e9" -> 8.486809501519e+22, which looked like a different job.
    if not (mol.HasProp("job_id") and mol.HasProp("smiles") and mol.HasProp("vina_score")):
        return "pose_mismatch", "pose lacks job_id, smiles or vina_score properties"
    if mol.GetProp("job_id") != jid:
        return "pose_mismatch", f"pose carries job_id {mol.GetProp('job_id')!r}"
    if mol.GetProp("smiles") != record["smiles"]:
        return "pose_mismatch", "pose smiles property differs from the record"
    try:
        if abs(float(mol.GetProp("vina_score")) - record["vina_score"]) > 6e-4:
            return "pose_mismatch", "pose vina_score property differs from the record"
    except ValueError:
        return "pose_mismatch", "pose has no readable vina_score property"
    same, reason = pose_identity(record["smiles"], mol)
    if not same:
        return "pose_mismatch", reason
    return None


def reconcile(run_dir):
    """
    Check every artifact in a run directory. Returns a dict with the classes listed in the module text (each a list of job ids or
    file names), `records` (id -> valid record) and `details` (id or name -> reason).
    """
    jobs = Path(run_dir) / "jobs"
    report = {"complete": [], "failed": [], "records": {}, "details": {}}
    for name in ANOMALY_CLASSES:
        report[name] = []
    seen = {}
    for path in sorted(jobs.glob("*.json")):
        stem = path.stem
        try:
            record = json.loads(path.read_text())
        except (json.JSONDecodeError, UnicodeDecodeError):
            report["incomplete_records"].append(stem)
            report["details"][stem] = "record is not valid JSON"
            continue
        problems = validate_record(record)
        if problems:
            report["incomplete_records"].append(stem)
            report["details"][stem] = "; ".join(problems)
            continue
        key = (record["smiles"], record["vina_seed"], record["prep_seed"])
        if job_id(*key) != stem or record["job_id"] != stem:
            report["id_mismatch"].append(stem)
            report["details"][stem] = "file name, record job_id and recomputed id disagree"
            continue
        if key in seen:
            report["duplicate_keys"].append(stem)
            continue
        seen[key] = stem
        report["records"][stem] = record
        if record["status"] != OK:
            report["failed"].append(stem)
            if (jobs / f"{stem}.sdf").exists():
                report["orphan_poses"].append(stem)          # a failed job must not keep a pose
                report["details"][stem] = "failed job has a pose file"
            continue
        problem = check_pose(jobs, stem, record)
        if problem is None:
            report["complete"].append(stem)
        else:
            report[problem[0]].append(stem)
            report["details"][stem] = problem[1]
    committed = {p.stem for p in jobs.glob("*.json")}
    for p in sorted(jobs.glob("*.sdf")):
        if p.stem not in committed and p.stem not in report["orphan_poses"]:
            report["orphan_poses"].append(p.stem)
            report["details"].setdefault(p.stem, "pose without a committed record")
    for p in sorted(jobs.iterdir()):
        if p.suffix not in (".json", ".sdf") and not p.name.endswith((".sdf.tmp", ".json.part")) and p.is_file():
            report["unattributable"].append(p.name)
            report["details"][p.name] = "unexpected file in the job directory"
    return report


def anomalies(report):
    """The reconcile findings that must never be silently accepted, as {class: ids}."""
    return {n: sorted(set(report[n])) for n in ANOMALY_CLASSES if report[n]}
