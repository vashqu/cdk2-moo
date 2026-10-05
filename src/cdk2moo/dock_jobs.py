"""
Durable, resumable docking jobs.

The first docking script (scripts/08d, results/08d_docking_scores.csv, data/structures/poses_seed*.sdf) kept scores in a CSV and poses
in an SDF file that it re-created on resume. Its saved artifacts were audited and are complete, so they are kept as historical "legacy"
results with no provenance manifest. This module is the scheme used for every new docking run; it never touches the legacy files.

Layout of a run directory:
  manifest.json          receptor hash, box, docking settings, ligand-preparation seed, versions
  maps/                  grid maps with their own hash manifest
  jobs/<job_id>.json     the committed record of one job (written LAST, atomically)
  jobs/<job_id>.sdf      the best pose of that job (written first, atomically; its SHA-256 is committed in the record)
  quarantine/<stamp>/    artifacts set aside (never deleted) when a job is repaired, with a reason file
  .lock                  held by the one process that is writing; a second writer is refused

A job is one molecule x one Vina seed x one ligand-preparation seed. Validation of records and poses is in dock_validate.py.

Repair and retry semantics (explicit):
  * run_jobs(jobs, repair=False): any inconsistency anywhere in the store stops the run, with a diagnostic.
  * repair=True repairs ONLY inconsistent jobs that belong to the requested `jobs`: their artifacts are moved to quarantine/ (kept as
    evidence) and the jobs are redone. An inconsistency belonging to a job that was NOT requested is never touched; it makes the run
    refuse, so a restricted request (for example --limit) cannot silently delete unrelated work.
  * Files whose job identity cannot be recovered are never attributed to a request; they are moved to quarantine only with
    repair_scope="store", which repairs every attributable job and quarantines unattributable files, and is meant to be run on its own.
  * A recorded failure (prep_failed, dock_failed, pose_export_failed) is a valid result. It is redone only with retry_failed=True,
    and only for requested jobs.
"""

import datetime
import fcntl
import json
import math
import os
import shutil
import time
from contextlib import contextmanager
from pathlib import Path

from rdkit import Chem

from cdk2moo.dock_validate import (OK, PREP_FAILED, DOCK_FAILED, POSE_FAILED, STATUSES, ANOMALY_CLASSES, job_id, sha256_of,
                                   reconcile, anomalies, validate_record, read_pose)


class ResumeRefused(Exception):
    """Raised when a run directory cannot be safely resumed or written; the message says why."""


def pose_is_readable(path):
    """True if the SDF file holds a first record that RDKit can parse (an empty or garbled file is not readable)."""
    mols, readable = read_pose(path)
    return readable and bool(mols) and mols[0] is not None


def atomic_write_text(path, text):
    """Write a file so that a reader never sees a half-written version."""
    path = Path(path)
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(text)
    os.replace(tmp, path)


@contextmanager
def store_lock(run_dir):
    """Exclusive writer lock on a run directory. A second writer raises ResumeRefused instead of waiting or interleaving."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    handle = open(run_dir / ".lock", "w")
    try:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        handle.close()
        raise ResumeRefused(f"{run_dir} is being written by another process (lock held); refusing a concurrent writer")
    try:
        yield
    finally:
        fcntl.flock(handle, fcntl.LOCK_UN)
        handle.close()


# ---- provenance ---------------------------------------------------------------------------------

def check_manifest(run_dir, current):
    """
    Create the manifest on first use; on resume require it to match `current` exactly.

    `current` holds everything that could change a result: receptor hash, box, Vina settings, ligand-preparation seed, package
    versions. A run directory with job files but no manifest is a legacy directory: refused, not adopted with invented provenance.
    """
    run_dir = Path(run_dir)
    manifest = run_dir / "manifest.json"
    if not manifest.exists():
        if (run_dir / "jobs").exists() and any((run_dir / "jobs").iterdir()):
            raise ResumeRefused(f"{run_dir} holds job files but no manifest.json (legacy artifacts): not resuming, no provenance "
                                "can be established. Use a new run directory.")
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "jobs").mkdir(exist_ok=True)
        atomic_write_text(manifest, json.dumps(current, indent=2, sort_keys=True))
        return "created"
    saved = json.loads(manifest.read_text())
    different = sorted(k for k in set(saved) | set(current) if saved.get(k) != current.get(k))
    if different:
        details = "; ".join(f"{k}: saved {saved.get(k)!r} vs now {current.get(k)!r}" for k in different)
        raise ResumeRefused(f"cannot resume {run_dir}: settings differ ({details})")
    return "matches"


# ---- one job --------------------------------------------------------------------------------------

def run_job(run_dir, smiles, vina_seed, prep_seed, prepare, dock_fn):
    """
    Prepare, dock and commit one job; returns its record. `prepare(smiles, prep_seed)` returns (pdbqt, mol_with_hydrogens,
    n_rotatable) or None; `dock_fn(pdbqt, vina_seed)` returns (score, pose_mol). Both are arguments so tests can supply stand-ins.
    The pose is written and verified first, its hash is committed in the record, and the record is written last.
    """
    jobs = Path(run_dir) / "jobs"
    jid = job_id(smiles, vina_seed, prep_seed)
    start = time.time()
    record = {"job_id": jid, "smiles": smiles, "vina_seed": int(vina_seed), "prep_seed": int(prep_seed), "status": OK,
              "vina_score": None, "n_heavy_atoms": None, "n_rotatable": None, "seconds": 0.0, "error": "", "pose_sha256": None}
    try:
        prepared = prepare(smiles, prep_seed)
    except Exception as exc:                              # preparation is allowed to fail per molecule
        prepared, record["error"] = None, f"{type(exc).__name__}: {exc}"
    if prepared is None:
        record["status"] = PREP_FAILED
    else:
        pdbqt, mol_h, n_rot = prepared
        record["n_heavy_atoms"], record["n_rotatable"] = int(mol_h.GetNumHeavyAtoms()), int(n_rot)
        pose_mol = None
        try:
            score, pose_mol = dock_fn(pdbqt, vina_seed)
            score = float(score)
        except Exception as exc:
            score, record["status"], record["error"] = float("nan"), DOCK_FAILED, f"{type(exc).__name__}: {exc}"
        if record["status"] == OK and not math.isfinite(score):
            record["status"], record["error"] = DOCK_FAILED, "non-finite score"
        if record["status"] == OK:
            record["vina_score"] = score
            record = _write_pose(jobs, jid, record, pose_mol)
    if record["status"] != OK:
        if record["status"] != POSE_FAILED:
            record["vina_score"] = None
        stale = jobs / f"{jid}.sdf"
        if stale.exists():
            stale.unlink()                                  # a failed redo must not leave an old pose behind
    record["seconds"] = float(time.time() - start)
    atomic_write_text(jobs / f"{jid}.json", json.dumps(record))   # record last: this commits the job
    return record


def _write_pose(jobs, jid, record, pose_mol):
    """Export the pose: one record, verified readable, hash committed. Export problems become status pose_export_failed."""
    if pose_mol is None:
        record["status"], record["error"] = POSE_FAILED, "no pose molecule returned"
        return record
    try:
        pose_mol.SetProp("smiles", record["smiles"])
        pose_mol.SetProp("vina_score", f"{record['vina_score']:.3f}")
        pose_mol.SetProp("job_id", jid)
        tmp = jobs / f"{jid}.sdf.tmp"
        writer = Chem.SDWriter(str(tmp))
        writer.write(pose_mol)
        writer.close()
        mols, readable = read_pose(tmp)
        if not readable or len(mols) != 1 or mols[0] is None:
            tmp.unlink(missing_ok=True)
            record["status"], record["error"] = POSE_FAILED, "pose could not be read back as exactly one record"
            return record
        record["pose_sha256"] = sha256_of(tmp)
        os.replace(tmp, jobs / f"{jid}.sdf")                # pose first ...
    except Exception as exc:
        record["status"], record["error"] = POSE_FAILED, f"{type(exc).__name__}: {exc}"
        record["pose_sha256"] = None
    return record


# ---- running with scoped, evidence-preserving repair ---------------------------------------------------

def quarantine(run_dir, names, reasons):
    """Move artifacts (job ids or file names) into quarantine/<stamp>/ with a reason file; nothing is deleted."""
    run_dir = Path(run_dir)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    target = run_dir / "quarantine" / stamp
    target.mkdir(parents=True, exist_ok=True)
    moved = {}
    for name in names:
        for suffix in ("", ".json", ".sdf"):
            path = run_dir / "jobs" / f"{name}{suffix}"
            if path.is_file():
                shutil.move(str(path), str(target / path.name))
                moved.setdefault(name, []).append(path.name)
    (target / "reasons.json").write_text(json.dumps({"reasons": {n: reasons.get(n, "") for n in names}, "moved": moved}, indent=2))
    return target


def run_jobs(run_dir, jobs, prepare, dock_fn, prep_seed, retry_failed=False, repair=False, repair_scope="requested", progress=None):
    """
    Run `jobs` (a list of (smiles, vina_seed)) in a run directory, resuming safely (see the module text for the semantics).

    Returns (records, repaired) where `records` are all valid records in the store and `repaired` maps anomaly class -> ids that were
    quarantined and redone. Raises ResumeRefused for any unresolved inconsistency.
    """
    if repair_scope not in ("requested", "store"):
        raise ValueError("repair_scope must be 'requested' or 'store'")
    run_dir = Path(run_dir)
    requested = {job_id(s, v, prep_seed) for s, v in jobs}
    with store_lock(run_dir):
        (run_dir / "jobs").mkdir(exist_ok=True)
        for pattern in ("*.tmp", "*.part"):               # unfinished writes by a crashed writer: never a result
            for stale in (run_dir / "jobs").glob(pattern):
                stale.unlink()
        report = reconcile(run_dir)
        in_scope, external = {}, {}
        for cls, ids in anomalies(report).items():
            for jid in ids:
                owned = jid in requested or (repair_scope == "store")
                (in_scope if owned else external).setdefault(cls, []).append(jid)
        if external:
            raise ResumeRefused("inconsistent jobs outside the requested scope were left untouched, and the run is refused: "
                                + "; ".join(f"{c}={v[:3]}{'...' if len(v) > 3 else ''}" for c, v in external.items())
                                + ". Repair them in a run that requests them, or use repair_scope='store'.")
        if in_scope and not repair:
            raise ResumeRefused("inconsistent requested jobs, not resuming: "
                                + "; ".join(f"{c}={v[:3]}{'...' if len(v) > 3 else ''}" for c, v in in_scope.items())
                                + ". Pass repair=True to quarantine and redo them.")
        repaired = {}
        if in_scope:
            names = sorted({n for ids in in_scope.values() for n in ids})
            quarantine(run_dir, names, report["details"])      # evidence kept before anything is replaced
            repaired = in_scope
        report = reconcile(run_dir)
        redo_failed = set(report["failed"]) & requested if retry_failed else set()
        done = (set(report["complete"]) | set(report["failed"])) - redo_failed
        todo = [(s, v) for s, v in jobs if job_id(s, v, prep_seed) not in done]
        if redo_failed:
            quarantine(run_dir, sorted(redo_failed), {j: "recorded failure retried on request" for j in redo_failed})
        for n, (smiles, vina_seed) in enumerate(todo, 1):
            run_job(run_dir, smiles, vina_seed, prep_seed, prepare, dock_fn)
            if progress:
                progress(n, len(todo))
        return reconcile(run_dir)["records"], repaired


def print_progress(done, total, every=25):
    """Progress callback for run_jobs: one line every `every` jobs."""
    if done % every == 0 or done == total:
        print(f"  {done}/{total} jobs done", flush=True)


# ---- export for downstream scripts -------------------------------------------------------------------

def export(run_dir, out_dir, only_smiles=None):
    """
    Write legacy-compatible files from a reconciled run: `docking_scores.csv` (the columns of results/08d_docking_scores.csv plus
    job_id, prep_seed and pose_sha256) and `poses_seed<N>.sdf` (same `smiles` and `vina_score` properties). With `only_smiles`, only
    those molecules are exported (a scope's own sets); jobs that failed stay in the table with their status, so nothing disappears
    from a denominator. Refuses to export an inconsistent store.
    """
    import pandas as pd
    report = reconcile(run_dir)
    found = anomalies(report)
    if found:
        raise ResumeRefused(f"cannot export an inconsistent run directory: {sorted(found)}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    keep = {j: r for j, r in report["records"].items() if only_smiles is None or r["smiles"] in only_smiles}
    cols = ["smiles", "vina_seed", "vina_score", "n_heavy_atoms", "n_rotatable", "seconds", "status", "job_id", "prep_seed", "pose_sha256"]
    table = pd.DataFrame([{k: r[k] for k in cols} for r in keep.values()], columns=cols)
    table = table.sort_values(["vina_seed", "smiles"]).reset_index(drop=True)
    table.to_csv(out_dir / "docking_scores.csv", index=False)
    for seed in sorted(table["vina_seed"].unique()):
        ids = sorted(j for j in report["complete"] if j in keep and keep[j]["vina_seed"] == seed)
        text = "".join((Path(run_dir) / "jobs" / f"{j}.sdf").read_text() for j in ids)
        atomic_write_text(out_dir / f"poses_seed{seed}.sdf", text)
    return table
