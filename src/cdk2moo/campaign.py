"""
Campaign bookkeeping: hashing, source identity, manifests, the execution ledger, and the historical inventory.

A CAMPAIGN is one complete, clearly identified execution of the repaired pipeline, kept in
campaigns/<campaign_id>/ and never mixed with the historical results in data/, results/ and figures/.

    campaigns/<id>/inputs/    frozen copies of the raw inputs, with hashes (cached ChEMBL file, 4KD1 structure, crystal ligand)
    campaigns/<id>/shared/    policy-independent stages: curation, fingerprints, splits, surrogate evaluation, holdout definition,
                              receptor preparation and redocking
    campaigns/<id>/legacy/    GA-dependent stages for the legacy candidate-preparation policy
    campaigns/<id>/corrected/ the same for the corrected policy (kept separate; the two are never pooled)
    campaigns/<id>/docking_store/   one durable docking job store shared by every scope (identical jobs are computed once)
    campaigns/<id>/manifests/ per-stage manifests, `ledger.jsonl` (planned/running/completed/failed), plan and source identity

Every scope has its own processed/, results/, figures/ and structures/ directories (see config.py). A script run with no campaign
selected refuses to start; historical reanalysis must be selected explicitly with CDK2_CAMPAIGN=historical.
"""

import datetime
import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def now_utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def write_json_atomic(path, obj):
    """Write JSON so that a reader never sees a half-written file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str))
    os.replace(tmp, path)


def hash_files(paths, root=None):
    """{relative path: {sha256, bytes}} for existing files; a missing file is recorded as null, never skipped."""
    out = {}
    for p in sorted(Path(x) for x in paths):
        key = str(p.relative_to(root)) if root and p.is_absolute() and str(p).startswith(str(root)) else str(p)
        out[key] = {"sha256": sha256_file(p), "bytes": p.stat().st_size} if p.is_file() else None
    return out


def source_identity(project_root):
    """
    What source code produced a run, independent of Git: a hash over the content of every file in src/, scripts/ and tests/,
    plus the Git commit and the working-tree status (which lists modified and untracked files). A commit alone does not identify the
    code when files are untracked or modified, which is the normal state of this repository.
    """
    project_root = Path(project_root)
    files = sorted(p for d in ("src", "scripts", "tests") for p in (project_root / d).rglob("*.py"))
    per_file = {str(p.relative_to(project_root)): sha256_file(p) for p in files}
    tree = hashlib.sha256(json.dumps(per_file, sort_keys=True).encode()).hexdigest()

    def git(*args):
        try:
            return subprocess.run(["git", *args], cwd=project_root, capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception:
            return ""
    status = git("status", "--porcelain").splitlines()
    return {"source_tree_sha256": tree, "n_source_files": len(per_file), "files": per_file,
            "git_head": git("rev-parse", "HEAD") or None, "git_status_lines": status,
            "working_tree_clean": not status}


def environment_identity():
    """What this machine and environment looked like when the campaign ran (a NEW capture, not evidence about older runs)."""
    import importlib.metadata as metadata
    packages = ["numpy", "pandas", "scikit-learn", "scipy", "rdkit", "vina", "meeko", "posebusters", "gemmi",
                "matplotlib", "openmm", "pdbfixer"]
    versions = {}
    for name in packages:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return {"captured_utc": now_utc(), "python": platform.python_version(), "platform": platform.platform(),
            "machine": platform.machine(), "cpu_count": os.cpu_count(), "versions": versions}


def ledger_append(ledger_path, **event):
    """Append one event (stage, condition, status, ...) to the campaign ledger, one JSON object per line."""
    ledger_path = Path(ledger_path)
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    with open(ledger_path, "a") as fh:
        fh.write(json.dumps({"time_utc": now_utc(), **event}, default=str) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def ledger_read(ledger_path):
    ledger_path = Path(ledger_path)
    if not ledger_path.exists():
        return []
    return [json.loads(line) for line in ledger_path.read_text().splitlines() if line.strip()]


def historical_inventory(project_root, cutoff="2026-10-03T00:00:00"):
    """
    Hash every historical artifact (data, results, figures, structures, decisions, documentation) so that preservation can
    be verified later. Files are tagged "historical" if last modified before `cutoff` (the start of the repair pass) and
    "repair_pass" if modified at or after it (new files written by that pass); both are preserved. Campaign directories
    are excluded.
    """
    project_root = Path(project_root)
    cutoff_ts = datetime.datetime.fromisoformat(cutoff).timestamp()
    roots = ["data", "results", "figures", "reproducibility", "decisions.md", "CLAUDE.md", "README.md", "environment.yml"]
    inventory = {}
    for item in roots:
        base = project_root / item
        files = [base] if base.is_file() else sorted(p for p in base.rglob("*") if p.is_file())
        for p in files:
            rel = str(p.relative_to(project_root))
            inventory[rel] = {"sha256": sha256_file(p), "bytes": p.stat().st_size,
                              "modified_utc": datetime.datetime.fromtimestamp(p.stat().st_mtime, datetime.timezone.utc).isoformat(timespec="seconds"),
                              "tag": "historical" if p.stat().st_mtime < cutoff_ts else "repair_pass"}
    return inventory


def verify_inventory(project_root, inventory):
    """Compare current files with a recorded inventory. Returns {changed, missing, added} lists of relative paths."""
    project_root = Path(project_root)
    changed, missing = [], []
    for rel, info in inventory.items():
        p = project_root / rel
        if not p.exists():
            missing.append(rel)
        elif sha256_file(p) != info["sha256"]:
            changed.append(rel)
    current = {str(p.relative_to(project_root)) for item in ("data", "results", "figures", "reproducibility")
               for p in (project_root / item).rglob("*") if p.is_file()}
    return {"changed": changed, "missing": missing, "added": sorted(current - set(inventory))}
