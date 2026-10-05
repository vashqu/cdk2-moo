"""
Where a run writes its results, and a small manifest saying exactly what produced them.

Rules (so a smoke test, a relaxed-window run or a corrected-policy run can never replace a primary result):
  * The primary file name is RESERVED for the primary protocol: legacy policy, the primary size window, the
    default arms, seeds and generations. A run with any other setting that reaches the primary name raises.
  * With no explicit --out, a non-primary run gets a derived name that spells out what differs, e.g.
    05_ga_relaxed_populations.csv or 05_ga_seeds-42_gen-3_populations.csv.
  * An explicit --out is honoured, but never silently replaces an existing file; pass overwrite=True.
"""

import datetime
import importlib.metadata as metadata
import json
import os
import platform
import subprocess
from pathlib import Path

PACKAGES = ["numpy", "pandas", "scikit-learn", "scipy", "rdkit", "vina", "meeko", "posebusters", "gemmi", "matplotlib"]


def choose_output(results_dir, primary_name, out=None, tokens=(), overwrite=False):
    """
    Resolve the output path for a run.

    `tokens` lists what makes this run non-primary (empty for the primary protocol), e.g.
    ["relaxed", "seeds-42", "gen-3"]. Returns a Path. Raises ValueError if a non-primary run would take the
    reserved primary name, and FileExistsError if the file exists and overwrite is False.
    """
    results_dir = Path(results_dir).resolve()
    if out is not None:
        path = Path(out)
        path = path if path.is_absolute() else results_dir / path
        path = path.resolve()
    elif not tokens:
        path = results_dir / primary_name
    else:
        stem, suffix = primary_name.rsplit(".", 1)
        root = stem[:-len("_populations")] if stem.endswith("_populations") else stem
        ending = "_populations" if stem.endswith("_populations") else ""
        path = results_dir / f"{root}_{'_'.join(tokens)}{ending}.{suffix}"
    path = path.resolve()                                    # symlinks and relative spellings cannot dodge the reserved name
    if tokens and path.name == primary_name and path.parent == results_dir:
        raise ValueError(f"{primary_name} is reserved for the primary protocol, but this run differs from it ({', '.join(tokens)}). "
                         "Choose another --out.")
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists; use another --out or pass --overwrite.")
    return path


def run_tokens(relaxed=False, policy="legacy", arms=None, default_arms=None, seeds=None, default_seeds=None,
               generations=None, default_generations=None, extra=()):
    """The list of setting differences from the primary protocol, for naming outputs."""
    tokens = []
    if relaxed:
        tokens.append("relaxed")
    if policy != "legacy":
        tokens.append(policy)
    if arms is not None and list(arms) != list(default_arms):
        tokens.append("arms-" + "+".join(arms))
    if seeds is not None and list(seeds) != list(default_seeds):
        tokens.append("seeds-" + "+".join(str(s) for s in seeds))
    if generations is not None and generations != default_generations:
        tokens.append(f"gen-{generations}")
    return tokens + list(extra)


def commit_csv(frame, path, **to_csv_kwargs):
    """
    Write a table so that a partial file is never mistaken for a finished one: write <path>.part, then rename. A reader should only trust a
    result whose <path>.manifest.json exists (written afterwards by write_manifest); a lone .part file is an unfinished run.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    frame.to_csv(part, **to_csv_kwargs)
    os.replace(part, path)
    return path


def _version(name):
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def write_manifest(output_path, settings, project_root=None, inputs=None, extra_outputs=None):
    """
    Write <output>.manifest.json AFTER the output exists: the output's hash, the hashes of the inputs it was computed from, the effective
    settings and parameters, the identity of the source code (a content hash of src/, scripts/ and tests/, the Git commit and the working-tree
    status; a commit alone is not enough because the tree is modified and partly untracked), package versions and platform. The presence of
    the manifest marks the run as complete. It describes the run that wrote that file, never any other run.
    """
    from cdk2moo import campaign, config
    output_path = Path(output_path)
    root = Path(project_root or config.PROJECT_ROOT)
    outputs = {output_path.name: {"sha256": campaign.sha256_file(output_path), "bytes": output_path.stat().st_size}}
    for extra in extra_outputs or []:
        outputs[Path(extra).name] = {"sha256": campaign.sha256_file(extra), "bytes": Path(extra).stat().st_size}
    identity = campaign.source_identity(root)
    manifest = {
        "complete": True, "output": output_path.name, "outputs": outputs, "created_utc": campaign.now_utc(),
        "campaign": config.CAMPAIGN, "scope": config.SCOPE if config.CAMPAIGN not in (None, "historical") else None,
        "settings": settings, "parameters": config.effective_parameters(),
        "inputs": {name: ({"path": str(p), "sha256": campaign.sha256_file(p), "bytes": Path(p).stat().st_size} if Path(p).is_file() else None)
                   for name, p in (inputs or {}).items()},
        "source": {k: identity[k] for k in ("source_tree_sha256", "n_source_files", "git_head", "working_tree_clean", "git_status_lines")},
        "versions": {name: _version(name) for name in PACKAGES}, "python": platform.python_version(), "platform": platform.platform(),
    }
    path = output_path.with_name(output_path.name + ".manifest.json")
    campaign.write_json_atomic(path, manifest)
    return path
