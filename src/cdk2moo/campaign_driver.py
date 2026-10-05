"""
Running a campaign: the stage table, the ledger, per-stage manifests, shared-file linking, and a small parallel scheduler.

The campaign mechanism (config.py) makes every pipeline script read and write inside one scope of one campaign. This module decides which
scripts run in which scope, in what order, and records exactly what happened:

  * STAGES (below) lists every stage: its scope, script, arguments, dependencies, declared outputs and the CPU cores it needs.
  * A stage runs only when every dependency has a COMPLETE manifest whose recorded output hashes still match the files on disk. A script is
    never launched on inputs that changed after their producer finished.
  * Each stage writes manifests/stages/<scope>__<stage>.json: command, exit code, duration, the hashes of the upstream outputs it consumed,
    the hashes of its own declared outputs, the source-tree identity, and complete=true. The ledger (manifests/ledger.jsonl) records
    planned / running / completed / failed events with timestamps.
  * Files shared between scopes (curated data, fingerprints, splits, holdout assignments, receptor) are COPIED into each policy scope by a
    "link" stage with their hashes checked, and the core inputs are validated for row/identity alignment before anything consumes them.
  * Existing outputs of a stage that is not marked complete are moved to superseded/ (never deleted, never silently reused).
"""

import fnmatch
import glob
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from cdk2moo import campaign, config

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
POLICY_SCOPES = ("legacy", "corrected")
ROBUSTNESS_CONTROLS = ["scaffold_split", "scaffold_split_s43", "scaffold_split_s44", "scaffold_split_s45", "scaffold_split_s46"]

# files copied from the shared scope into each policy scope (relative to a scope root)
LINKED_FILES = ["processed/cdk2_ic50_curated.csv", "processed/cdk2_kikd_external.csv", "processed/ecfp4.npy", "processed/cdk2_splits.csv",
                "processed/eligible_starts.csv", "results/02_curation_funnel.json", "results/03_splits.json", "results/04_surrogate_metrics.csv",
                "results/04_test_predictions.csv", "results/05e_holdout_sets.csv", "results/05e_heldout_predictions.csv",
                "results/05e_setup_report.json", "results/08b_redocking.json", "structures/receptor.pdbqt", "structures/receptor_H.pdb",
                "structures/receptor_box.txt", "structures/4KD1.pdb", "structures/4kd1_1QK_crystal.sdf"]


def stage(sid, scope, script, args=(), deps=(), outputs=(), cores=1, note=""):
    return {"id": sid, "scope": scope, "script": script, "args": list(args), "deps": list(deps), "outputs": list(outputs), "cores": cores, "note": note}


def build_stages():
    """The full stage table. Dependencies are (scope, id) pairs written as 'scope/id'."""
    st = [
        stage("curate", "shared", "02_curate.py", outputs=["processed/cdk2_ic50_curated.csv", "processed/cdk2_kikd_external.csv", "results/02_curation_funnel.json"]),
        stage("features_splits", "shared", "03_features_and_splits.py", deps=["shared/curate"],
              outputs=["processed/ecfp4.npy", "processed/cdk2_splits.csv", "results/03_splits.json"]),
        stage("eligible_starts", "shared", "03b_eligible_starts.py", deps=["shared/curate"], outputs=["processed/eligible_starts.csv"]),
        stage("surrogate_eval", "shared", "04_train_surrogate.py", deps=["shared/features_splits"],
              outputs=["results/04_surrogate_metrics.csv", "results/04_test_predictions.csv"]),
        stage("holdout_setup", "shared", "05e_positive_control_setup.py", deps=["shared/features_splits"],
              outputs=["results/05e_holdout_sets.csv", "results/05e_heldout_predictions.csv", "results/05e_setup_report.json"]),
        stage("receptor", "shared", "08_prepare_receptor.py", outputs=["structures/receptor.pdbqt", "structures/receptor_H.pdb", "structures/receptor_box.txt"]),
        stage("redock", "shared", "08b_redock_native_ligand.py", deps=["shared/receptor"], outputs=["results/08b_redocking.json"]),
    ]
    shared_done = [f"shared/{s}" for s in ("curate", "features_splits", "eligible_starts", "surrogate_eval", "holdout_setup", "receptor", "redock")]
    previous_dock = None
    for scope in POLICY_SCOPES:
        P = lambda sid: f"{scope}/{sid}"                         # noqa: E731  (a tiny local name helper)
        st.append(stage("link", scope, "@link", deps=shared_done, note="copy shared inputs into this scope with hash checks and alignment validation"))
        st += [
            stage("ga_primary", scope, "05_run_ga.py", deps=[P("link")], outputs=["results/05_ga_populations.csv"]),
            stage("size_percentile", scope, "05b_size_conditioned_similarity.py", deps=[P("ga_primary")],
                  outputs=["results/05_ga_populations_scored.csv", "results/05_training_loo_similarity.csv"]),
            stage("summarize_ga", scope, "05c_summarize_ga.py", deps=[P("size_percentile")], outputs=["results/05_ga_per_run.csv", "results/05_ga_final_summary.csv"]),
            stage("descriptive", scope, "05d_ga_descriptive_checks.py", deps=[P("size_percentile")], outputs=["@log"]),
            stage("holdout_molecule", scope, "05f_positive_control_ga.py", ["--controls", "molecule"], [P("link")], ["results/05f_molecule_populations.csv"]),
            stage("holdout_scaffold", scope, "05f_positive_control_ga.py", ["--controls", "scaffold"], [P("link")], ["results/05f_scaffold_populations.csv"]),
        ]
        for control in ROBUSTNESS_CONTROLS:
            st.append(stage(f"robust_{control}", scope, "05f_positive_control_ga.py", ["--controls", control], [P("link")],
                            [f"results/05f_{control}_populations.csv"]))
        st += [
            stage("robustness_eval", scope, "05h_scaffold_surrogate_robustness.py", deps=[P("ga_primary")] + [P(f"robust_{c}") for c in ROBUSTNESS_CONTROLS],
                  outputs=["results/05h_robustness_summary.csv"]),
            stage("holdout_eval", scope, "05g_positive_control_eval.py", deps=[P("holdout_molecule"), P("holdout_scaffold")],
                  outputs=["results/05g_positive_control_summary_v2.csv"]),
            stage("pareto", scope, "06_pareto.py", deps=[P("size_percentile")], outputs=["results/06_pareto_front.csv"]),
            stage("ad_reliability", scope, "07_ad_reliability.py", deps=[P("link")], outputs=["results/07_stage4_ad.csv", "results/07_reliability_bins.csv"]),
            stage("ad_overlay", scope, "07b_ad_generated_overlay.py", deps=[P("size_percentile"), P("pareto"), P("ad_reliability")],
                  outputs=["results/07_generated_regions.csv"]),
            stage("alerts", scope, "07c_structural_alerts.py", deps=[P("ad_overlay")], outputs=["results/07_generated_alerts.csv", "results/07_front_alerts.csv"]),
        ]
        for experiment in ("primary", "molecule", "scaffold"):
            st.append(stage(f"constrained_{experiment}", scope, "10_run_constrained_ga.py", ["--experiments", experiment], [P("link")],
                            [f"results/10_constrained_{experiment}_populations.csv"]))
        st += [
            stage("tradeoff", scope, "10b_analyze_tradeoff.py", deps=[P("ga_primary"), P("holdout_molecule"), P("holdout_scaffold")]
                  + [P(f"constrained_{e}") for e in ("primary", "molecule", "scaffold")], outputs=["results/10_tradeoff_v2.csv"]),
            stage("endpoints", scope, "13_hypothesis_endpoints.py", deps=[P("size_percentile"), P("holdout_molecule"), P("holdout_scaffold")]
                  + [P(f"constrained_{e}") for e in ("primary", "molecule", "scaffold")],
                  outputs=["results/13_h1_summary.csv", "results/13_h2_summary.csv", "results/13_h4_summary.csv", "results/13_h5_summary.csv"]),
            stage("dock_sets", scope, "08c_select_docking_sets.py", deps=[P("size_percentile"), P("pareto"), P("holdout_molecule"), P("holdout_scaffold")],
                  outputs=["results/08c_docking_sets.csv", "results/08c_selection_report.json"]),
            stage("dock", scope, "08d_dock_sets.py", deps=[P("dock_sets")] + ([previous_dock] if previous_dock else []),
                  outputs=["results/08d_export/docking_scores.csv"], cores=8, note="all docking stages need the full machine; the store is shared and has one writer"),
            stage("dock_analysis", scope, "08e_docking_analysis.py", deps=[P("dock"), P("posebusters")], outputs=["results/08e_docking_by_molecule.csv"]),
            stage("posebusters", scope, "08f_posebusters.py", deps=[P("dock")], outputs=["results/posebusters_audit/posebusters_checks.csv"]),
            stage("hinge", scope, "08g_hinge_contacts.py", deps=[P("dock")], outputs=["results/08g_hinge_contacts.csv", "results/08g_hinge_by_set.csv"]),
            stage("figures", scope, "11_make_figures.py", deps=[P("dock_analysis"), P("posebusters"), P("hinge"), P("summarize_ga"), P("pareto"), P("ad_reliability")],
                  outputs=["figures/final/fig1_surrogate.png", "figures/final/fig2_trajectories.png", "figures/final/fig3_pareto_fronts.png", "figures/final/fig3b_pareto_structures.png", "figures/final/fig4_similarity_vs_prediction.png",
                           "figures/final/fig5_docking_P.png", "figures/final/fig5_docking_M.png", "figures/final/fig5_docking_S.png", "figures/final/fig6_pose_checks.png"]),
            stage("table1", scope, "11b_table1.py", deps=[P("dock_analysis"), P("posebusters"), P("hinge")], outputs=["results/table1.csv", "figures/final/table1.md"]),
            stage("audit", scope, "12_audit_leakage_and_artifacts.py", deps=[P("dock"), P("posebusters"), P("hinge")], outputs=["results/audit/leakage_summary.csv"]),
        ]
        previous_dock = P("dock")
    return st


STAGES = build_stages()
BY_KEY = {f"{s['scope']}/{s['id']}": s for s in STAGES}


# ---- paths ---------------------------------------------------------------------------------------------

def campaign_root(project_root, cid):
    return Path(project_root) / "campaigns" / cid


def scope_root(project_root, cid, scope):
    return campaign_root(project_root, cid) / scope


def stage_manifest_path(project_root, cid, key):
    scope, sid = key.split("/")
    return campaign_root(project_root, cid) / "manifests" / "stages" / f"{scope}__{sid}.json"


def ledger_path(project_root, cid):
    return campaign_root(project_root, cid) / "manifests" / "ledger.jsonl"


# ---- initialisation ---------------------------------------------------------------------------------------

def init_campaign(project_root, cid, historical_record_file):
    """
    Create the campaign tree, copy and hash the FROZEN inputs (the cached ChEMBL download and the original structural inputs from the historical
    data directory), capture the source identity and a new environment snapshot, and record every stage as planned. Idempotent; refuses to
    change an input already frozen.
    """
    project_root = Path(project_root)
    root = campaign_root(project_root, cid)
    for sub in ("inputs/raw", "inputs/structures", "manifests/stages", "manifests/logs", "shared", "legacy", "corrected", "docking_store", "plan"):
        (root / sub).mkdir(parents=True, exist_ok=True)
    frozen = {"inputs/raw/chembl_CHEMBL301_activities.csv": project_root / "data/raw/chembl_CHEMBL301_activities.csv",
              "inputs/structures/4KD1.pdb": project_root / "data/structures/4KD1.pdb",
              "inputs/structures/4kd1_1QK_crystal.sdf": project_root / "data/structures/4kd1_1QK_crystal.sdf"}
    inputs = {}
    for rel, source in frozen.items():
        target = root / rel
        if not source.exists():
            raise FileNotFoundError(f"frozen input missing: {source}")
        if target.exists():
            if campaign.sha256_file(target) != campaign.sha256_file(source):
                raise RuntimeError(f"{target} differs from its source; frozen inputs are never changed")
        else:
            shutil.copy2(source, target)
        inputs[rel] = {"sha256": campaign.sha256_file(target), "bytes": target.stat().st_size, "copied_from": str(source.relative_to(project_root))}
    historical = json.loads(Path(historical_record_file).read_text())["files"]
    for rel, info in inputs.items():
        original = "data/" + rel.split("/", 1)[1].replace("raw/", "raw/").replace("structures/", "structures/")
        recorded = historical.get(original)
        info["matches_historical_record"] = bool(recorded and recorded["sha256"] == info["sha256"])
    note = {"inputs": inputs,
            "provenance": "The ChEMBL activity table is the CACHED download made in the historical run (5,367 records for CHEMBL301); no newer ChEMBL "
                          "release was downloaded. 4KD1.pdb and the crystal ligand SDF were downloaded from RCSB in the historical run and are the "
                          "structural inputs; the campaign never downloads in their place."}
    if not (root / "manifests/inputs.json").exists():
        campaign.write_json_atomic(root / "manifests/inputs.json", note)
    if not (root / "manifests/source_identity_at_init.json").exists():
        campaign.write_json_atomic(root / "manifests/source_identity_at_init.json", campaign.source_identity(project_root))
        campaign.write_json_atomic(root / "manifests/environment_at_init.json", campaign.environment_identity())
    done = {e["stage"] for e in campaign.ledger_read(ledger_path(project_root, cid)) if e.get("status") == "planned"}
    for s in STAGES:
        key = f"{s['scope']}/{s['id']}"
        if key not in done:
            campaign.ledger_append(ledger_path(project_root, cid), stage=key, status="planned", script=s["script"], args=s["args"], deps=s["deps"])
    return inputs


# ---- linking and validation ---------------------------------------------------------------------------------------

def validate_core_inputs(processed_dir, results_dir):
    """
    Check at a stage boundary that the curated molecules, fingerprints, splits, eligible-start list and holdout assignments describe the same
    molecules in the same order. Raises ValueError with the first problem found.
    """
    from cdk2moo.audit import check_alignment
    from cdk2moo.features import ecfp4
    processed_dir, results_dir = Path(processed_dir), Path(results_dir)
    curated = pd.read_csv(processed_dir / "cdk2_ic50_curated.csv")
    fps = np.load(processed_dir / "ecfp4.npy")
    if len(fps) != len(curated):
        raise ValueError(f"{len(fps)} fingerprints for {len(curated)} curated molecules")
    sample = np.random.default_rng(0).choice(len(curated), min(25, len(curated)), replace=False)
    if not (ecfp4(curated["std_smiles"].iloc[sample].tolist(), config.FP_RADIUS, config.FP_BITS) == fps[sample]).all():
        raise ValueError("stored fingerprints do not match the curated structures at the sampled rows")
    splits = pd.read_csv(processed_dir / "cdk2_splits.csv")
    for seed in sorted(splits["seed"].unique()):
        check_alignment(splits[splits["seed"] == seed], curated)
    for name in ("eligible_starts.csv",):
        table = pd.read_csv(processed_dir / name)
        if len(table) != len(curated) or not (table["inchikey"].to_numpy() == curated["inchikey"].to_numpy()).all():
            raise ValueError(f"{name} does not line up with the curated molecules")
    roles = pd.read_csv(results_dir / "05e_holdout_sets.csv")
    if len(roles) != len(curated) or not (roles["inchikey"].to_numpy() == curated["inchikey"].to_numpy()).all():
        raise ValueError("05e_holdout_sets.csv does not line up with the curated molecules")
    return {"curated": len(curated), "fingerprint_bits": int(fps.shape[1]), "split_seeds": int(splits["seed"].nunique())}


def link_shared(project_root, cid, scope):
    """Copy the shared inputs into a policy scope (hash-checked), validate alignment, and write link_manifest.json in the scope."""
    shared, target = scope_root(project_root, cid, "shared"), scope_root(project_root, cid, scope)
    linked = {}
    for rel in LINKED_FILES:
        source, dest = shared / rel, target / rel
        if not source.exists():
            raise FileNotFoundError(f"shared input {rel} does not exist; run the shared stages first")
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            if campaign.sha256_file(dest) != campaign.sha256_file(source):
                raise RuntimeError(f"{dest} differs from the shared file; linked inputs are never modified")
        else:
            shutil.copy2(source, dest)
        linked[rel] = {"sha256": campaign.sha256_file(dest), "bytes": dest.stat().st_size}
    checks = validate_core_inputs(target / "processed", target / "results")
    campaign.write_json_atomic(target / "link_manifest.json", {"linked_from": str(shared), "files": linked, "validated": checks,
                                                              "linked_utc": campaign.now_utc()})
    return checks


# ---- running one stage ---------------------------------------------------------------------------------------------

def _declared_outputs(project_root, cid, s):
    root = scope_root(project_root, cid, s["scope"])
    found, missing = [], []
    for pattern in s["outputs"]:
        if pattern == "@log":
            continue
        matches = sorted(Path(p) for p in glob.glob(str(root / pattern)))
        (found.extend(matches) if matches else missing.append(pattern))
    return found, missing


def stage_complete(project_root, cid, key):
    """True if the stage's manifest says complete and every recorded output still has the recorded hash."""
    path = stage_manifest_path(project_root, cid, key)
    if not path.exists():
        return False
    manifest = json.loads(path.read_text())
    if not manifest.get("complete"):
        return False
    for rel, info in manifest["outputs"].items():
        file = Path(project_root) / rel
        if not file.is_file() or campaign.sha256_file(file) != info["sha256"]:
            return False
    return True


def run_stage(project_root, cid, key, log_dir):
    """Run one stage as a subprocess in its scope. Returns (exit_code, seconds). Writes the stage manifest only on success."""
    project_root = Path(project_root)
    s = BY_KEY[key]
    ledger = ledger_path(project_root, cid)
    root = scope_root(project_root, cid, s["scope"])
    started = time.time()
    # outputs of an unfinished earlier attempt are set aside as evidence, never reused or deleted
    stale, _ = _declared_outputs(project_root, cid, s)
    if stale:
        stamp = time.strftime("%Y%m%dT%H%M%S")
        for f in stale:
            dest = campaign_root(project_root, cid) / "superseded" / stamp / s["scope"] / f.relative_to(root)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(f), str(dest))
    log_path = Path(log_dir) / f"{s['scope']}__{s['id']}.log"
    campaign.ledger_append(ledger, stage=key, status="running", log=str(log_path.relative_to(project_root)))
    if s["script"] == "@link":
        try:
            checks = link_shared(project_root, cid, s["scope"])
            log_path.write_text(json.dumps(checks))
            code = 0
        except Exception as exc:
            log_path.write_text(f"{type(exc).__name__}: {exc}")
            code = 1
    else:
        env = {**os.environ, "CDK2_CAMPAIGN": cid, "CDK2_SCOPE": s["scope"], "CDK2_PROJECT_ROOT": str(project_root), "PYTHONUNBUFFERED": "1"}
        command = [sys.executable, str(SCRIPTS / s["script"]), *s["args"]]
        with open(log_path, "w") as fh:
            fh.write("$ " + " ".join(command) + f"\n# CDK2_CAMPAIGN={cid} CDK2_SCOPE={s['scope']}\n")
            fh.flush()
            code = subprocess.run(command, env=env, stdout=fh, stderr=subprocess.STDOUT, cwd=project_root).returncode
    seconds = time.time() - started
    outputs, missing = _declared_outputs(project_root, cid, s)
    if code == 0 and missing:
        code = 2
        with open(log_path, "a") as fh:
            fh.write(f"\n# declared outputs missing after a clean exit: {missing}\n")
    if code != 0:
        tail = log_path.read_text().splitlines()[-6:]
        campaign.ledger_append(ledger, stage=key, status="failed", exit_code=code, seconds=round(seconds, 1), log_tail=tail)
        return code, seconds
    upstream = {dep: json.loads(stage_manifest_path(project_root, cid, dep).read_text())["outputs"] for dep in s["deps"]
                if stage_manifest_path(project_root, cid, dep).exists()}
    identity = campaign.source_identity(project_root)
    manifest = {"complete": True, "stage": key, "script": s["script"], "args": s["args"], "exit_code": 0, "seconds": round(seconds, 1),
                "finished_utc": campaign.now_utc(), "log": str(log_path.relative_to(project_root)), "campaign": cid, "scope": s["scope"],
                "outputs": campaign.hash_files(outputs, project_root), "log_sha256": campaign.sha256_file(log_path),
                "consumed_upstream_outputs": upstream,
                "source": {k: identity[k] for k in ("source_tree_sha256", "git_head", "working_tree_clean", "git_status_lines")}}
    campaign.write_json_atomic(stage_manifest_path(project_root, cid, key), manifest)
    campaign.ledger_append(ledger, stage=key, status="completed", seconds=round(seconds, 1), n_outputs=len(outputs))
    return 0, seconds


# ---- the scheduler -------------------------------------------------------------------------------------------------

def run_campaign(project_root, cid, cores=8, only=None, poll=5.0):
    """
    Run every stage (or those matching the `only` patterns, e.g. ['shared/*', 'legacy/ga_*']) whose dependencies are complete, keeping the
    total cores of running stages within `cores`. Stops scheduling after a failure but lets running stages finish. Returns the failed keys.
    """
    project_root = Path(project_root)
    log_dir = campaign_root(project_root, cid) / "manifests" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    wanted = [s for s in STAGES if only is None or any(fnmatch.fnmatch(f"{s['scope']}/{s['id']}", p) for p in only)]
    keys = [f"{s['scope']}/{s['id']}" for s in wanted]
    done = {k for k in keys if stage_complete(project_root, cid, k)}
    external_done = lambda dep: stage_complete(project_root, cid, dep)           # noqa: E731
    running, failed = {}, []
    from concurrent.futures import ThreadPoolExecutor
    pool = ThreadPoolExecutor(max_workers=max(1, cores))
    try:
        while True:
            for key in [k for k, f in running.items() if f.done()]:
                code, _ = running.pop(key).result()
                (done.add(key) if code == 0 else failed.append(key))
            used = sum(BY_KEY[k]["cores"] for k in running)
            if not failed:
                for key in keys:
                    s = BY_KEY[key]
                    if key in done or key in running or key in failed:
                        continue
                    if not all(d in done or external_done(d) for d in s["deps"]):
                        continue
                    if used + s["cores"] > cores and running:
                        continue
                    running[key] = pool.submit(run_stage, project_root, cid, key, log_dir)
                    used += s["cores"]
            if not running:
                break
            time.sleep(poll)
    finally:
        pool.shutdown(wait=True)
    return failed


def status(project_root, cid):
    """Latest ledger status per stage: {key: status}, plus counts."""
    latest = {}
    for event in campaign.ledger_read(ledger_path(project_root, cid)):
        latest[event["stage"]] = event["status"]
    for key in BY_KEY:
        if key in latest and latest[key] == "completed" and not stage_complete(project_root, cid, key):
            latest[key] = "completed_but_outputs_changed"
    counts = {}
    for v in latest.values():
        counts[v] = counts.get(v, 0) + 1
    return latest, counts
