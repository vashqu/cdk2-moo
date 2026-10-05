"""Integration tests for the campaign mechanism (Part 2D): downstream stages follow the selected campaign and cannot silently read historical
data; the driver verifies dependencies by hash. Fixtures are deliberately different between 'historical' and 'campaign'.
Run: python -m unittest discover tests"""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

from cdk2moo import campaign, campaign_driver as driver

ROOT = Path(__file__).resolve().parents[1]
N = 150


def write_fixture(base, tag, atom_range, seed):
    """processed/ + results/ files for 05b under `base`, with every molecule marked by `tag` so reads can be told apart."""
    rng = np.random.default_rng(seed)
    (base / "processed").mkdir(parents=True, exist_ok=True)
    (base / "results").mkdir(parents=True, exist_ok=True)
    curated = pd.DataFrame({"inchikey": [f"{tag}{i:04d}" for i in range(N)], "std_smiles": ["C" * (1 + i % 5) for i in range(N)],
                            "n_heavy_atoms": rng.integers(*atom_range, N)})
    curated.to_csv(base / "processed" / "cdk2_ic50_curated.csv", index=False)
    np.save(base / "processed" / "ecfp4.npy", (rng.random((N, 64)) < 0.3).astype(np.uint8))
    pop = pd.DataFrame({"smiles": [f"{tag}_mol{i}" for i in range(20)], "max_tanimoto": rng.random(20),
                        "n_heavy_atoms": rng.integers(*atom_range, 20), "generation": 1, "birth_generation": 1, "arm": "multi_real", "seed": 1})
    pop.to_csv(base / "results" / "05_ga_populations.csv", index=False)


def run_script(project, script, env_extra, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith("CDK2_")}
    env.update({"CDK2_PROJECT_ROOT": str(project), **env_extra})
    return subprocess.run([sys.executable, str(ROOT / "scripts" / script), *args], env=env, capture_output=True, text=True)


class DownstreamFollowsTheCampaign(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.project = Path(tempfile.mkdtemp())
        write_fixture(cls.project / "data", "HIST", (20, 30), 1)                                  # historical layout: data/processed, results
        (cls.project / "results").mkdir(exist_ok=True)
        pd.read_csv(cls.project / "data/processed/cdk2_ic50_curated.csv")                          # sanity: exists
        os.replace(cls.project / "data/results/05_ga_populations.csv", cls.project / "results/05_ga_populations.csv")
        write_fixture(cls.project / "campaigns" / "new" / "legacy", "NEW", (40, 50), 2)           # the selected campaign scope, different data

    def test_a_script_with_no_campaign_selected_refuses_to_start(self):
        result = run_script(self.project, "05b_size_conditioned_similarity.py", {})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No campaign selected", result.stderr + result.stdout)
        self.assertFalse((self.project / "NO_CAMPAIGN_SELECTED").exists())            # nothing was created either

    def test_the_selected_campaign_is_read_and_historical_files_are_not(self):
        result = run_script(self.project, "05b_size_conditioned_similarity.py", {"CDK2_CAMPAIGN": "new", "CDK2_SCOPE": "legacy"})
        self.assertEqual(result.returncode, 0, result.stderr[-600:])
        out = pd.read_csv(self.project / "campaigns/new/legacy/results/05_ga_populations_scored.csv")
        self.assertTrue(out["smiles"].str.startswith("NEW").all())
        self.assertFalse(out["smiles"].str.contains("HIST").any())
        ref = pd.read_csv(self.project / "campaigns/new/legacy/results/05_training_loo_similarity.csv")
        self.assertTrue(ref["inchikey"].str.startswith("NEW").all())
        self.assertFalse((self.project / "results/05_ga_populations_scored.csv").exists())   # nothing written to the historical results

    def test_historical_reanalysis_is_possible_only_by_explicit_selection(self):
        result = run_script(self.project, "05b_size_conditioned_similarity.py", {"CDK2_CAMPAIGN": "historical"}, "--out", "05_hist_scored.csv")
        self.assertEqual(result.returncode, 0, result.stderr[-600:])
        out = pd.read_csv(self.project / "results/05_hist_scored.csv")
        self.assertTrue(out["smiles"].str.startswith("HIST").all())

    def test_a_campaign_scope_missing_its_upstream_file_fails_instead_of_falling_back(self):
        empty = self.project / "campaigns" / "empty" / "legacy"
        (empty / "processed").mkdir(parents=True)
        result = run_script(self.project, "05b_size_conditioned_similarity.py", {"CDK2_CAMPAIGN": "empty", "CDK2_SCOPE": "legacy"})
        self.assertNotEqual(result.returncode, 0)                                      # does not read data/processed or results/
        self.assertIn("No such file", result.stderr)

    def test_policy_scopes_refuse_the_other_policy(self):
        result = run_script(self.project, "05_run_ga.py", {"CDK2_CAMPAIGN": "new", "CDK2_SCOPE": "legacy"}, "--policy", "corrected")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("holds policy 'legacy'", result.stderr + result.stdout)


class StageTableTests(unittest.TestCase):
    def test_graph_is_complete_and_acyclic_and_every_script_exists(self):
        keys = set(driver.BY_KEY)
        for key, s in driver.BY_KEY.items():
            self.assertTrue(set(s["deps"]) <= keys, key)
            if s["script"] != "@link":
                self.assertTrue((ROOT / "scripts" / s["script"]).exists(), s["script"])
        order, seen = [], set()

        def visit(key, trail=()):
            self.assertNotIn(key, trail, f"cycle through {key}")
            if key in seen:
                return
            for dep in driver.BY_KEY[key]["deps"]:
                visit(dep, trail + (key,))
            seen.add(key)
            order.append(key)
        for key in keys:
            visit(key)
        self.assertEqual(len(order), len(keys))

    def test_both_policy_scopes_have_the_same_stages_and_never_share_outputs(self):
        ids = {scope: {s["id"] for s in driver.STAGES if s["scope"] == scope} for scope in driver.POLICY_SCOPES}
        self.assertEqual(ids["legacy"], ids["corrected"])
        for s in driver.STAGES:
            for dep in s["deps"]:
                dep_scope = dep.split("/")[0]
                self.assertIn(dep_scope, ("shared", s["scope"]) if dep_scope != s["scope"] else (s["scope"],)) if not (
                    s["id"] == "dock" and dep.endswith("/dock")) else None

    def test_the_two_docking_stages_are_serialized_through_the_shared_store(self):
        self.assertIn("legacy/dock", driver.BY_KEY["corrected/dock"]["deps"])
        self.assertEqual(driver.BY_KEY["legacy/dock"]["cores"], 8)


class DriverTests(unittest.TestCase):
    """The scheduler, ledger and manifests, on a toy stage table with tiny scripts."""

    def setUp(self):
        self.project = Path(tempfile.mkdtemp())
        self.scripts = self.project / "scripts"
        self.scripts.mkdir()
        (self.scripts / "make_a.py").write_text("import os,pathlib;p=pathlib.Path(os.environ['CDK2_PROJECT_ROOT'])/'campaigns'/os.environ['CDK2_CAMPAIGN']/os.environ['CDK2_SCOPE']/'results';p.mkdir(parents=True,exist_ok=True);(p/'a.txt').write_text('A')")
        (self.scripts / "use_a.py").write_text("import os,pathlib;b=pathlib.Path(os.environ['CDK2_PROJECT_ROOT'])/'campaigns'/os.environ['CDK2_CAMPAIGN']/os.environ['CDK2_SCOPE']/'results';(b/'b.txt').write_text((b/'a.txt').read_text()+'B')")
        (self.scripts / "fail.py").write_text("raise SystemExit(3)")
        toy = [driver.stage("a", "shared", "make_a.py", outputs=["results/a.txt"]),
               driver.stage("b", "shared", "use_a.py", deps=["shared/a"], outputs=["results/b.txt"]),
               driver.stage("c", "shared", "fail.py", deps=["shared/b"], outputs=["results/c.txt"]),
               driver.stage("d", "shared", "use_a.py", deps=["shared/c"], outputs=["results/b.txt"])]
        self.patches = [mock.patch.object(driver, "STAGES", toy), mock.patch.object(driver, "BY_KEY", {f"{s['scope']}/{s['id']}": s for s in toy}),
                        mock.patch.object(driver, "SCRIPTS", self.scripts)]
        for p in self.patches:
            p.start()
        (campaign_root := driver.campaign_root(self.project, "t")).mkdir(parents=True)

    def tearDown(self):
        for p in self.patches:
            p.stop()

    def test_stages_run_in_dependency_order_with_hashes_ledger_and_resume(self):
        failed = driver.run_campaign(self.project, "t", cores=2, only=["shared/a", "shared/b"], poll=0.05)
        self.assertEqual(failed, [])
        manifest = json.loads(driver.stage_manifest_path(self.project, "t", "shared/b").read_text())
        self.assertTrue(manifest["complete"])
        self.assertIn("shared/a", manifest["consumed_upstream_outputs"])          # records the upstream hashes it consumed
        events = [(e["stage"], e["status"]) for e in campaign.ledger_read(driver.ledger_path(self.project, "t"))]
        self.assertLess(events.index(("shared/a", "completed")), events.index(("shared/b", "running")))
        before = len(events)
        driver.run_campaign(self.project, "t", cores=2, only=["shared/a", "shared/b"], poll=0.05)   # resume: nothing reruns
        self.assertEqual(len(campaign.ledger_read(driver.ledger_path(self.project, "t"))), before)

    def test_a_failure_is_recorded_and_blocks_its_dependents(self):
        failed = driver.run_campaign(self.project, "t", cores=2, poll=0.05)
        self.assertEqual(failed, ["shared/c"])
        latest, _ = driver.status(self.project, "t")
        self.assertEqual(latest["shared/c"], "failed")
        self.assertNotIn("shared/d", latest)                                       # never launched

    def test_changed_upstream_output_invalidates_the_stage_so_it_is_rerun_not_trusted(self):
        driver.run_campaign(self.project, "t", cores=2, only=["shared/a"], poll=0.05)
        out = driver.scope_root(self.project, "t", "shared") / "results" / "a.txt"
        out.write_text("TAMPERED")
        self.assertFalse(driver.stage_complete(self.project, "t", "shared/a"))
        driver.run_campaign(self.project, "t", cores=2, only=["shared/a"], poll=0.05)            # redone; old output kept as evidence
        self.assertEqual(out.read_text(), "A")
        self.assertTrue(any((driver.campaign_root(self.project, "t") / "superseded").rglob("a.txt")))

    def test_linking_refuses_a_modified_shared_file(self):
        shared = driver.scope_root(self.project, "t", "shared")
        scope = driver.scope_root(self.project, "t", "legacy")
        for rel in driver.LINKED_FILES:
            (shared / rel).parent.mkdir(parents=True, exist_ok=True)
            (shared / rel).write_text("x")
            (scope / rel).parent.mkdir(parents=True, exist_ok=True)
        (scope / driver.LINKED_FILES[0]).write_text("different")
        with mock.patch.object(driver, "validate_core_inputs", return_value={}):
            with self.assertRaisesRegex(RuntimeError, "differs from the shared file"):
                driver.link_shared(self.project, "t", "legacy")


class CoreInputValidationTests(unittest.TestCase):
    def test_misaligned_core_inputs_are_refused_at_the_boundary(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "processed").mkdir()
        (tmp / "results").mkdir()
        curated = pd.DataFrame({"inchikey": [f"K{i}" for i in range(6)], "std_smiles": ["CCO", "CCN", "CCC", "CCCl", "CCF", "CCBr"]})
        curated.to_csv(tmp / "processed" / "cdk2_ic50_curated.csv", index=False)
        from cdk2moo.features import ecfp4
        fps = ecfp4(curated["std_smiles"].tolist())
        np.save(tmp / "processed" / "ecfp4.npy", fps[::-1])                               # rows in the wrong order
        with self.assertRaisesRegex(ValueError, "fingerprints do not match"):
            driver.validate_core_inputs(tmp / "processed", tmp / "results")
        np.save(tmp / "processed" / "ecfp4.npy", fps[:4])
        with self.assertRaisesRegex(ValueError, "4 fingerprints for 6"):
            driver.validate_core_inputs(tmp / "processed", tmp / "results")


if __name__ == "__main__":
    unittest.main()
