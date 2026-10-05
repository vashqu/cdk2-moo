"""Tests for output-path safety and manifests (Task 3). Run: python -m unittest discover tests"""

import json
import tempfile
import unittest
from pathlib import Path

from cdk2moo.run_outputs import choose_output, run_tokens, write_manifest

PRIMARY = "05_ga_populations.csv"
ARMS = ["multi_real", "multi_scrambled", "activity_only", "druglike_only"]
SEEDS = [42, 43, 44, 45, 46]


def tokens(**kw):
    base = dict(relaxed=False, policy="legacy", arms=ARMS, default_arms=ARMS, seeds=SEEDS, default_seeds=SEEDS,
                generations=50, default_generations=50)
    base.update(kw)
    return run_tokens(**base)


class OutputTests(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())

    def test_primary_protocol_gets_the_primary_name(self):
        self.assertEqual(choose_output(self.dir, PRIMARY, tokens=tokens()).name, PRIMARY)

    def test_relaxed_without_out_cannot_target_the_primary_filename(self):
        path = choose_output(self.dir, PRIMARY, tokens=tokens(relaxed=True))
        self.assertNotEqual(path.name, PRIMARY)
        self.assertIn("relaxed", path.name)

    def test_shortened_smoke_run_cannot_overwrite_the_primary_result(self):
        (self.dir / PRIMARY).write_text("historical")
        path = choose_output(self.dir, PRIMARY, tokens=tokens(seeds=[42], generations=3, arms=["multi_real"]), overwrite=True)
        self.assertNotEqual(path.name, PRIMARY)
        for part in ("seeds-42", "gen-3", "arms-multi_real"):
            self.assertIn(part, path.name)
        self.assertEqual((self.dir / PRIMARY).read_text(), "historical")

    def test_explicit_primary_name_for_a_nonprimary_run_is_refused_even_with_overwrite(self):
        with self.assertRaisesRegex(ValueError, "reserved"):
            choose_output(self.dir, PRIMARY, out=PRIMARY, tokens=tokens(generations=3), overwrite=True)

    def test_corrected_policy_has_its_own_name(self):
        path = choose_output(self.dir, PRIMARY, tokens=tokens(policy="corrected"))
        self.assertIn("corrected", path.name)

    def test_explicit_output_path_is_honoured(self):
        target = self.dir / "mine.csv"
        self.assertEqual(choose_output(self.dir, PRIMARY, out=str(target), tokens=tokens(relaxed=True)), target.resolve())
        self.assertEqual(choose_output(self.dir, PRIMARY, out="mine.csv", tokens=tokens(relaxed=True)), target.resolve())

    def test_existing_output_is_protected_unless_overwrite_is_requested(self):
        target = self.dir / "mine.csv"
        target.write_text("x")
        with self.assertRaises(FileExistsError):
            choose_output(self.dir, PRIMARY, out=str(target))
        self.assertEqual(choose_output(self.dir, PRIMARY, out=str(target), overwrite=True), target.resolve())

    def test_existing_primary_file_is_protected_by_default(self):
        (self.dir / PRIMARY).write_text("historical")
        with self.assertRaises(FileExistsError):
            choose_output(self.dir, PRIMARY, tokens=tokens())

    def test_a_symlinked_or_relative_spelling_cannot_dodge_the_reserved_name(self):
        link_dir = self.dir.parent / (self.dir.name + "_link")
        link_dir.symlink_to(self.dir)
        with self.assertRaisesRegex(ValueError, "reserved"):
            choose_output(self.dir, PRIMARY, out=str(link_dir / PRIMARY), tokens=tokens(generations=3), overwrite=True)
        with self.assertRaisesRegex(ValueError, "reserved"):
            choose_output(self.dir, PRIMARY, out=str(self.dir / "sub" / ".." / PRIMARY), tokens=tokens(generations=3), overwrite=True)

    def test_commit_writes_the_output_atomically_and_the_manifest_marks_completion(self):
        import pandas as pd
        from cdk2moo.run_outputs import commit_csv
        out = self.dir / "run.csv"
        commit_csv(pd.DataFrame({"a": [1, 2]}), out, index=False)
        self.assertTrue(out.exists())
        self.assertFalse((self.dir / "run.csv.part").exists())
        self.assertFalse((self.dir / "run.csv.manifest.json").exists())           # not complete until a manifest exists
        inp = self.dir / "in.txt"
        inp.write_text("input")
        data = json.loads(write_manifest(out, {"k": 1}, inputs={"in": inp, "missing": self.dir / "nope"}).read_text())
        self.assertTrue(data["complete"])
        self.assertEqual(data["outputs"]["run.csv"]["bytes"], out.stat().st_size)
        self.assertEqual(len(data["inputs"]["in"]["sha256"]), 64)
        self.assertIsNone(data["inputs"]["missing"])                              # a missing input is recorded as missing
        self.assertIn("source_tree_sha256", data["source"])
        self.assertIn("working_tree_clean", data["source"])
        self.assertIn("FP_RADIUS", data["parameters"])                            # effective parameters are recorded

    def test_manifest_records_effective_settings_and_versions(self):
        out = self.dir / "run.csv"
        out.write_text("x")
        path = write_manifest(out, {"seeds": [42], "generations": 3, "policy": "corrected"})
        data = json.loads(path.read_text())
        self.assertEqual(data["settings"]["generations"], 3)
        self.assertIn("rdkit", data["versions"])
        self.assertEqual(data["output"], "run.csv")
        self.assertTrue(path.name.endswith(".manifest.json"))


if __name__ == "__main__":
    unittest.main()
