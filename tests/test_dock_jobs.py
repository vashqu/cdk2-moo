"""Regression tests for the docking job store (Part 2A and 2B), with docking mocked and temporary stores.
Run: python -m unittest discover tests"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rdkit import Chem
from rdkit.Chem import AllChem

from cdk2moo import dock_jobs as dj
from cdk2moo import dock_validate as dv
from cdk2moo import docking


def embed(smiles, seed=1):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


class Counter:
    """Fake preparation and docking that count their calls; poses are real 3-D embeddings of the job's own molecule."""

    def __init__(self, fail_prep=(), pose_none=()):
        self.prepared, self.docked = 0, 0
        self.fail_prep, self.pose_none = set(fail_prep), set(pose_none)

    def prepare(self, smiles, prep_seed):
        self.prepared += 1
        if smiles in self.fail_prep:
            return None
        return "PDBQT:" + smiles, Chem.AddHs(Chem.MolFromSmiles(smiles)), 2

    def dock(self, pdbqt, vina_seed):
        self.docked += 1
        smiles = pdbqt.split(":", 1)[1]
        return (-7.0, None) if smiles in self.pose_none else (-7.0 - 0.1 * len(smiles), embed(smiles))


A, B, C = ("CCO", 42), ("CCN", 42), ("CCCl", 42)
JOBS = [A, B, C]


def jid(job, prep=42):
    return dj.job_id(job[0], job[1], prep)


def run(tmp, fake, jobs, **kw):
    return dj.run_jobs(tmp, jobs, fake.prepare, fake.dock, prep_seed=42, **kw)


def rewrite_record(tmp, job, **changes):
    path = Path(tmp) / "jobs" / f"{jid(job)}.json"
    record = json.loads(path.read_text())
    record.update(changes)
    path.write_text(json.dumps(record))


def recommit_pose(tmp, job, text):
    """Replace a pose file's text AND commit the matching hash, so only the content/metadata check can object."""
    sdf = Path(tmp) / "jobs" / f"{jid(job)}.sdf"
    sdf.write_text(text)
    rewrite_record(tmp, job, pose_sha256=dv.sha256_of(sdf))


class StoreCase(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        dj.check_manifest(self.tmp, {"exhaustiveness": 8})

    def sdf(self, job):
        return (self.tmp / "jobs" / f"{jid(job)}.sdf")


class ResumeTests(StoreCase):
    def test_two_jobs_stop_resume_third_everything_present_once(self):
        fake = Counter()
        run(self.tmp, fake, JOBS[:2])
        run(self.tmp, fake, JOBS)
        self.assertEqual(fake.docked, 3)
        self.assertEqual(len(dj.reconcile(self.tmp)["complete"]), 3)
        table = dj.export(self.tmp, self.tmp / "export")
        self.assertEqual(sorted(table["smiles"]), ["CCCl", "CCN", "CCO"])
        poses = [m.GetProp("smiles") for m in Chem.SDMolSupplier(str(self.tmp / "export" / "poses_seed42.sdf"))]
        self.assertEqual(sorted(poses), ["CCCl", "CCN", "CCO"])

    def test_resuming_a_completed_run_docks_nothing(self):
        fake = Counter()
        run(self.tmp, fake, JOBS)
        before = (fake.prepared, fake.docked)
        run(self.tmp, fake, JOBS)
        self.assertEqual((fake.prepared, fake.docked), before)

    def test_interruption_between_pose_and_record_leaves_a_detected_orphan(self):
        with mock.patch.object(dj, "atomic_write_text", side_effect=RuntimeError("killed before commit")):
            with self.assertRaises(RuntimeError):
                run(self.tmp, Counter(), [A])
        self.assertTrue(self.sdf(A).exists())
        self.assertFalse((self.tmp / "jobs" / f"{jid(A)}.json").exists())
        self.assertEqual(dj.reconcile(self.tmp)["orphan_poses"], [jid(A)])
        with self.assertRaises(dj.ResumeRefused):
            run(self.tmp, Counter(), [A])
        _, repaired = run(self.tmp, Counter(), [A], repair=True)
        self.assertEqual(repaired, {"orphan_poses": [jid(A)]})
        self.assertEqual(dj.reconcile(self.tmp)["complete"], [jid(A)])
        self.assertTrue(any((self.tmp / "quarantine").rglob(f"{jid(A)}.sdf")))     # the orphan was kept as evidence

    def test_interruption_during_pose_write_leaves_only_a_temporary_file_that_is_never_a_result(self):
        (self.tmp / "jobs" / f"{jid(A)}.sdf.tmp").write_text("half written")
        fake = Counter()
        run(self.tmp, fake, [A])
        self.assertEqual(fake.docked, 1)
        self.assertFalse((self.tmp / "jobs" / f"{jid(A)}.sdf.tmp").exists())
        self.assertEqual(dj.reconcile(self.tmp)["complete"], [jid(A)])

    def test_pose_none_never_produces_ok_and_the_export_error_is_recorded(self):
        records, _ = run(self.tmp, Counter(pose_none={"CCO"}), [A])
        record = records[jid(A)]
        self.assertEqual(record["status"], dv.POSE_FAILED)
        self.assertIn("no pose", record["error"])
        self.assertEqual(dj.reconcile(self.tmp)["complete"], [])
        self.assertFalse(self.sdf(A).exists())

    def test_failures_are_distinguished_and_retried_only_when_asked(self):
        fake = Counter(fail_prep={"CCN"})
        records, _ = run(self.tmp, fake, [A, B])
        self.assertEqual(records[jid(B)]["status"], dv.PREP_FAILED)
        before = fake.prepared
        run(self.tmp, fake, [A, B])
        self.assertEqual(fake.prepared, before)
        run(self.tmp, fake, [A, B], retry_failed=True)
        self.assertEqual(fake.prepared, before + 1)
        self.assertTrue(any((self.tmp / "quarantine").rglob(f"{jid(B)}.json")))

    def test_dock_exception_is_recorded_as_dock_failed(self):
        def boom(pdbqt, vina_seed):
            raise RuntimeError("vina crashed")
        records, _ = dj.run_jobs(self.tmp, [A], Counter().prepare, boom, prep_seed=42)
        self.assertEqual(records[jid(A)]["status"], dv.DOCK_FAILED)

    def test_job_identity_includes_preparation_seed(self):
        self.assertNotEqual(dj.job_id("CCO", 42, 1), dj.job_id("CCO", 42, 2))


class ScopedRepairTests(StoreCase):
    """A restricted request must never silently delete an unrelated job."""

    def test_repair_of_a_requested_job_does_not_touch_an_inconsistent_unrequested_job(self):
        run(self.tmp, Counter(), [A, B])
        self.sdf(B).unlink()                                                  # B becomes inconsistent (score without pose)
        before = (self.tmp / "jobs" / f"{jid(B)}.json").read_text()
        with self.assertRaisesRegex(dj.ResumeRefused, "outside the requested scope"):
            run(self.tmp, Counter(), [A], repair=True)                        # request A only
        self.assertEqual((self.tmp / "jobs" / f"{jid(B)}.json").read_text(), before)    # B's record is still there, unchanged
        self.assertFalse((self.tmp / "quarantine").exists())

    def test_requesting_the_inconsistent_job_repairs_it_and_keeps_the_evidence(self):
        run(self.tmp, Counter(), [A, B])
        self.sdf(B).unlink()
        fake = Counter()
        _, repaired = run(self.tmp, fake, [B], repair=True)
        self.assertEqual(repaired, {"missing_poses": [jid(B)]})
        self.assertEqual(fake.docked, 1)
        self.assertEqual(sorted(dj.reconcile(self.tmp)["complete"]), sorted([jid(A), jid(B)]))
        reasons = json.loads(next((self.tmp / "quarantine").rglob("reasons.json")).read_text())
        self.assertIn(jid(B), reasons["reasons"])
        self.assertTrue(any((self.tmp / "quarantine").rglob(f"{jid(B)}.json")))

    def test_without_repair_even_a_requested_inconsistency_stops_the_run(self):
        run(self.tmp, Counter(), [A])
        self.sdf(A).unlink()
        with self.assertRaisesRegex(dj.ResumeRefused, "inconsistent requested jobs"):
            run(self.tmp, Counter(), [A])

    def test_unattributable_files_are_never_deleted_by_a_restricted_repair(self):
        run(self.tmp, Counter(), [A])
        (self.tmp / "jobs" / "garbage.json").write_text("{broken")           # identity cannot be recovered
        (self.tmp / "jobs" / "notes.txt").write_text("unexpected")
        with self.assertRaises(dj.ResumeRefused):
            run(self.tmp, Counter(), [A], repair=True)
        self.assertTrue((self.tmp / "jobs" / "garbage.json").exists())
        self.assertTrue((self.tmp / "jobs" / "notes.txt").exists())
        run(self.tmp, Counter(), [A], repair=True, repair_scope="store")      # explicit store-wide repair: quarantined, not deleted
        self.assertFalse((self.tmp / "jobs" / "garbage.json").exists())
        self.assertTrue(any((self.tmp / "quarantine").rglob("garbage.json")))
        self.assertTrue(any((self.tmp / "quarantine").rglob("notes.txt")))

    def test_classes_of_inconsistency_stay_distinct(self):
        run(self.tmp, Counter(fail_prep={"CCCl"}), JOBS)
        self.sdf(A).unlink()                                                   # missing pose
        self.sdf(B).write_text("not an sdf")                                   # unreadable pose
        (self.tmp / "jobs" / f"{dv.job_id('CCCC', 42, 42)}.sdf").write_text("x")    # orphan
        report = dj.reconcile(self.tmp)
        self.assertEqual(report["missing_poses"], [jid(A)])
        self.assertEqual(report["unreadable_poses"], [jid(B)])
        self.assertEqual(report["orphan_poses"], [dv.job_id("CCCC", 42, 42)])
        self.assertEqual(report["failed"], [jid(C)])                           # a recorded failure is not an anomaly
        self.assertNotIn("failed", dv.anomalies(report))


class AssociationTests(StoreCase):
    """A readable SDF is not enough: the pose must belong to the job whose record points at it."""

    def setUp(self):
        super().setUp()
        run(self.tmp, Counter(), [A, B])

    def test_swapped_pose_files_are_detected_for_both_jobs(self):
        a, b = self.sdf(A), self.sdf(B)
        text_a, text_b = a.read_text(), b.read_text()
        a.write_text(text_b)
        b.write_text(text_a)
        report = dj.reconcile(self.tmp)
        self.assertEqual(sorted(report["pose_mismatch"]), sorted([jid(A), jid(B)]))
        self.assertIn("hash", report["details"][jid(A)])

    def test_swapped_files_with_recommitted_hashes_are_still_caught_by_the_embedded_identity(self):
        text_a, text_b = self.sdf(A).read_text(), self.sdf(B).read_text()
        recommit_pose(self.tmp, A, text_b)
        recommit_pose(self.tmp, B, text_a)
        report = dj.reconcile(self.tmp)
        self.assertEqual(sorted(report["pose_mismatch"]), sorted([jid(A), jid(B)]))
        self.assertIn("job_id", report["details"][jid(A)])

    def test_changing_only_the_record_job_id_is_detected(self):
        rewrite_record(self.tmp, A, job_id="0123456789abcdef")
        self.assertEqual(dj.reconcile(self.tmp)["id_mismatch"], [jid(A)])

    def test_changing_only_pose_metadata_is_detected(self):
        text = self.sdf(A).read_text()
        recommit_pose(self.tmp, A, text.replace(">  <smiles>  (1) \nCCO", ">  <smiles>  (1) \nCCC") if "CCO" in text else text)
        # same hash committed, so only the metadata check can object
        self.assertIn(jid(A), dj.reconcile(self.tmp)["pose_mismatch"])
        self.assertIn("smiles property", dj.reconcile(self.tmp)["details"][jid(A)])

    def test_changing_the_score_property_in_the_pose_is_detected(self):
        text = self.sdf(A).read_text()
        score = json.loads((self.tmp / "jobs" / f"{jid(A)}.json").read_text())["vina_score"]
        recommit_pose(self.tmp, A, text.replace(f"{score:.3f}", f"{score + 1:.3f}"))
        self.assertIn("vina_score property", dj.reconcile(self.tmp)["details"][jid(A)])

    def test_changing_pose_contents_while_keeping_metadata_is_detected_by_the_hash(self):
        text = self.sdf(A).read_text()
        self.sdf(A).write_text(text.replace("0.", "1.", 3))                   # perturb coordinates, keep every property
        report = dj.reconcile(self.tmp)
        self.assertIn(jid(A), report["pose_mismatch"])
        self.assertIn("hash", report["details"][jid(A)])

    def test_a_different_molecule_with_spoofed_metadata_is_detected_by_chemical_identity(self):
        other = Chem.AddHs(Chem.MolFromSmiles("CCCCCC"))
        AllChem.EmbedMolecule(other, randomSeed=1)
        other.SetProp("smiles", "CCO")
        other.SetProp("vina_score", f"{json.loads((self.tmp / 'jobs' / (jid(A) + '.json')).read_text())['vina_score']:.3f}")
        other.SetProp("job_id", jid(A))
        writer = Chem.SDWriter(str(self.sdf(A)))
        writer.write(other)
        writer.close()
        recommit_pose(self.tmp, A, self.sdf(A).read_text())
        self.assertIn("different connectivity", dj.reconcile(self.tmp)["details"][jid(A)])

    def test_an_sdf_with_two_records_is_refused(self):
        text = self.sdf(A).read_text()
        recommit_pose(self.tmp, A, text + text)
        self.assertIn("2 pose records", dj.reconcile(self.tmp)["details"][jid(A)])

    def test_invalid_status_and_malformed_field_types_are_not_silent_completions(self):
        for change, expect in [({"status": "weird"}, "status"), ({"vina_score": "-7.0"}, "vina_score"),
                               ({"vina_seed": True}, "vina_seed"), ({"vina_score": float("nan")}, "not finite"),
                               ({"pose_sha256": None}, "pose hash"), ({"extra": 1}, "unexpected field")]:
            with self.subTest(change=change):
                rewrite_record(self.tmp, A, **change)
                report = dj.reconcile(self.tmp)
                self.assertIn(jid(A), report["incomplete_records"])
                self.assertIn(expect, report["details"][jid(A)])
                self.assertNotIn(jid(A), report["complete"] + report["failed"])        # never a completed failure either
                run(self.tmp, Counter(), [A], repair=True)                            # restore a valid record for the next case

    def test_a_committed_failure_cannot_carry_a_pose(self):
        run(self.tmp, Counter(fail_prep={"CCCl"}), [C])
        (self.tmp / "jobs" / f"{jid(C)}.sdf").write_text(self.sdf(A).read_text())
        self.assertIn(jid(C), dj.reconcile(self.tmp)["orphan_poses"])


class IdentityPolicyTests(unittest.TestCase):
    def check(self, smiles, pose_smiles):
        return dv.pose_identity(smiles, embed(pose_smiles))

    def test_hydrogens_and_coordinates_do_not_change_identity(self):
        self.assertTrue(self.check("CCO", "CCO")[0])

    def test_unspecified_stereo_in_the_input_is_not_compared(self):
        self.assertTrue(self.check("CC(O)CC", "C[C@H](O)CC")[0])
        self.assertTrue(self.check("CC(O)CC", "C[C@@H](O)CC")[0])

    def test_specified_stereo_in_the_input_must_match(self):
        self.assertTrue(self.check("C[C@H](O)CC", "C[C@H](O)CC")[0])
        ok, reason = self.check("C[C@H](O)CC", "C[C@@H](O)CC")
        self.assertFalse(ok)
        self.assertIn("stereochemistry", reason)

    def test_different_molecules_are_never_identical_even_with_matching_fingerprint_features(self):
        self.assertFalse(self.check("CCO", "CCN")[0])


class LockAndManifestTests(unittest.TestCase):
    def test_a_second_writer_is_refused(self):
        tmp = Path(tempfile.mkdtemp())
        with dj.store_lock(tmp):
            with self.assertRaisesRegex(dj.ResumeRefused, "another process"):
                with dj.store_lock(tmp):
                    pass
        with dj.store_lock(tmp):                                               # released after the first writer finishes
            pass

    def test_run_jobs_refuses_while_another_writer_holds_the_store(self):
        tmp = Path(tempfile.mkdtemp())
        dj.check_manifest(tmp, {"x": 1})
        with dj.store_lock(tmp):
            with self.assertRaises(dj.ResumeRefused):
                dj.run_jobs(tmp, [A], Counter().prepare, Counter().dock, prep_seed=42)

    def test_changed_settings_refuse_resume_with_the_differing_key(self):
        tmp = Path(tempfile.mkdtemp())
        dj.check_manifest(tmp, {"receptor_sha256": "aaa", "exhaustiveness": 8, "ligand_prep_seed": 42})
        self.assertEqual(dj.check_manifest(tmp, {"receptor_sha256": "aaa", "exhaustiveness": 8, "ligand_prep_seed": 42}), "matches")
        with self.assertRaisesRegex(dj.ResumeRefused, "receptor_sha256"):
            dj.check_manifest(tmp, {"receptor_sha256": "bbb", "exhaustiveness": 8, "ligand_prep_seed": 42})

    def test_legacy_directory_without_manifest_is_refused_not_adopted(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "jobs").mkdir()
        (tmp / "jobs" / "abc.json").write_text("{}")
        with self.assertRaisesRegex(dj.ResumeRefused, "legacy"):
            dj.check_manifest(tmp, {"x": 1})
        self.assertFalse((tmp / "manifest.json").exists())


class ExportTests(StoreCase):
    def test_export_keeps_failures_visible_and_can_be_limited_to_a_scope(self):
        run(self.tmp, Counter(fail_prep={"CCN"}), JOBS)
        table = dj.export(self.tmp, self.tmp / "export", only_smiles={"CCO", "CCN"})
        self.assertEqual(sorted(table["smiles"]), ["CCN", "CCO"])
        self.assertEqual(table.set_index("smiles").loc["CCN", "status"], dv.PREP_FAILED)     # stays in the denominator
        poses = [m.GetProp("smiles") for m in Chem.SDMolSupplier(str(self.tmp / "export" / "poses_seed42.sdf"))]
        self.assertEqual(poses, ["CCO"])

    def test_export_refuses_an_inconsistent_store(self):
        run(self.tmp, Counter(), [A])
        self.sdf(A).unlink()
        with self.assertRaises(dj.ResumeRefused):
            dj.export(self.tmp, self.tmp / "export")


class GridMapTests(unittest.TestCase):
    def make_maps(self, tmp, receptor):
        centre, size = self.box
        files = {}
        for name in ["maps.A.map", "maps.C.map", "maps.e.map"]:
            (tmp / name).write_text(name * 10)
            files[name] = {"bytes": (tmp / name).stat().st_size, "sha256": dv.sha256_of(tmp / name)}
        manifest = {"receptor_sha256": dv.sha256_of(receptor), "box_center": centre, "box_size": size, "files": files}
        (tmp / docking.MAPS_MANIFEST).write_text(json.dumps(manifest))

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.box = ([1.0, 2.0, 3.0], [22.0, 22.0, 22.0])
        self.receptor = self.tmp / "receptor.pdbqt"
        self.receptor.write_text("ATOM")
        self.maps = self.tmp / "maps"
        self.maps.mkdir()
        self.make_maps(self.maps, self.receptor)

    def test_complete_matching_maps_are_accepted(self):
        self.assertTrue(docking.validate_grid_maps(self.maps, self.receptor, box=self.box))
        with self.assertRaisesRegex(dj.ResumeRefused, "different search box"):
            docking.validate_grid_maps(self.maps, self.receptor, box=([9.0, 9.0, 9.0], [22.0, 22.0, 22.0]))

    def test_incomplete_or_legacy_or_tampered_maps_are_refused(self):
        (self.maps / "maps.C.map").unlink()
        with self.assertRaisesRegex(dj.ResumeRefused, "missing"):
            docking.validate_grid_maps(self.maps, self.receptor, box=self.box)
        self.make_maps(self.maps, self.receptor)
        (self.maps / "maps.A.map").write_text("tampered")
        with self.assertRaisesRegex(dj.ResumeRefused, "hash"):
            docking.validate_grid_maps(self.maps, self.receptor, box=self.box)
        (self.maps / docking.MAPS_MANIFEST).unlink()
        with self.assertRaisesRegex(dj.ResumeRefused, "legacy"):
            docking.validate_grid_maps(self.maps, self.receptor, box=self.box)


if __name__ == "__main__":
    unittest.main()


class NumericLookingJobIds(unittest.TestCase):
    """A job id made only of digits, or digits with one 'e', must still match the id stored in its pose file."""

    def test_pose_with_numeric_looking_job_ids_is_consistent(self):
        smiles = "CCO"
        for jid in ("4445476191389595", "84868095015190e9"):
            with tempfile.TemporaryDirectory() as tmp:
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                AllChem.EmbedMolecule(mol, randomSeed=1)
                mol.SetProp("job_id", jid)
                mol.SetProp("smiles", smiles)
                mol.SetProp("vina_score", "-5.000")
                writer = Chem.SDWriter(str(Path(tmp) / f"{jid}.sdf"))
                writer.write(mol)
                writer.close()
                record = {"smiles": smiles, "vina_score": -5.0, "pose_sha256": dv.sha256_of(Path(tmp) / f"{jid}.sdf")}
                self.assertIsNone(dv.check_pose(tmp, jid, record), jid)
