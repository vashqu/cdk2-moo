"""Regression tests for the PoseBusters audit (Task 1). Run: python -m unittest discover tests"""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from cdk2moo import posebusters_audit as pa

REQ = pa.REQUIRED_CHECKS


def make_report(ids, overrides=None, drop=()):
    """A PoseBusters-like report: every required check True, with chosen values replaced."""
    data = {c: [True] * len(ids) for c in REQ if c not in drop}
    for name, values in (overrides or {}).items():
        data[name] = values
    return pd.DataFrame(data, index=pd.Index(ids, name="pose_id"))


def records_for(ids, smiles=None, unparsed=()):
    return [{"position": i, "pose_id": pid, "smiles": (smiles or ids)[i],
             "mol": None if pid in unparsed else object()} for i, pid in enumerate(ids)]


class AggregationTests(unittest.TestCase):
    def test_check_with_true_false_missing_stays_in_report(self):
        ids = ["a", "b", "c"]
        report = make_report(ids, {"internal_energy": [True, False, np.nan]})
        self.assertEqual(report["internal_energy"].dtype, object)   # the old dtype==bool filter would drop it
        out = pa.aggregate_status(report, records_for(ids))
        self.assertIn("internal_energy", out.columns)
        self.assertEqual(out.set_index("pose_id")["status"].to_dict(), {"a": "pass", "b": "fail", "c": "unevaluable"})

    def test_all_missing_column_remains_required(self):
        ids = ["a", "b"]
        report = make_report(ids, {"internal_energy": [np.nan, np.nan]})
        out = pa.aggregate_status(report, records_for(ids))
        self.assertTrue((out["status"] == "unevaluable").all())
        self.assertTrue((out["n_missing"] == 1).all())

    def test_missing_required_column_raises_clear_diagnostic(self):
        report = make_report(["a"], drop=("internal_energy",))
        with self.assertRaisesRegex(ValueError, "internal_energy"):
            pa.aggregate_status(report, records_for(["a"]))

    def test_missing_internal_energy_prevents_complete_pass(self):
        report = make_report(["a"], {"internal_energy": [np.nan]})
        row = pa.aggregate_status(report, records_for(["a"])).iloc[0]
        self.assertFalse(row["pb_pass"])
        self.assertEqual(row["status"], "unevaluable")
        self.assertEqual(row["missing_checks"], "internal_energy")

    def test_failure_outranks_missing_and_missing_counts_are_kept(self):
        report = make_report(["a"], {"bond_angles": [False], "internal_energy": [np.nan]})
        row = pa.aggregate_status(report, records_for(["a"])).iloc[0]
        self.assertEqual(row["status"], "fail")
        self.assertEqual((row["n_failed"], row["n_missing"]), (1, 1))

    def test_reordered_report_rows_map_to_correct_molecules(self):
        ids = ["p0", "p1", "p2"]
        report = make_report(ids, {"bond_angles": [True, False, True]}).iloc[[2, 0, 1]]
        out = pa.aggregate_status(report, records_for(ids, smiles=["CCO", "CCN", "CCC"])).set_index("smiles")
        self.assertEqual(out.loc["CCN", "status"], "fail")
        self.assertEqual(out.loc["CCO", "status"], "pass")

    def test_unparsed_record_keeps_its_position_and_is_unevaluable(self):
        ids = ["p0", "p1", "p2"]
        report = make_report(["p0", "p2"])
        out = pa.aggregate_status(report, records_for(ids, smiles=["A", "B", "C"], unparsed=("p1",)))
        self.assertEqual(out["record_position"].tolist(), [0, 1, 2])
        self.assertEqual(out.set_index("smiles").loc["B", "status"], "unevaluable")
        self.assertEqual(out.set_index("smiles").loc["C", "status"], "pass")     # later identity not shifted

    def test_duplicate_pose_ids_detected(self):
        with self.assertRaisesRegex(ValueError, "duplicate pose ids"):
            pa.aggregate_status(make_report(["a"]), records_for(["a", "a"]))

    def test_summary_counts_and_denominator(self):
        report = make_report(["a", "b", "c"], {"bond_angles": [True, False, True], "internal_energy": [True, True, np.nan]})
        table = pa.aggregate_status(report, records_for(["a", "b", "c"]))
        s = pa.summarize(table, expected=3)
        self.assertEqual((s["passed"], s["failed"], s["unevaluable"], s["denominator_for_rates"]), (1, 1, 1, 3))


class JoinTests(unittest.TestCase):
    def test_join_keeps_every_molecule_and_never_passes_a_missing_audit(self):
        table = pa.aggregate_status(make_report(["a"]), records_for(["a"], smiles=["CCO"]))
        molecules = pd.DataFrame({"smiles": ["CCO", "CCO", "CCN"], "set": ["x", "y", "x"]})
        merged = pa.attach_to_molecules(molecules, table)
        self.assertEqual(len(merged), 3)                                   # cardinality preserved
        self.assertEqual(merged.set_index(["set", "smiles"]).loc[("x", "CCN"), "status"], "no_audit_result")
        self.assertFalse(merged.loc[merged["smiles"] == "CCN", "pb_pass"].iloc[0])

    def test_duplicate_audit_keys_raise(self):
        table = pd.DataFrame({"smiles": ["CCO", "CCO"], "status": ["pass", "pass"], "pb_pass": [True, True],
                              "n_failed": [0, 0], "n_missing": [0, 0], "failed_checks": ["", ""], "missing_checks": ["", ""]})
        with self.assertRaisesRegex(ValueError, "not unique"):
            pa.attach_to_molecules(pd.DataFrame({"smiles": ["CCO"]}), table)


class SdfTests(unittest.TestCase):
    GOOD = ("\n  RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    2.5000    1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0\n  2  3  1  0\nM  END\n> <smiles>\n{smi}\n\n$$$$\n")

    def test_malformed_record_does_not_shift_later_identities(self):
        text = (self.GOOD.format(smi="CCO") + "garbage record\nnot a molfile\n> <smiles>\nBROKEN\n\n$$$$\n"
                + self.GOOD.format(smi="CCN"))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "poses.sdf"
            path.write_text(text)
            records = pa.read_sdf_records(path)
        self.assertEqual([r["position"] for r in records], [0, 1, 2])
        self.assertEqual([r["smiles"] for r in records], ["CCO", "BROKEN", "CCN"])
        self.assertIsNone(records[1]["mol"])
        self.assertIsNotNone(records[2]["mol"])
        self.assertEqual(records[2]["mol"].GetProp("_Name"), "pose000002")

    def test_identity_read_from_rdkit_style_property_header(self):
        text = self.GOOD.format(smi="CCO").replace("> <smiles>", ">  <smiles>  (1) ")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "poses.sdf"
            path.write_text(text)
            records = pa.read_sdf_records(path)
        self.assertEqual(records[0]["smiles"], "CCO")

    def test_unexpected_row_loss_detected(self):
        class DropsOne:
            def bust(self, mols, true, protein, full_report=False):
                names = [m.GetProp("_Name") for m in mols][:-1]
                return pd.DataFrame({c: [True] * len(names) for c in REQ},
                                    index=pd.MultiIndex.from_arrays([["f"] * len(names), names, [0] * len(names)],
                                                                    names=["file", "molecule", "position"]))
        text = self.GOOD.format(smi="CCO") + self.GOOD.format(smi="CCN")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "poses.sdf"
            path.write_text(text)
            records = pa.read_sdf_records(path)
        with self.assertRaisesRegex(ValueError, "lacks 1 parsed poses"):
            pa.run_posebusters(records, "protein.pdb", DropsOne())

    def test_duplicate_ids_from_posebusters_detected(self):
        class Repeats:
            def bust(self, mols, true, protein, full_report=False):
                names = ["pose000000", "pose000000"]
                return pd.DataFrame({c: [True, True] for c in REQ},
                                    index=pd.MultiIndex.from_arrays([["f", "f"], names, [0, 0]],
                                                                    names=["file", "molecule", "position"]))
        text = self.GOOD.format(smi="CCO") + self.GOOD.format(smi="CCN")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "poses.sdf"
            path.write_text(text)
            records = pa.read_sdf_records(path)
        with self.assertRaisesRegex(ValueError, "repeated pose id"):
            pa.run_posebusters(records, "protein.pdb", Repeats())


class InstalledApiTests(unittest.TestCase):
    """Compares the frozen check list with the installed PoseBusters on three real poses."""

    def test_required_checks_match_installed_dock_config_output(self):
        from posebusters import PoseBusters
        root = Path(__file__).resolve().parents[1]
        poses = root / "data" / "structures" / "redock_poses_seed42.sdf"
        protein = root / "data" / "structures" / "receptor_H.pdb"
        if not poses.exists() or not protein.exists():
            self.skipTest("structure files not present")
        records = pa.read_sdf_records(poses)
        report = pa.run_posebusters(records, protein, PoseBusters(config="dock"))
        self.assertEqual(set(report.columns), set(REQ))
        for name in REQ:                                                   # inapplicable checks are True, never NaN
            if "cofactor" in name or "water" in name:
                self.assertTrue(report[name].notna().all())


if __name__ == "__main__":
    unittest.main()
