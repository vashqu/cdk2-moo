"""Tests for the matched, deterministic initialization (Part 2F). Run: python -m unittest discover tests"""

import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
from rdkit import Chem

from cdk2moo import config, starts

PYRIDONE, PYRIDINOL = "O=c1cccc[nH]1", "Oc1ccccn1"


def atoms(smiles_list):
    return [Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s in smiles_list]


class EligibilityTests(unittest.TestCase):
    def test_reasons_are_recorded_for_each_exclusion(self):
        smiles = [PYRIDONE, PYRIDINOL, "CCCC[Sn](CCCC)CCCC", "CCO", PYRIDONE, "CC(C)Cc1ccccc1"]
        rows = starts.eligible_rows(smiles, atoms(smiles), (5, 50))
        by_row = {r["row"]: r for r in rows}
        self.assertTrue(by_row[0]["eligible"])
        self.assertEqual(by_row[1]["reason"], "not_idempotent")                    # the stored structure is not the standardized one
        self.assertIn("rejected_standardize_failed", by_row[2]["reason"])
        self.assertEqual(by_row[3]["reason"], "size")                              # 3 heavy atoms, below the window
        self.assertEqual(by_row[4]["reason"], "duplicate_after_standardization")
        self.assertTrue(by_row[5]["eligible"])
        self.assertEqual(by_row[0]["tautomer_status"], "Completed")

    def test_eligibility_is_independent_of_the_cache_state(self):
        smiles = [PYRIDONE, "CC(C)Cc1ccccc1"]
        cold = starts.eligible_rows(smiles, atoms(smiles), (5, 50))
        warm_cache = {}
        starts.eligible_rows(smiles, atoms(smiles), (5, 50), warm_cache)
        warm = starts.eligible_rows(smiles, atoms(smiles), (5, 50), warm_cache)
        self.assertEqual(cold, warm)


class SelectionTests(unittest.TestCase):
    SMILES = [f"C{'C' * i}O" for i in range(30)]

    def test_selection_is_deterministic_and_depends_only_on_seed_eligible_and_allowed(self):
        a = starts.select_start_population(range(30), range(30), self.SMILES, 10, seed=3)
        b = starts.select_start_population(list(reversed(range(30))), list(reversed(range(30))), self.SMILES, 10, seed=3)
        self.assertEqual(list(a[0]), list(b[0]))                                   # input order does not matter
        self.assertEqual(a[1], b[1])
        c = starts.select_start_population(range(30), range(30), self.SMILES, 10, seed=4)
        self.assertNotEqual(list(a[0]), list(c[0]))

    def test_allowed_set_restricts_the_draw(self):
        rows, _ = starts.select_start_population(range(30), range(10, 20), self.SMILES, 5, seed=1)
        self.assertTrue(all(10 <= r < 20 for r in rows))

    def test_never_silently_smaller(self):
        with self.assertRaisesRegex(ValueError, "only 4 eligible"):
            starts.select_start_population([1, 2, 3, 4, 5], [1, 2, 3, 4], self.SMILES, 5, seed=1)

    def test_paired_runs_get_identical_starts(self):
        legacy = starts.select_start_population(range(30), range(30), self.SMILES, 8, seed=42)
        corrected = starts.select_start_population(range(30), range(30), self.SMILES, 8, seed=42)
        self.assertEqual(legacy[1], corrected[1])


class LoadTests(unittest.TestCase):
    def table(self, keys):
        return pd.DataFrame({"row": range(len(keys)), "inchikey": keys, "eligible_primary": [True, False, True],
                             "eligible_relaxed": [True, True, True]})

    def test_alignment_is_verified_and_the_window_column_is_selected(self):
        curated = pd.DataFrame({"inchikey": ["a", "b", "c"]})
        tmp = Path(tempfile.mkdtemp())
        self.table(["a", "b", "c"]).to_csv(tmp / "eligible_starts.csv", index=False)
        with mock.patch.object(config, "PROCESSED_DIR", tmp):
            self.assertEqual(list(starts.load_eligible(curated)), [0, 2])
            self.assertEqual(list(starts.load_eligible(curated, relaxed=True)), [0, 1, 2])
            self.table(["a", "x", "c"]).to_csv(tmp / "eligible_starts.csv", index=False)
            with self.assertRaisesRegex(ValueError, "does not line up"):
                starts.load_eligible(curated)


if __name__ == "__main__":
    unittest.main()
