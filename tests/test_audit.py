"""Tests for the read-only leakage audits (Task 8). Run: python -m unittest discover tests"""

import unittest

import numpy as np
import pandas as pd

from cdk2moo import audit

CURATED = pd.DataFrame({"inchikey": ["AAAAAAAAAAAAAA-BBBBBBBBBB-N", "AAAAAAAAAAAAAA-CCCCCCCCCC-N", "DDDDDDDDDDDDDD-EEEEEEEEEE-N",
                                     "FFFFFFFFFFFFFF-GGGGGGGGGG-N", "HHHHHHHHHHHHHH-IIIIIIIIII-N", "JJJJJJJJJJJJJJ-KKKKKKKKKK-N"]})


class AlignmentTests(unittest.TestCase):
    def splits(self):
        return pd.DataFrame({"row": np.arange(6), "inchikey": CURATED["inchikey"], "seed": 1})

    def test_reordered_split_rows_are_realigned_by_row(self):
        shuffled = self.splits().sample(frac=1, random_state=0)
        ordered = audit.check_alignment(shuffled, CURATED)
        self.assertEqual(ordered["row"].tolist(), list(range(6)))
        self.assertEqual(ordered["inchikey"].tolist(), CURATED["inchikey"].tolist())

    def test_misaligned_identity_is_detected(self):
        bad = self.splits()
        bad["inchikey"] = bad["inchikey"].iloc[::-1].to_numpy()          # keys reordered independently of rows
        with self.assertRaisesRegex(ValueError, "misaligned"):
            audit.check_alignment(bad, CURATED)

    def test_missing_identity_keys_give_a_clear_diagnostic(self):
        with self.assertRaisesRegex(ValueError, "lacks the identity column 'inchikey'"):
            audit.check_alignment(self.splits().drop(columns="inchikey"), CURATED)
        broken = self.splits()
        broken.loc[2, "inchikey"] = np.nan
        with self.assertRaisesRegex(ValueError, "missing inchikey"):
            audit.check_alignment(broken, CURATED)

    def test_incomplete_row_coverage_is_detected(self):
        with self.assertRaisesRegex(ValueError, "cover 0..n-1"):
            audit.check_alignment(self.splits().iloc[:5], CURATED)


class OverlapTests(unittest.TestCase):
    def test_identity_overlap_full_key_versus_connectivity_block(self):
        keys = CURATED["inchikey"].to_numpy()
        # row 0 and row 1 are stereoisomers: same connectivity block, different full key
        self.assertEqual(audit.key_overlap([0, 2], [1, 3], keys), (0, 0))                   # full identity: none shared
        self.assertEqual(audit.key_overlap([0, 2], [1, 3], audit.connectivity_keys(keys)), (1, 1))   # connectivity: row 1 matches row 0
        self.assertEqual(audit.key_overlap([0, 2], [0, 3], keys), (1, 1))                   # a true duplicate

    def test_fingerprint_overlap_counts_identical_bit_vectors(self):
        fps = np.zeros((4, 16), dtype=np.uint8)
        fps[0, [1, 5]] = 1
        fps[1, [1, 5]] = 1            # identical to row 0
        fps[2, [2, 3]] = 1
        fps[3, [1, 5, 9]] = 1
        keys = audit.fingerprint_keys(fps)
        self.assertEqual(audit.key_overlap([0], [1, 2, 3], keys), (1, 1))

    def test_scaffold_overlap(self):
        scaffolds = np.array(["c1ccccc1", "c1ccccc1", "C1CCCCC1", "c1ccncc1"], dtype=object)
        self.assertEqual(audit.key_overlap([0, 2], [1, 3], scaffolds), (1, 1))

    def test_any_document_overlap_is_stricter_than_primary_document(self):
        documents = ["D1", "D2|D3", "D3", "D4"]
        primary = ["D1", "D2", "D3", "D4"]
        out = audit.document_overlap([0, 1], [2, 3], documents, primary)
        self.assertEqual(out["any_document"], 1)           # molecule 2 shares D3 with molecule 1's second document
        self.assertEqual(out["primary_document"], 0)       # its primary document D3 is not a primary document in training

    def test_input_order_does_not_change_overlap_counts(self):
        keys = np.array(["a", "b", "a", "c", "b"], dtype=object)
        train, test = [0, 1], [2, 3, 4]
        a = audit.key_overlap(train, test, keys)
        b = audit.key_overlap(train[::-1], test[::-1], keys)
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
