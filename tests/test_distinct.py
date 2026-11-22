"""The distinct-standardized-molecules guard (src/cdk2moo/distinct.py)."""
import unittest

import numpy as np
import pandas as pd

from cdk2moo.distinct import require_distinct_standardized


def table(**override):
    rows = pd.DataFrame({"smiles": ["CCO", "CCN", "CCC", "CCO"], "arm": "a", "seed": 42, "generation": [0, 0, 0, 1], "policy": "corrected"})
    for k, v in override.items():
        rows[k] = v
    return rows


class DistinctTests(unittest.TestCase):
    def test_a_corrected_table_with_unique_populations_passes_and_the_same_smiles_in_different_generations_is_fine(self):
        self.assertEqual(require_distinct_standardized(table())["populations"], 2)

    def test_a_repeated_smiles_within_one_population_is_refused(self):
        with self.assertRaisesRegex(ValueError, "repeat a SMILES"):
            require_distinct_standardized(table(generation=[0, 0, 0, 0]))

    def test_legacy_rows_are_refused(self):
        with self.assertRaisesRegex(ValueError, "legacy"):
            require_distinct_standardized(table(policy="legacy"))

    def test_a_missing_policy_column_is_refused(self):
        with self.assertRaisesRegex(ValueError, "no `policy` column"):
            require_distinct_standardized(table().drop(columns="policy"))

    def test_constrained_floors_are_separate_populations(self):
        t = pd.concat([table(min_similarity=0.4), table(min_similarity=0.6)], ignore_index=True)
        self.assertEqual(require_distinct_standardized(t)["populations"], 4)
        t = table(min_similarity=np.nan)                      # unconstrained rows carry NaN
        self.assertEqual(require_distinct_standardized(t)["populations"], 2)


if __name__ == "__main__":
    unittest.main()
