"""Tests for the per-run trajectory statistics (campaign plan section 3). Run: python -m unittest discover tests"""

import unittest

import numpy as np
import pandas as pd

from cdk2moo import endpoints as ep


def run_table(medians, n_generated=20):
    """A run whose generation g has `n_generated` generated molecules all equal to medians[g] (generation 0 = starting molecules)."""
    rows = []
    for g, m in enumerate(medians):
        births = np.zeros(n_generated, int) if g == 0 else np.full(n_generated, g)
        rows.append(pd.DataFrame({"generation": g, "birth_generation": births, "x": m}))
    return pd.concat(rows, ignore_index=True)


class TrendTests(unittest.TestCase):
    def test_monotone_series(self):
        stats = ep.trend_statistics(ep.generation_table(run_table([1, 2, 3, 4, 5]), "x"))
        self.assertEqual((stats["net_rise"], stats["fraction_nondecreasing"]), (4, 1.0))
        self.assertAlmostEqual(stats["spearman_generation"], 1.0)

    def test_a_single_dip_breaks_literal_monotonicity_but_not_the_trend(self):
        stats = ep.trend_statistics(ep.generation_table(run_table([1, 2, 3, 2.9, 4, 5, 6, 7]), "x"))
        self.assertLess(stats["fraction_nondecreasing"], 1.0)
        self.assertGreater(stats["spearman_generation"], 0.9)

    def test_generations_with_too_few_generated_molecules_do_not_qualify(self):
        run = run_table([1, 2, 3, 4])
        run = run[~((run["generation"] == 2) & (run.index % 20 >= 5))]                 # leaves 5 generated molecules in generation 2
        table = ep.generation_table(run, "x")
        self.assertEqual([r["qualified"] for r in table], [True, True, False, True])
        self.assertEqual(ep.trend_statistics(table)["n_qualified"], 3)                # reported, not silently dropped

    def test_survivors_from_generation_zero_are_not_counted_as_generated(self):
        run = run_table([1, 2, 3])
        run.loc[(run["generation"] == 2) & (run.index % 20 < 15), "birth_generation"] = 0   # 15 of 20 are starting molecules
        table = ep.generation_table(run, "x")
        self.assertEqual(table[2]["n_generated"], 5)
        self.assertFalse(table[2]["qualified"])


class AssociationTests(unittest.TestCase):
    def test_opposite_series_correlate_negatively(self):
        a = ep.generation_table(run_table([0, 1, 2, 3, 4, 5]), "x")
        b = ep.generation_table(run_table([9, 8, 7, 6, 5, 4]), "x")
        self.assertAlmostEqual(ep.association_across_generations(a, b)[0], -1.0)

    def test_generation_zero_is_excluded_by_default(self):
        a = ep.generation_table(run_table([100, 1, 2, 3, 4]), "x")
        b = ep.generation_table(run_table([0, 1, 2, 3, 4]), "x")
        rho, n = ep.association_across_generations(a, b)
        self.assertAlmostEqual(rho, 1.0)
        self.assertEqual(n, 4)


class ExchangeRateTests(unittest.TestCase):
    def test_defined_only_when_something_is_given_up(self):
        point, low, high, share = ep.bootstrap_ratio([10, 12, 8, 11, 9], [-0.5, -0.4, -0.6, -0.5, -0.5])
        self.assertAlmostEqual(point, 20.0, places=6)
        self.assertLessEqual(low, point)
        self.assertGreaterEqual(high, point)
        point, _, _, share = ep.bootstrap_ratio([10, 12, 8, 11, 9], [0.1, 0.2, 0.0, 0.1, 0.3])
        self.assertTrue(np.isnan(point))                                              # nothing given up: undefined, not infinite
        self.assertEqual(share, 0.0)


if __name__ == "__main__":
    unittest.main()
