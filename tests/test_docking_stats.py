"""Tests for the docking statistics (Part 3 definitions). Run: python -m unittest discover tests"""

import unittest

import numpy as np

from cdk2moo import docking_stats as ds


class PSupTests(unittest.TestCase):
    def test_direction_ties_and_empty(self):
        self.assertEqual(ds.p_better([-9, -9], [-5, -5]), 1.0)        # more negative = better
        self.assertEqual(ds.p_better([-5, -5], [-9, -9]), 0.0)
        self.assertEqual(ds.p_better([-7], [-7]), 0.5)                 # a tie counts half
        self.assertTrue(np.isnan(ds.p_better([], [-7])))

    def test_stratification_removes_a_pure_size_effect(self):
        rng = np.random.default_rng(0)
        atoms_a = np.r_[np.full(40, 22), np.full(40, 32)]               # sample a is bigger on average than b
        atoms_b = np.r_[np.full(60, 22), np.full(20, 32)]
        score = lambda atoms: -7 - 0.1 * (atoms - 20) + rng.normal(0, 0.05, len(atoms))   # score depends on size only
        a, b = score(atoms_a), score(atoms_b)
        raw = ds.p_better(a, b)
        stratified, used, total = ds.stratified_p_better(a, atoms_a, b, atoms_b)
        self.assertGreater(raw, 0.6)                                    # raw comparison is fooled by size
        self.assertAlmostEqual(stratified, 0.5, delta=0.08)             # within strata there is no difference
        self.assertEqual(total, 6400)
        self.assertLessEqual(used, total)

    def test_strata_with_too_few_molecules_do_not_contribute(self):
        value, used, _ = ds.stratified_p_better([-8] * 3, [22] * 3, [-7] * 3, [22] * 3)
        self.assertTrue(np.isnan(value))
        self.assertEqual(used, 0)


class BootstrapTests(unittest.TestCase):
    def test_cluster_bootstrap_is_wider_than_naive_when_clusters_carry_the_effect(self):
        rng = np.random.default_rng(1)
        # five clusters with strong cluster-level offsets: effectively 5 replicates, not 100 molecules
        offsets = rng.normal(0, 1.0, 5)
        a = np.concatenate([o + rng.normal(0, 0.1, 20) for o in offsets])
        ca = np.repeat([f"s{i}" for i in range(5)], 20)
        b = rng.normal(0, 0.1, 100)
        cb = np.arange(100)                                              # reference molecules as their own clusters
        out = ds.compare_samples(a, np.full(100, 25), ca, b, np.full(100, 25), cb, n_boot=300)
        self.assertGreater(out["p_sup_high"] - out["p_sup_low"], 0.25)   # honest interval is wide
        self.assertEqual((out["clusters_a"], out["clusters_b"]), (5, 100))

    def test_bootstrap_is_reproducible(self):
        a, b = np.linspace(-9, -7, 40), np.linspace(-8, -6, 40)
        args = (a, np.full(40, 25), np.arange(40) % 5, b, np.full(40, 25), np.arange(40) % 4)
        self.assertEqual(ds.compare_samples(*args, n_boot=200), ds.compare_samples(*args, n_boot=200))


class GateTests(unittest.TestCase):
    def test_gate_logic(self):
        self.assertEqual(ds.validity_gate(0.55, 0.52, 0.66, 0.70), (True, "moderate"))
        self.assertEqual(ds.validity_gate(0.55, 0.52, 0.60, 0.70), (True, "weak"))
        self.assertEqual(ds.validity_gate(0.45, 0.52, 0.55, 0.70), (False, "failed"))     # a lower bound at or below 0.5 fails the gate
        self.assertEqual(ds.validity_gate(float("nan"), 0.6, 0.7, 0.7), (False, "undefined"))


if __name__ == "__main__":
    unittest.main()
