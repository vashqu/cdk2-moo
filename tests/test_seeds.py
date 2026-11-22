"""Seed design of the campaign: 10 GA seeds, one label permutation per GA seed, 80%-of-runs decision rule."""
import unittest

import numpy as np

from cdk2moo import config
from cdk2moo.surrogate import scramble_labels


class SeedDesignTests(unittest.TestCase):
    def test_ten_ga_seeds_starting_at_42(self):
        self.assertEqual(config.GA_SEEDS, list(range(42, 52)))

    def test_first_seed_keeps_the_permutation_of_the_first_campaign(self):
        self.assertEqual(config.scramble_seed_for(42), 20260930)

    def test_every_ga_seed_has_its_own_permutation(self):
        seeds = {config.scramble_seed_for(s) for s in config.GA_SEEDS}
        self.assertEqual(len(seeds), len(config.GA_SEEDS))
        y = np.arange(200.0)
        permutations = {tuple(scramble_labels(y, s)[:10]) for s in seeds}
        self.assertEqual(len(permutations), len(seeds))

    def test_majority_rule_is_four_of_five_and_eight_of_ten(self):
        self.assertEqual(config.runs_needed(5), 4)
        self.assertEqual(config.runs_needed(10), 8)
        self.assertEqual(config.runs_needed(20), 16)


if __name__ == "__main__":
    unittest.main()
