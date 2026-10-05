"""The repaired GA, under the legacy policy, must reproduce the saved HISTORICAL populations (generations 0-3, seed 42).
It reads the historical files by explicit path, since pipeline paths are campaign-scoped and refuse to resolve without a selection."""

import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.ga import run_ga
from cdk2moo.surrogate import fit_forest, scramble_labels

ROOT = Path(__file__).resolve().parents[1]
MAIN = ROOT / "results" / "05_ga_populations.csv"
CONSTRAINED = ROOT / "results" / "10_constrained_primary_populations.csv"


@unittest.skipUnless(MAIN.exists() and CONSTRAINED.exists(), "historical populations not present")
class LegacyEquivalence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        df = pd.read_csv(ROOT / "data" / "processed" / "cdk2_ic50_curated.csv")      # the HISTORICAL files, selected by explicit path
        cls.fps = np.load(ROOT / "data" / "processed" / "ecfp4.npy")
        y = df["pactivity"].to_numpy()
        cls.real = fit_forest(cls.fps, y, config.SURROGATE_SEED)
        cls.scrambled = fit_forest(cls.fps, scramble_labels(y, config.SCRAMBLE_SEED), config.SCRAMBLE_SEED)
        sizes = df["n_heavy_atoms"].to_numpy()
        pool = np.where((sizes >= config.GA_MIN_HEAVY_ATOMS) & (sizes <= config.GA_MAX_HEAVY_ATOMS))[0]
        rng = np.random.default_rng(42)
        cls.start = df["std_smiles"].iloc[rng.choice(pool, config.GA_POP_SIZE, replace=False)]
        cls.size_range = (config.GA_MIN_HEAVY_ATOMS, config.GA_MAX_HEAVY_ATOMS)

    def compare(self, new, old):
        for generation in range(0, 4):
            a = new[new["generation"] == generation].sort_values("smiles")
            b = old[old["generation"] == generation].sort_values("smiles")
            self.assertEqual(a["smiles"].tolist(), b["smiles"].tolist(), f"generation {generation}")
            np.testing.assert_allclose(a["fitness"].to_numpy(), b["fitness"].to_numpy(), atol=2e-5)

    def test_unconstrained_run_matches_saved_population(self):
        old = pd.read_csv(MAIN, usecols=["smiles", "fitness", "generation", "arm", "seed"])
        old = old[(old["arm"] == "multi_real") & (old["seed"] == 42)]
        new = run_ga(self.start, "multi_real", self.real, self.scrambled, self.fps, 42, self.size_range, n_generations=3)
        self.compare(new, old)

    def test_constrained_run_matches_saved_population(self):
        old = pd.read_csv(CONSTRAINED, usecols=["smiles", "fitness", "generation", "arm", "seed", "min_similarity"])
        old = old[(old["arm"] == "multi_real") & (old["seed"] == 42) & np.isclose(old["min_similarity"], 0.6)]
        new = run_ga(self.start, "multi_real", self.real, self.scrambled, self.fps, 42, self.size_range,
                     n_generations=3, min_similarity=0.6)
        self.compare(new, old)
        self.assertTrue((new["max_tanimoto"] >= 0.6).all())


if __name__ == "__main__":
    unittest.main()
