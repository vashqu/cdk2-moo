"""Tests for the similarity floor as a feasibility constraint (Task 4). Scoring and proposals are stubbed."""

import unittest
from unittest import mock

import numpy as np
import pandas as pd

from cdk2moo import config, ga

SIZE = (1, 50)


def stub_scorer(table):
    """score_molecules replacement: smiles -> predicted activity and similarity from a lookup."""
    def score(smiles_list, real, scrambled, fps):
        return pd.DataFrame([{"smiles": s, "pred_real": table[s]["act"], "pred_scrambled": 6.0, "qed": 0.5, "sa": 3.0,
                              "max_tanimoto": table[s]["sim"], "n_heavy_atoms": 5} for s in smiles_list])
    return score


def child(smiles, parent="CCO"):
    return {"smiles": smiles, "proposal_smiles": smiles, "origin": "mutation", "parent_a": parent, "parent_b": "",
            "tautomer_status": "not_applicable"}


def run(table, starts, batches, floor, generations=None):
    """Run the GA with population size 3 and the given child batches (one list per generation)."""
    proposals = iter([[child(s) for s in batch] for batch in batches])
    with mock.patch.object(config, "GA_POP_SIZE", 3), \
         mock.patch.object(ga, "score_molecules", stub_scorer(table)), \
         mock.patch.object(ga, "_propose_children", side_effect=lambda *a, **k: next(proposals)):
        return ga.run_ga(starts, "multi_real", None, None, None, 1, SIZE,
                         n_generations=generations or len(batches), min_similarity=floor)


STARTS = ["CCO", "CCN", "CCC"]


class FloorTests(unittest.TestCase):
    def table(self, **extra):
        base = {s: {"act": 6.0, "sim": 0.9} for s in STARTS}
        base.update(extra)
        return base

    def test_no_survivor_lies_below_the_threshold(self):
        table = self.table(CCCC={"act": 9.0, "sim": 0.3}, CCF={"act": 7.0, "sim": 0.8})
        history = run(table, STARTS, [["CCCC", "CCF"]], floor=0.6)
        self.assertTrue((history["max_tanimoto"] >= 0.6).all())

    def test_infeasible_molecule_with_excellent_activity_is_excluded(self):
        table = self.table(CCCC={"act": 9.0, "sim": 0.3})
        history = run(table, STARTS, [["CCCC"]], floor=0.6)
        self.assertNotIn("CCCC", set(history["smiles"]))
        self.assertEqual(history.attrs["rejects"]["infeasible_similarity"], 1)

    def test_floor_off_keeps_legacy_behaviour(self):
        table = self.table(CCCC={"act": 9.0, "sim": 0.3})
        history = run(table, STARTS, [["CCCC"]], floor=None)
        self.assertIn("CCCC", set(history[history["generation"] == 1]["smiles"]))

    def test_feasible_zero_fitness_parents_survive_when_children_are_infeasible(self):
        table = {s: {"act": 3.0, "sim": 0.9} for s in STARTS}               # activity below the scale: fitness exactly 0
        table.update({"CCCC": {"act": 9.0, "sim": 0.2}, "CCF": {"act": 9.0, "sim": 0.1}, "CCCl": {"act": 9.0, "sim": 0.3}})
        history = run(table, STARTS, [["CCCC", "CCF", "CCCl"]], floor=0.6)
        final = history[history["generation"] == 1]
        self.assertEqual(sorted(final["smiles"]), sorted(STARTS))             # zero-fitness feasible parents kept
        self.assertTrue((final["fitness"] == 0).all())

    def test_all_infeasible_child_batch_leaves_valid_parents(self):
        table = self.table(CCCC={"act": 9.0, "sim": 0.2}, CCF={"act": 9.0, "sim": 0.2})
        history = run(table, STARTS, [["CCCC", "CCF"], ["CCCC", "CCF"]], floor=0.6)
        for generation in (1, 2):
            self.assertEqual(sorted(history[history["generation"] == generation]["smiles"]), sorted(STARTS))

    def test_nonfinite_child_similarity_is_infeasible_under_a_floor(self):
        table = self.table(CCCC={"act": 9.0, "sim": float("nan")})
        history = run(table, STARTS, [["CCCC"]], floor=0.6)
        self.assertNotIn("CCCC", set(history["smiles"]))

    def test_invalid_initial_population_raises_an_informative_error(self):
        table = self.table(CCF={"act": 6.0, "sim": 0.3})
        with self.assertRaisesRegex(ValueError, "1 of 3 starting molecules violate the similarity floor 0.6"):
            run(table, ["CCO", "CCN", "CCF"], [[]], floor=0.6)

    def test_nan_similarity_in_the_initial_population_raises(self):
        table = self.table(CCF={"act": 6.0, "sim": float("nan")})
        with self.assertRaisesRegex(ValueError, "violate the similarity floor"):
            run(table, ["CCO", "CCN", "CCF"], [[]], floor=0.6)

    def test_invalid_threshold_values_are_rejected(self):
        for bad in (-0.1, 1.5, float("nan"), float("inf"), "0.5", True):
            with self.subTest(bad=bad), self.assertRaisesRegex(ValueError, "min_similarity"):
                run(self.table(), STARTS, [[]], floor=bad)

    def test_threshold_boundaries_are_allowed(self):
        for ok in (0.0, 1.0):                               # molecules at similarity exactly 1.0 satisfy a floor of 1.0
            run({s: {"act": 6.0, "sim": 1.0} for s in STARTS}, STARTS, [[]], floor=ok)


if __name__ == "__main__":
    unittest.main()
