"""Tests for the preparation policies and the closed-shell rule (Tasks 5 and 6)."""

import itertools
import unittest
from unittest import mock

import numpy as np
import pandas as pd
from rdkit import Chem

from cdk2moo import config, ga
from cdk2moo.candidates import prepare_candidate
from cdk2moo.mutations import count_radicals, is_valid, rejection_reason
from cdk2moo.standardize import standardize_with_status

SIZE = (1, 50)
PYRIDINOL, PYRIDONE = "Oc1ccccn1", "O=c1cccc[nH]1"        # a tautomer pair: same compound, different SMILES


def mol(smiles):
    return Chem.MolFromSmiles(smiles)


class StandardizationPolicyTests(unittest.TestCase):
    def test_tautomer_pair_maps_to_the_same_scored_input_under_corrected_policy(self):
        a = prepare_candidate(mol(PYRIDINOL), "corrected", SIZE)
        b = prepare_candidate(mol(PYRIDONE), "corrected", SIZE)
        self.assertIsNone(a[1])
        self.assertEqual(a[0], b[0])

    def test_legacy_and_corrected_policies_are_distinguishable(self):
        legacy = {prepare_candidate(mol(s), "legacy", SIZE)[0] for s in (PYRIDINOL, PYRIDONE)}
        corrected = {prepare_candidate(mol(s), "corrected", SIZE)[0] for s in (PYRIDINOL, PYRIDONE)}
        self.assertEqual(len(legacy), 2)                  # legacy scores the two tautomers as different molecules
        self.assertEqual(len(corrected), 1)

    def test_tautomer_enumeration_status_is_reported(self):
        out, reason, status = standardize_with_status(mol(PYRIDINOL))
        self.assertIsNone(reason)
        self.assertEqual(status, "Completed")

    def test_incomplete_enumeration_is_counted_not_hidden(self):
        rejects = {}
        with mock.patch("cdk2moo.candidates.standardize_with_status", return_value=(mol("CCO"), None, "MaxTransformsReached")):
            smiles, reason, status = prepare_candidate(mol("CCO"), "corrected", SIZE, rejects=rejects)
        self.assertEqual((smiles, status), ("CCO", "MaxTransformsReached"))
        self.assertEqual(rejects["tautomer_incomplete"], 1)

    def test_post_standardization_size_is_enforced(self):
        # the sodium salt's organic part, after neutralising, is acetic acid: 4 heavy atoms
        smiles, reason, _ = prepare_candidate(mol("CC(=O)[O-].[Na+]"), "corrected", (5, 50))
        self.assertIsNone(smiles)
        self.assertEqual(reason, "after_standardization_size")

    def test_invalid_structures_are_rejected_with_a_reason(self):
        smiles, reason, _ = prepare_candidate(mol("C[Sn](C)C"), "corrected", SIZE)
        self.assertIsNone(smiles)
        self.assertIn("inorganic", reason)

    def test_final_string_rebuilds_the_same_molecule(self):
        smiles, _, _ = prepare_candidate(mol(PYRIDONE), "corrected", SIZE)
        self.assertEqual(Chem.MolToSmiles(mol(smiles)), smiles)
        again, _, _ = prepare_candidate(mol(smiles), "corrected", SIZE)
        self.assertEqual(again, smiles)                   # idempotent

    def test_unknown_policy_is_refused(self):
        with self.assertRaises(ValueError):
            prepare_candidate(mol("CCO"), "other", SIZE)


class ProposalPipelineTests(unittest.TestCase):
    def population(self):
        return pd.DataFrame({"smiles": ["CCO", "CCN", "CCC"], "fitness": [0.5, 0.5, 0.5]})

    def test_two_proposals_that_standardize_to_one_structure_are_deduplicated(self):
        proposals = itertools.cycle([mol(PYRIDINOL), mol(PYRIDONE), mol("CCCCC")])
        with mock.patch.object(config, "GA_N_CHILDREN", 2), \
             mock.patch.object(ga, "crossover", return_value=None), \
             mock.patch.object(ga, "mutate", side_effect=lambda *a, **k: next(proposals)):
            children = ga._propose_children(self.population(), self.population()["smiles"], np.random.default_rng(0),
                                            SIZE, policy="corrected", cache={}, rejects={})
        self.assertEqual(len(children), 2)
        self.assertEqual(len({c["smiles"] for c in children}), 2)
        taut = [c for c in children if c["proposal_smiles"] in (Chem.MolToSmiles(mol(PYRIDINOL)), Chem.MolToSmiles(mol(PYRIDONE)))]
        self.assertEqual(len(taut), 1)                    # one tautomer kept; its proposal string is preserved for tracing

    def test_stored_smiles_and_scored_smiles_agree(self):
        scored = []

        def scorer(smiles_list, real, scrambled, fps):
            scored.extend(smiles_list)
            return pd.DataFrame({"smiles": smiles_list, "pred_real": 6.0, "pred_scrambled": 6.0, "qed": 0.5, "sa": 3.0,
                                 "max_tanimoto": 0.9, "n_heavy_atoms": 5})

        proposals = itertools.cycle([mol(PYRIDINOL), mol("CCCCC"), mol("CCCCCC")])
        with mock.patch.object(config, "GA_POP_SIZE", 3), mock.patch.object(config, "GA_N_CHILDREN", 2), \
             mock.patch.object(ga, "score_molecules", scorer), mock.patch.object(ga, "crossover", return_value=None), \
             mock.patch.object(ga, "mutate", side_effect=lambda *a, **k: next(proposals)):
            history = ga.run_ga(["CCO", "CCN", "CCC"], "multi_real", None, None, None, 0, SIZE, n_generations=2, policy="corrected")
        self.assertTrue(set(history["smiles"]) <= set(scored))      # every stored string is exactly a string that was scored
        for stored in set(history["smiles"]):                        # and it is already in its standardized form
            self.assertEqual(prepare_candidate(mol(stored), "corrected", SIZE)[0], stored)
        self.assertTrue((history["policy"] == "corrected").all())
        self.assertIn("proposal_smiles", history.columns)          # provenance of the proposal is kept

    def test_canonicalization_cannot_silently_shrink_the_initial_population(self):
        with self.assertRaisesRegex(ValueError, "merged 1 starting molecules into duplicates"):
            ga._prepare_start([PYRIDINOL, PYRIDONE, "CCC"], "corrected", SIZE, {}, {})

    def test_failing_start_molecule_raises(self):
        with self.assertRaisesRegex(ValueError, "fail the 'corrected' policy"):
            ga._prepare_start(["CCO", "C[Sn](C)C"], "corrected", SIZE, {}, {})

    def test_legacy_start_is_plain_canonicalization(self):
        self.assertEqual(ga._prepare_start([PYRIDINOL, PYRIDONE], "legacy", SIZE, {}, {}),
                         [Chem.MolToSmiles(mol(PYRIDINOL)), Chem.MolToSmiles(mol(PYRIDONE))])


class RadicalTests(unittest.TestCase):
    def test_sanitizable_radical_is_identified_and_rejected_by_closed_shell_policy(self):
        radical = mol("C[CH2]")                             # ethyl radical: sanitizes, one atom with an unpaired electron
        self.assertIsNotNone(radical)
        self.assertEqual(count_radicals(radical), 1)
        self.assertEqual(rejection_reason(radical, SIZE, closed_shell=True), "radical")

    def test_legacy_policy_still_admits_the_radical_and_is_identifiable(self):
        radical = mol("C[CH2]")
        self.assertIsNone(rejection_reason(radical, SIZE, closed_shell=False))
        self.assertTrue(is_valid(radical, SIZE))             # default is legacy
        self.assertFalse(is_valid(radical, SIZE, closed_shell=True))

    def test_ordinary_closed_shell_molecules_stay_accepted(self):
        for smiles in ("CCO", "c1ccccc1O", "CC(=O)Nc1ccc(O)cc1"):
            with self.subTest(smiles=smiles):
                self.assertEqual(count_radicals(mol(smiles)), 0)
                self.assertIsNone(rejection_reason(mol(smiles), SIZE, closed_shell=True))

    def test_rejection_reasons_are_distinct(self):
        self.assertEqual(rejection_reason(mol("CCO.CCN"), SIZE), "disconnected")
        self.assertEqual(rejection_reason(mol("CCO"), (5, 50)), "size")
        self.assertEqual(rejection_reason(None, SIZE), "no_molecule")


if __name__ == "__main__":
    unittest.main()
