"""Tests for positive-control evaluation semantics (Task 7). Run: python -m unittest discover tests"""

import unittest

import numpy as np

from cdk2moo import recovery as rc

L_ALA, D_ALA = "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"     # enantiomers; the tautomer step strips the alpha-carbon stereo
R_BUT, S_BUT = "C[C@H](O)CC", "C[C@@H](O)CC"              # enantiomers whose stereo the standardization keeps
PYRIDINOL, PYRIDONE = "Oc1ccccn1", "O=c1cccc[nH]1"        # tautomers of one compound
KW = dict(thresholds=[0.5, 0.6], annotation_cutoff=8.0)


def reference(smiles, activity=None):
    return rc.build_reference(smiles, activity if activity is not None else [9.0] * len(smiles))


class IdentityTests(unittest.TestCase):
    def test_stereoisomers_with_identical_fingerprints_are_not_automatic_exact_matches(self):
        out = rc.evaluate_group([S_BUT], reference([R_BUT, "CCCCCCO"]), held_idx=[0], **KW)
        self.assertEqual(out["fingerprint_identical_weighted"], 1.0)      # ECFP4 ignores chirality: fingerprints collide ...
        self.assertEqual(out["exact_std_isomeric_weighted"], 0.0)         # ... but these are different molecules
        self.assertEqual(out["exact_strict_stereo_weighted"], 0.0)
        self.assertEqual(out["exact_std_flat_weighted"], 1.0)             # same connectivity only

    def test_tautomer_step_merges_alpha_carbon_stereo_and_the_strict_key_does_not(self):
        out = rc.evaluate_group([D_ALA], reference([L_ALA, "CCCCCCO"]), held_idx=[0], **KW)
        self.assertEqual(out["exact_std_isomeric_weighted"], 1.0)         # documented RDKit behaviour: stereo removed there
        self.assertEqual(out["exact_strict_stereo_weighted"], 0.0)        # stereo-strict key tells them apart
        self.assertEqual(out["fingerprint_identical_weighted"], 1.0)

    def test_same_stereoisomer_is_an_exact_match_under_every_key(self):
        out = rc.evaluate_group([S_BUT], reference([S_BUT, "CCCCCCO"]), held_idx=[0], **KW)
        for key in ("exact_std_isomeric_weighted", "exact_std_flat_weighted", "exact_strict_stereo_weighted"):
            self.assertEqual(out[key], 1.0)

    def test_known_exact_standardized_match_is_counted(self):
        out = rc.evaluate_group([PYRIDINOL], reference([PYRIDONE, "CCCCO"]), held_idx=[0], **KW)
        self.assertEqual(out["exact_std_isomeric_weighted"], 1.0)         # tautomers standardize to one structure
        self.assertEqual(out["exact_strict_stereo_weighted"], 0.0)        # the strict key does not merge tautomers
        self.assertLess(out["fingerprint_identical_weighted"], 1.0)       # raw fingerprints differ: only the keys match

    def test_unrelated_molecule_is_not_exact(self):
        out = rc.evaluate_group(["CCCCCC"], reference([PYRIDONE, "CCCCO"]), held_idx=[0], **KW)
        self.assertEqual((out["exact_std_isomeric_weighted"], out["exact_std_flat_weighted"]), (0.0, 0.0))


class NearestNeighbourTests(unittest.TestCase):
    # reference: [held-out active H, other training molecule X]; X is also the query (a starting molecule)
    H, X = "c1ccc2[nH]ccc2c1", "c1ccc2occc2c1"

    def test_starting_molecule_does_not_win_by_matching_itself(self):
        ref = reference([self.H, self.X])
        with_self = rc.evaluate_group([self.X], ref, held_idx=[0], **KW)
        without_self = rc.evaluate_group([self.X], ref, held_idx=[0], exclude_self=True, **KW)
        self.assertEqual(with_self["nearest_is_heldout_weighted"], 0.0)        # matches itself
        self.assertEqual(without_self["nearest_is_heldout_weighted"], 1.0)     # own identity excluded: the only neighbour left

    def test_generated_molecule_that_reproduces_a_reference_is_not_excluded(self):
        out = rc.evaluate_group([self.H], reference([self.H, self.X]), held_idx=[0], **KW)   # default exclude_self=False
        self.assertEqual(out["exact_std_isomeric_weighted"], 1.0)
        self.assertEqual(out["nearest_is_heldout_weighted"], 1.0)

    def test_query_order_does_not_change_identity_mapping(self):
        ref = reference([L_ALA, PYRIDONE, "CCCCO", self.X])
        queries = [D_ALA, PYRIDINOL, "CCCCCCC", self.X]
        a = rc.evaluate_group(queries, ref, held_idx=[0, 1], **KW)
        b = rc.evaluate_group(queries[::-1], ref, held_idx=[0, 1], **KW)
        for key in ("exact_std_isomeric_weighted", "exact_std_flat_weighted", "nearest_is_heldout_weighted", "proximity_ge_0.6_weighted"):
            self.assertAlmostEqual(a[key], b[key])
        sim = rc.tanimoto_matrix(rc.ecfp4(queries).astype(np.float32), ref["fps"])
        idx_forward = rc.nearest_reference(sim)[0]
        idx_reversed = rc.nearest_reference(sim[::-1])[0]
        self.assertEqual(list(idx_forward), list(idx_reversed[::-1]))

    def test_duplicates_have_explicit_weighted_and_distinct_semantics(self):
        ref = reference([L_ALA, "CCCCO"])
        out = rc.evaluate_group([L_ALA, L_ALA, "CCCCCCC"], ref, held_idx=[0], **KW)
        self.assertEqual((out["n_weighted"], out["n_distinct"]), (3, 2))
        self.assertAlmostEqual(out["exact_std_isomeric_weighted"], 2 / 3)       # per occurrence
        self.assertAlmostEqual(out["exact_std_isomeric_distinct"], 1 / 2)       # per distinct structure

    def test_ties_use_the_lowest_reference_index_and_are_reported(self):
        # two reference rows with identical fingerprints (a stereoisomer pair); the lower index is not held-out
        ref = reference(["CCCCO", L_ALA, D_ALA])
        out = rc.evaluate_group([L_ALA], ref, held_idx=[2], exclude_self=False, thresholds=[0.6], annotation_cutoff=8.0)
        sim = rc.tanimoto_matrix(rc.ecfp4([L_ALA]).astype(np.float32), ref["fps"])
        index, _, tied, _ = rc.nearest_reference(sim)
        self.assertEqual(int(index[0]), 1)                                  # lowest of the tied rows 1 and 2
        self.assertEqual(int(tied.sum()), 2)
        self.assertEqual(out["nearest_is_heldout_weighted"], 0.0)           # tie-break picks row 1: not held-out
        self.assertEqual(out["any_tied_nearest_is_heldout_weighted"], 1.0)  # but a tied neighbour is
        self.assertEqual(out["nearest_tie_share_weighted"], 1.0)

    def test_reference_prevalence_is_a_descriptive_property_not_a_null(self):
        out = rc.evaluate_group(["CCCCO"], reference([L_ALA, "CCCCO", "CCN", "CCC"]), held_idx=[0], **KW)
        self.assertAlmostEqual(out["reference_prevalence_of_heldout"], 0.25)
        self.assertNotIn("chance", " ".join(out))

    def test_recall_counts_distinct_generated_molecules_once(self):
        ref = reference([L_ALA, "CCCCO"])
        share, hits, n_distinct = rc.recall_of_held_actives([L_ALA, L_ALA, "CCCCCCC"], ref, [0], 0.6)
        self.assertEqual((share, hits, n_distinct), (1.0, 1, 2))


class DistinctMetricOrderInvarianceTests(unittest.TestCase):
    ACETONE, ENOL = "CC(C)=O", "CC(O)=C"                    # a ketone/enol pair that standardizes to one structure
    REF = ["CCC(C)=O", "CCCCCCO", "c1ccccc1O", "CCN"]

    def ref(self):
        return rc.build_reference(self.REF, [9.0, 5.0, 6.0, 5.0])

    def test_the_pair_is_one_standardized_identity(self):
        keys = rc.standardized_keys([self.ACETONE, self.ENOL])
        self.assertEqual(keys[0][0], keys[1][0])

    def test_distinct_metrics_are_identical_in_either_query_order(self):
        for first, second in [(self.ACETONE, self.ENOL), (self.ENOL, self.ACETONE)]:
            out = rc.evaluate_group([first, second, "CCCCCCC"], self.ref(), [0], thresholds=[0.5, 0.6], annotation_cutoff=8.0)
            self.assertEqual(out["n_distinct"], 2)
        a = rc.evaluate_group([self.ACETONE, self.ENOL, "CCCCCCC"], self.ref(), [0], thresholds=[0.3, 0.6], annotation_cutoff=8.0)
        b = rc.evaluate_group([self.ENOL, self.ACETONE, "CCCCCCC"], self.ref(), [0], thresholds=[0.3, 0.6], annotation_cutoff=8.0)
        for key in a:
            if key.endswith("_distinct") or key.startswith("n_distinct") or "_distinct_" in key:
                np.testing.assert_equal(a[key], b[key], err_msg=key)

    def test_distinct_metrics_use_the_standardized_representation_not_the_first_raw_form(self):
        ref = self.ref()
        enol_first = rc.evaluate_group([self.ENOL], ref, [0], thresholds=[0.6], annotation_cutoff=8.0)
        ketone_first = rc.evaluate_group([self.ACETONE], ref, [0], thresholds=[0.6], annotation_cutoff=8.0)
        # the raw forms have different fingerprints (weighted metric differs) but one identity (distinct metric identical)
        self.assertNotEqual(enol_first["median_proximity_to_held_actives_weighted"], ketone_first["median_proximity_to_held_actives_weighted"])
        self.assertEqual(enol_first["median_proximity_to_held_actives_distinct"], ketone_first["median_proximity_to_held_actives_distinct"])

    def test_shuffled_queries_give_identical_distinct_output(self):
        queries = [self.ACETONE, self.ENOL, L_ALA, D_ALA, "CCCCCCC", self.ACETONE, "CCN"]
        base = rc.evaluate_group(queries, self.ref(), [0, 1], thresholds=[0.5], annotation_cutoff=8.0)
        for seed in range(4):
            shuffled = list(np.random.default_rng(seed).permutation(queries))
            out = rc.evaluate_group(shuffled, self.ref(), [0, 1], thresholds=[0.5], annotation_cutoff=8.0)
            for key in base:
                if "weighted" not in key:
                    np.testing.assert_equal(out[key], base[key], err_msg=key)

    def test_duplicates_change_only_the_weighted_metrics(self):
        ref = self.ref()
        once = rc.evaluate_group([self.ACETONE, "c1ccccc1"], ref, [0], thresholds=[0.1], annotation_cutoff=8.0)
        thrice = rc.evaluate_group([self.ACETONE] * 3 + ["c1ccccc1"], ref, [0], thresholds=[0.1], annotation_cutoff=8.0)
        for key in once:
            if key.endswith("_distinct"):
                np.testing.assert_equal(once[key], thrice[key], err_msg=key)
        self.assertEqual((once["proximity_ge_0.1_weighted"], thrice["proximity_ge_0.1_weighted"]), (0.5, 0.75))   # weighting counts occurrences

    def test_specified_stereoisomers_stay_distinct_identities(self):
        out = rc.evaluate_group([R_BUT, S_BUT], self.ref(), [0], thresholds=[0.6], annotation_cutoff=8.0)
        self.assertEqual(out["n_distinct"], 2)

    def test_empty_group_is_defined_not_an_error(self):
        out = rc.evaluate_group([], self.ref(), [0], thresholds=[0.6], annotation_cutoff=8.0)
        self.assertEqual((out["status"], out["n_weighted"], out["n_distinct"]), ("empty", 0, 0))
        self.assertEqual(rc.recall_of_held_actives([], self.ref(), [0], 0.6)[1:], (0, 0))
        self.assertEqual(rc.proximity_to_held([], self.ref(), [0]).size, 0)

    def test_exhausted_reference_after_self_exclusion_is_reported_not_hidden(self):
        ref = rc.build_reference(["CCCCO"], [9.0])
        out = rc.evaluate_group(["CCCCO"], ref, [0], thresholds=[0.6], annotation_cutoff=8.0, exclude_self=True)
        self.assertEqual(out["no_reference_neighbour_weighted"], 1)
        self.assertEqual(out["n_with_reference_weighted"], 0)
        self.assertTrue(np.isnan(out["nearest_is_heldout_weighted"]))

    def test_standardization_failures_and_incomplete_enumeration_are_reported(self):
        out = rc.evaluate_group(["CCCCCCC", "C[Sn](C)C"], self.ref(), [0], thresholds=[0.6], annotation_cutoff=8.0)
        self.assertEqual((out["n_standardization_failed"], out["n_distinct_standardization_failed"]), (1, 1))
        self.assertIn("n_tautomer_incomplete", out)
        self.assertEqual(out["n_distinct"], 2)                             # the failure is its own identity, still counted

    def test_recall_does_not_depend_on_generated_order(self):
        ref = self.ref()
        a = rc.recall_of_held_actives([self.ACETONE, self.ENOL, "CCCCCCC"], ref, [0, 1], 0.5)
        b = rc.recall_of_held_actives(["CCCCCCC", self.ENOL, self.ACETONE], ref, [0, 1], 0.5)
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
