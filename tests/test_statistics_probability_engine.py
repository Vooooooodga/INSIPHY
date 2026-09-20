import math

import tempfile

import unittest

from pathlib import Path

from types import SimpleNamespace

from unittest.mock import patch

import numpy as np

from scipy.optimize import brentq as scipy_brentq

from intraphy.io import (
    read_structural_site_matrix,
    read_tsv,
    structural_site_observed_state,
    write_structural_site_matrix,
)

from intraphy.parsimony import (
    _event_type,
    _parsimony_tables,
    _placement_status,
    _structural_relation,
    infer_single_copy_parsimony,
)

from intraphy.structural_phylogeny import (
    PROFILE_DROP_95,
    _ascertainment_log_probability,
    _dataset_log_likelihood,
    _expected_transition_count,
    _fit_valid_for_posterior,
    _pattern_log_likelihood,
    _posterior_messages,
    _profile_contrast_interval,
    _profile_interval,
    _selected_for_ascertainment,
    _transition_matrix,
    _transition_event_type,
    _transition_structural_relation,
    _validated_tree_rows,
    fit_model,
    infer_single_copy_phylogeny,
)

from intraphy.tree import SpeciesTree

from intraphy.structural_sites import structural_matrix_annotation_view

from support_statistics_probability_engine import ProbabilityEngineTestsSupport

class ProbabilityEngineTests(ProbabilityEngineTestsSupport, unittest.TestCase):
    def test_pruning_and_posteriors_match_explicit_enumeration(self):
        tree = self._small_tree()
        observations = {"A": 0, "B": 1, "C": "unknown"}
        gain = 0.3
        loss = 0.8
        multiplier = 4.0
        foreground_children = frozenset({"b"})
        root_presence = 0.35

        expected = self._enumerate_probability(
            tree, observations, gain, loss, multiplier, foreground_children, root_presence
        )
        observed_log_likelihood = _pattern_log_likelihood(
            tree, observations, gain, loss, multiplier, foreground_children, root_presence
        )
        node, edge = _posterior_messages(
            tree, observations, gain, loss, multiplier, foreground_children, root_presence
        )

        self.assertAlmostEqual(math.exp(observed_log_likelihood), expected["likelihood"], places=12)
        for node_id, probabilities in node.items():
            self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=12)
            np.testing.assert_allclose(probabilities, expected["node"][node_id], rtol=1e-11, atol=1e-11)
        for edge_id, probabilities in edge.items():
            self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=12)
            np.testing.assert_allclose(probabilities, expected["edge"][edge_id], rtol=1e-11, atol=1e-11)

    def test_ascertainment_matches_enumerated_observed_tip_patterns(self):
        tree = self._small_tree()
        observations = {"A": 0, "B": 1, "C": "unknown"}
        args = (0.3, 0.8, 4.0, frozenset({"b"}), 0.35)
        probability = self._enumerate_probability(tree, observations, *args)["likelihood"]
        zero = self._enumerate_probability(tree, {"A": 0, "B": 0, "C": "unknown"}, *args)["likelihood"]
        one = self._enumerate_probability(tree, {"A": 1, "B": 1, "C": "unknown"}, *args)["likelihood"]
        theta = np.log([0.3, 0.8, 4.0])
        for mode, selected_probability in (
            ("complete-universe", 1.0),
            ("observed-at-least-one", 1.0 - zero),
            ("variable-only", 1.0 - zero - one),
        ):
            with self.subTest(mode=mode):
                actual = _dataset_log_likelihood(
                    tree, [(observations, 3)], "ARD_FOREGROUND", theta,
                    frozenset({"b"}), mode, "fixed", 0.35,
                )
                self.assertAlmostEqual(actual, 3.0 * math.log(probability / selected_probability), places=12)

    def test_zero_branch_impossible_pattern_has_zero_likelihood(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "0"},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": "0"},
            ]
        )
        value = _pattern_log_likelihood(tree, {"A": 0, "B": 1}, 0.2, 0.5, 1.0, frozenset(), 0.5)
        self.assertEqual(value, -math.inf)
        with self.assertRaisesRegex(ValueError, "positive-probability"):
            _posterior_messages(tree, {"A": 0, "B": 1}, 0.2, 0.5, 1.0, frozenset(), 0.5)
        for mode in ("complete-universe", "observed-at-least-one", "variable-only"):
            with self.subTest(mode=mode):
                value = _dataset_log_likelihood(
                    tree, [({"A": 0, "B": 1}, 1)], "ARD", np.log([0.2, 0.5]),
                    frozenset(), mode, "fixed", 0.5,
                )
                self.assertEqual(value, -math.inf)

    def test_expected_transition_count_matches_closed_form_unconditional_gain_count(self):
        gain = 0.25
        loss = 0.75
        branch_length = 1.4
        matrix = _transition_matrix(gain, loss, branch_length)
        joint = np.zeros((2, 2), dtype=float)
        joint[0, 0] = matrix[0, 0]
        joint[0, 1] = matrix[0, 1]
        expected = _expected_transition_count(joint, gain, loss, branch_length, 1.0, 0, 1)

        total = gain + loss
        stationary_zero = loss / total
        closed = gain * (
            stationary_zero * branch_length
            + (1.0 - stationary_zero) * (1.0 - math.exp(-total * branch_length)) / total
        )
        self.assertAlmostEqual(expected, closed, places=11)
        self.assertEqual(_expected_transition_count(joint, gain, loss, 0.0, 1.0, 0, 1), 0.0)

    def test_symmetric_bridge_expected_counts_match_even_and_odd_jump_formulas(self):
        rate, length, multiplier = 0.4, 1.7, 2.0
        mean_jumps = rate * multiplier * length
        for end, total_jumps in (
            (0, mean_jumps * math.tanh(mean_jumps)),
            (1, mean_jumps / math.tanh(mean_jumps)),
        ):
            joint = np.zeros((2, 2), dtype=float)
            joint[0, end] = 1.0
            with self.subTest(end=end):
                gain = _expected_transition_count(joint, rate, rate, length, multiplier, 0, 1)
                loss = _expected_transition_count(joint, rate, rate, length, multiplier, 1, 0)
                self.assertAlmostEqual(gain, 0.5 * (total_jumps + end), places=11)
                self.assertAlmostEqual(loss, 0.5 * (total_jumps - end), places=11)

    def test_large_polytomy_known_patterns_remain_finite_and_normalized(self):
        rows = [{"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"}]
        observations = {}
        for index in range(400):
            node = f"t{index}"
            label = f"T{index}"
            rows.append({"node_id": node, "parent_id": "root", "label": label, "branch_length": "1"})
            observations[label] = 0 if index < 200 else 1
        tree = SpeciesTree(rows)
        rate = 0.001
        changed = -0.5 * math.expm1(-2.0 * rate)
        expected_log_likelihood = 200 * (math.log1p(-changed) + math.log(changed))
        value = _pattern_log_likelihood(tree, observations, rate, rate, 1.0, frozenset(), 0.5)
        self.assertTrue(math.isfinite(value))
        self.assertLess(expected_log_likelihood, -1000.0)
        self.assertAlmostEqual(value, expected_log_likelihood, places=9)
        node, edge = _posterior_messages(tree, observations, rate, rate, 1.0, frozenset(), 0.5)
        np.testing.assert_allclose(node["root"], [0.5, 0.5], rtol=0.0, atol=1e-10)
        for probabilities in node.values():
            self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=10)
        for (_parent, child), probabilities in edge.items():
            self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=10)
            expected = np.zeros((2, 2), dtype=float)
            expected[:, observations[tree.label[child]]] = 0.5
            np.testing.assert_allclose(probabilities, expected, rtol=0.0, atol=1e-10)
