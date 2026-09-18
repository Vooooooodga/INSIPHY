import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.optimize import brentq as scipy_brentq

from insiphy.io import read_tsv
from insiphy.structural_phylogeny import (
    PROFILE_DROP_95,
    _ascertainment_log_probability,
    _dataset_log_likelihood,
    _expected_transition_count,
    _fit_valid_for_posterior,
    _pattern_log_likelihood,
    _posterior_messages,
    _profile_interval,
    _selected_for_ascertainment,
    _transition_matrix,
    _validated_tree_rows,
    fit_model,
    infer_single_copy_phylogeny,
)
from insiphy.tree import SpeciesTree


class TreeValidationTests(unittest.TestCase):
    def test_species_tree_records_missing_lengths_and_accepts_zero_branches(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "0"},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": ""},
            ]
        )
        self.assertIsNone(tree.length["root"])
        self.assertEqual(tree.branch_length("a"), 0.0)
        with self.assertRaises(SystemExit):
            tree.branch_length("b")

    def test_species_tree_rejects_invalid_topology(self):
        with self.assertRaises(SystemExit):
            SpeciesTree(
                [
                    {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "0"},
                    {"node_id": "root", "parent_id": "", "label": "root2", "branch_length": "0"},
                ]
            )
        with self.assertRaises(SystemExit):
            SpeciesTree(
                [
                    {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "0"},
                    {"node_id": "a", "parent_id": "missing", "label": "A", "branch_length": "1"},
                ]
            )
        with self.assertRaises(SystemExit):
            SpeciesTree(
                [
                    {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "0"},
                    {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "1"},
                    {"node_id": "b", "parent_id": "root", "label": "A", "branch_length": "1"},
                ]
            )
        with self.assertRaises(SystemExit):
            SpeciesTree(
                [
                    {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "0"},
                    {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "-1"},
                ]
            )
        with self.assertRaises(SystemExit):
            SpeciesTree(
                [
                    {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "0"},
                    {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "inf"},
                ]
            )
        with self.assertRaises(SystemExit):
            SpeciesTree(
                [
                    {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "0"},
                    {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "nan"},
                ]
            )
        with self.assertRaises(SystemExit):
            SpeciesTree(
                [
                    {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "0"},
                    {"node_id": "a", "parent_id": "b", "label": "A", "branch_length": "1"},
                    {"node_id": "b", "parent_id": "a", "label": "B", "branch_length": "1"},
                ]
            )

    def test_likelihood_branch_length_modes(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree_file = Path(tmp) / "species_tree.tsv"
            tree_file.write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\n"
                "a\troot\tA\t\n"
                "b\troot\tB\t0\n"
            )
            with self.assertRaises(SystemExit):
                _validated_tree_rows(tree_file, "supplied")
            rows = _validated_tree_rows(tree_file, "unit")
            tree = SpeciesTree(rows)
            self.assertEqual(tree.branch_length("a"), 1.0)
            self.assertEqual(tree.branch_length("b"), 1.0)


class AscertainmentTests(unittest.TestCase):
    def test_small_discovery_probability_avoids_subtraction_cancellation(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "1"},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": "1"},
            ]
        )
        observations = {"A": 0, "B": 1}
        with patch("insiphy.structural_phylogeny._pattern_log_likelihood", return_value=-1e-20):
            selected = _ascertainment_log_probability(
                tree, observations, 0.2, 0.3, 1.0, frozenset(), 0.5, "observed-at-least-one"
            )
        self.assertAlmostEqual(selected, math.log(1e-20), places=12)
        with patch(
            "insiphy.structural_phylogeny._pattern_log_likelihood",
            side_effect=[-2e-20, math.log(1e-20)],
        ):
            selected = _ascertainment_log_probability(
                tree, observations, 0.2, 0.3, 1.0, frozenset(), 0.5, "variable-only"
            )
        self.assertAlmostEqual(selected, math.log(1e-20), places=12)

    def test_three_ascertainment_modes_select_and_condition_differently(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "1"},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": "1"},
            ]
        )
        all_zero = {"A": 0, "B": 0}
        all_one = {"A": 1, "B": 1}
        variable = {"A": 0, "B": 1}

        self.assertFalse(_selected_for_ascertainment(all_zero, "observed-at-least-one"))
        self.assertTrue(_selected_for_ascertainment(all_one, "observed-at-least-one"))
        self.assertFalse(_selected_for_ascertainment(all_one, "variable-only"))
        self.assertTrue(_selected_for_ascertainment(variable, "variable-only"))
        self.assertTrue(_selected_for_ascertainment(all_zero, "complete-universe"))

        theta = np.array([math.log(0.2)], dtype=float)
        complete = _dataset_log_likelihood(
            tree, [(variable, 1)], "ER", theta, frozenset(), "complete-universe", "fixed", 0.5
        )
        at_least_one = _dataset_log_likelihood(
            tree, [(variable, 1)], "ER", theta, frozenset(), "observed-at-least-one", "fixed", 0.5
        )
        variable_only = _dataset_log_likelihood(
            tree, [(variable, 1)], "ER", theta, frozenset(), "variable-only", "fixed", 0.5
        )
        self.assertGreater(at_least_one, complete)
        self.assertGreater(variable_only, at_least_one)

    def test_zero_length_transition_is_identity(self):
        matrix = _transition_matrix(0.2, 0.4, 0.0)
        self.assertAlmostEqual(float(matrix[0, 0]), 1.0)
        self.assertAlmostEqual(float(matrix[1, 1]), 1.0)
        self.assertAlmostEqual(float(matrix[0, 1]), 0.0)
        self.assertAlmostEqual(float(matrix[1, 0]), 0.0)


class ProbabilityEngineTests(unittest.TestCase):
    def _small_tree(self):
        return SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "0.4"},
                {"node_id": "inner", "parent_id": "root", "label": "inner", "branch_length": "0.7"},
                {"node_id": "b", "parent_id": "inner", "label": "B", "branch_length": "0.5"},
                {"node_id": "c", "parent_id": "inner", "label": "C", "branch_length": "0.3"},
            ]
        )

    def _enumerate_probability(self, tree, observations, gain, loss, multiplier, foreground_children, root_presence):
        nodes = tree.preorder()
        prior = {0: 1.0 - root_presence, 1: root_presence}
        totals = {
            "likelihood": 0.0,
            "node": {node: np.zeros(2, dtype=float) for node in nodes},
            "edge": {edge: np.zeros((2, 2), dtype=float) for edge in tree.edges()},
        }
        for values in range(1 << len(nodes)):
            states = {node: (values >> idx) & 1 for idx, node in enumerate(nodes)}
            probability = prior[states[tree.root]]
            for parent, child in tree.edges():
                branch_multiplier = multiplier if child in foreground_children else 1.0
                matrix = _transition_matrix(gain, loss, tree.branch_length(child), branch_multiplier)
                probability *= float(matrix[states[parent], states[child]])
            for leaf in tree.leaves:
                observed = observations.get(tree.label[leaf], "unknown")
                if observed in {0, 1} and states[leaf] != observed:
                    probability = 0.0
                    break
            totals["likelihood"] += probability
            if probability == 0.0:
                continue
            for node in nodes:
                totals["node"][node][states[node]] += probability
            for parent, child in tree.edges():
                totals["edge"][(parent, child)][states[parent], states[child]] += probability
        for node in nodes:
            totals["node"][node] /= totals["likelihood"]
        for edge in list(totals["edge"]):
            totals["edge"][edge] /= totals["likelihood"]
        return totals

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


class ProfileIntervalTests(unittest.TestCase):
    def test_nuisance_optimizer_failure_inside_brentq_has_no_finite_endpoint(self):
        root_evaluations = 0
        failed_inside_root = False

        def optimizer(objective, start, **kwargs):
            nonlocal failed_inside_root
            if root_evaluations == 3 and not failed_inside_root:
                failed_inside_root = True
                return SimpleNamespace(success=False, fun=math.inf)
            return SimpleNamespace(success=True, fun=objective(np.zeros_like(start)))

        def root_finder(function, lower, upper):
            def evaluate(value):
                nonlocal root_evaluations
                root_evaluations += 1
                return function(value)

            return scipy_brentq(evaluate, lower, upper)

        with patch("insiphy.structural_phylogeny.minimize", side_effect=optimizer), patch(
            "insiphy.structural_phylogeny.brentq", side_effect=root_finder
        ):
            interval = _profile_interval(
                lambda theta: 0.5 * float(theta @ theta), np.zeros(2), 0,
                [(-4.0, 4.0), (-4.0, 4.0)], 0.0, transform=float,
            )
        self.assertTrue(failed_inside_root)
        self.assertIsNone(interval["low"])
        self.assertIsNone(interval["raw_low"])
        self.assertAlmostEqual(interval["high"], math.sqrt(2.0 * PROFILE_DROP_95), places=10)
        self.assertEqual(interval["status"], "profile_optimization_failed")

    def test_failed_or_nonfinite_profile_optimization_has_no_bounds(self):
        for success in (False, True):
            with self.subTest(success=success), patch(
                "insiphy.structural_phylogeny.minimize",
                return_value=SimpleNamespace(success=success, fun=math.inf),
            ):
                interval = _profile_interval(
                    lambda theta: float(theta @ theta), np.zeros(2), 0,
                    [(-4.0, 4.0), (-4.0, 4.0)], 0.0,
                )
            self.assertIsNone(interval["low"])
            self.assertIsNone(interval["high"])
            self.assertIsNone(interval["raw_low"])
            self.assertIsNone(interval["raw_high"])
            self.assertEqual(interval["status"], "profile_optimization_failed")

    def test_failed_root_search_has_no_bounds(self):
        with patch("insiphy.structural_phylogeny.brentq", side_effect=ValueError("no root")):
            interval = _profile_interval(
                lambda theta: float(theta[0] ** 2), [0.0], 0, [(-4.0, 4.0)], 0.0,
            )
        self.assertIsNone(interval["low"])
        self.assertIsNone(interval["high"])
        self.assertEqual(interval["status"], "profile_root_failed")

    def test_flat_profile_labels_numerical_search_limits(self):
        interval = _profile_interval(lambda theta: 0.0, [0.0], 0, [(-4.0, 4.0)], 0.0, transform=float)
        self.assertEqual(interval["status"], "range_limited_both")
        self.assertEqual((interval["raw_low"], interval["raw_high"]), (-4.0, 4.0))

    def test_failed_endpoint_keeps_opposite_range_limit_label(self):
        interval = _profile_interval(
            lambda theta: math.nan if theta[0] < 0.0 else 0.0,
            [0.0], 0, [(-4.0, 4.0)], 0.0, transform=float,
        )
        self.assertIsNone(interval["low"])
        self.assertEqual(interval["high"], 4.0)
        self.assertEqual(interval["status"], "profile_optimization_failed;upper_range_limited")


class EstimatedRateInferenceTests(unittest.TestCase):
    def setUp(self):
        self.tree = SpeciesTree(
            [{"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"}]
            + [
                {"node_id": label.lower(), "parent_id": "root", "label": label, "branch_length": "1"}
                for label in "ABCD"
            ]
        )
        self.patterns = [
            ({"A": 1, "B": 1, "C": "unknown", "D": "unknown"}, 1),
            ({"A": "unknown", "B": "unknown", "C": 1, "D": 1}, 1),
            ({"A": 1, "B": 1, "C": 1, "D": 1}, 1),
        ]
        self.optimizer = patch(
            "insiphy.structural_phylogeny.minimize",
            side_effect=lambda objective, start, **kwargs: SimpleNamespace(
                x=np.asarray(start, dtype=float), fun=1e-6, success=True, message="numerical convergence",
            ),
        )
        self.information = patch(
            "insiphy.structural_phylogeny._observed_information",
            side_effect=lambda objective, theta: (True, np.ones(len(theta)), 1.0),
        )
        self.profiles = patch(
            "insiphy.structural_phylogeny._profile_interval", side_effect=self._flat_profile,
        )

    @staticmethod
    def _flat_profile(objective, optimum, index, bounds, max_log_likelihood, transform=math.exp):
        lower, upper = bounds[index]
        return {
            "low": transform(lower), "high": transform(upper), "status": "range_limited_both",
            "raw_low": lower, "raw_high": upper,
        }

    def test_all_present_mixed_missing_invalidates_estimated_rates_despite_positive_hessian(self):
        with self.optimizer, self.information, self.profiles:
            for model in ("ER", "ARD", "ARD_FOREGROUND"):
                with self.subTest(model=model):
                    fit = fit_model(self.tree, self.patterns, model, foreground_children=frozenset({"a"}))
                    self.assertTrue(fit["converged"])
                    self.assertFalse(fit["boundary"])
                    self.assertEqual(fit["information_condition"], 1.0)
                    self.assertEqual(fit["intervals"][-1]["status"], "range_limited_both")
                    self.assertEqual(fit["fit_status"], "not_estimable")
                    self.assertEqual(fit["inference_status"], "no_observed_contrast")
                    self.assertFalse(fit["identifiable"])
                    self.assertFalse(_fit_valid_for_posterior(fit))

    def test_open_profiles_alone_do_not_invalidate_a_fit_with_observed_contrast(self):
        with self.optimizer, self.information, self.profiles:
            fit = fit_model(self.tree, [({"A": 0, "B": 1, "C": "unknown", "D": 1}, 1)], "ER")
        self.assertEqual(fit["intervals"][-1]["status"], "range_limited_both")
        self.assertEqual(fit["fit_status"], "success")
        self.assertEqual(fit["inference_status"], "success")
        self.assertTrue(_fit_valid_for_posterior(fit))

    def test_fixed_parameter_conditional_posterior_allows_no_observed_contrast(self):
        rate, root_presence = 0.2, 0.4
        nodes, _edges = _posterior_messages(
            self.tree, self.patterns[0][0], rate, rate, 1.0, frozenset(), root_presence,
        )
        changed = -0.5 * math.expm1(-2.0 * rate)
        present_weight = root_presence * (1.0 - changed) ** 2
        absent_weight = (1.0 - root_presence) * changed ** 2
        self.assertAlmostEqual(nodes["root"][1], present_weight / (present_weight + absent_weight), places=12)

    def test_no_contrast_exports_diagnostics_without_formal_posteriors(self):
        site_rows = [
            {
                "family_id": "fam", "layer": "exon_presence", "site_id": f"S{index}",
                "species": species, "state": "present" if state == 1 else "unknown",
                "state_0": "absent", "state_1": "present", "evidence": "fixture",
            }
            for index, (pattern, _weight) in enumerate(self.patterns)
            for species, state in pattern.items()
        ]
        with self.optimizer, self.information, self.profiles:
            for model in ("er-ard", "foreground"):
                with self.subTest(model=model), tempfile.TemporaryDirectory() as tmp:
                    input_dir = Path(tmp) / "input"
                    output_dir = Path(tmp) / "output"
                    input_dir.mkdir()
                    output_dir.mkdir()
                    (input_dir / "species_tree.tsv").write_text(
                        "node_id\tparent_id\tlabel\tbranch_length\nroot\t\troot\tNA\n"
                        + "".join(f"{label.lower()}\troot\t{label}\t1\n" for label in "ABCD")
                    )
                    foreground = input_dir / "foreground.tsv"
                    foreground.write_text("parent_id\tchild_id\nroot\ta\n")
                    with patch(
                        "insiphy.structural_phylogeny.build_structural_site_matrix", return_value=(site_rows, []),
                    ), patch("insiphy.structural_phylogeny._posterior_messages") as posterior:
                        infer_single_copy_phylogeny(
                            input_dir, output_dir, model=model,
                            foreground_branches=foreground if model == "foreground" else None,
                        )
                    posterior.assert_not_called()
                    fits = read_tsv(output_dir / "model_fits.tsv")
                    self.assertEqual(len(fits), 2)
                    for fit in fits:
                        self.assertEqual(fit["fit_status"], "not_estimable")
                        self.assertEqual(fit["inference_status"], "no_observed_contrast")
                        self.assertEqual(fit["identifiable"], "false")
                        self.assertEqual(fit["converged"], "true")
                        self.assertNotEqual(fit["gain_rate"], "NA")
                        self.assertEqual(fit["root_presence_ci_status"], "range_limited_both")
                    tests = read_tsv(output_dir / "model_tests.tsv")
                    self.assertEqual(tests[0]["p_value"], "NA")
                    self.assertEqual(tests[0]["q_value"], "NA")
                    self.assertEqual(tests[0]["reference_distribution"], "not_available")
                    self.assertEqual(tests[0]["inference_status"], "no_observed_contrast")
                    changes = read_tsv(output_dir / "structural_changes.tsv")
                    self.assertEqual(len(changes), len(self.patterns))
                    for change in changes:
                        self.assertEqual(change["structural_change_type"], "posterior_not_reported")
                        self.assertEqual(change["endpoint_change_probability"], "NA")
                        self.assertIn("inference_status=no_observed_contrast", change["conditioning"])
                    self.assertEqual(read_tsv(output_dir / "node_state_posteriors.tsv"), [])
                    self.assertEqual(read_tsv(output_dir / "branch_transition_posteriors.tsv"), [])


class PosteriorStatusTests(unittest.TestCase):
    def _fit(self, model, status, identifiable):
        intervals = [{"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None}]
        if model in {"ARD", "ARD_FOREGROUND"}:
            intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        if model == "ARD_FOREGROUND":
            intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        return {
            "model": model,
            "theta": np.zeros(2),
            "gain_rate": 0.1,
            "loss_rate": 0.1,
            "foreground_multiplier": 1.0,
            "root_presence": 0.5,
            "log_likelihood": -2.0 if model == "ER" else -1.9,
            "parameter_count": 2 if model == "ER" else 3,
            "aic": 8.0 if model == "ER" else 9.8,
            "converged": True,
            "fit_status": status,
            "optimizer_message": status,
            "boundary": False,
            "intervals": intervals,
            "start_count": 1,
            "identifiable": identifiable,
            "information_condition": math.inf,
            "root_frequency": "estimated",
        }

    def test_inference_selects_patterns_for_each_ascertainment_mode(self):
        patterns = {
            "zero": (0, 0),
            "one": (1, 1),
            "variable": (0, 1),
            "variable_repeat": (0, 1),
            "outside_catalogue": (1, 0),
        }
        site_rows = [
            {
                "family_id": "fam", "layer": "exon_presence", "site_id": site,
                "species": species, "state": "present" if state else "absent",
                "state_0": "absent", "state_1": "present", "evidence": "fixture",
            }
            for site, values in patterns.items()
            for species, state in zip(("A", "B"), values)
        ]
        selected = {
            "observed-at-least-one": {"one", "variable", "variable_repeat", "outside_catalogue"},
            "variable-only": {"variable", "variable_repeat", "outside_catalogue"},
            "complete-universe": {"zero", "one", "variable", "variable_repeat"},
        }
        for mode, expected_sites in selected.items():
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                input_dir = Path(tmp) / "input"
                output_dir = Path(tmp) / "output"
                input_dir.mkdir()
                output_dir.mkdir()
                (input_dir / "species_tree.tsv").write_text(
                    "node_id\tparent_id\tlabel\tbranch_length\n"
                    "root\t\troot\tNA\na\troot\tA\t1\nb\troot\tB\t1\n"
                )
                (input_dir / "structural_site_universe.tsv").write_text(
                    "family_id\tlayer\tsite_id\n"
                    + "".join(f"fam\texon_presence\t{site}\n" for site in patterns if site != "outside_catalogue")
                )
                with patch(
                    "insiphy.structural_phylogeny.build_structural_site_matrix",
                    return_value=(site_rows, []),
                ), patch(
                    "insiphy.structural_phylogeny.fit_model",
                    side_effect=[self._fit("ER", "nonidentifiable", False), self._fit("ARD", "nonidentifiable", False)],
                ) as fit:
                    infer_single_copy_phylogeny(input_dir, output_dir, ascertainment=mode)
                expected_patterns = {}
                for site in expected_sites:
                    values = patterns[site]
                    expected_patterns[values] = expected_patterns.get(values, 0) + 1
                self.assertEqual(fit.call_count, 2)
                for call in fit.call_args_list:
                    actual_patterns = {
                        (observation["A"], observation["B"]): weight
                        for observation, weight in call.args[1]
                    }
                    self.assertEqual(actual_patterns, expected_patterns)
                    self.assertEqual(call.kwargs["ascertainment"], mode)
                changes = read_tsv(output_dir / "structural_changes.tsv")
                self.assertEqual({row["site_id"] for row in changes}, expected_sites)

    def test_invalid_likelihood_fit_produces_diagnostics_without_posteriors(self):
        site_rows = [
            {
                "family_id": "fam",
                "layer": "exon_presence",
                "site_id": "S1",
                "species": "A",
                "state": "present",
                "state_0": "absent",
                "state_1": "present",
                "evidence": "fixture",
            },
            {
                "family_id": "fam",
                "layer": "exon_presence",
                "site_id": "S1",
                "species": "B",
                "state": "absent",
                "state_0": "absent",
                "state_1": "present",
                "evidence": "fixture",
            },
        ]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\n"
                "a\troot\tA\t1\n"
                "b\troot\tB\t1\n"
            )
            with patch(
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "insiphy.structural_phylogeny.fit_model",
                side_effect=[
                    self._fit("ER", "nonidentifiable", False),
                    self._fit("ARD", "nonidentifiable", False),
                ],
            ):
                infer_single_copy_phylogeny(input_dir, output_dir)

            tests = read_tsv(output_dir / "model_tests.tsv")
            changes = read_tsv(output_dir / "structural_changes.tsv")
            nodes = read_tsv(output_dir / "node_state_posteriors.tsv")
            branches = read_tsv(output_dir / "branch_transition_posteriors.tsv")

            self.assertEqual(tests[0]["p_value"], "NA")
            self.assertEqual(tests[0]["reference_distribution"], "not_available")
            self.assertEqual(changes[0]["structural_change_type"], "posterior_not_reported")
            self.assertEqual(changes[0]["endpoint_change_probability"], "NA")
            self.assertFalse(nodes)
            self.assertFalse(branches)

    def test_complete_universe_requires_catalog(self):
        site_rows = [
            {
                "family_id": "fam",
                "layer": "exon_presence",
                "site_id": "S1",
                "species": "A",
                "state": "present",
                "state_0": "absent",
                "state_1": "present",
                "evidence": "fixture",
            },
            {
                "family_id": "fam",
                "layer": "exon_presence",
                "site_id": "S1",
                "species": "B",
                "state": "absent",
                "state_0": "absent",
                "state_1": "present",
                "evidence": "fixture",
            },
        ]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\n"
                "a\troot\tA\t1\n"
                "b\troot\tB\t1\n"
            )
            with patch(
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), self.assertRaises(SystemExit):
                infer_single_copy_phylogeny(
                    input_dir,
                    output_dir,
                    ascertainment="complete-universe",
                )

    def test_zero_probability_site_gets_diagnostic_not_nan_posterior(self):
        site_rows = [
            {
                "family_id": "fam",
                "layer": "exon_presence",
                "site_id": "S1",
                "species": "A",
                "state": "absent",
                "state_0": "absent",
                "state_1": "present",
                "evidence": "fixture",
            },
            {
                "family_id": "fam",
                "layer": "exon_presence",
                "site_id": "S1",
                "species": "B",
                "state": "present",
                "state_0": "absent",
                "state_1": "present",
                "evidence": "fixture",
            },
        ]
        valid_fit = self._fit("ER", "success", True)
        valid_fit["aic"] = 1.0
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\n"
                "a\troot\tA\t0\n"
                "b\troot\tB\t0\n"
            )
            with patch(
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "insiphy.structural_phylogeny.fit_model",
                side_effect=[
                    valid_fit,
                    self._fit("ARD", "nonidentifiable", False),
                ],
            ):
                infer_single_copy_phylogeny(input_dir, output_dir)

            changes = read_tsv(output_dir / "structural_changes.tsv")
            nodes = read_tsv(output_dir / "node_state_posteriors.tsv")
            branches = read_tsv(output_dir / "branch_transition_posteriors.tsv")
            self.assertEqual(changes[0]["structural_change_type"], "posterior_not_reported")
            self.assertIn("site_likelihood_zero", changes[0]["conditioning"])
            self.assertFalse(nodes)
            self.assertFalse(branches)


if __name__ == "__main__":
    unittest.main()
