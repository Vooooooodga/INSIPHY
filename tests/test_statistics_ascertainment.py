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

from support_statistics_ascertainment import AscertainmentTestsSupport

class AscertainmentTests(AscertainmentTestsSupport, unittest.TestCase):
    def test_small_discovery_probability_avoids_subtraction_cancellation(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "1"},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": "1"},
            ]
        )
        observations = {"A": 0, "B": 1}
        with patch("intraphy.inference.ctmc._pattern_log_likelihood", return_value=-1e-20):
            selected = _ascertainment_log_probability(
                tree, observations, 0.2, 0.3, 1.0, frozenset(), 0.5, "observed-at-least-one"
            )
        self.assertAlmostEqual(selected, math.log(1e-20), places=12)
        with patch(
            "intraphy.inference.ctmc._pattern_log_likelihood",
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
