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

from support_statistics_profile_interval import ProfileIntervalTestsSupport

class ProfileIntervalTests(ProfileIntervalTestsSupport, unittest.TestCase):
    def test_nuisance_optimizer_failure_inside_brentq_has_no_finite_endpoint(self):
        root_evaluations = 0
        failed_inside_root = False

        def optimizer(objective, start, **kwargs):
            nonlocal failed_inside_root
            if root_evaluations == 3 and not failed_inside_root:
                failed_inside_root = True
                return SimpleNamespace(success=False, fun=math.inf)
            values = np.zeros_like(start)
            return SimpleNamespace(success=True, fun=objective(values), x=values)

        def root_finder(function, lower, upper):
            def evaluate(value):
                nonlocal root_evaluations
                root_evaluations += 1
                return function(value)

            return scipy_brentq(evaluate, lower, upper)

        with patch("intraphy.inference.fitting.minimize", side_effect=optimizer), patch(
            "intraphy.inference.fitting.brentq", side_effect=root_finder
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
                "intraphy.inference.fitting.minimize",
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
        with patch("intraphy.inference.fitting.brentq", side_effect=ValueError("no root")):
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

    def test_rate_ratio_profile_saves_constrained_endpoint_parameters(self):
        interval = _profile_contrast_interval(
            lambda theta: 0.5 * float(theta @ theta),
            np.zeros(2),
            [(-4.0, 4.0), (-4.0, 4.0)],
            0.0,
            [1.0, -1.0],
            "gain_loss_rate_ratio",
        )
        self.assertEqual(interval["status"], "two_sided")
        self.assertIsNotNone(interval["theta_low"])
        self.assertIsNotNone(interval["theta_high"])
        self.assertAlmostEqual(
            interval["theta_low"][0] - interval["theta_low"][1],
            interval["raw_low"],
            places=6,
        )
