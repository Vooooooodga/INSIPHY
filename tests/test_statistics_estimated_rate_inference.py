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

from support_statistics_estimated_rate_inference import EstimatedRateInferenceTestsSupport

class EstimatedRateInferenceTests(EstimatedRateInferenceTestsSupport, unittest.TestCase):
    def test_all_present_mixed_missing_invalidates_estimated_rates_despite_positive_hessian(self):
        with self.optimizer, self.information, self.profiles:
            for model in ("ER", "ARD", "ARD_FOREGROUND"):
                with self.subTest(model=model):
                    fit = fit_model(self.tree, self.patterns, model, foreground_children=frozenset({"a"}))
                    self.assertFalse(fit["converged"])
                    self.assertFalse(fit["boundary"])
                    self.assertTrue(math.isinf(fit["information_condition"]))
                    self.assertEqual(fit["intervals"][-1]["status"], "no_observed_contrast")
                    self.assertEqual(fit["fit_status"], "not_fitted")
                    self.assertEqual(fit["inference_status"], "no_observed_contrast")
                    self.assertFalse(fit["identifiable"])
                    self.assertIsNone(fit["gain_rate"])
                    self.assertIsNone(fit["loss_rate"])
                    self.assertIsNone(fit["root_presence"])
                    self.assertIsNone(fit["aic"])
                    self.assertFalse(fit["posterior_available"])
                    self.assertEqual(fit["posterior_unavailable_reason"], "no_observed_contrast")
                    self.assertFalse(_fit_valid_for_posterior(fit))

    def test_open_profiles_alone_do_not_invalidate_a_fit_with_observed_contrast(self):
        with self.optimizer, self.information, self.profiles:
            fit = fit_model(self.tree, [({"A": 0, "B": 1, "C": "unknown", "D": 1}, 1)], "ER")
        self.assertEqual(fit["intervals"][-1]["status"], "range_limited_both")
        self.assertEqual(fit["fit_status"], "success")
        self.assertEqual(fit["inference_status"], "success")
        self.assertTrue(fit["posterior_available"])
        self.assertEqual(fit["posterior_unavailable_reason"], "NA")
        self.assertEqual(fit["posterior_sensitivity_status"], "unavailable_open_profile_interval")
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
                        "intraphy.structural_sites.build_structural_site_matrix", return_value=(site_rows, []),
                    ), patch("intraphy.inference.posterior_export._posterior_messages") as posterior:
                        infer_single_copy_phylogeny(
                            input_dir, output_dir, model=model,
                            foreground_branches=foreground if model == "foreground" else None,
                        )
                    posterior.assert_not_called()
                    fits = read_tsv(output_dir / "model_fits.tsv")
                    self.assertEqual(len(fits), 2)
                    for fit in fits:
                        self.assertEqual(fit["fit_status"], "not_fitted")
                        self.assertEqual(fit["inference_status"], "no_observed_contrast")
                        self.assertEqual(fit["identifiable"], "false")
                        self.assertEqual(fit["converged"], "false")
                        self.assertEqual(fit["posterior_available"], "false")
                        self.assertEqual(fit["posterior_unavailable_reason"], "no_observed_contrast")
                        self.assertEqual(fit["gain_rate"], "NA")
                        self.assertEqual(fit["loss_rate"], "NA")
                        self.assertEqual(fit["root_presence"], "NA")
                        self.assertEqual(fit["aic"], "NA")
                        self.assertEqual(fit["root_presence_ci_status"], "no_observed_contrast")
                    tests = read_tsv(output_dir / "model_tests.tsv")
                    self.assertEqual(tests[0]["p_value"], "NA")
                    self.assertEqual(tests[0]["q_value"], "NA")
                    self.assertEqual(tests[0]["lrt_available"], "false")
                    self.assertEqual(tests[0]["lrt_unavailable_reason"], "no_observed_contrast")
                    self.assertEqual(tests[0]["reference_distribution"], "asymptotic_chi_square")
                    self.assertEqual(tests[0]["small_sample_accuracy"], "unassessed")
                    self.assertEqual(tests[0]["posterior_available"], "false")
                    self.assertEqual(tests[0]["inference_status"], "no_observed_contrast")
                    changes = read_tsv(output_dir / "structural_changes.tsv")
                    self.assertEqual(len(changes), len(self.patterns))
                    for change in changes:
                        self.assertEqual(change["structural_change_type"], "posterior_not_reported")
                        self.assertEqual(change["endpoint_change_probability"], "NA")
                        self.assertIn("inference_status=no_observed_contrast", change["conditioning"])
                    self.assertEqual(read_tsv(output_dir / "node_state_posteriors.tsv"), [])
                    self.assertEqual(read_tsv(output_dir / "branch_transition_posteriors.tsv"), [])
