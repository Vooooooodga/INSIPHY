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

from support_statistics_posterior_status import PosteriorStatusTestsSupport

class PosteriorStatusTests(PosteriorStatusTestsSupport, unittest.TestCase):
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
                "discovery_rule": "curated_complete_universe", "observation_mask": "observed",
            }
            for site, values in patterns.items()
            for species, state in zip(("A", "B"), values)
        ]
        selected = {
            "observed-at-least-one": {"one", "variable", "variable_repeat", "outside_catalogue"},
            "variable-only": {"variable", "variable_repeat", "outside_catalogue"},
            "complete-universe": set(patterns),
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
                with patch(
                    "intraphy.structural_sites.build_structural_site_matrix",
                    return_value=(site_rows, []),
                ), patch(
                    "intraphy.inference.layer_fitting.fit_model",
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
                "intraphy.structural_sites.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "intraphy.inference.layer_fitting.fit_model",
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
            self.assertEqual(tests[0]["lrt_available"], "false")
            self.assertEqual(tests[0]["reference_distribution"], "asymptotic_chi_square")
            self.assertEqual(tests[0]["small_sample_accuracy"], "unassessed")
            self.assertEqual(changes[0]["structural_change_type"], "posterior_not_reported")
            self.assertEqual(changes[0]["endpoint_change_probability"], "NA")
            self.assertFalse(nodes)
            self.assertFalse(branches)

    def test_all_unknown_layer_is_reported_without_calling_optimizer(self):
        site_rows = [
            {
                "family_id": "fam",
                "layer": "exon_presence",
                "site_id": "S1",
                "species": species,
                "state": "unknown",
                "state_0": "absent",
                "state_1": "present",
                "evidence": "fixture",
            }
            for species in ("A", "B")
        ]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\na\troot\tA\t1\nb\troot\tB\t1\n"
            )
            with patch(
                "intraphy.structural_sites.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch("intraphy.inference.layer_fitting.fit_model") as fit:
                infer_single_copy_phylogeny(input_dir, output_dir)
            fit.assert_not_called()
            fits = read_tsv(output_dir / "model_fits.tsv")
            self.assertEqual(len(fits), 2)
            self.assertTrue(all(row["fit_status"] == "not_fitted" for row in fits))
            self.assertTrue(all(row["inference_status"] == "no_observed_states" for row in fits))
            self.assertTrue(all(row["posterior_available"] == "false" for row in fits))
            tests = read_tsv(output_dir / "model_tests.tsv")
            self.assertEqual(tests[0]["lrt_available"], "false")
            self.assertEqual(tests[0]["lrt_unavailable_reason"], "no_observed_states")
            self.assertEqual(tests[0]["small_sample_accuracy"], "unassessed")
            changes = read_tsv(output_dir / "structural_changes.tsv")
            self.assertEqual(len(changes), 1)
            self.assertEqual(changes[0]["posterior_unavailable_reason"], "no_observed_states")

    def test_likelihood_reads_frozen_matrix_without_rebuilding(self):
        rows = [
            {
                "family_id": "fam", "layer": "exon_presence", "site_id": "S1",
                "species": species, "state": "unknown", "state_0": "absent",
                "state_1": "present", "evidence": "fixture",
            }
            for species in ("A", "B")
        ]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            matrix = input_dir / "frozen.tsv"
            write_structural_site_matrix(matrix, rows)
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\na\troot\tA\t1\nb\troot\tB\t1\n"
            )
            with patch("intraphy.structural_sites.build_structural_site_matrix") as build:
                infer_single_copy_phylogeny(
                    input_dir, output_dir, structural_site_matrix_path=matrix
                )
            build.assert_not_called()
            scope = read_tsv(output_dir / "phylogeny_scope.tsv")[0]
            self.assertEqual(scope["matrix_mode"], "frozen_matrix")
            self.assertEqual(scope["matrix_schema_version"], "3")

    def test_open_profiles_leave_posterior_sensitivity_ranges_unavailable(self):
        site_rows = [
            {
                "family_id": "fam",
                "layer": "exon_presence",
                "site_id": "S1",
                "species": species,
                "state": state,
                "state_0": "absent",
                "state_1": "present",
                "evidence": "fixture",
            }
            for species, state in (("A", "absent"), ("B", "present"))
        ]
        null_fit = self._fit("ER", "success", True)
        alternative_fit = self._fit("ARD", "success", True)
        for fit in (null_fit, alternative_fit):
            for interval in fit["intervals"]:
                interval["status"] = "range_limited_both"
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\na\troot\tA\t1\nb\troot\tB\t1\n"
            )
            with patch(
                "intraphy.structural_sites.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "intraphy.inference.layer_fitting.fit_model",
                side_effect=[null_fit, alternative_fit],
            ):
                infer_single_copy_phylogeny(input_dir, output_dir)
            tests = read_tsv(output_dir / "model_tests.tsv")
            self.assertEqual(tests[0]["lrt_available"], "true")
            self.assertEqual(tests[0]["reference_distribution"], "asymptotic_chi_square")
            nodes = read_tsv(output_dir / "node_state_posteriors.tsv")
            self.assertTrue(nodes)
            self.assertTrue(all(row["profile_probability_low"] == "NA" for row in nodes))
            self.assertTrue(all(row["profile_probability_high"] == "NA" for row in nodes))
            self.assertTrue(
                all(row["uncertainty_status"] == "unavailable_open_profile_interval" for row in nodes)
            )

    def test_open_tested_contrast_disables_asymptotic_lrt(self):
        site_rows = [
            {
                "family_id": "fam", "layer": "exon_presence", "site_id": "S1",
                "species": species, "state": state, "state_0": "absent", "state_1": "present",
                "evidence": "fixture",
            }
            for species, state in (("A", "absent"), ("B", "present"))
        ]
        null_fit = self._fit("ER", "success", True)
        alternative_fit = self._fit("ARD", "success", True)
        alternative_fit["contrast_profile"]["status"] = "upper_range_limited"
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\na\troot\tA\t1\nb\troot\tB\t1\n"
            )
            with patch(
                "intraphy.structural_sites.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "intraphy.inference.layer_fitting.fit_model",
                side_effect=[null_fit, alternative_fit],
            ):
                infer_single_copy_phylogeny(input_dir, output_dir)
            test = read_tsv(output_dir / "model_tests.tsv")[0]
            self.assertEqual(test["lrt_available"], "false")
            self.assertEqual(
                test["lrt_unavailable_reason"],
                "tested_contrast_profile_upper_range_limited",
            )
            self.assertEqual(test["reference_distribution"], "asymptotic_chi_square")
            self.assertEqual(test["small_sample_accuracy"], "unassessed")

    def test_finite_profile_endpoints_produce_posterior_envelopes(self):
        site_rows = [
            {
                "family_id": "fam", "layer": "exon_presence", "site_id": "S1",
                "species": species, "state": state, "state_0": "absent", "state_1": "present",
                "evidence": "fixture",
            }
            for species, state in (("A", "absent"), ("B", "present"))
        ]
        null_fit = self._fit("ER", "success", True)
        null_fit["aic"] = 1.0
        support = [np.array([math.log(0.05), 0.0]), np.array([math.log(0.2), 0.0])]
        for index, interval in enumerate(null_fit["intervals"]):
            interval.update(
                {
                    "status": "two_sided", "theta_low": support[0], "theta_high": support[1],
                    "low": 0.05 if index == 0 else 0.4,
                    "high": 0.2 if index == 0 else 0.6,
                }
            )
        null_fit["profile_support_thetas"] = support
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\na\troot\tA\t1\nb\troot\tB\t1\n"
            )
            with patch(
                "intraphy.structural_sites.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "intraphy.inference.layer_fitting.fit_model",
                side_effect=[null_fit, self._fit("ARD", "nonidentifiable", False)],
            ):
                infer_single_copy_phylogeny(input_dir, output_dir)
            nodes = read_tsv(output_dir / "node_state_posteriors.tsv")
            branches = read_tsv(output_dir / "branch_transition_posteriors.tsv")
            self.assertTrue(nodes)
            self.assertTrue(branches)
            self.assertTrue(all(row["profile_probability_low"] != "NA" for row in nodes))
            self.assertTrue(all(row["profile_probability_high"] != "NA" for row in nodes))
            self.assertTrue(all(row["uncertainty_status"] == "finite_profile_envelope" for row in nodes))
            self.assertTrue(
                all(row["profile_transition_probability_low"] != "NA" for row in branches)
            )

    def test_linked_sites_pause_ctmc_lrt(self):
        site_rows = [
            {
                "family_id": "fam", "layer": "splice_junction", "site_id": site,
                "species": species, "state": state, "state_0": "absent", "state_1": "present",
                "evidence": "fixture", "linked_group_id": "split_group_1",
                "annotation_view": "repertoire",
                "transcript_scope": "annotated_transcript_repertoire",
                "site_kind": "within_element_junction",
            }
            for site in ("J1", "J2")
            for species, state in (("A", "absent"), ("B", "present"))
        ]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\na\troot\tA\t1\nb\troot\tB\t1\n"
            )
            with patch(
                "intraphy.structural_sites.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "intraphy.inference.layer_fitting.fit_model",
                side_effect=[self._fit("ER", "success", True), self._fit("ARD", "success", True)],
            ):
                infer_single_copy_phylogeny(input_dir, output_dir)
            test = read_tsv(output_dir / "model_tests.tsv")[0]
            self.assertEqual(test["lrt_available"], "false")
            self.assertEqual(test["lrt_unavailable_reason"], "correlated_linked_sites_not_modelled")
            self.assertEqual(test["correlated_linked_group_count"], "1")

    def test_complete_universe_rejects_id_only_catalogue_without_species_observations(self):
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
                "discovery_rule": "independent_catalogue",
                "observation_mask": "observed",
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
            (input_dir / "structural_site_universe.tsv").write_text(
                "family_id\tlayer\tsite_id\nfam\texon_presence\tS1\n"
            )
            with patch(
                "intraphy.structural_sites.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), self.assertRaises(SystemExit):
                infer_single_copy_phylogeny(
                    input_dir,
                    output_dir,
                    ascertainment="complete-universe",
                )
