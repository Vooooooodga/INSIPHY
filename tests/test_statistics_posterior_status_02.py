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
    def test_complete_universe_rejects_default_derived_discovery_rule(self):
        site_rows = [
            {
                "family_id": "fam", "layer": "exon_presence", "site_id": "S1",
                "species": species, "state": state, "state_0": "absent", "state_1": "present",
                "evidence": "fixture", "discovery_rule": "derived_from_element_correspondence",
                "observation_mask": "observed",
            }
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
            ), self.assertRaises(SystemExit):
                infer_single_copy_phylogeny(
                    input_dir, output_dir, ascertainment="complete-universe"
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
                "intraphy.structural_sites.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "intraphy.inference.layer_fitting.fit_model",
                side_effect=[
                    valid_fit,
                    self._fit("ARD", "nonidentifiable", False),
                ],
            ):
                infer_single_copy_phylogeny(input_dir, output_dir)

            changes = read_tsv(output_dir / "structural_changes.tsv")
            nodes = read_tsv(output_dir / "node_state_posteriors.tsv")
            branches = read_tsv(output_dir / "branch_transition_posteriors.tsv")
            tests = read_tsv(output_dir / "model_tests.tsv")
            self.assertEqual(changes[0]["structural_change_type"], "posterior_not_reported")
            self.assertIn("site_likelihood_zero", changes[0]["conditioning"])
            self.assertEqual(tests[0]["lrt_available"], "false")
            self.assertEqual(tests[0]["posterior_available"], "true")
            self.assertEqual(tests[0]["posterior_model"], "ER")
            self.assertFalse(nodes)
            self.assertFalse(branches)
