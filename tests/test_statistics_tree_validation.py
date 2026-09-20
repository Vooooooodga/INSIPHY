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

from support_statistics_tree_validation import TreeValidationTestsSupport

class TreeValidationTests(TreeValidationTestsSupport, unittest.TestCase):
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
