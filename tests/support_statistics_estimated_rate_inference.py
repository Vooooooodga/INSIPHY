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


class EstimatedRateInferenceTestsSupport:
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
            "intraphy.inference.fitting.minimize",
            side_effect=lambda objective, start, **kwargs: SimpleNamespace(
                x=np.asarray(start, dtype=float), fun=1e-6, success=True, message="numerical convergence",
            ),
        )
        self.information = patch(
            "intraphy.inference.fitting._observed_information",
            side_effect=lambda objective, theta: (True, np.ones(len(theta)), 1.0),
        )
        self.profiles = patch(
            "intraphy.inference.fitting._profile_interval", side_effect=self._flat_profile,
        )

    @staticmethod
    def _flat_profile(objective, optimum, index, bounds, max_log_likelihood, transform=math.exp):
        lower, upper = bounds[index]
        return {
            "low": transform(lower), "high": transform(upper), "status": "range_limited_both",
            "raw_low": lower, "raw_high": upper,
        }
