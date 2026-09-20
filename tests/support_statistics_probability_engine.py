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


class ProbabilityEngineTestsSupport:
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
