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


class ParsimonySemanticsTestsSupport:
    @staticmethod
    def _run_linked_junctions(site_states, linked_group_ids=None):
        all_species = sorted({species for states in site_states for species in states})
        rows = []
        for index, states in enumerate(site_states, start=1):
            linked_group_id = (
                linked_group_ids[index - 1]
                if linked_group_ids is not None
                else "LG_fam_element_reference_SP_S_GC_g_TX_tx"
            )
            for species in all_species:
                rows.append(
                    {
                        "family_id": "fam",
                        "layer": "splice_junction",
                        "site_id": f"J{index}",
                        "species": species,
                        "state": states[species],
                        "state_0": "absent",
                        "state_1": "present",
                        "evidence": "fixture",
                        "annotation_view": "repertoire",
                        "transcript_scope": "annotated_transcript_repertoire",
                        "linked_group_id": linked_group_id,
                        "site_kind": "within_element_junction",
                    }
                )
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            matrix = input_dir / "matrix.tsv"
            write_structural_site_matrix(matrix, rows)
            tree_rows = ["node_id\tparent_id\tlabel", "root\t\troot"]
            tree_rows.extend(
                f"{species.lower()}\troot\t{species}" for species in all_species
            )
            (input_dir / "species_tree.tsv").write_text("\n".join(tree_rows) + "\n")
            infer_single_copy_parsimony(
                input_dir, output_dir, structural_site_matrix_path=matrix
            )
            return {
                "compound": read_tsv(output_dir / "compound_structural_events.tsv", optional=True),
                "compound_file_exists": (output_dir / "compound_structural_events.tsv").exists(),
                "gene_summary": read_tsv(output_dir / "gene_change_summary.tsv"),
                "branch": read_tsv(output_dir / "branch_structural_events.tsv"),
                "node": read_tsv(output_dir / "node_structural_states.tsv"),
                "site": read_tsv(output_dir / "structural_site_summary.tsv"),
            }

    @staticmethod
    def _enumerate_optima(tree, observations):
        nodes = tree.preorder()
        minimum = math.inf
        histories = []
        for encoded in range(1 << len(nodes)):
            states = {node: (encoded >> index) & 1 for index, node in enumerate(nodes)}
            if any(
                observations.get(tree.label[leaf], "unknown") in {0, 1}
                and states[leaf] != observations[tree.label[leaf]]
                for leaf in tree.leaves
            ):
                continue
            score = sum(states[parent] != states[child] for parent, child in tree.edges())
            if score < minimum:
                minimum = score
                histories = [states]
            elif score == minimum:
                histories.append(states)
        node_states = {
            node: {history[node] for history in histories}
            for node in nodes
        }
        branch_pairs = {
            edge: {(history[edge[0]], history[edge[1]]) for history in histories}
            for edge in tree.edges()
        }
        return float(minimum), node_states, branch_pairs
