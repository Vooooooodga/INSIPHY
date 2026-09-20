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

from support_statistics_parsimony_semantics import ParsimonySemanticsTestsSupport

class ParsimonySemanticsTests(ParsimonySemanticsTestsSupport, unittest.TestCase):
    def test_schema_two_and_explicit_legacy_presence_rows_are_supported(self):
        with tempfile.TemporaryDirectory() as tmp:
            for version in ("2", "1-legacy", ""):
                with self.subTest(version=version):
                    matrix = Path(tmp) / f"schema-{version or 'missing'}.tsv"
                    columns = (
                        "family_id\tlayer\tsite_id\tspecies\tstate\tstate_0\tstate_1"
                    )
                    if version:
                        columns += "\tschema_version"
                    matrix.write_text(
                        columns + "\n"
                        + "fam\texon_presence\tS1\tA\tpresent\tabsent\tpresent"
                        + (f"\t{version}" if version else "") + "\n"
                    )
                    observed = read_structural_site_matrix(matrix)[0]
                    self.assertEqual(
                        observed["schema_version"], version or "1-legacy"
                    )

    def test_annotation_view_conflicting_with_transcript_scope_is_rejected(self):
        row = {
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "exonic", "state_0": "not_exonic",
            "state_1": "exonic", "annotation_view": "repertoire",
            "transcript_scope": "explicit_canonical_transcript",
        }
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "conflicts with transcript_scope"):
                write_structural_site_matrix(Path(tmp) / "conflict.tsv", [row])

    def test_observation_masks_are_bounded_and_legacy_missing_masks_mask_state(self):
        row = {
            "family_id": "fam", "layer": "exon_presence", "site_id": "S1",
            "species": "A", "state": "present", "state_0": "absent",
            "state_1": "present", "observation_mask": "ambiguous",
        }
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "invalid observation_mask"):
                write_structural_site_matrix(Path(tmp) / "ambiguous.tsv", [row])

            for legacy_mask in ("masked", "excluded", "unobserved", "unknown"):
                with self.subTest(mask=legacy_mask):
                    row["observation_mask"] = legacy_mask
                    matrix = Path(tmp) / f"{legacy_mask}.tsv"
                    write_structural_site_matrix(matrix, [row])
                    observed = read_structural_site_matrix(matrix)[0]
                    self.assertEqual(observed["observation_mask"], "missing")
                    self.assertEqual(structural_site_observed_state(observed), "unknown")

    def test_legacy_scope_aliases_match_matrix_view_inference(self):
        rows = [{
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "exonic", "state_0": "not_exonic",
            "state_1": "exonic", "schema_version": "2",
            "transcript_scope": "annotated_transcript_repertoire",
        }]
        self.assertEqual(structural_matrix_annotation_view(rows), "repertoire")

    def test_structural_matrix_rejects_cross_field_contradictions(self):
        base = {
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "exonic", "state_0": "not_exonic",
            "state_1": "exonic", "evidence": "fixture",
            "annotation_view": "repertoire",
            "transcript_scope": "annotated_transcript_repertoire",
            "applicability": "applicable", "observation_mask": "observed",
        }
        contradictory = [
            {**base, "applicability": "inapplicable"},
            {**base, "state": "unknown"},
            {**base, "observation_mask": "generated_all_zero"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            for index, row in enumerate(contradictory):
                with self.subTest(index=index):
                    with self.assertRaises(ValueError):
                        write_structural_site_matrix(Path(tmp) / f"bad-{index}.tsv", [row])

    def test_within_element_junction_matrix_drives_parsimony_event_name(self):
        rows = [{
            "family_id": "fam", "layer": "splice_junction", "site_id": "J1",
            "species": species, "state": state, "state_0": "absent",
            "state_1": "present", "evidence": "fixture",
            "annotation_view": "repertoire",
            "transcript_scope": "annotated_transcript_repertoire",
            "site_kind": "within_element_junction",
        } for species, state in (("A", "absent"), ("B", "present"))]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            matrix = input_dir / "matrix.tsv"
            write_structural_site_matrix(matrix, rows)
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\nroot\t\troot\na\troot\tA\nb\troot\tB\n"
            )
            infer_single_copy_parsimony(
                input_dir, output_dir, structural_site_matrix_path=matrix
            )
            events = read_tsv(output_dir / "branch_structural_events.tsv")
        changed = {row["event_type"] for row in events if row["event_type"] != "no_change"}
        self.assertTrue(changed)
        self.assertLessEqual(changed, {"intron_gain", "intron_loss"})
        relations = {
            row["structural_relation"]
            for row in events
            if row["event_type"] != "no_change"
        }
        self.assertLessEqual(relations, {"exon_split", "exon_fusion"})
