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
    def test_forward_backward_messages_retain_every_global_optimum(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root"},
                {"node_id": "a", "parent_id": "root", "label": "A"},
                {"node_id": "inner", "parent_id": "root", "label": "inner"},
                {"node_id": "b", "parent_id": "inner", "label": "B"},
                {"node_id": "c", "parent_id": "inner", "label": "C"},
            ]
        )
        observations = {"A": 0, "B": 1, "C": "unknown"}
        expected = self._enumerate_optima(tree, observations)
        observed = _parsimony_tables(tree, observations)
        self.assertEqual(observed, expected)
        self.assertTrue(math.isfinite(observed[0]))

    def test_required_and_possible_use_all_optimal_endpoint_pairs(self):
        tree = SpeciesTree(
            [{"node_id": "root", "parent_id": "", "label": "root"}]
            + [
                {"node_id": label.lower(), "parent_id": "root", "label": label}
                for label in "ABC"
            ]
        )
        _minimum, _nodes, pairs = _parsimony_tables(tree, {"A": 0, "B": 0, "C": 1})
        self.assertEqual(pairs[("root", "c")], {(0, 1)})
        self.assertEqual(
            _placement_status("observed_contrast", (0, 1), pairs[("root", "c")]),
            "required",
        )
        self.assertIn((0, 0), pairs[("root", "a")])
        self.assertEqual(
            _placement_status("observed_contrast", (0, 0), pairs[("root", "a")]),
            "no_change",
        )
        ambiguous_tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root"},
                {"node_id": "a", "parent_id": "root", "label": "A"},
                {"node_id": "b", "parent_id": "root", "label": "B"},
            ]
        )
        _minimum, _nodes, ambiguous_pairs = _parsimony_tables(
            ambiguous_tree, {"A": 0, "B": 1}
        )
        self.assertEqual(
            _placement_status(
                "observed_contrast", (0, 1), ambiguous_pairs[("root", "b")]
            ),
            "possible",
        )

    def test_sequence_presence_event_names_describe_dna(self):
        self.assertEqual(_event_type("exon_presence", ("fam", "exon_presence", "S"), 0, 1, set()), "sequence_gain")
        self.assertEqual(_event_type("exon_presence", ("fam", "exon_presence", "S"), 1, 0, set()), "sequence_loss")

    def test_role_and_junction_event_names_match_biological_layers(self):
        key = ("fam", "splice_junction", "J1")
        self.assertEqual(_event_type("splice_junction", key, 0, 1, {key}), "intron_gain")
        self.assertEqual(_event_type("splice_junction", key, 1, 0, {key}), "intron_loss")
        self.assertEqual(
            _structural_relation("splice_junction", key, 0, 1, {key}), "exon_split"
        )
        self.assertEqual(
            _structural_relation("splice_junction", key, 1, 0, {key}), "exon_fusion"
        )
        self.assertEqual(_transition_event_type("splice_junction", 0, 1), "intron_gain")
        self.assertEqual(
            _transition_structural_relation(
                "splice_junction", "within_element_junction", 0, 1
            ),
            "exon_split",
        )
        self.assertEqual(_event_type("exon_role", key, 0, 1, set()), "exon_role_gain")
        self.assertEqual(_event_type("exon_role", key, 1, 0, set()), "exon_role_loss")

    def test_all_unknown_and_no_contrast_have_explicit_reasons(self):
        rows = []
        for site_id, state in (("unknown", "unknown"), ("constant", "absent")):
            for species in ("A", "B"):
                rows.append(
                    {
                        "family_id": "fam",
                        "layer": "exon_presence",
                        "site_id": site_id,
                        "species": species,
                        "state": state,
                        "state_0": "absent",
                        "state_1": "present",
                        "evidence": "fixture",
                    }
                )
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\nroot\t\troot\na\troot\tA\nb\troot\tB\n"
            )
            with patch("intraphy.structural_sites.build_structural_site_matrix", return_value=(rows, [])):
                infer_single_copy_parsimony(input_dir, output_dir)
            summaries = {row["site_id"]: row for row in read_tsv(output_dir / "structural_site_summary.tsv")}
            self.assertEqual(summaries["unknown"]["unavailable_reason"], "no_observed_states")
            self.assertEqual(summaries["constant"]["unavailable_reason"], "event_direction_not_applicable")
            self.assertEqual(summaries["constant"]["inference_status"], "observed_conserved")
            self.assertEqual(summaries["constant"]["min_changes"], "0")
            self.assertEqual(summaries["unknown"]["inference_status"], "not_analyzable")
            for filename in (
                "node_structural_states.tsv",
                "branch_structural_events.tsv",
                "structural_site_summary.tsv",
            ):
                for row in read_tsv(output_dir / filename):
                    self.assertEqual(
                        row["inference_method"], "equal_cost_maximum_parsimony"
                    )
                    self.assertTrue(row["support_kind"])
                    self.assertEqual(
                        row["conditional_scope"], "homologous_sequence_unit_presence"
                    )
                    self.assertTrue(row["unavailable_reason"])

    def test_multiple_cutpoints_are_separate_character_changes(self):
        for count in (2, 3):
            for states, event in (({"A": "absent", "B": "absent", "C": "present"}, "intron_gain"),
                                  ({"A": "present", "B": "present", "C": "absent"}, "intron_loss")):
                outputs = self._run_linked_junctions([states] * count)
                self.assertFalse(outputs["compound_file_exists"])
                changes = [r for r in outputs["branch"] if r["event_type"] != "no_change"]
                self.assertEqual(len(changes), count)
                self.assertEqual({r["event_type"] for r in changes}, {event})
                self.assertEqual({r["placement_status"] for r in changes}, {"required"})
                self.assertEqual(len({r["event_id"] for r in changes}), count)
                self.assertEqual(outputs["gene_summary"][0]["minimum_character_changes"], str(count))
            outputs = self._run_linked_junctions([{"A": "absent", "B": "present"}] * count)
            self.assertFalse(outputs["compound_file_exists"])
            self.assertEqual(outputs["gene_summary"][0]["minimum_character_changes"], str(count))
            changes = [r for r in outputs["branch"] if r["event_type"] != "no_change"]
            self.assertEqual(len(changes), 2 * count)
            self.assertEqual({r["support_kind"] for r in changes}, {"possible"})

    def test_compound_pattern_requires_same_branch_and_direction(self):
        scenarios = {
            "different_branches": [
                {"A": "absent", "B": "absent", "C": "absent", "D": "present"},
                {"A": "absent", "B": "absent", "C": "present", "D": "absent"},
            ],
            "different_directions": [
                {"A": "absent", "B": "absent", "C": "present"},
                {"A": "present", "B": "present", "C": "absent"},
            ],
        }
        for scenario, site_states in scenarios.items():
            with self.subTest(scenario=scenario):
                outputs = self._run_linked_junctions(site_states)
                self.assertEqual(outputs["compound"], [])
                self.assertTrue(
                    any(row["event_type"] != "no_change" for row in outputs["branch"])
                )

    def test_distinct_transcript_path_groups_do_not_form_a_compound_event(self):
        outputs = self._run_linked_junctions(
            [
                {"A": "absent", "B": "present", "C": "present"},
                {"A": "absent", "B": "present", "C": "present"},
            ],
            linked_group_ids=[
                "LG_fam_E1_ref_SP_B_GC_gB_TX_tx1",
                "LG_fam_E1_ref_SP_B_GC_gB_TX_tx2",
            ],
        )
        self.assertEqual(outputs["compound"], [])
        changed = [row for row in outputs["branch"] if row["event_type"] != "no_change"]
        self.assertEqual({row["support_kind"] for row in changed}, {"required"})

    def test_parsimony_reads_frozen_matrix_without_rebuilding(self):
        rows = [
            {
                "family_id": "fam", "layer": "exon_presence", "site_id": "S1",
                "species": species, "state": state, "state_0": "absent", "state_1": "present",
                "evidence": "fixture", "discovery_rule": "complete_catalogue",
            }
            for species, state in (("A", "absent"), ("B", "present"))
        ]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            matrix = input_dir / "frozen.tsv"
            write_structural_site_matrix(matrix, rows)
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\nroot\t\troot\na\troot\tA\nb\troot\tB\n"
            )
            with patch("intraphy.structural_sites.build_structural_site_matrix") as build:
                infer_single_copy_parsimony(
                    input_dir, output_dir, structural_site_matrix_path=matrix
                )
            build.assert_not_called()
            scope = read_tsv(output_dir / "phylogeny_scope.tsv")[0]
            self.assertEqual(scope["matrix_mode"], "frozen_matrix")
            self.assertEqual(scope["matrix_schema_version"], "3")
            self.assertEqual(scope["annotation_view"], "repertoire")

    def test_frozen_matrix_rejects_a_different_annotation_view(self):
        rows = [
            {
                "family_id": "fam", "layer": "exon_role", "site_id": "S1",
                "species": species, "state": state,
                "state_0": "not_exonic", "state_1": "exonic",
                "evidence": "fixture",
                "transcript_scope": "explicit_canonical_transcript",
            }
            for species, state in (("A", "not_exonic"), ("B", "exonic"))
        ]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            output_dir = Path(tmp) / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            matrix = input_dir / "frozen.tsv"
            write_structural_site_matrix(matrix, rows)
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\nroot\t\troot\na\troot\tA\nb\troot\tB\n"
            )
            with self.assertRaisesRegex(SystemExit, "annotation view"):
                infer_single_copy_parsimony(
                    input_dir,
                    output_dir,
                    structural_site_matrix_path=matrix,
                    annotation_view="repertoire",
                )

    def test_schema_three_records_and_checks_annotation_view(self):
        rows = [{
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "exonic", "state_0": "not_exonic",
            "state_1": "exonic", "evidence": "fixture",
            "annotation_view": "canonical",
            "transcript_scope": "explicit_canonical_transcript",
        }]
        with tempfile.TemporaryDirectory() as tmp:
            matrix = Path(tmp) / "matrix.tsv"
            write_structural_site_matrix(matrix, rows)
            observed = read_structural_site_matrix(matrix)
        self.assertEqual(observed[0]["schema_version"], "3")
        self.assertEqual(observed[0]["annotation_view"], "canonical")
        self.assertEqual(structural_matrix_annotation_view(observed), "canonical")
        mixed = observed + [{
            **observed[0], "species": "B", "annotation_view": "repertoire",
            "transcript_scope": "annotated_transcript_repertoire",
        }]
        with self.assertRaisesRegex(ValueError, "mixes canonical and repertoire"):
            structural_matrix_annotation_view(mixed)

    def test_matrix_enforces_canonical_labels_and_site_metadata(self):
        base = {
            "family_id": "fam", "layer": "splice_junction", "site_id": "J1",
            "species": "A", "state": "present", "state_0": "absent",
            "state_1": "present", "evidence": "fixture",
            "annotation_view": "repertoire",
            "transcript_scope": "annotated_transcript_repertoire",
            "linked_group_id": "LG_path", "site_kind": "within_element_junction",
            "discovery_rule": "unique_exact_position_correspondence",
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "matrix.tsv"
            with self.assertRaisesRegex(ValueError, "noncanonical structural-state labels"):
                write_structural_site_matrix(path, [{
                    **base, "state": "absent", "state_0": "present", "state_1": "absent",
                }])
            for field, value in (
                ("linked_group_id", "LG_other"),
                ("site_kind", "between_element_junction"),
                ("discovery_rule", "other_rule"),
                ("annotation_view", "canonical"),
            ):
                with self.subTest(field=field), self.assertRaisesRegex(
                    ValueError, "inconsistent structural-site metadata"
                ):
                    second = {**base, "species": "B", field: value}
                    if field == "annotation_view":
                        second["transcript_scope"] = "explicit_canonical_transcript"
                    write_structural_site_matrix(path, [base, second])

    def test_matrix_writer_has_stable_site_order_and_missing_rows_require_a_reason(self):
        rows = [
            {
                "family_id": "fam", "layer": "exon_presence", "site_id": site,
                "species": species, "state": "present", "state_0": "absent",
                "state_1": "present", "evidence": "fixture",
            }
            for site, species in (("S2", "B"), ("S1", "B"), ("S1", "A"))
        ]
        with tempfile.TemporaryDirectory() as tmp:
            matrix = Path(tmp) / "matrix.tsv"
            write_structural_site_matrix(matrix, rows)
            observed = read_structural_site_matrix(matrix)
            self.assertEqual(
                [(row["site_id"], row["species"]) for row in observed],
                [("S1", "A"), ("S1", "B"), ("S2", "B")],
            )
            with self.assertRaisesRegex(ValueError, "requires observation_reason"):
                write_structural_site_matrix(matrix, [{
                    **rows[0], "state": "unknown", "evidence": "NA",
                    "observation_mask": "missing", "observation_reason": "NA",
                }])

    def test_frozen_matrix_requires_an_explicit_row_for_every_tree_tip(self):
        rows = [{
            "family_id": "fam", "layer": "exon_presence", "site_id": "S1",
            "species": "A", "state": "present", "state_0": "absent",
            "state_1": "present", "evidence": "fixture",
        }]
        with tempfile.TemporaryDirectory() as tmp:
            input_dir = Path(tmp) / "input"
            input_dir.mkdir()
            matrix = input_dir / "matrix.tsv"
            write_structural_site_matrix(matrix, rows)
            (input_dir / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\tNA\na\troot\tA\t1\nb\troot\tB\t1\n"
            )
            for method in (infer_single_copy_parsimony, infer_single_copy_phylogeny):
                with self.subTest(method=method.__name__), self.assertRaisesRegex(
                    SystemExit, "explicit row for every selected tree tip"
                ):
                    method(
                        input_dir,
                        Path(tmp) / method.__name__,
                        structural_site_matrix_path=matrix,
                    )

    def test_schema_three_rejects_unknown_annotation_view(self):
        rows = [{
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "exonic", "state_0": "not_exonic",
            "state_1": "exonic", "evidence": "fixture",
            "annotation_view": "unknown",
            "transcript_scope": "explicit_canonical_transcript",
        }]
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "annotation_view"):
                write_structural_site_matrix(Path(tmp) / "matrix.tsv", rows)

    def test_legacy_v2_role_without_inferable_scope_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            matrix = Path(tmp) / "legacy.tsv"
            matrix.write_text(
                "family_id\tlayer\tsite_id\tspecies\tstate\tstate_0\tstate_1\t"
                "schema_version\ttranscript_scope\n"
                "fam\texon_role\tS1\tA\texonic\tnot_exonic\texonic\t2\tcatalogue_defined\n"
            )
            with self.assertRaisesRegex(ValueError, "inferable transcript_scope"):
                read_structural_site_matrix(matrix)

    def test_legacy_v2_role_without_scope_is_not_defaulted_to_repertoire(self):
        with tempfile.TemporaryDirectory() as tmp:
            matrix = Path(tmp) / "legacy.tsv"
            matrix.write_text(
                "family_id\tlayer\tsite_id\tspecies\tstate\tstate_0\tstate_1\tschema_version\n"
                "fam\texon_role\tS1\tA\texonic\tnot_exonic\texonic\t2\n"
            )
            with self.assertRaisesRegex(ValueError, "inferable transcript_scope"):
                read_structural_site_matrix(matrix)

    def test_schema_v3_requires_view_and_supported_version_is_enforced(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing_view = Path(tmp) / "missing-view.tsv"
            missing_view.write_text(
                "family_id\tlayer\tsite_id\tspecies\tstate\tstate_0\tstate_1\t"
                "schema_version\ttranscript_scope\n"
                "fam\texon_role\tS1\tA\texonic\tnot_exonic\texonic\t3\t"
                "explicit_canonical_transcript\n"
            )
            with self.assertRaisesRegex(ValueError, "schema-v3.*annotation_view"):
                read_structural_site_matrix(missing_view)

            unsupported = [{
                "family_id": "fam", "layer": "exon_presence", "site_id": "S1",
                "species": "A", "state": "present", "state_0": "absent",
                "state_1": "present", "schema_version": "4",
                "annotation_view": "view_independent",
            }]
            with self.assertRaisesRegex(ValueError, "unsupported.*schema_version"):
                write_structural_site_matrix(Path(tmp) / "unsupported.tsv", unsupported)

            unsupported_read = Path(tmp) / "unsupported-read.tsv"
            unsupported_read.write_text(
                "family_id\tlayer\tsite_id\tspecies\tstate\tstate_0\tstate_1\t"
                "schema_version\tannotation_view\n"
                "fam\texon_presence\tS1\tA\tpresent\tabsent\tpresent\t4\t"
                "view_independent\n"
            )
            with self.assertRaisesRegex(ValueError, "unsupported.*schema_version"):
                read_structural_site_matrix(unsupported_read)
