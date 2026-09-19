import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.optimize import brentq as scipy_brentq

from insiphy.io import (
    read_structural_site_matrix,
    read_tsv,
    structural_site_observed_state,
    write_structural_site_matrix,
)
from insiphy.parsimony import (
    _event_type,
    _parsimony_tables,
    _placement_status,
    _structural_relation,
    infer_single_copy_parsimony,
)
from insiphy.structural_phylogeny import (
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
from insiphy.tree import SpeciesTree
from insiphy.structural_sites import structural_matrix_annotation_view


class ParsimonySemanticsTests(unittest.TestCase):
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
                "compound": read_tsv(output_dir / "compound_structural_events.tsv"),
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
            with patch("insiphy.parsimony.build_structural_site_matrix", return_value=(rows, [])):
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

    def test_required_and_possible_one_to_n_compound_patterns(self):
        for junction_count in (2, 3):
            segment_count = junction_count + 1
            with self.subTest(segment_count=segment_count, support="required_split"):
                outputs = self._run_linked_junctions(
                    [{"A": "absent", "B": "absent", "C": "present"}]
                    * junction_count
                )
                compounds = outputs["compound"]
                self.assertEqual(len(compounds), 1)
                compound = compounds[0]
                self.assertEqual(compound["pattern_status"], "supported_compound_pattern")
                self.assertEqual(compound["support_kind"], "required")
                self.assertEqual(
                    compound["compound_event_type"], f"one_to_{segment_count}_split"
                )
                self.assertEqual(compound["junction_count"], str(junction_count))
                self.assertEqual(compound["segment_count"], str(segment_count))
                self.assertEqual(
                    compound["joint_support"],
                    "all_components_required_on_same_branch_and_direction",
                )
                self.assertEqual(
                    compound["inference_method"], "equal_cost_maximum_parsimony"
                )
                self.assertEqual(
                    compound["conditional_scope"],
                    "homologous_sequence_present;supplied_transcript_view=repertoire",
                )
                self.assertEqual(compound["unavailable_reason"], "NA")
                changed = [
                    row for row in outputs["branch"]
                    if row["event_type"] != "no_change"
                ]
                self.assertEqual(len(changed), junction_count)
                self.assertEqual({row["event_type"] for row in changed}, {"intron_gain"})

            with self.subTest(segment_count=segment_count, support="required_fusion"):
                outputs = self._run_linked_junctions(
                    [{"A": "present", "B": "present", "C": "absent"}]
                    * junction_count
                )
                compounds = outputs["compound"]
                self.assertEqual(len(compounds), 1)
                self.assertEqual(
                    compounds[0]["compound_event_type"],
                    f"{segment_count}_to_one_fusion",
                )
                self.assertEqual(compounds[0]["support_kind"], "required")

            with self.subTest(segment_count=segment_count, support="possible"):
                outputs = self._run_linked_junctions(
                    [{"A": "absent", "B": "present"}] * junction_count
                )
                compounds = outputs["compound"]
                self.assertEqual(len(compounds), 2)
                self.assertEqual(
                    {row["compound_event_type"] for row in compounds},
                    {f"one_to_{segment_count}_split", f"{segment_count}_to_one_fusion"},
                )
                for compound in compounds:
                    self.assertEqual(compound["pattern_status"], "candidate_compound_pattern")
                    self.assertEqual(compound["support_kind"], "possible_non_joint")
                    self.assertEqual(
                        compound["joint_support"],
                        "not_established_from_marginal_site_support",
                    )
                    self.assertEqual(
                        compound["unavailable_reason"],
                        "joint_history_not_inferred_from_possible_components",
                    )
                    self.assertFalse(
                        any("probability" in field for field in compound)
                    )
                changed = [
                    row for row in outputs["branch"]
                    if row["event_type"] != "no_change"
                ]
                unchanged = [
                    row for row in outputs["branch"]
                    if row["event_type"] == "no_change"
                ]
                self.assertEqual(
                    {row["support_kind"] for row in changed}, {"possible_non_joint"}
                )
                self.assertEqual(
                    {row["support_kind"] for row in unchanged}, {"no_change"}
                )

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
            with patch("insiphy.parsimony.build_structural_site_matrix") as build:
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


class TreeValidationTests(unittest.TestCase):
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


class AscertainmentTests(unittest.TestCase):
    def test_small_discovery_probability_avoids_subtraction_cancellation(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "1"},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": "1"},
            ]
        )
        observations = {"A": 0, "B": 1}
        with patch("insiphy.structural_phylogeny._pattern_log_likelihood", return_value=-1e-20):
            selected = _ascertainment_log_probability(
                tree, observations, 0.2, 0.3, 1.0, frozenset(), 0.5, "observed-at-least-one"
            )
        self.assertAlmostEqual(selected, math.log(1e-20), places=12)
        with patch(
            "insiphy.structural_phylogeny._pattern_log_likelihood",
            side_effect=[-2e-20, math.log(1e-20)],
        ):
            selected = _ascertainment_log_probability(
                tree, observations, 0.2, 0.3, 1.0, frozenset(), 0.5, "variable-only"
            )
        self.assertAlmostEqual(selected, math.log(1e-20), places=12)

    def test_three_ascertainment_modes_select_and_condition_differently(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "1"},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": "1"},
            ]
        )
        all_zero = {"A": 0, "B": 0}
        all_one = {"A": 1, "B": 1}
        variable = {"A": 0, "B": 1}

        self.assertFalse(_selected_for_ascertainment(all_zero, "observed-at-least-one"))
        self.assertTrue(_selected_for_ascertainment(all_one, "observed-at-least-one"))
        self.assertFalse(_selected_for_ascertainment(all_one, "variable-only"))
        self.assertTrue(_selected_for_ascertainment(variable, "variable-only"))
        self.assertTrue(_selected_for_ascertainment(all_zero, "complete-universe"))

        theta = np.array([math.log(0.2)], dtype=float)
        complete = _dataset_log_likelihood(
            tree, [(variable, 1)], "ER", theta, frozenset(), "complete-universe", "fixed", 0.5
        )
        at_least_one = _dataset_log_likelihood(
            tree, [(variable, 1)], "ER", theta, frozenset(), "observed-at-least-one", "fixed", 0.5
        )
        variable_only = _dataset_log_likelihood(
            tree, [(variable, 1)], "ER", theta, frozenset(), "variable-only", "fixed", 0.5
        )
        self.assertGreater(at_least_one, complete)
        self.assertGreater(variable_only, at_least_one)

    def test_zero_length_transition_is_identity(self):
        matrix = _transition_matrix(0.2, 0.4, 0.0)
        self.assertAlmostEqual(float(matrix[0, 0]), 1.0)
        self.assertAlmostEqual(float(matrix[1, 1]), 1.0)
        self.assertAlmostEqual(float(matrix[0, 1]), 0.0)
        self.assertAlmostEqual(float(matrix[1, 0]), 0.0)


class ProbabilityEngineTests(unittest.TestCase):
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

    def test_pruning_and_posteriors_match_explicit_enumeration(self):
        tree = self._small_tree()
        observations = {"A": 0, "B": 1, "C": "unknown"}
        gain = 0.3
        loss = 0.8
        multiplier = 4.0
        foreground_children = frozenset({"b"})
        root_presence = 0.35

        expected = self._enumerate_probability(
            tree, observations, gain, loss, multiplier, foreground_children, root_presence
        )
        observed_log_likelihood = _pattern_log_likelihood(
            tree, observations, gain, loss, multiplier, foreground_children, root_presence
        )
        node, edge = _posterior_messages(
            tree, observations, gain, loss, multiplier, foreground_children, root_presence
        )

        self.assertAlmostEqual(math.exp(observed_log_likelihood), expected["likelihood"], places=12)
        for node_id, probabilities in node.items():
            self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=12)
            np.testing.assert_allclose(probabilities, expected["node"][node_id], rtol=1e-11, atol=1e-11)
        for edge_id, probabilities in edge.items():
            self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=12)
            np.testing.assert_allclose(probabilities, expected["edge"][edge_id], rtol=1e-11, atol=1e-11)

    def test_ascertainment_matches_enumerated_observed_tip_patterns(self):
        tree = self._small_tree()
        observations = {"A": 0, "B": 1, "C": "unknown"}
        args = (0.3, 0.8, 4.0, frozenset({"b"}), 0.35)
        probability = self._enumerate_probability(tree, observations, *args)["likelihood"]
        zero = self._enumerate_probability(tree, {"A": 0, "B": 0, "C": "unknown"}, *args)["likelihood"]
        one = self._enumerate_probability(tree, {"A": 1, "B": 1, "C": "unknown"}, *args)["likelihood"]
        theta = np.log([0.3, 0.8, 4.0])
        for mode, selected_probability in (
            ("complete-universe", 1.0),
            ("observed-at-least-one", 1.0 - zero),
            ("variable-only", 1.0 - zero - one),
        ):
            with self.subTest(mode=mode):
                actual = _dataset_log_likelihood(
                    tree, [(observations, 3)], "ARD_FOREGROUND", theta,
                    frozenset({"b"}), mode, "fixed", 0.35,
                )
                self.assertAlmostEqual(actual, 3.0 * math.log(probability / selected_probability), places=12)

    def test_zero_branch_impossible_pattern_has_zero_likelihood(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": "0"},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": "0"},
            ]
        )
        value = _pattern_log_likelihood(tree, {"A": 0, "B": 1}, 0.2, 0.5, 1.0, frozenset(), 0.5)
        self.assertEqual(value, -math.inf)
        with self.assertRaisesRegex(ValueError, "positive-probability"):
            _posterior_messages(tree, {"A": 0, "B": 1}, 0.2, 0.5, 1.0, frozenset(), 0.5)
        for mode in ("complete-universe", "observed-at-least-one", "variable-only"):
            with self.subTest(mode=mode):
                value = _dataset_log_likelihood(
                    tree, [({"A": 0, "B": 1}, 1)], "ARD", np.log([0.2, 0.5]),
                    frozenset(), mode, "fixed", 0.5,
                )
                self.assertEqual(value, -math.inf)

    def test_expected_transition_count_matches_closed_form_unconditional_gain_count(self):
        gain = 0.25
        loss = 0.75
        branch_length = 1.4
        matrix = _transition_matrix(gain, loss, branch_length)
        joint = np.zeros((2, 2), dtype=float)
        joint[0, 0] = matrix[0, 0]
        joint[0, 1] = matrix[0, 1]
        expected = _expected_transition_count(joint, gain, loss, branch_length, 1.0, 0, 1)

        total = gain + loss
        stationary_zero = loss / total
        closed = gain * (
            stationary_zero * branch_length
            + (1.0 - stationary_zero) * (1.0 - math.exp(-total * branch_length)) / total
        )
        self.assertAlmostEqual(expected, closed, places=11)
        self.assertEqual(_expected_transition_count(joint, gain, loss, 0.0, 1.0, 0, 1), 0.0)

    def test_symmetric_bridge_expected_counts_match_even_and_odd_jump_formulas(self):
        rate, length, multiplier = 0.4, 1.7, 2.0
        mean_jumps = rate * multiplier * length
        for end, total_jumps in (
            (0, mean_jumps * math.tanh(mean_jumps)),
            (1, mean_jumps / math.tanh(mean_jumps)),
        ):
            joint = np.zeros((2, 2), dtype=float)
            joint[0, end] = 1.0
            with self.subTest(end=end):
                gain = _expected_transition_count(joint, rate, rate, length, multiplier, 0, 1)
                loss = _expected_transition_count(joint, rate, rate, length, multiplier, 1, 0)
                self.assertAlmostEqual(gain, 0.5 * (total_jumps + end), places=11)
                self.assertAlmostEqual(loss, 0.5 * (total_jumps - end), places=11)

    def test_large_polytomy_known_patterns_remain_finite_and_normalized(self):
        rows = [{"node_id": "root", "parent_id": "", "label": "root", "branch_length": "NA"}]
        observations = {}
        for index in range(400):
            node = f"t{index}"
            label = f"T{index}"
            rows.append({"node_id": node, "parent_id": "root", "label": label, "branch_length": "1"})
            observations[label] = 0 if index < 200 else 1
        tree = SpeciesTree(rows)
        rate = 0.001
        changed = -0.5 * math.expm1(-2.0 * rate)
        expected_log_likelihood = 200 * (math.log1p(-changed) + math.log(changed))
        value = _pattern_log_likelihood(tree, observations, rate, rate, 1.0, frozenset(), 0.5)
        self.assertTrue(math.isfinite(value))
        self.assertLess(expected_log_likelihood, -1000.0)
        self.assertAlmostEqual(value, expected_log_likelihood, places=9)
        node, edge = _posterior_messages(tree, observations, rate, rate, 1.0, frozenset(), 0.5)
        np.testing.assert_allclose(node["root"], [0.5, 0.5], rtol=0.0, atol=1e-10)
        for probabilities in node.values():
            self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=10)
        for (_parent, child), probabilities in edge.items():
            self.assertAlmostEqual(float(probabilities.sum()), 1.0, places=10)
            expected = np.zeros((2, 2), dtype=float)
            expected[:, observations[tree.label[child]]] = 0.5
            np.testing.assert_allclose(probabilities, expected, rtol=0.0, atol=1e-10)


class ProfileIntervalTests(unittest.TestCase):
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

        with patch("insiphy.structural_phylogeny.minimize", side_effect=optimizer), patch(
            "insiphy.structural_phylogeny.brentq", side_effect=root_finder
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
                "insiphy.structural_phylogeny.minimize",
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
        with patch("insiphy.structural_phylogeny.brentq", side_effect=ValueError("no root")):
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


class EstimatedRateInferenceTests(unittest.TestCase):
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
            "insiphy.structural_phylogeny.minimize",
            side_effect=lambda objective, start, **kwargs: SimpleNamespace(
                x=np.asarray(start, dtype=float), fun=1e-6, success=True, message="numerical convergence",
            ),
        )
        self.information = patch(
            "insiphy.structural_phylogeny._observed_information",
            side_effect=lambda objective, theta: (True, np.ones(len(theta)), 1.0),
        )
        self.profiles = patch(
            "insiphy.structural_phylogeny._profile_interval", side_effect=self._flat_profile,
        )

    @staticmethod
    def _flat_profile(objective, optimum, index, bounds, max_log_likelihood, transform=math.exp):
        lower, upper = bounds[index]
        return {
            "low": transform(lower), "high": transform(upper), "status": "range_limited_both",
            "raw_low": lower, "raw_high": upper,
        }

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
                        "insiphy.structural_phylogeny.build_structural_site_matrix", return_value=(site_rows, []),
                    ), patch("insiphy.structural_phylogeny._posterior_messages") as posterior:
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


class PosteriorStatusTests(unittest.TestCase):
    def _fit(self, model, status, identifiable):
        intervals = [{"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None}]
        if model in {"ARD", "ARD_FOREGROUND"}:
            intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        if model == "ARD_FOREGROUND":
            intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        theta_length = 2 if model == "ER" else 3 if model == "ARD" else 4
        if model == "ARD":
            contrast = {
                "name": "gain_loss_rate_ratio", "low": 0.5, "high": 2.0,
                "status": "two_sided", "raw_low": math.log(0.5), "raw_high": math.log(2.0),
                "theta_low": np.zeros(theta_length), "theta_high": np.zeros(theta_length),
            }
        elif model == "ARD_FOREGROUND":
            contrast = {
                "name": "foreground_multiplier", "low": 0.5, "high": 2.0,
                "status": "two_sided", "raw_low": math.log(0.5), "raw_high": math.log(2.0),
                "theta_low": np.zeros(theta_length), "theta_high": np.zeros(theta_length),
            }
        else:
            contrast = {"name": "NA", "status": "not_applicable"}
        return {
            "model": model,
            "theta": np.zeros(theta_length),
            "gain_rate": 0.1,
            "loss_rate": 0.1,
            "foreground_multiplier": 1.0,
            "root_presence": 0.5,
            "log_likelihood": -2.0 if model == "ER" else -1.9,
            "parameter_count": theta_length,
            "aic": 8.0 if model == "ER" else 9.8,
            "converged": True,
            "fit_status": status,
            "optimizer_message": status,
            "boundary": False,
            "intervals": intervals,
            "start_count": 1,
            "identifiable": identifiable,
            "information_condition": math.inf,
            "root_frequency": "estimated",
            "profile_support_thetas": [],
            "contrast_profile": contrast,
        }

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
                    "insiphy.structural_phylogeny.build_structural_site_matrix",
                    return_value=(site_rows, []),
                ), patch(
                    "insiphy.structural_phylogeny.fit_model",
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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "insiphy.structural_phylogeny.fit_model",
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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch("insiphy.structural_phylogeny.fit_model") as fit:
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
            with patch("insiphy.structural_phylogeny.build_structural_site_matrix") as build:
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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "insiphy.structural_phylogeny.fit_model",
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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "insiphy.structural_phylogeny.fit_model",
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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "insiphy.structural_phylogeny.fit_model",
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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "insiphy.structural_phylogeny.fit_model",
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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), self.assertRaises(SystemExit):
                infer_single_copy_phylogeny(
                    input_dir,
                    output_dir,
                    ascertainment="complete-universe",
                )

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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
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
                "insiphy.structural_phylogeny.build_structural_site_matrix",
                return_value=(site_rows, []),
            ), patch(
                "insiphy.structural_phylogeny.fit_model",
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


if __name__ == "__main__":
    unittest.main()
