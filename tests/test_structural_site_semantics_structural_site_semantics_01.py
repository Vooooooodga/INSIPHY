import json

import tempfile

import unittest

from pathlib import Path

from intraphy.io import normalize_structural_site_row, read_tsv

from intraphy.structural_sites import (
    _complete_tree_tip_observations,
    _continuous_reference_block_covers,
    _element_site_rows,
    _junction_site_rows,
    _mapped_reference_coordinate,
    _merge_catalogue_observations,
    _parse_projected_reference_blocks,
)

from support_structural_site_semantics_structural_site_semantics import StructuralSiteSemanticsTestsSupport

class StructuralSiteSemanticsTests(StructuralSiteSemanticsTestsSupport, unittest.TestCase):
    def test_predicted_exon_supports_dna_presence_and_leaves_role_unknown(self):
        occurrences, elements = self._base_element()
        completion = [{
            "evidence_id": "evB", "homology_id": "H1", "family_id": "fam",
            "species": "B", "gene_copy_id": "gB", "homologous_dna_presence": "present",
            "homologous_dna_evidence": "protein_coding_projection",
            "candidate_resolution_status": "resolved", "predicted_role": "CDS",
            "predicted_role_blocks": json.dumps([{
                "block_id": "pb1", "target_contig": "chrB", "target_start": 100,
                "target_end": 129, "target_strand": "+",
            }]),
        }]
        rows = _element_site_rows(occurrences, elements, completion, {"fam"}, {"fam": {"A", "B"}})
        by_layer = self._by_layer(rows, "B")
        self.assertEqual(by_layer["exon_presence"]["state"], "present")
        self.assertEqual(by_layer["exon_role"]["state"], "unknown")
        self.assertEqual(by_layer["exon_presence"]["observation_mask"], "observed")
        self.assertEqual(by_layer["exon_role"]["observation_mask"], "missing")
        self.assertIn("predicted_exonic_role_candidate", by_layer["exon_role"]["evidence"])

    def test_completion_assignment_does_not_silently_choose_one_local_element(self):
        occurrences, elements = self._base_element()
        elements[0]["element_id"] = "EG1.local_001"
        elements.append({**elements[0], "element_id": "EG1.local_002"})
        completion = [{
            "evidence_id": "evB", "homology_id": "H1", "family_id": "fam",
            "species": "B", "gene_copy_id": "gB", "homologous_dna_presence": "present",
            "homologous_dna_evidence": "protein_coding_projection",
            "candidate_resolution_status": "resolved", "correspondence_status": "resolved",
            "source_parent_occurrence_id": "A_ex", "predicted_role": "CDS",
        }]
        rows = _element_site_rows(
            occurrences, elements, completion, {"fam"}, {"fam": {"A", "B"}}
        )
        b_presence = [
            row for row in rows
            if row["species"] == "B" and row["layer"] == "exon_presence"
        ]
        self.assertEqual(len(b_presence), 2)
        self.assertEqual({row["state"] for row in b_presence}, {"unknown"})
        self.assertTrue(all(
            "completion_element_assignment_ambiguous" in row["evidence"]
            for row in b_presence
        ))

    def test_known_nonexonic_homologous_dna_is_retained(self):
        occurrences, elements = self._base_element()
        occurrences.append({
            "occurrence_id": "B_in", "family_id": "fam", "species": "B",
            "gene_copy_id": "gB", "role": "intron", "presence_status": "present",
            "contig": "chrB", "start": "100", "end": "120", "strand": "+",
        })
        elements.append({
            "element_id": "EG1", "homology_id": "H1", "occurrence_id": "B_in",
            "element_class": "candidate_source", "membership_call": "core_member",
            "membership_edge_eligible": "1", "position_edge_eligible": "1",
            "correspondence_status": "resolved", "genomic_matched_blocks": "chrB:100-120:+",
        })
        rows = _element_site_rows(occurrences, elements, [], {"fam"}, {"fam": {"A", "B"}})
        by_layer = self._by_layer(rows, "B")
        self.assertEqual(by_layer["exon_presence"]["state"], "present")
        self.assertEqual(by_layer["exon_role"]["state"], "not_exonic")

    def test_ordered_flank_deletion_makes_role_inapplicable(self):
        occurrences, elements = self._base_element()
        completion = [{
            "evidence_id": "delB", "homology_id": "H1", "family_id": "fam",
            "species": "B", "gene_copy_id": "gB", "homologous_dna_presence": "absent",
            "homologous_dna_evidence": "ordered_flank_deletion",
            "candidate_resolution_status": "resolved",
        }]
        rows = _element_site_rows(occurrences, elements, completion, {"fam"}, {"fam": {"A", "B"}})
        by_layer = self._by_layer(rows, "B")
        self.assertEqual(by_layer["exon_presence"]["state"], "absent")
        self.assertEqual(by_layer["exon_role"]["state"], "unknown")
        self.assertEqual(by_layer["exon_role"]["applicability"], "inapplicable")

    def test_numeric_zero_presence_makes_role_inapplicable(self):
        occurrences, elements = self._base_element()
        occurrences[0]["presence_status"] = "0"
        rows = _element_site_rows(
            occurrences, elements, [], {"fam"}, {"fam": {"A"}},
        )
        by_layer = self._by_layer(rows, "A")

        self.assertEqual(by_layer["exon_presence"]["state"], "absent")
        self.assertEqual(by_layer["exon_role"]["state"], "unknown")
        self.assertEqual(by_layer["exon_role"]["applicability"], "inapplicable")
        self.assertEqual(by_layer["exon_role"]["observation_mask"], "missing")

    def test_prediction_conflict_uses_local_block_intersection(self):
        overlapping = self._prediction_overlap_case(110, 115)
        separate = self._prediction_overlap_case(200, 210)
        self.assertEqual(overlapping["state"], "unknown")
        self.assertEqual(overlapping["observation_reason"], "local_annotation_conflict")
        self.assertEqual(separate["state"], "not_exonic")

    def test_local_overlap_requires_resolved_same_strand(self):
        antisense = self._prediction_overlap_case(110, 115, "-")
        missing = self._prediction_overlap_case(110, 115, "NA")
        self.assertEqual(antisense["state"], "not_exonic")
        self.assertEqual(missing["state"], "unknown")
        self.assertIn("local_overlap_strand_unavailable", missing["evidence"])

    def test_repertoire_any_exonic_path_wins_over_intronic_path(self):
        occurrences, elements = self._base_element()
        occurrences[0].update(contig="chr", start="100", end="110", transcript_id="tx_ex")
        elements[0]["genomic_matched_blocks"] = "chr:100-110:+"
        occurrences.extend([
            {"occurrence_id": "left", "family_id": "fam", "species": "A", "gene_copy_id": "gA",
             "role": "exon", "presence_status": "present", "contig": "chr", "start": "50", "end": "90", "strand": "+"},
            {"occurrence_id": "right", "family_id": "fam", "species": "A", "gene_copy_id": "gA",
             "role": "exon", "presence_status": "present", "contig": "chr", "start": "120", "end": "150", "strand": "+"},
        ])
        paths = [
            {"family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "tx_ex",
             "occurrence_id": "A_ex", "path_status": "annotated_transcript_path"},
            {"family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "tx_skip",
             "occurrence_id": "left", "path_status": "annotated_transcript_path"},
            {"family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "tx_skip",
             "occurrence_id": "right", "path_status": "annotated_transcript_path"},
        ]
        rows = _element_site_rows(occurrences, elements, [], {"fam"}, {"fam": {"A"}}, paths)
        role = self._by_layer(rows, "A")["exon_role"]
        self.assertEqual(role["state"], "exonic")
        self.assertIn("path_role:tx_ex:exonic", role["evidence"])
        self.assertIn("path_role:tx_skip:not_exonic", role["evidence"])
        self.assertIn("alternative_usage", role["evidence"])

    def test_path_local_conflict_stays_unknown_with_presence_observed(self):
        occurrences, elements = self._base_element()
        occurrences[0].update(contig="chr", start="100", end="120", transcript_id="tx_conflict")
        elements[0]["genomic_matched_blocks"] = "chr:90-110:+"
        paths = [{
            "family_id": "fam", "species": "A", "gene_copy_id": "gA",
            "transcript_id": "tx_conflict", "occurrence_id": "A_ex",
            "path_role": "exonic", "path_status": "annotated_transcript_path",
        }]
        rows = _element_site_rows(
            occurrences, elements, [], {"fam"}, {"fam": {"A"}}, paths,
        )
        by_layer = self._by_layer(rows, "A")

        self.assertEqual(by_layer["exon_presence"]["state"], "present")
        self.assertEqual(by_layer["exon_role"]["state"], "unknown")
        self.assertIn("local_path_role_conflict", by_layer["exon_role"]["evidence"])

    def test_repertoire_exonic_path_wins_over_a_conflicting_alternative_path(self):
        occurrences, elements = self._base_element()
        occurrences[0].update(contig="chr", start="100", end="120")
        elements[0]["genomic_matched_blocks"] = "chr:100-110:+"
        occurrences.append({
            "occurrence_id": "conflict_ex", "family_id": "fam", "species": "A",
            "gene_copy_id": "gA", "role": "exon", "presence_status": "present",
            "contig": "chr", "start": "105", "end": "130", "strand": "+",
        })
        paths = [
            {"family_id": "fam", "species": "A", "gene_copy_id": "gA",
             "transcript_id": "tx_exonic", "occurrence_id": "A_ex",
             "path_status": "annotated_transcript_path"},
            {"family_id": "fam", "species": "A", "gene_copy_id": "gA",
             "transcript_id": "tx_conflict", "occurrence_id": "conflict_ex",
             "path_status": "annotated_transcript_path"},
        ]
        rows = _element_site_rows(
            occurrences, elements, [], {"fam"}, {"fam": {"A"}}, paths,
        )
        role = self._by_layer(rows, "A")["exon_role"]
        self.assertEqual(role["state"], "exonic")
        self.assertIn("path_role:tx_conflict:conflicting", role["evidence"])

    def test_repertoire_exonic_path_wins_over_completion_conflict_candidate(self):
        occurrences, elements = self._base_element()
        occurrences[0].update(contig="chr", start="100", end="120")
        elements[0]["genomic_matched_blocks"] = "chr:100-120:+"
        paths = [{
            "family_id": "fam", "species": "A", "gene_copy_id": "gA",
            "transcript_id": "tx_exonic", "occurrence_id": "A_ex",
            "path_status": "annotated_transcript_path",
        }]
        completion = [{
            "evidence_id": "candidate_conflict", "homology_id": "H1",
            "family_id": "fam", "species": "A", "gene_copy_id": "gA",
            "homologous_dna_presence": "present",
            "homologous_dna_evidence": "local_dna_alignment",
            "candidate_resolution_status": "resolved",
            "dna_aligned_blocks": json.dumps([{
                "block_id": "dna", "target_contig": "chr",
                "target_start": 100, "target_end": 120, "target_strand": "+",
            }]),
            "annotation_conflict_blocks": json.dumps([{
                "block_id": "conflict", "target_contig": "chr",
                "target_start": 105, "target_end": 110, "target_strand": "+",
            }]),
        }]
        rows = _element_site_rows(
            occurrences, elements, completion, {"fam"}, {"fam": {"A"}}, paths,
        )
        role = self._by_layer(rows, "A")["exon_role"]
        self.assertEqual(role["state"], "exonic")
        self.assertEqual(role["observation_reason"], "repertoire_any_exonic_path")
        self.assertIn("local_annotation_prediction_conflict", role["evidence"])

    def test_canonical_view_requires_explicit_canonical_path(self):
        occurrences, elements = self._base_element()
        paths = [{
            "family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "tx1",
            "occurrence_id": "A_ex", "path_status": "annotated_transcript_path",
        }]
        rows = _element_site_rows(occurrences, elements, [], {"fam"}, {"fam": {"A"}},
                                  transcript_paths=paths, annotation_view="canonical")
        role = self._by_layer(rows, "A")["exon_role"]
        self.assertEqual(role["state"], "unknown")
        self.assertIn("canonical_transcript_not_explicitly_recorded", role["evidence"])

    def test_canonical_view_without_transcript_paths_keeps_role_unknown(self):
        occurrences, elements = self._base_element()
        rows = _element_site_rows(
            occurrences,
            elements,
            [],
            {"fam"},
            {"fam": {"A"}},
            transcript_paths=[],
            annotation_view="canonical",
        )
        role = self._by_layer(rows, "A")["exon_role"]
        self.assertEqual(role["state"], "unknown")
        self.assertEqual(role["observation_mask"], "missing")
        self.assertIn("canonical_transcript_path_unavailable", role["evidence"])

    def test_noncanonical_completion_overlap_is_candidate_in_canonical_view(self):
        occurrences, elements = self._base_element()
        completion = [{
            "evidence_id": "evB", "homology_id": "H1", "family_id": "fam",
            "species": "B", "gene_copy_id": "gB", "homologous_dna_presence": "present",
            "homologous_dna_evidence": "local_dna_alignment",
            "candidate_resolution_status": "resolved",
            "dna_aligned_blocks": json.dumps([{
                "block_id": "db1", "target_contig": "chrB", "target_start": 100,
                "target_end": 129, "target_strand": "+",
            }]),
            "supplied_annotation_overlaps": json.dumps([{
                "role": "CDS", "strand_relation": "sense", "contains_block": True,
            }]),
        }]
        paths = [{
            "family_id": "fam", "species": "B", "gene_copy_id": "gB",
            "transcript_id": "tx_noncanonical", "occurrence_id": "B_ex",
            "path_status": "annotated_transcript_path",
        }]
        rows = _element_site_rows(
            occurrences,
            elements,
            completion,
            {"fam"},
            {"fam": {"A", "B"}},
            transcript_paths=paths,
            annotation_view="canonical",
        )
        role = self._by_layer(rows, "B")["exon_role"]
        self.assertEqual(role["state"], "unknown")
        self.assertEqual(role["annotation_view"], "canonical")
        self.assertNotIn("supplied_repertoire_exonic_overlap", role["evidence"])

    def test_legacy_membership_without_blocks_cannot_create_local_conflict(self):
        occurrences, elements = self._base_element()
        elements[0].pop("genomic_matched_blocks")
        occurrences[0]["role"] = "intron"
        completion = [{
            "evidence_id": "pred", "homology_id": "H1", "family_id": "fam", "species": "A",
            "homologous_dna_presence": "present", "homologous_dna_evidence": "protein_coding_projection",
            "candidate_resolution_status": "resolved", "predicted_role": "CDS",
            "predicted_role_blocks": json.dumps([{
                "target_contig": "chrA", "target_start": 1, "target_end": 5,
            }]),
        }]
        rows = _element_site_rows(occurrences, elements, completion, {"fam"}, {"fam": {"A"}})
        role = self._by_layer(rows, "A")["exon_role"]
        self.assertEqual(role["state"], "unknown")
        self.assertEqual(role["observation_reason"], "actual_local_blocks_unavailable")

    def test_ambiguous_position_does_not_create_junction(self):
        fixture = self._junction_fixture(ambiguous_middle=True)
        rows = _junction_site_rows(*fixture, {"fam"}, {"fam": {"R", "S", "C", "P", "D"}})
        self.assertEqual(rows, [])
        diagnostics = read_tsv(fixture[1] / "splice_boundary_correspondence.tsv")
        self.assertTrue(diagnostics)
        self.assertTrue(all(row["position_edge_eligible"] == "0" for row in diagnostics if row["species"] == "S"))

    def test_one_to_three_junctions_share_linked_group_and_keep_distinct_cutpoints(self):
        fixture = self._junction_fixture()
        rows = _junction_site_rows(*fixture, {"fam"}, {"fam": {"R", "S", "C", "P", "D"}})
        split_rows = [row for row in rows if row["species"] == "S"]
        self.assertEqual(len(split_rows), 2)
        self.assertEqual({row["state"] for row in split_rows}, {"present"})
        self.assertEqual(len({row["linked_group_id"] for row in split_rows}), 1)
        self.assertTrue(all("_TX_split" in row["linked_group_id"] for row in split_rows))
        self.assertEqual(len({row["site_id"] for row in split_rows}), 2)

    def test_missing_tree_tip_is_materialized_as_unknown_with_a_reason(self):
        row = normalize_structural_site_row({
            "family_id": "fam", "layer": "exon_presence", "site_id": "E1",
            "species": "A", "state": "present", "state_0": "absent",
            "state_1": "present", "evidence": "fixture",
        })
        completed = _complete_tree_tip_observations([row], {"A", "B"})
        by_species = {item["species"]: item for item in completed}
        self.assertEqual(by_species["B"]["state"], "unknown")
        self.assertEqual(by_species["B"]["observation_mask"], "missing")
        self.assertEqual(
            by_species["B"]["observation_reason"],
            "tree_tip_missing_from_structural_evidence",
        )

    def test_junction_absence_requires_one_block_spanning_both_sides(self):
        fixture = self._junction_fixture()
        rows = _junction_site_rows(*fixture, {"fam"}, {"fam": {"R", "S", "C", "P", "D"}})
        by_species_site = {(row["species"], row["site_id"]): row for row in rows}
        site_ids = sorted({row["site_id"] for row in rows})
        self.assertEqual(by_species_site[("C", site_ids[0])]["state"], "absent")
        self.assertEqual(by_species_site[("C", site_ids[1])]["state"], "absent")
        self.assertEqual(by_species_site[("P", site_ids[0])]["state"], "unknown")
        self.assertEqual(by_species_site[("P", site_ids[1])]["state"], "unknown")

    def test_spanning_block_is_evaluated_per_acceptable_alignment_candidate(self):
        spanning = {
            "query_occurrence_id": "query", "subject_occurrence_id": "reference",
            "match_status": "mapped", "position_edge_eligible": "1",
            "candidate_resolution": "resolved", "alignment_strand": "+",
            "matched_blocks": "1-90:1-90;1-5:1-5",
        }
        self.assertTrue(
            _continuous_reference_block_covers(
                [spanning], "query", "reference", 30, 31
            )
        )
        discordant = {**spanning, "matched_blocks": "1-5:1-5"}
        self.assertFalse(
            _continuous_reference_block_covers(
                [spanning, discordant], "query", "reference", 30, 31
            )
        )

        retained_alternatives = {
            **spanning,
            "matched_blocks": "1-5:1-5",
            "retained_candidate_ids": "candidate_1;candidate_2",
            "candidate_assessments": json.dumps([
                {
                    "candidate_id": "candidate_1",
                    "aligned_blocks": [[1, 90, 1, 90], [1, 5, 1, 5]],
                },
                {
                    "candidate_id": "candidate_2",
                    "aligned_blocks": [[10, 80, 10, 80]],
                },
            ]),
        }
        self.assertTrue(
            _continuous_reference_block_covers(
                [retained_alternatives], "query", "reference", 30, 31
            )
        )
        retained_alternatives["candidate_assessments"] = json.dumps([
            {
                "candidate_id": "candidate_1",
                "aligned_blocks": [[1, 90, 1, 90]],
            },
            {
                "candidate_id": "candidate_2",
                "aligned_blocks": [[1, 5, 1, 5]],
            },
        ])
        self.assertFalse(
            _continuous_reference_block_covers(
                [retained_alternatives], "query", "reference", 30, 31
            )
        )

    def test_position_projection_uses_the_retained_candidate(self):
        row = {
            "query_occurrence_id": "query", "subject_occurrence_id": "reference",
            "match_status": "mapped", "position_edge_eligible": "1",
            "candidate_resolution": "resolved", "alignment_strand": "+",
            "matched_blocks": "1-5:1-5",
            "retained_candidate_ids": "candidate_2",
            "candidate_assessments": json.dumps([
                {"candidate_id": "candidate_1", "aligned_blocks": [[1, 5, 1, 5]]},
                {"candidate_id": "candidate_2", "aligned_blocks": [[1, 90, 1, 90]]},
            ]),
        }
        self.assertEqual(
            _parse_projected_reference_blocks(row),
            [{
                "query_start": 1, "query_end": 90,
                "target_start": 1, "target_end": 90, "strand": "+",
            }],
        )

    def test_reverse_relative_alignment_projects_both_directions(self):
        query_to_reference = {
            "query_occurrence_id": "query", "subject_occurrence_id": "reference",
            "match_status": "mapped", "position_edge_eligible": "1",
            "candidate_resolution": "resolved", "alignment_strand": "-",
            "matched_blocks": "1-5:11-15",
        }
        self.assertEqual(
            _mapped_reference_coordinate(
                [query_to_reference], "query", "reference", 1
            ),
            15,
        )
        reference_to_query = {
            **query_to_reference,
            "query_occurrence_id": "reference",
            "subject_occurrence_id": "query",
        }
        self.assertEqual(
            _mapped_reference_coordinate(
                [reference_to_query], "query", "reference", 11
            ),
            5,
        )

    def test_repeated_overlap_stays_unknown(self):
        fixture = self._junction_fixture()
        rows = _junction_site_rows(*fixture, {"fam"}, {"fam": {"R", "S", "C", "P", "D"}})
        repeated = [row for row in rows if row["species"] == "D"]
        self.assertEqual({row["state"] for row in repeated}, {"unknown"})
        self.assertTrue(all("repeated_overlap_position_ambiguous" in row["evidence"] for row in repeated))
