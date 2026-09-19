import json
import tempfile
import unittest
from pathlib import Path

from insiphy.io import normalize_structural_site_row, read_tsv
from insiphy.structural_sites import (
    _complete_tree_tip_observations,
    _continuous_reference_block_covers,
    _element_site_rows,
    _junction_site_rows,
    _mapped_reference_coordinate,
    _merge_catalogue_observations,
    _parse_projected_reference_blocks,
)


class StructuralSiteSemanticsTests(unittest.TestCase):
    def _base_element(self):
        occurrences = [{
            "occurrence_id": "A_ex", "family_id": "fam", "species": "A",
            "gene_copy_id": "gA", "role": "CDS", "presence_status": "present",
            "contig": "chrA", "start": "1", "end": "30", "strand": "+",
            "source_feature_id": "exA",
        }]
        elements = [{
            "element_id": "EG1", "homology_id": "H1", "occurrence_id": "A_ex",
            "family_id": "fam", "species": "A", "gene_copy_id": "gA",
            "element_class": "exon_like", "membership_call": "core_member",
            "membership_edge_eligible": "1", "position_edge_eligible": "1",
            "correspondence_status": "resolved",
            "genomic_matched_blocks": "chrA:1-30:+", "parent_feature_ids": "exA",
        }]
        return occurrences, elements

    @staticmethod
    def _by_layer(rows, species):
        return {row["layer"]: row for row in rows if row["species"] == species}

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

    def _prediction_overlap_case(self, prediction_start, prediction_end, prediction_strand="+"):
        occurrences, elements = self._base_element()
        occurrences.append({
            "occurrence_id": "B_in", "family_id": "fam", "species": "B",
            "gene_copy_id": "gB", "role": "intron", "presence_status": "present",
            "contig": "chrB", "start": "50", "end": "250", "strand": "+",
        })
        elements.append({
            "element_id": "EG1", "homology_id": "H1", "occurrence_id": "B_in",
            "element_class": "candidate_source", "membership_call": "core_member",
            "membership_edge_eligible": "1", "position_edge_eligible": "1",
            "correspondence_status": "resolved", "genomic_matched_blocks": "chrB:100-120:+",
        })
        completion = [{
            "evidence_id": "predB", "homology_id": "H1", "family_id": "fam",
            "species": "B", "gene_copy_id": "gB", "homologous_dna_presence": "present",
            "homologous_dna_evidence": "protein_coding_projection",
            "candidate_resolution_status": "resolved", "predicted_role": "CDS",
            "predicted_role_blocks": json.dumps([{
                "block_id": "p", "target_contig": "chrB", "target_start": prediction_start,
                "target_end": prediction_end, "target_strand": prediction_strand,
            }]),
        }]
        rows = _element_site_rows(occurrences, elements, completion, {"fam"}, {"fam": {"A", "B"}})
        return self._by_layer(rows, "B")["exon_role"]

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

    def _junction_fixture(self, ambiguous_middle=False, missing_phase=False):
        root = Path(self.tempdir.name)
        input_dir, output_dir = root / "input", root / "output"
        input_dir.mkdir(); output_dir.mkdir()
        phase = "." if missing_phase else "0"
        paths = [
            ("R", "gR", "ref_tx", 1, "ref", "CDS", "0"),
            ("S", "gS", "split", 1, "left", "CDS", phase),
            ("S", "gS", "split", 2, "i1", "intron", "."),
            ("S", "gS", "split", 3, "middle", "CDS", phase),
            ("S", "gS", "split", 4, "i2", "intron", "."),
            ("S", "gS", "split", 5, "right", "CDS", phase),
            ("C", "gC", "continuous", 1, "continuous", "CDS", "0"),
            ("P", "gP", "partial", 1, "partial", "CDS", "0"),
            ("D", "gD", "repeat", 1, "repeat", "CDS", "0"),
        ]
        header = "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\tphase\tpath_status\n"
        lines = [header] + [f"p{index}\tfam\t{sp}\t{gene}\t{tx}\t{rank}\t{occ}\t{role}\t{ph}\tannotated_transcript_path\n"
                            for index, (sp, gene, tx, rank, occ, role, ph) in enumerate(paths)]
        (input_dir / "transcript_paths.tsv").write_text("".join(lines))
        match_rows = [
            ("m1", "left", "ref", "1-30:1-30", "1", "resolved"),
            ("m2", "middle", "ref", "1-30:31-60", "0" if ambiguous_middle else "1",
             "ambiguous" if ambiguous_middle else "resolved"),
            ("m3", "right", "ref", "1-30:61-90", "1", "resolved"),
            ("m4", "continuous", "ref", "1-90:1-90;1-5:1-5", "1", "resolved"),
            ("m5", "partial", "ref", "1-30:1-30", "1", "resolved"),
            ("m6", "repeat", "ref", "1-90:1-90", "1", "resolved"),
        ]
        (output_dir / "segment_matches.tsv").write_text(
            "match_id\tquery_occurrence_id\tsubject_occurrence_id\tmatch_status\tmatched_blocks\talignment_strand\tposition_edge_eligible\tcandidate_resolution\n" +
            "".join(f"{mid}\t{query}\t{subject}\tmapped\t{blocks}\t+\t{eligible}\t{resolution}\n"
                    for mid, query, subject, blocks, eligible, resolution in match_rows))
        occurrences = []
        coordinates = {"ref": (1, 90), "left": (1, 30), "i1": (31, 40),
                       "middle": (41, 70), "i2": (71, 80), "right": (81, 110),
                       "continuous": (1, 90), "partial": (1, 30), "repeat": (1, 90)}
        species = {"ref": ("R", "gR"), "left": ("S", "gS"), "i1": ("S", "gS"),
                   "middle": ("S", "gS"), "i2": ("S", "gS"), "right": ("S", "gS"),
                   "continuous": ("C", "gC"), "partial": ("P", "gP"), "repeat": ("D", "gD")}
        for occurrence_id, (start, end) in coordinates.items():
            sp, gene = species[occurrence_id]
            occurrences.append({
                "occurrence_id": occurrence_id, "family_id": "fam", "species": sp,
                "gene_copy_id": gene, "role": "intron" if occurrence_id.startswith("i") else "CDS",
                "presence_status": "present", "contig": f"chr{sp}", "start": str(start),
                "end": str(end), "strand": "+", "phase": phase if sp == "S" else "0",
            })
        elements = []
        for occurrence_id in ("ref", "left", "middle", "right", "continuous", "partial", "repeat"):
            relation = "repeated_overlap" if occurrence_id == "repeat" else (
                "complementary" if occurrence_id in {"left", "middle", "right"} else "single")
            elements.append({
                "element_id": "EG1", "occurrence_id": occurrence_id,
                "element_class": "exon_like", "membership_call": "core_member",
                "position_edge_eligible": "0" if occurrence_id == "middle" and ambiguous_middle else "1",
                "correspondence_status": "ambiguous" if occurrence_id == "middle" and ambiguous_middle else "resolved",
                "reference_occurrence_id": "ref", "reference_coverage_relation": relation,
            })
        return input_dir, output_dir, occurrences, elements

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

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

    def test_missing_coding_phase_is_recorded_without_erasing_exact_junction(self):
        fixture = self._junction_fixture(missing_phase=True)
        rows = _junction_site_rows(*fixture, {"fam"}, {"fam": {"R", "S", "C", "P", "D"}})
        split_rows = [row for row in rows if row["species"] == "S"]
        self.assertEqual({row["state"] for row in split_rows}, {"present"})
        self.assertTrue(all("junction_phase_unknown" in row["evidence"] for row in split_rows))

    def test_canonical_view_excludes_noncanonical_junction_paths(self):
        fixture = self._junction_fixture()
        path = fixture[0] / "transcript_paths.tsv"
        lines = path.read_text().splitlines()
        header = lines[0].split("\t")
        status_index = header.index("path_status")
        transcript_index = header.index("transcript_id")
        rewritten = [lines[0]]
        for line in lines[1:]:
            fields = line.split("\t")
            if fields[transcript_index] == "split":
                fields[status_index] = "canonical_transcript_path"
            rewritten.append("\t".join(fields))
        path.write_text("\n".join(rewritten) + "\n")

        rows = _junction_site_rows(
            *fixture,
            {"fam"},
            {"fam": {"R", "S", "C", "P", "D"}},
            annotation_view="canonical",
        )
        self.assertTrue(rows)
        self.assertEqual(
            {row["state"] for row in rows if row["species"] == "S"}, {"present"}
        )
        self.assertEqual(
            {row["state"] for row in rows if row["species"] == "C"}, {"unknown"}
        )
        self.assertEqual(
            {row["transcript_scope"] for row in rows},
            {"explicit_canonical_transcript"},
        )

    def test_site_only_complete_universe_entry_does_not_create_all_zero_state(self):
        rows = _merge_catalogue_observations(
            [],
            [{"family_id": "fam", "layer": "exon_presence", "site_id": "S_only"},
             {"family_id": "fam", "layer": "exon_presence", "site_id": "S_explicit",
              "species": "A", "state": "absent", "discovery_rule": "independent_catalogue",
              "observation_mask": "explicit_all_zero"},
             {"family_id": "fam", "layer": "exon_presence", "site_id": "S_generated",
              "species": "B", "state": "absent", "observation_mask": "generated_all_zero"}],
            {"fam"}, {"fam": {"A", "B"}}, annotation_view="repertoire",
        )
        by_site = {row["site_id"]: row for row in rows}
        self.assertEqual(set(by_site), {"S_explicit", "S_generated"})
        self.assertEqual(by_site["S_explicit"]["discovery_rule"], "independent_catalogue")
        self.assertEqual(by_site["S_generated"]["discovery_rule"], "independent_catalogue")
        self.assertEqual(by_site["S_explicit"]["observation_mask"], "explicit_all_zero")
        self.assertEqual(by_site["S_generated"]["observation_mask"], "generated_all_zero")

    def test_catalogue_role_view_must_match_analysis_and_is_preserved(self):
        base = {
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "not_exonic", "state_0": "not_exonic",
            "state_1": "exonic", "annotation_view": "canonical",
            "transcript_scope": "explicit_canonical_transcript",
        }
        with self.assertRaisesRegex(ValueError, "conflicts with this analysis"):
            _merge_catalogue_observations(
                [], [base], {"fam"}, {"fam": {"A"}}, annotation_view="repertoire"
            )

        repertoire = {
            **base,
            "annotation_view": "repertoire",
            "transcript_scope": "annotated_transcript_repertoire",
        }
        with self.assertRaisesRegex(ValueError, "conflicts with this analysis"):
            _merge_catalogue_observations(
                [], [repertoire], {"fam"}, {"fam": {"A"}}, annotation_view="canonical"
            )

        rows = _merge_catalogue_observations(
            [], [base], {"fam"}, {"fam": {"A"}}, annotation_view="canonical"
        )
        self.assertEqual(rows[0]["annotation_view"], "canonical")

        conflicting_scope = {
            **base,
            "annotation_view": "repertoire",
        }
        with self.assertRaisesRegex(ValueError, "conflicts with transcript_scope"):
            _merge_catalogue_observations(
                [], [conflicting_scope], {"fam"}, {"fam": {"A"}},
                annotation_view="repertoire",
            )

    def test_catalogue_role_without_explicit_or_inferable_view_is_rejected(self):
        row = {
            "family_id": "fam", "layer": "splice_junction", "site_id": "J1",
            "state": "present", "state_0": "absent",
            "state_1": "present", "transcript_scope": "catalogue_defined",
        }
        with self.assertRaisesRegex(ValueError, "require an explicit or inferable"):
            _merge_catalogue_observations(
                [], [row], {"fam"}, {"fam": {"A"}}, annotation_view="canonical"
            )

    def test_catalogue_unknown_mask_is_not_replaced_with_observed(self):
        row = {
            "family_id": "fam", "layer": "exon_role", "site_id": "S1",
            "species": "A", "state": "exonic", "state_0": "not_exonic",
            "state_1": "exonic", "annotation_view": "repertoire",
            "transcript_scope": "annotated_transcript_repertoire",
            "observation_mask": "ambiguous",
        }
        merged = _merge_catalogue_observations(
            [], [row], {"fam"}, {"fam": {"A"}}, annotation_view="repertoire"
        )
        with self.assertRaisesRegex(ValueError, "invalid observation_mask"):
            normalize_structural_site_row(merged[0])


if __name__ == "__main__":
    unittest.main()
