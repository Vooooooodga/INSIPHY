import json

import tempfile

import unittest

from pathlib import Path

from unittest.mock import patch

from intraphy.correspondence import _membership_match_details, infer_correspondence, match_total, tree_distances

from intraphy.coordinates import parse_legacy_blocks

from intraphy.elements import collect_element_profiles

from intraphy.io import read_tsv, write_tsv

from intraphy.cli import main as cli_main

from intraphy.preprocess import MATCH_FIELDS, _apply_ordered_candidate_chains, _species_tree_distances, assess_short_candidate_thresholds, cheap_match_evidence, cluster_segments, copy_order_context, extract_gene, graph_components, introns_from_path, match_evidence, read_annotation_for_gene, transcript_cds_length

from intraphy.alignment import AlignmentStats, phase_compatibility

from intraphy.structural_sites import (
    _completion_presence,
    _element_site_rows,
    _genomically_contiguous,
    _junction_site_rows,
)

from support_observations_observation_semantics import ObservationSemanticsTestsSupport

class ObservationSemanticsTests(ObservationSemanticsTestsSupport, unittest.TestCase):
    def test_same_locus_alternative_exons_share_element_without_merging_tandem_repeat(self):
        occurrences = [
            {
                "occurrence_id": "alt_short",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx1",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "3535611",
                "end": "3536164",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
            {
                "occurrence_id": "alt_long",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx2",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "3535611",
                "end": "3536182",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
            {
                "occurrence_id": "nearby_repeat",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx3",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "3537000",
                "end": "3537553",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
        ]

        homology, matches = cluster_segments(occurrences, {}, threads=1)
        by_occ = {row["occurrence_id"]: row["homology_id"] for row in homology}
        overlap_match = [row for row in matches if {row["query_occurrence_id"], row["subject_occurrence_id"]} == {"alt_short", "alt_long"}]

        self.assertEqual(by_occ["alt_short"], by_occ["alt_long"])
        self.assertNotEqual(by_occ["alt_short"], by_occ["nearby_repeat"])
        self.assertEqual(overlap_match[0]["alignment_backend"], "genomic_overlap")
        self.assertNotEqual(overlap_match[0]["projected_reference_blocks"], "NA")

    def test_same_locus_overlap_with_shared_transcript_membership_is_not_alternative_evidence(self):
        occurrences = [
            {
                "occurrence_id": "left",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx1;tx_shared",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "100",
                "end": "180",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
            {
                "occurrence_id": "right",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx2;tx_shared",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "120",
                "end": "200",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
        ]

        homology, matches = cluster_segments(occurrences, {}, threads=1)
        by_occ = {row["occurrence_id"]: row["homology_id"] for row in homology}

        self.assertNotEqual(by_occ["left"], by_occ["right"])
        self.assertFalse(matches)

    def test_bounded_short_context_retains_candidates_without_claiming_absence(self):
        common = {
            "family_id": "fam",
            "presence_status": "present",
            "contig": "chr1",
            "strand": "+",
            "phase": ".",
            "splice_motif_score": "0.5",
        }
        left = {
            **common,
            "occurrence_id": "left",
            "species": "A",
            "gene_copy_id": "gA",
            "transcript_id": "txA",
            "source_feature_id": "featureA",
            "role": "CDS",
            "start": "101",
            "end": "124",
        }
        right = {
            **common,
            "occurrence_id": "right",
            "species": "B",
            "gene_copy_id": "gB",
            "transcript_id": "txB",
            "source_feature_id": "featureB",
            "role": "intron",
            "start": "501",
            "end": "524",
        }
        sequence = "ACGTTGCAAGTCCTGATCGTACGA"
        evidence = match_evidence(
            left,
            right,
            {"left": sequence, "right": sequence},
            {"left": {}, "right": {}},
            context_aligner="minimap2",
        )

        self.assertEqual(evidence["alignment_requested_backend"], "anchored_short_alignment")
        self.assertEqual(evidence["score_scheme"], "nt_blastn_v1")
        self.assertEqual(evidence["short_context_route"], "bounded_local")
        self.assertTrue(evidence["candidate_records"])
        self.assertNotEqual(evidence["matched_blocks"], "NA")
        self.assertEqual(evidence["query_parent_feature_ids"], "featureA")
        self.assertEqual(evidence["subject_transcript_ids"], "txB")
        self.assertEqual(evidence["true_absence_eligible"], 0)

    def test_candidate_chain_fields_keep_structural_scores_diagnostic(self):
        sequence = "ACGTTGCAAGTCCTGATCGTACGA"
        occurrences = [
            {
                "occurrence_id": "left",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "txA",
                "source_feature_id": "featureA",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chrA",
                "start": "101",
                "end": "124",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
            {
                "occurrence_id": "right",
                "family_id": "fam",
                "species": "B",
                "gene_copy_id": "gB",
                "transcript_id": "txB",
                "source_feature_id": "featureB",
                "role": "intron",
                "presence_status": "present",
                "contig": "chrB",
                "start": "501",
                "end": "524",
                "strand": "+",
                "phase": ".",
                "splice_motif_score": "0.5",
            },
        ]
        emitted = []
        _homology, matches = cluster_segments(
            occurrences,
            {"left": sequence, "right": sequence},
            context_aligner="minimap2",
            threads=1,
            match_writer=emitted.append,
        )

        self.assertFalse(matches)
        self.assertEqual(len(emitted), 1)
        match = emitted[0]
        self.assertEqual(match["match_status"], "candidate_unanchored")
        self.assertEqual(match["candidate_resolution"], "candidate")
        self.assertEqual(match["score_scheme"], "nt_blastn_v1")
        self.assertEqual(match["correspondence_score"], match["sequence_score"])
        self.assertNotEqual(match["candidate_ids"], "NA")
        self.assertNotEqual(match["query_genomic_matched_blocks"], "NA")
        self.assertEqual(match["flanking_anchor_status"], "no_independent_homologous_flanks_on_same_path")
        self.assertEqual(match["membership_edge_eligible"], 0)
        self.assertEqual(match["position_edge_eligible"], 0)
        self.assertEqual(match["true_absence_eligible"], 0)

    def test_short_query_uses_explicit_long_target_but_skips_unbounded_locus(self):
        query_sequence = "ACGTTGCAAGTCCTGATCGTACGA"
        target_sequence = "T" * 180 + query_sequence + "T" * 196
        query = {
            "occurrence_id": "query", "role": "CDS", "contig": "chrQ",
            "start": "1", "end": str(len(query_sequence)), "strand": "+",
            "source_feature_id": "query_feature",
        }
        bounded = {
            "occurrence_id": "bounded", "role": "intron", "contig": "chrT",
            "start": "1001", "end": str(1000 + len(target_sequence)), "strand": "+",
            "boundary_class": "annotated_segment",
        }
        bounded_evidence = match_evidence(
            query,
            bounded,
            {"query": query_sequence, "bounded": target_sequence},
            {"query": {}, "bounded": {}},
            context_aligner="internal",
        )
        right_short = match_evidence(
            {**bounded, "occurrence_id": "bounded_left"},
            {**query, "occurrence_id": "short_right"},
            {"bounded_left": target_sequence, "short_right": query_sequence},
            {"bounded_left": {}, "short_right": {}},
            context_aligner="internal",
        )
        unbounded = {**bounded, "occurrence_id": "locus", "boundary_class": "whole_locus_search_interval"}
        unbounded_evidence = match_evidence(
            query,
            unbounded,
            {"query": query_sequence, "locus": target_sequence},
            {"query": {}, "locus": {}},
            context_aligner="internal",
        )

        self.assertEqual(bounded_evidence["alignment_requested_backend"], "anchored_short_alignment")
        self.assertEqual(bounded_evidence["short_context_route"], "bounded_local")
        self.assertEqual(bounded_evidence["query_length"], len(query_sequence))
        self.assertEqual(bounded_evidence["target_length"], len(target_sequence))
        self.assertEqual(right_short["alignment_requested_backend"], "anchored_short_alignment")
        self.assertEqual(right_short["alignment_input_transposed"], 1)
        self.assertEqual(right_short["query_length"], len(target_sequence))
        self.assertEqual(right_short["target_length"], len(query_sequence))
        self.assertEqual(right_short["query_mapped_contig"], "chrT")
        self.assertEqual(right_short["subject_mapped_contig"], "chrQ")
        right_candidate = right_short["candidate_records"][0]
        self.assertEqual(right_candidate["query_length"], len(target_sequence))
        self.assertEqual(right_candidate["target_length"], len(query_sequence))
        self.assertEqual(right_candidate["query_occurrence_id"], "bounded_left")
        self.assertEqual(right_candidate["target_occurrence_id"], "short_right")
        self.assertEqual(right_candidate["short_sequence_coverage"], 1.0)
        self.assertEqual(unbounded_evidence["alignment_requested_backend"], "internal")
        self.assertEqual(unbounded_evidence["short_context_route"], "not_used")
        self.assertEqual(unbounded_evidence["search_interval"], "NA")
        self.assertEqual(unbounded_evidence["candidate_records"][0]["search_interval"], "NA")

    def test_match_evidence_contract_survives_primary_and_candidate_rows(self):
        sequence = "ACGTTGCAAGTCCTGATCGTACGA"
        common = {
            "family_id": "fam", "role": "CDS", "presence_status": "present",
            "strand": "+", "phase": "0", "splice_motif_score": "0.5",
        }
        occurrences = [
            {**common, "occurrence_id": "q", "species": "A", "gene_copy_id": "gA",
             "transcript_id": "qtx", "source_feature_id": "qf", "contig": "chrQ",
             "start": "1", "end": str(len(sequence))},
            {**common, "occurrence_id": "t", "species": "B", "gene_copy_id": "gB",
             "transcript_id": "ttx", "source_feature_id": "tf", "contig": "chrT",
             "start": "101", "end": str(100 + len(sequence))},
        ]
        emitted = []
        cluster_segments(
            occurrences,
            {"q": sequence, "t": sequence},
            context_aligner="internal",
            match_writer=emitted.append,
        )
        row = emitted[0]
        candidate = json.loads(row["candidate_assessments"])[0]

        for field in (
            "candidate_id", "alternative_candidate_ids", "aligned_blocks", "gap_blocks",
            "sequence_kind", "backend", "backend_version", "raw_score", "nt_identity",
            "aa_identity", "known_aligned_pairs", "unknown_aligned_pairs",
            "query_covered_bases", "target_covered_bases", "query_length", "target_length",
            "relative_strand", "mapq", "left_anchor_id", "right_anchor_id",
            "search_interval", "enumeration_complete", "incomplete_reason",
            "chain_configuration", "chain_score_delta",
        ):
            self.assertIn(field, MATCH_FIELDS)
            self.assertIn(field, row)
        self.assertEqual(row["candidate_id"], candidate["candidate_id"])
        self.assertEqual(candidate["sequence_kind"], "nucleotide")
        self.assertEqual(candidate["query_length"], len(sequence))
        self.assertEqual(candidate["target_length"], len(sequence))
        self.assertEqual(candidate["known_aligned_pairs"] + candidate["unknown_aligned_pairs"], candidate["aligned_pairs"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "segment_matches.tsv"
            write_tsv(path, [row], MATCH_FIELDS)
            round_tripped = read_tsv(path)[0]
        for field in (
            "candidate_id", "alternative_candidate_ids", "aligned_blocks", "gap_blocks",
            "sequence_kind", "backend", "backend_version", "score_scheme", "raw_score",
            "nt_identity", "aa_identity", "known_aligned_pairs", "unknown_aligned_pairs",
            "query_covered_bases", "target_covered_bases", "query_length", "target_length",
            "relative_strand", "mapq", "left_anchor_id", "right_anchor_id",
            "search_interval", "enumeration_complete", "incomplete_reason",
            "chain_configuration", "chain_score_delta",
        ):
            self.assertEqual(round_tripped[field], str(row[field]))

    def test_multiple_retained_positions_keep_membership_without_hard_coordinates(self):
        occurrences = [
            {"occurrence_id": "query", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "contig": "chrA", "start": "1", "end": "100", "strand": "+"},
            {"occurrence_id": "subject", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "contig": "chrB", "start": "1", "end": "200", "strand": "+"},
        ]
        row = {
            "match_id": "match_00001",
            "query_occurrence_id": "query",
            "subject_occurrence_id": "subject",
            "match_status": "mapped",
            "enumeration_complete": 1,
            "candidate_enumeration_status": "complete",
            "short_context_route": "not_used",
            "protein_hard_observation_eligible": 0,
            "_candidate_records": [
                {"candidate_id": "match_00001.candidate_001", "accepted": 1, "query_start": 1, "query_end": 40, "target_start": 1, "target_end": 40, "score": 80, "score_scheme": "nt_blastn_v1", "strand": "+", "aligned_blocks": [(1, 40, 1, 40)]},
                {"candidate_id": "match_00001.candidate_002", "accepted": 1, "query_start": 1, "query_end": 40, "target_start": 101, "target_end": 140, "score": 80, "score_scheme": "nt_blastn_v1", "strand": "+", "aligned_blocks": [(1, 40, 101, 140)]},
            ],
        }

        _apply_ordered_candidate_chains([row], {item["occurrence_id"]: item for item in occurrences}, occurrences)

        self.assertEqual(row["membership_edge_eligible"], 1)
        self.assertEqual(row["position_edge_eligible"], 0)
        self.assertEqual(row["candidate_resolution"], "ambiguous")
        self.assertEqual(row["position_edge_reason"], "multiple_accepted_retained_coordinate_placements")
        self.assertEqual(row["true_absence_eligible"], 0)

    def test_external_overlap_score_and_candidate_enumeration_are_explicit(self):
        sequence = "A" * 301
        left = {"occurrence_id": "left", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "contig": "chrA", "start": "1", "end": "301", "strand": "+"}
        right = {"occurrence_id": "right", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "role": "CDS", "contig": "chrB", "start": "1", "end": "301", "strand": "+"}
        stats = AlignmentStats(
            identity=1.0,
            coverage=1.0,
            score=0.0,
            query_start=1,
            query_end=301,
            target_start=1,
            target_end=301,
            cigar="301=",
            backend="mafft_overlap",
            query_coverage=1.0,
            target_coverage=1.0,
            aligned_pairs=301,
            aligned_blocks=[(1, 301, 1, 301)],
            matches=301,
            alignment_mode="overlap_projection",
        )
        with patch("intraphy.mapping.pairwise_matches.overlap_alignment_stats", return_value=stats):
            evidence = match_evidence(left, right, {"left": sequence, "right": sequence}, {"left": {}, "right": {}})

        self.assertEqual(evidence["raw_alignment_score"], 602.0)
        self.assertEqual(evidence["score_scheme"], "nt_blastn_v1/sum_of_column_scores")
        self.assertFalse(evidence["enumeration_complete"])
        self.assertEqual(evidence["candidate_enumeration_status"], "unassessed")

    def test_ambiguous_membership_does_not_publish_element_coordinates(self):
        occurrences = [
            {"occurrence_id": "query", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "source_feature_id": "qf", "transcript_id": "qt", "start": "1", "end": "40"},
            {"occurrence_id": "subject", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "source_feature_id": "sf", "transcript_id": "st", "start": "1", "end": "100"},
        ]
        element_rows = [
            {"element_id": "EG1", "occurrence_id": row["occurrence_id"], "species": row["species"], "gene_copy_id": row["gene_copy_id"]}
            for row in occurrences
        ]
        scored = [{
            "match_id": "match_00001",
            "query_occurrence_id": "query",
            "subject_occurrence_id": "subject",
            "match_status": "mapped",
            "candidate_resolution": "ambiguous",
            "membership_edge_eligible": 1,
            "membership_edge_reason": "retained_sequence_membership",
            "position_edge_eligible": 0,
            "position_edge_reason": "multiple_accepted_retained_coordinate_placements",
            "candidate_ids": "match_00001.candidate_001;match_00001.candidate_002",
            "retained_candidate_ids": "match_00001.candidate_001;match_00001.candidate_002",
            "matched_blocks": "1-40:1-40",
            "query_genomic_matched_blocks": "chrA:1-40:+",
            "subject_genomic_matched_blocks": "chrB:1-40:+",
        }]

        _membership_match_details(element_rows, occurrences, scored)
        query_row = next(row for row in element_rows if row["occurrence_id"] == "query")

        self.assertEqual(query_row["membership_edge_eligible"], 1)
        self.assertEqual(query_row["position_edge_eligible"], 0)
        self.assertEqual(query_row["matched_blocks"], "NA")
        self.assertEqual(query_row["genomic_matched_blocks"], "NA")
        self.assertEqual(query_row["correspondence_status"], "ambiguous")
