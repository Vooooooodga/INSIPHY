import json

import tempfile

import unittest

from pathlib import Path

from types import SimpleNamespace

from unittest.mock import patch

from intraphy.alignment import AlignmentStats, local_alignment_stats

from intraphy.annotation import (
    _cached_protein_projection,
    _flank_support_from_spanning_blocks,
    _ordered_anchor_deletion_support,
    _protein_projection_from_rows,
    _reference_protein_context,
    complete_annotation,
    completion_call,
    generate_sequence_evidence,
)

from intraphy.io import read_tsv, write_tsv as write_rows

from intraphy.candidate_chain import ChainCandidate, ChainPathMembership, ordered_candidate_chain

from intraphy.coordinates import CoordinateBlock, Interval0

from intraphy.correspondence import _membership_match_details

from intraphy.preprocess import _rerun_anchor_bounded_short_candidates

def write_tsv(path, header, rows):
    path.write_text("\n".join(["\t".join(header), *["\t".join(row) for row in rows]]) + "\n")

from support_evidence_local_correspondence_contract import LocalCorrespondenceContractTestsSupport

class LocalCorrespondenceContractTests(LocalCorrespondenceContractTestsSupport, unittest.TestCase):
    def test_identical_path_in_two_transcript_contexts_is_one_mapping(self):
        context_a = ChainPathMembership(
            "q_tx_a", "t_tx_a", 1, 1, "chrQ", "chrT", "+", "+",
        )
        context_b = ChainPathMembership(
            "q_tx_b", "t_tx_b", 1, 1, "chrQ", "chrT", "+", "+",
        )
        candidates = (
            ChainCandidate(
                "left", Interval0(0, 4), Interval0(0, 4), 5.0,
                "nt_blastn_v1",
                path_memberships=(context_a, context_b),
            ),
            ChainCandidate(
                "right", Interval0(4, 8), Interval0(4, 8), 5.0,
                "nt_blastn_v1",
                path_memberships=(
                    ChainPathMembership(
                        "q_tx_a", "t_tx_a", 2, 2,
                        "chrQ", "chrT", "+", "+",
                    ),
                    ChainPathMembership(
                        "q_tx_b", "t_tx_b", 2, 2,
                        "chrQ", "chrT", "+", "+",
                    ),
                ),
            ),
        )
        result = ordered_candidate_chain(
            candidates,
            score_delta=0.0,
            context_anchor_ids={
                context_a.context: ({"left"}, {"right"}),
                context_b.context: ({"left"}, {"right"}),
            },
        )
        self.assertEqual(len(result.context_summaries), 2)
        self.assertEqual(result.best_path_count_capped, 1)
        self.assertEqual(result.near_optimal_path_count_capped, 1)
        # Each context has one path, but the global candidate set has two
        # alternatives. No node occurs in every retained context.
        self.assertEqual(result.ambiguous_ids, frozenset())
        self.assertEqual(result.ambiguity_status, "unique_within_reported_candidates")

    def test_disconnected_equal_complete_chains_are_ambiguous(self):
        context = ChainPathMembership(
            "q_tx", "t_tx", 1, 1, "chrQ", "chrT", "+", "+",
        ).context

        def candidate(name, query_start0, query_end0, target_start0, target_end0, order):
            return ChainCandidate(
                name,
                Interval0(query_start0, query_end0),
                Interval0(target_start0, target_end0),
                5.0,
                "nt_blastn_v1",
                path_memberships=(ChainPathMembership(
                    "q_tx", "t_tx", order, order,
                    "chrQ", "chrT", "+", "+",
                ),),
            )

        candidates = (
            candidate("left_a", 0, 4, 20, 24, 1),
            candidate("right_a", 4, 8, 24, 28, 2),
            candidate("left_b", 20, 24, 0, 4, 1),
            candidate("right_b", 24, 28, 4, 8, 2),
        )
        result = ordered_candidate_chain(
            candidates,
            score_delta=0.0,
            context_anchor_ids={
                context: ({"left_a", "left_b"}, {"right_a", "right_b"}),
            },
        )
        self.assertEqual(result.best_score, 10.0)
        self.assertEqual(result.best_path_count_capped, 2)
        self.assertEqual(result.near_optimal_path_count_capped, 2)
        self.assertEqual(
            result.ambiguous_ids,
            frozenset({"left_a", "right_a", "left_b", "right_b"}),
        )
        self.assertEqual(result.ambiguity_status, "multiple_near_optimal_chains")

    def test_distinct_transcript_paths_are_scored_in_separate_contexts(self):
        context_a = ChainPathMembership(
            "q_tx_a", "t_tx_a", 1, 1, "chrQ", "chrT", "+", "+",
        ).context
        context_b = ChainPathMembership(
            "q_tx_b", "t_tx_b", 1, 1, "chrQ", "chrT", "+", "+",
        ).context
        candidates = (
            ChainCandidate(
                "a_left", Interval0(0, 2), Interval0(0, 2), 5.0,
                "nt_blastn_v1", path_memberships=(ChainPathMembership(
                    "q_tx_a", "t_tx_a", 1, 1,
                    "chrQ", "chrT", "+", "+",
                ),),
            ),
            ChainCandidate(
                "a_right", Interval0(2, 4), Interval0(2, 4), 5.0,
                "nt_blastn_v1", path_memberships=(ChainPathMembership(
                    "q_tx_a", "t_tx_a", 2, 2,
                    "chrQ", "chrT", "+", "+",
                ),),
            ),
            ChainCandidate(
                "b_left", Interval0(10, 12), Interval0(20, 22), 5.0,
                "nt_blastn_v1", path_memberships=(ChainPathMembership(
                    "q_tx_b", "t_tx_b", 1, 1,
                    "chrQ", "chrT", "+", "+",
                ),),
            ),
            ChainCandidate(
                "b_right", Interval0(12, 14), Interval0(22, 24), 5.0,
                "nt_blastn_v1", path_memberships=(ChainPathMembership(
                    "q_tx_b", "t_tx_b", 2, 2,
                    "chrQ", "chrT", "+", "+",
                ),),
            ),
        )
        result = ordered_candidate_chain(
            candidates,
            score_delta=0.0,
            context_anchor_ids={
                context_a: ({"a_left"}, {"a_right"}),
                context_b: ({"b_left"}, {"b_right"}),
            },
        )
        self.assertEqual(len(result.context_summaries), 2)
        self.assertEqual(result.best_path_count_capped, 1)
        self.assertEqual(result.ambiguous_ids, frozenset({"a_left", "a_right", "b_left", "b_right"}))
        self.assertEqual(result.ambiguity_status, "multiple_near_optimal_chains")

    def test_shared_flanks_remain_resolved_across_two_tied_internal_paths(self):
        context = ChainPathMembership(
            "q_tx", "t_tx", 1, 1, "chrQ", "chrT", "+", "+",
        ).context

        def candidate(name, query, target, order):
            return ChainCandidate(
                name, Interval0(*query), Interval0(*target), 5.0,
                "nt_blastn_v1",
                path_memberships=(ChainPathMembership(
                    "q_tx", "t_tx", order, order,
                    "chrQ", "chrT", "+", "+",
                ),),
            )

        result = ordered_candidate_chain(
            (
                candidate("left", (0, 2), (0, 2), 1),
                candidate("branch_a", (2, 4), (2, 4), 2),
                candidate("branch_b", (2, 4), (4, 6), 2),
                candidate("right", (4, 6), (6, 8), 3),
            ),
            score_delta=0.0,
            context_anchor_ids={context: ({"left"}, {"right"})},
        )
        self.assertEqual(result.best_path_count_capped, 2)
        self.assertEqual(
            result.ambiguous_ids,
            frozenset({"branch_a", "branch_b"}),
        )

    def test_disjoint_subintervals_make_two_memberships_and_localize_ambiguity(self):
        occurrences = [
            {"occurrence_id": "left", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "source_feature_id": "left_exon", "transcript_id": "txA"},
            {"occurrence_id": "right", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "source_feature_id": "right_exon", "transcript_id": "txA"},
            {"occurrence_id": "fused", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "source_feature_id": "fused_exon", "transcript_id": "txB"},
        ]
        element_rows = [
            {"element_id": "EG1", "homology_id": "HC1", "occurrence_id": row["occurrence_id"], "family_id": "fam", "species": row["species"], "gene_copy_id": row["gene_copy_id"], "confidence": "1", "membership_call": "core_member"}
            for row in occurrences
        ]

        def hard(match_id, query, subject, query_blocks, subject_blocks):
            return {
                "match_id": match_id,
                "query_occurrence_id": query,
                "subject_occurrence_id": subject,
                "match_status": "mapped",
                "candidate_resolution": "resolved",
                "membership_edge_eligible": 1,
                "position_edge_eligible": 1,
                "candidate_ids": f"{match_id}.candidate_001",
                "retained_candidate_ids": f"{match_id}.candidate_001",
                "matched_blocks": "1-10:1-10",
                "query_genomic_matched_blocks": query_blocks,
                "subject_genomic_matched_blocks": subject_blocks,
            }

        scored = [
            hard("m_left", "left", "fused", "chrA:1-10:+", "chrB:101-110:+"),
            hard("m_right", "right", "fused", "chrA:21-30:+", "chrB:121-130:+"),
            {
                "match_id": "m_competing",
                "query_occurrence_id": "fused",
                "subject_occurrence_id": "left",
                "match_status": "candidate_ambiguous",
                "candidate_resolution": "candidate",
                "membership_edge_eligible": 0,
                "position_edge_eligible": 0,
                "candidate_ids": "m_competing.candidate_001",
                "query_genomic_matched_blocks": "NA",
                "subject_genomic_matched_blocks": "NA",
                "candidate_assessments": json.dumps([{
                    "candidate_id": "m_competing.candidate_001",
                    "accepted": 1,
                    "query_genomic_matched_blocks": "chrB:105-108:+",
                    "subject_genomic_matched_blocks": "chrA:5-8:+",
                }]),
            },
        ]
        _membership_match_details(element_rows, occurrences, scored)

        fused_rows = sorted(
            (row for row in element_rows if row["occurrence_id"] == "fused"),
            key=lambda row: json.loads(row["actual_matched_blocks"])[0]["target_start"],
        )
        self.assertEqual(len(fused_rows), 2)
        self.assertNotEqual(fused_rows[0]["element_id"], fused_rows[1]["element_id"])
        self.assertEqual(
            [json.loads(row["actual_matched_blocks"])[0]["target_start"] for row in fused_rows],
            [101, 121],
        )
        self.assertEqual([row["parent_feature_ids"] for row in fused_rows], ["fused_exon", "fused_exon"])
        self.assertEqual(
            [row["matched_blocks"] for row in fused_rows],
            [
                "m_left:subject:1-10:1-10",
                "m_right:subject:1-10:1-10",
            ],
        )
        self.assertEqual(fused_rows[0]["correspondence_status"], "ambiguous")
        self.assertEqual(fused_rows[0]["membership_edge_eligible"], 0)
        self.assertEqual(fused_rows[1]["correspondence_status"], "resolved")
        self.assertEqual(fused_rows[1]["membership_edge_eligible"], 1)

    def test_unanchored_protein_projection_stays_candidate_evidence(self):
        occurrences = [
            {"occurrence_id": "a", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "source_feature_id": "a_exon"},
            {"occurrence_id": "b", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "source_feature_id": "b_exon"},
        ]
        element_rows = [
            {"element_id": "EG1", "homology_id": "HC1", "occurrence_id": row["occurrence_id"], "family_id": "fam", "species": row["species"], "gene_copy_id": row["gene_copy_id"], "confidence": "1", "membership_call": "core_member"}
            for row in occurrences
        ]
        scored = [{
            "match_id": "protein_candidate",
            "query_occurrence_id": "a",
            "subject_occurrence_id": "b",
            "match_status": "mapped",
            "candidate_resolution": "ambiguous",
            "correspondence_basis": "annotated_CDS_protein",
            "membership_edge_eligible": 1,
            "position_edge_eligible": 1,
            "protein_candidate_evidence_available": 1,
            "protein_hard_observation_eligible": 0,
            "candidate_ids": "protein:a>b",
            "query_genomic_matched_blocks": "chrA:11-20:+",
            "subject_genomic_matched_blocks": "chrB:31-40:+",
        }]
        _membership_match_details(element_rows, occurrences, scored)

        self.assertEqual(len(element_rows), 2)
        self.assertTrue(all(row["actual_matched_blocks"] == "[]" for row in element_rows))
        self.assertTrue(all(row["membership_edge_eligible"] == 0 for row in element_rows))
        self.assertTrue(all(row["correspondence_status"] == "ambiguous" for row in element_rows))

    def test_transcript_aliases_share_one_repeat_instance(self):
        occurrences = [
            {"occurrence_id": "alias_1", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "source_feature_id": "exon_a", "transcript_id": "tx1"},
            {"occurrence_id": "alias_2", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "source_feature_id": "exon_a", "transcript_id": "tx2"},
            {"occurrence_id": "duplicate", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "source_feature_id": "exon_b", "transcript_id": "tx1"},
            {"occurrence_id": "reference", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "source_feature_id": "exon_ref", "transcript_id": "txB"},
        ]
        element_rows = [
            {"element_id": "EG1", "homology_id": "HC1", "occurrence_id": row["occurrence_id"], "family_id": "fam", "species": row["species"], "gene_copy_id": row["gene_copy_id"], "confidence": "1", "membership_call": "core_member"}
            for row in occurrences
        ]
        scored = []
        for index, (occurrence_id, query_blocks) in enumerate((
            ("alias_1", "chrA:101-110:+"),
            ("alias_2", "chrA:101-110:+"),
            ("duplicate", "chrA:201-210:+"),
        ), start=1):
            scored.append({
                "match_id": f"m{index}",
                "query_occurrence_id": occurrence_id,
                "subject_occurrence_id": "reference",
                "match_status": "mapped",
                "candidate_resolution": "resolved",
                "membership_edge_eligible": 1,
                "position_edge_eligible": 1,
                "candidate_ids": f"c{index}",
                "retained_candidate_ids": f"c{index}",
                "query_genomic_matched_blocks": query_blocks,
                "subject_genomic_matched_blocks": "chrB:301-310:+",
            })
        _membership_match_details(element_rows, occurrences, scored)

        by_occurrence = {row["occurrence_id"]: row for row in element_rows}
        self.assertEqual(
            by_occurrence["alias_1"]["repeat_instance_id"],
            by_occurrence["alias_2"]["repeat_instance_id"],
        )
        self.assertNotEqual(
            by_occurrence["alias_1"]["repeat_instance_id"],
            by_occurrence["duplicate"]["repeat_instance_id"],
        )
        self.assertTrue(all(row["membership_edge_eligible"] == 1 for row in element_rows))

    def test_short_candidate_is_rerun_in_the_actual_interflank_interval(self):
        occurrences = {
            "short": {"occurrence_id": "short", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "txA", "contig": "chrA", "start": "101", "end": "112", "strand": "+"},
            "target": {"occurrence_id": "target", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "contig": "chrB", "start": "5", "end": "20", "strand": "+"},
            "left_q": {"occurrence_id": "left_q", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "txA", "contig": "chrA", "start": "1", "end": "4", "strand": "+"},
            "left_t": {"occurrence_id": "left_t", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "contig": "chrB", "start": "1", "end": "4", "strand": "+"},
            "right_q": {"occurrence_id": "right_q", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "txA", "contig": "chrA", "start": "201", "end": "204", "strand": "+"},
            "right_t": {"occurrence_id": "right_t", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "contig": "chrB", "start": "21", "end": "24", "strand": "+"},
        }
        left_record = {
            "candidate_id": "left_anchor", "accepted": 1,
            "aligned_blocks": (CoordinateBlock(Interval0(0, 4), Interval0(0, 4)),),
            "query_genomic_blocks0": (Interval0(0, 4),),
            "target_genomic_blocks0": (Interval0(0, 4),),
        }
        right_record = {
            "candidate_id": "right_anchor", "accepted": 1,
            "aligned_blocks": (CoordinateBlock(Interval0(0, 4), Interval0(0, 4)),),
            "query_genomic_blocks0": (Interval0(200, 204),),
            "target_genomic_blocks0": (Interval0(20, 24),),
        }
        focus_record = {
            "candidate_id": "focus_feature_candidate", "accepted": 1,
            "aligned_blocks": (CoordinateBlock(Interval0(0, 12), Interval0(0, 12)),),
        }
        rows = [
            {"match_id": "left_match", "query_occurrence_id": "left_q", "subject_occurrence_id": "left_t", "retained_candidate_ids": "left_anchor", "_position_edge_eligible": True, "_candidate_records": [left_record]},
            {"match_id": "right_match", "query_occurrence_id": "right_q", "subject_occurrence_id": "right_t", "retained_candidate_ids": "right_anchor", "_position_edge_eligible": True, "_candidate_records": [right_record]},
            {
                "match_id": "focus_match", "query_occurrence_id": "short", "subject_occurrence_id": "target",
                "alignment_input_transposed": 0, "short_context_route": "feature_bounded_candidate",
                "retained_candidate_ids": "NA", "threshold": "0.7",
                "_candidate_records": [focus_record],
            },
        ]
        query = "ACGTTGCAAGTC"
        gene_loci = {
            ("B", "gB"): {
                "sequence": "AAAA" + query + "TTTT" + "CCCC",
                "contig": "chrB", "strand": "+", "interval": Interval0(0, 24),
            },
        }
        transcript_paths = [
            {**occurrences[occurrence_id], "path_rank": rank}
            for occurrence_id, rank in (
                ("left_q", 1), ("short", 2), ("right_q", 3),
                ("left_t", 1), ("target", 2), ("right_t", 3),
            )
        ]
        changed = _rerun_anchor_bounded_short_candidates(
            rows,
            occurrences,
            {"short": query},
            gene_loci,
            transcript_paths=transcript_paths,
        )
        focus = rows[-1]
        self.assertTrue(changed)
        self.assertEqual(focus["short_context_route"], "anchor_bounded_local")
        self.assertEqual(focus["left_anchor_id"], "left_anchor")
        self.assertEqual(focus["right_anchor_id"], "right_anchor")
        self.assertEqual(json.loads(focus["search_interval"])["start"], 5)
        self.assertEqual(json.loads(focus["search_interval"])["end"], 20)
        self.assertEqual(focus["subject_genomic_matched_blocks"], "chrB:5-16:+")
        self.assertEqual(focus["membership_edge_eligible"], 0)
