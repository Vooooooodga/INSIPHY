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

from support_evidence_evidence_semantics import EvidenceSemanticsTestsSupport

class EvidenceSemanticsTests(EvidenceSemanticsTestsSupport, unittest.TestCase):
    def test_no_hit_at_search_limit_remains_unknown_not_sequence_absence(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_minimal_case(Path(tmp))
            write_tsv(
                input_dir / "gene_loci.tsv",
                [
                    "species", "gene_copy_id", "contig", "strand",
                    "annotation_start", "annotation_end", "linked_start", "linked_end",
                    "search_start", "search_end", "contig_length", "genome_fasta",
                    "max_extension", "range_status",
                ],
                [
                    ["B", "B_gene", "chrB", "+", "20", "23", "20", "23",
                     "1", "40", "40", "unused.fa", "10", "complete"],
                ],
            )
            no_hit = AlignmentStats(
                identity=0.0, coverage=0.0, score=0.0,
                query_end=4, target_end=40, backend="minimap2",
                enumeration_complete=True,
            )
            with patch("intraphy.evidence.nucleotide_search.local_alignment_stats", return_value=no_hit):
                evidence = generate_sequence_evidence(
                    input_dir, result_dir, aligner="minimap2"
                )
            row = evidence[0]
            self.assertEqual(row["evidence_status"], "ambiguous")
            self.assertEqual(row["homologous_dna_presence"], "unknown")
            self.assertNotEqual(row["homologous_dna_presence"], "absent")
            self.assertEqual(row["candidate_resolution_status"], "candidate_search_incomplete")
            self.assertEqual(row["candidate_search_complete"], "false")
            self.assertEqual(row["original_annotation_start"], 20)
            self.assertEqual(row["linked_feature_end"], 23)
            self.assertEqual(row["search_bound_start"], 1)
            self.assertEqual(
                row["hit_search_limit_status"],
                "no_hit_with_left_search_limit_reached",
            )
            completed = complete_annotation(input_dir, result_dir)
            self.assertEqual(completed[0]["completion_call"], "ambiguous_evidence")
            self.assertEqual(completed[0]["homologous_dna_presence"], "unknown")

    def test_ambiguous_protein_containers_retain_individual_query_intervals(self):
        context = {"protein_id": "A|g|tx", "query_start": 3, "query_end": 5}
        first = {
            "protein_id": "A|g|tx", "query_start": 1, "query_end": 6,
            "target_start": 10, "target_end": 27, "strand": "+", "identity": 1.0,
            "parent_id": "p1", "projection_status": "cds_target",
        }
        second = {**first, "query_start": 2, "target_start": 40, "target_end": 54, "parent_id": "p2"}
        candidates = []
        projected, reason = _protein_projection_from_rows(
            context, [first, second], "B|g|chrB:101-200:-", 0.7, 0.6, candidates
        )
        self.assertIsNone(projected)
        self.assertEqual(reason, "protein_projection_ambiguous_multiple_possible_mappings")
        self.assertEqual([(row["start"], row["end"]) for row in candidates], [(174, 191), (147, 161)])
        self.assertEqual([row["protein_cds_query_start"] for row in candidates], [1, 2])
        for row in candidates:
            self.assertEqual(row["interval_scope"], "projected_cds_container")
            self.assertEqual((row["protein_overlap_query_start"], row["protein_overlap_query_end"]), (3, 5))
            self.assertEqual(row["strand"], "-")

    def test_reference_protein_context_uses_real_transcript_not_joined_id(self):
        representative = {
            "occurrence_id": "A_e1",
            "species": "A",
            "gene_copy_id": "A_gene",
            "transcript_id": "tx1;tx2",
            "start": "1",
            "end": "12",
            "cds_length": "12",
            "cds_phase": "0",
        }
        transcript_paths = [
            {
                "species": "A",
                "gene_copy_id": "A_gene",
                "transcript_id": "tx1",
                "path_rank": "1",
                "occurrence_id": "A_e1",
                "coding_status": "coding",
                "cds_length": "9",
                "cds_phase": "0",
            },
            {
                "species": "A",
                "gene_copy_id": "A_gene",
                "transcript_id": "tx2",
                "path_rank": "1",
                "occurrence_id": "A_e1",
                "coding_status": "coding",
                "cds_length": "12",
                "cds_phase": "0",
            },
        ]
        proteins = {"A|A_gene|tx2": "MMMM"}
        context = _reference_protein_context(representative, transcript_paths, proteins, {"A_e1": representative})
        self.assertEqual(context["protein_id"], "A|A_gene|tx2")
        self.assertEqual(context["transcript_id"], "tx2")
        self.assertEqual(context["query_start"], 1)
        self.assertEqual(context["query_end"], 4)
        self.assertEqual(context["cds_length"], 12)

    def test_protein_context_uses_transcript_cds_intervals_and_initial_phase(self):
        first = {"occurrence_id": "first", "coding_status": "coding", "cds_length": "6", "cds_phase": "0", "start": "1", "end": "6"}
        exon = {"occurrence_id": "shared", "species": "A", "gene_copy_id": "g", "transcript_id": "tx;other", "coding_status": "coding", "cds_length": "6", "start": "10", "end": "15"}
        common = {"species": "A", "gene_copy_id": "g", "transcript_id": "tx", "coding_status": "coding"}
        paths = [
            {**common, "occurrence_id": "first", "path_rank": "1", "cds_intervals": "2-6", "cds_phase": "2"},
            {**common, "occurrence_id": "shared", "path_rank": "2", "cds_intervals": "10-15", "cds_phase": "0"},
        ]
        proteins = {"A|g|tx": "MMM"}
        occurrences = {"first": first, "shared": exon}
        context = _reference_protein_context(exon, paths, proteins, occurrences)
        self.assertEqual((context["query_start"], context["query_end"]), (2, 3))
        self.assertEqual(context["cds_phase"], "0")
        for overrides in ({"cds_length": "0"}, {"cds_length": "6", "coding_status": "noncoding"}):
            with self.subTest(**overrides):
                modified_paths = [paths[0], {**paths[1], **overrides}]
                self.assertIsNone(_reference_protein_context(exon, modified_paths, proteins, occurrences))

    def test_deletion_support_uses_flanks_without_penalizing_deleted_interval(self):
        source_occurrences = [
            {"occurrence_id": "left_s", "contig": "chrA", "start": "1", "end": "4", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "mid_s", "contig": "chrA", "start": "5", "end": "8", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "right_s", "contig": "chrA", "start": "9", "end": "12", "strand": "+", "presence_status": "present"},
        ]
        target_occurrences = [
            {"occurrence_id": "left_t", "contig": "chrB", "start": "1", "end": "4", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "right_t", "contig": "chrB", "start": "5", "end": "8", "strand": "+", "presence_status": "present"},
        ]
        element_by_occurrence = {
            "left_s": "Eleft",
            "mid_s": "Emid",
            "right_s": "Eright",
            "left_t": "Eleft",
            "right_t": "Eright",
        }
        spanning = AlignmentStats(
            identity=8 / 12,
            coverage=8 / 12,
            score=8.0,
            query_start=1,
            query_end=12,
            target_start=1,
            target_end=8,
            cigar="4M4I4M",
            backend="minimap2",
            query_coverage=8 / 12,
            target_coverage=1.0,
            aligned_pairs=8,
            strand="+",
            aligned_blocks=[(1, 4, 1, 4), (9, 12, 5, 8)],
            matches=8,
            gap_bases=4,
        )
        with patch("intraphy.evidence.deletion_support.local_alignment_stats", return_value=spanning) as align:
            supported, reason, provenance = _ordered_anchor_deletion_support(
                source_occurrences,
                target_occurrences,
                element_by_occurrence,
                source_occurrences[1],
                "Eleft",
                "Eright",
                "A|gene|chrA:1-12:+",
                "AAAACCCCGGGG",
                "B|gene|chrB:1-8:+",
                "AAAAGGGG",
                0.9,
                0.9,
            )
        align.assert_called_once_with("AAAACCCCGGGG", "AAAAGGGG", backend="minimap2", threads=1)
        self.assertTrue(supported)
        self.assertEqual(reason, "source_expected_exon_deleted_in_target_spanning_alignment")
        self.assertEqual(provenance["left_flank_identity"], "1")
        self.assertEqual(provenance["right_flank_identity"], "1")
        self.assertEqual(provenance["left_flank_paired_bases"], 4)
        self.assertEqual(provenance["right_flank_paired_bases"], 4)
        self.assertEqual(provenance["left_flank_matches"], 4)
        self.assertEqual(provenance["right_flank_mismatches"], 0)

    def test_deletion_support_keeps_iupac_ambiguous_intervals_unknown(self):
        source_occurrences = [
            {"occurrence_id": "left_s", "contig": "chrA", "start": "1", "end": "4", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "mid_s", "contig": "chrA", "start": "5", "end": "8", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "right_s", "contig": "chrA", "start": "9", "end": "12", "strand": "+", "presence_status": "present"},
        ]
        target_occurrences = [
            {"occurrence_id": "left_t", "contig": "chrB", "start": "1", "end": "4", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "right_t", "contig": "chrB", "start": "6", "end": "9", "strand": "+", "presence_status": "present"},
        ]
        element_by_occurrence = {
            "left_s": "Eleft", "mid_s": "Emid", "right_s": "Eright",
            "left_t": "Eleft", "right_t": "Eright",
        }
        ambiguous_codes = "RYSWKMBDHVN"
        for base in ambiguous_codes:
            with self.subTest(location="source_deletion", base=base):
                source_sequence = f"AAAA{base}CCCGGGG"
                with patch("intraphy.evidence.deletion_support.local_alignment_stats") as align:
                    supported, reason, _provenance = _ordered_anchor_deletion_support(
                        source_occurrences, target_occurrences, element_by_occurrence,
                        source_occurrences[1], "Eleft", "Eright",
                        "A|gene|chrA:1-12:+", source_sequence,
                        "B|gene|chrB:1-9:+", "AAAAGGGG", 0.9, 0.9,
                    )
                self.assertFalse(supported)
                self.assertEqual(reason, "deletion_interval_contains_non_acgt_bases")
                align.assert_not_called()

            with self.subTest(location="target_between_flanks", base=base):
                target_sequence = f"AAAA{base}GGGG"
                with patch("intraphy.evidence.deletion_support.local_alignment_stats") as align:
                    supported, reason, _provenance = _ordered_anchor_deletion_support(
                        source_occurrences, target_occurrences, element_by_occurrence,
                        source_occurrences[1], "Eleft", "Eright",
                        "A|gene|chrA:1-12:+", "AAAACCCCGGGG",
                        "B|gene|chrB:1-9:+", target_sequence, 0.9, 0.9,
                    )
                self.assertFalse(supported)
                self.assertEqual(
                    reason,
                    "target_sequence_between_flanking_homologs_contains_non_acgt_bases",
                )
                align.assert_not_called()

    def test_acceptable_alternative_without_deletion_blocks_keeps_absence_unknown(self):
        source_occurrences = [
            {"occurrence_id": "left_s", "contig": "chrA", "start": "1", "end": "4", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "mid_s", "contig": "chrA", "start": "5", "end": "8", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "right_s", "contig": "chrA", "start": "9", "end": "12", "strand": "+", "presence_status": "present"},
        ]
        target_occurrences = [
            {"occurrence_id": "left_t", "contig": "chrB", "start": "1", "end": "4", "strand": "+", "presence_status": "present"},
            {"occurrence_id": "right_t", "contig": "chrB", "start": "5", "end": "8", "strand": "+", "presence_status": "present"},
        ]
        element_by_occurrence = {
            "left_s": "Eleft", "mid_s": "Emid", "right_s": "Eright",
            "left_t": "Eleft", "right_t": "Eright",
        }
        spanning = AlignmentStats(
            identity=1.0, coverage=8 / 12, score=8.0,
            query_start=1, query_end=12, target_start=1, target_end=8,
            cigar="4M4I4M", backend="minimap2", query_coverage=8 / 12,
            target_coverage=1.0, aligned_pairs=8, strand="+",
            aligned_blocks=[(1, 4, 1, 4), (9, 12, 5, 8)], matches=8,
            gap_bases=4, hit_count=2, ambiguous_hit_count=1,
            alternative_hits=[{"identity": "0.95", "coverage": "0.95"}],
        )
        with patch("intraphy.evidence.deletion_support.local_alignment_stats", return_value=spanning):
            supported, reason, _provenance = _ordered_anchor_deletion_support(
                source_occurrences, target_occurrences, element_by_occurrence,
                source_occurrences[1], "Eleft", "Eright",
                "A|gene|chrA:1-12:+", "AAAACCCCGGGG",
                "B|gene|chrB:1-8:+", "AAAAGGGG", 0.9, 0.9,
            )
        self.assertFalse(supported)
        self.assertEqual(reason, "alternative_alignment_concordance_unverifiable")

    def test_flank_identity_counts_only_bases_on_spanning_path(self):
        blocks = [{"query_start": 1, "query_end": 4, "target_start": 1, "target_end": 4}]
        support = _flank_support_from_spanning_blocks("AAAA", "TTTTAAAA", blocks, (1, 4), (1, 4))
        self.assertEqual(support["identity"], 0.0)
        self.assertEqual(support["coverage"], 1.0)
        self.assertEqual(support["mismatches"], 4)
        support = _flank_support_from_spanning_blocks("AAAA", "AANA", blocks, (1, 4), (1, 4))
        self.assertEqual(support["identity"], 1.0)
        self.assertEqual(support["coverage"], 0.75)
        self.assertEqual(support["paired_bases"], 3)

    def test_unanchored_repeated_mapping_remains_descriptive_and_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_minimal_case(Path(tmp), target_role="noncoding_exon", repeated_hit=True)
            stat = AlignmentStats(
                identity=1.0,
                coverage=1.0,
                score=4.0,
                query_start=1,
                query_end=4,
                target_start=20,
                target_end=23,
                cigar="4M",
                backend="minimap2",
                query_coverage=1.0,
                target_coverage=0.1,
                aligned_pairs=4,
                strand="+",
                aligned_blocks=[(1, 4, 20, 23)],
                matches=4,
                alignment_mode="local",
                mapping_quality=0,
                hit_count=2,
                ambiguous_hit_count=1,
                alternative_hits=[{"rank": 2, "target_start": 30, "target_end": 33, "strand": "+", "identity": "1", "coverage": "1", "mapping_quality": 0, "is_secondary": 1}],
            )
            def align(query, target, **kwargs):
                return stat if len(target) == 40 else AlignmentStats(identity=0.0, coverage=0.0, score=0.0)

            with patch("intraphy.evidence.nucleotide_search.local_alignment_stats", side_effect=align):
                rows = generate_sequence_evidence(input_dir, result_dir, aligner="minimap2")
            row = next(row for row in rows if row["homology_id"] == "H1")
            self.assertEqual(row["evidence_status"], "ambiguous")
            self.assertEqual(row["inferred_role"], "unknown")
            self.assertEqual(row["primary_mapping_status"], "whole_locus_descriptive_candidate")
            self.assertEqual(row["correspondence_status"], "unknown")
            self.assertEqual(row["homologous_dna_presence"], "unknown")
            self.assertEqual((row["start"], row["end"]), ("NA", "NA"))
            self.assertEqual([candidate["start"] for candidate in json.loads(row["interval_candidates"])], [20, 30])
            completed = complete_annotation(input_dir, result_dir)
            row = next(row for row in completed if row["homology_id"] == "H1")
            self.assertEqual(row["completion_call"], "ambiguous_evidence")
            self.assertEqual(row["correspondence_status"], "unknown")
            self.assertEqual(len(json.loads(row["interval_candidates"])), 2)
