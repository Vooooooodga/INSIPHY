import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from insiphy.alignment import AlignmentStats, local_alignment_stats
from insiphy.annotation import (
    _cached_protein_projection,
    _flank_support_from_spanning_blocks,
    _ordered_anchor_deletion_support,
    _protein_projection_from_rows,
    _reference_protein_context,
    complete_annotation,
    completion_call,
    generate_sequence_evidence,
)
from insiphy.io import read_tsv, write_tsv as write_rows
from insiphy.candidate_chain import ChainCandidate, ChainPathMembership, ordered_candidate_chain
from insiphy.coordinates import CoordinateBlock, Interval0
from insiphy.correspondence import _membership_match_details
from insiphy.preprocess import _rerun_anchor_bounded_short_candidates


def write_tsv(path, header, rows):
    path.write_text("\n".join(["\t".join(header), *["\t".join(row) for row in rows]]) + "\n")


class EvidenceSemanticsTests(unittest.TestCase):
    def write_minimal_case(self, root, target_role="intron", repeated_hit=False):
        input_dir = root / "input"
        result_dir = root / "result"
        input_dir.mkdir()
        result_dir.mkdir()
        occurrence_header = [
            "occurrence_id",
            "family_id",
            "species",
            "gene_copy_id",
            "transcript_id",
            "segment_id",
            "contig",
            "start",
            "end",
            "strand",
            "role",
            "presence_status",
            "coding_status",
            "cds_length",
            "cds_phase",
        ]
        write_tsv(
            input_dir / "segment_occurrences.tsv",
            occurrence_header,
            [
                ["A_e1", "fam", "A", "A_gene", "tx1", "s1", "chrA", "1", "4", "+", "noncoding_exon", "present", "noncoding", "0", "."],
                ["B_i1", "fam", "B", "B_gene", "tx1", "s0", "chrB", "20", "23", "+", target_role, "present", "noncoding", "0", "."],
            ],
        )
        write_tsv(
            input_dir / "transcript_paths.tsv",
            [
                "path_id",
                "family_id",
                "species",
                "gene_copy_id",
                "transcript_id",
                "path_rank",
                "occurrence_id",
                "role",
                "contig",
                "start",
                "end",
                "strand",
                "phase",
                "path_status",
                "coding_status",
                "cds_phase",
            ],
            [
                ["A_p1", "fam", "A", "A_gene", "tx1", "1", "A_e1", "noncoding_exon", "chrA", "1", "4", "+", ".", "canonical", "noncoding", "."],
                ["B_p1", "fam", "B", "B_gene", "tx1", "1", "B_i1", target_role, "chrB", "20", "23", "+", ".", "canonical", "noncoding", "."],
            ],
        )
        write_tsv(
            result_dir / "element_correspondence.tsv",
            [
                "element_id",
                "family_id",
                "homology_id",
                "occurrence_id",
                "species",
                "gene_copy_id",
                "element_class",
                "display_role",
                "source_label",
                "support_type",
                "confidence",
                "membership_score",
                "membership_call",
            ],
            [
                ["EG1", "fam", "H1", "A_e1", "A", "A_gene", "exon_like", "noncoding_exon", "fam", "annotated", "medium", "0.4", "ambiguous_member"],
                ["EG0", "fam", "H0", "B_i1", "B", "B_gene", "candidate_source" if target_role == "intron" else "exon_like", target_role, "fam", "annotated", "high", "1", "core_member"],
            ],
        )
        (input_dir / "segment_sequences.fasta").write_text(">A_e1\nACGT\n>B_i1\nACGT\n")
        target = "T" * 19 + "ACGT" + ("T" * 6 + "ACGT" + "T" * 7 if repeated_hit else "T" * 17)
        (input_dir / "gene_loci.fasta").write_text(f">A|A_gene|chrA:1-4:+\nACGT\n>B|B_gene|chrB:1-40:+\n{target}\n")
        (input_dir / "protein_sequences.fasta").write_text("")
        return input_dir, result_dir

    def write_anchored_case(self, root):
        input_dir = root / "input"
        result_dir = root / "result"
        input_dir.mkdir()
        result_dir.mkdir()
        occurrence_header = [
            "occurrence_id", "family_id", "species", "gene_copy_id", "transcript_id",
            "segment_id", "contig", "start", "end", "strand", "role",
            "presence_status", "coding_status", "cds_length", "cds_phase",
        ]
        bounded_target = "T" * 8 + "ACGT" + "T" * 8
        right_start = 5 + len(bounded_target)
        write_tsv(
            input_dir / "segment_occurrences.tsv",
            occurrence_header,
            [
                ["A_left", "fam", "A", "A_gene", "tx1", "sl", "chrA", "1", "4", "+", "noncoding_exon", "present", "noncoding", "0", "."],
                ["A_mid", "fam", "A", "A_gene", "tx1", "sm", "chrA", "5", "8", "+", "noncoding_exon", "present", "noncoding", "0", "."],
                ["A_right", "fam", "A", "A_gene", "tx1", "sr", "chrA", "9", "12", "+", "noncoding_exon", "present", "noncoding", "0", "."],
                ["B_left", "fam", "B", "B_gene", "tx1", "tl", "chrB", "1", "4", "+", "noncoding_exon", "present", "noncoding", "0", "."],
                ["B_gap", "fam", "B", "B_gene", "tx1", "tg", "chrB", "5", str(right_start - 1), "+", "intron", "present", "noncoding", "0", "."],
                ["B_right", "fam", "B", "B_gene", "tx1", "tr", "chrB", str(right_start), str(right_start + 3), "+", "noncoding_exon", "present", "noncoding", "0", "."],
            ],
        )
        path_header = [
            "path_id", "family_id", "species", "gene_copy_id", "transcript_id",
            "path_rank", "occurrence_id", "role", "contig", "start", "end", "strand",
            "phase", "path_status", "coding_status", "cds_phase",
        ]
        write_tsv(
            input_dir / "transcript_paths.tsv",
            path_header,
            [
                ["A_p1", "fam", "A", "A_gene", "tx1", "1", "A_left", "noncoding_exon", "chrA", "1", "4", "+", ".", "canonical", "noncoding", "."],
                ["A_p2", "fam", "A", "A_gene", "tx1", "2", "A_mid", "noncoding_exon", "chrA", "5", "8", "+", ".", "canonical", "noncoding", "."],
                ["A_p3", "fam", "A", "A_gene", "tx1", "3", "A_right", "noncoding_exon", "chrA", "9", "12", "+", ".", "canonical", "noncoding", "."],
                ["B_p1", "fam", "B", "B_gene", "tx1", "1", "B_left", "noncoding_exon", "chrB", "1", "4", "+", ".", "canonical", "noncoding", "."],
                ["B_p2", "fam", "B", "B_gene", "tx1", "2", "B_gap", "intron", "chrB", "5", str(right_start - 1), "+", ".", "canonical", "noncoding", "."],
                ["B_p3", "fam", "B", "B_gene", "tx1", "3", "B_right", "noncoding_exon", "chrB", str(right_start), str(right_start + 3), "+", ".", "canonical", "noncoding", "."],
            ],
        )
        element_header = [
            "element_id", "family_id", "homology_id", "occurrence_id", "species",
            "gene_copy_id", "element_class", "display_role", "source_label",
            "support_type", "confidence", "membership_score", "membership_call",
        ]
        write_tsv(
            result_dir / "element_correspondence.tsv",
            element_header,
            [
                ["Eleft", "fam", "Hleft", "A_left", "A", "A_gene", "exon_like", "noncoding_exon", "fam", "annotated", "high", "1", "core_member"],
                ["Eleft", "fam", "Hleft", "B_left", "B", "B_gene", "exon_like", "noncoding_exon", "fam", "annotated", "high", "1", "core_member"],
                ["Emid", "fam", "Hmid", "A_mid", "A", "A_gene", "exon_like", "noncoding_exon", "fam", "annotated", "high", "1", "core_member"],
                ["Eright", "fam", "Hright", "A_right", "A", "A_gene", "exon_like", "noncoding_exon", "fam", "annotated", "high", "1", "core_member"],
                ["Eright", "fam", "Hright", "B_right", "B", "B_gene", "exon_like", "noncoding_exon", "fam", "annotated", "high", "1", "core_member"],
            ],
        )
        (input_dir / "segment_sequences.fasta").write_text(
            ">A_left\nGGGG\n>A_mid\nACGT\n>A_right\nCCCC\n"
            ">B_left\nAAAA\n>B_gap\n" + bounded_target + "\n>B_right\nCCCC\n"
        )
        (input_dir / "gene_loci.fasta").write_text(
            ">A|A_gene|chrA:1-12:+\nGGGGACGTCCCC\n"
            f">B|B_gene|chrB:1-{right_start + 3}:+\nAAAA{bounded_target}CCCC\n"
        )
        (input_dir / "protein_sequences.fasta").write_text("")
        return input_dir, result_dir, bounded_target

    def test_whole_locus_hit_without_double_flank_remains_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_minimal_case(Path(tmp))
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
            )
            with patch("insiphy.annotation.local_alignment_stats", return_value=stat) as align:
                rows = generate_sequence_evidence(input_dir, result_dir, aligner="miniprot")
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["evidence_status"], "ambiguous")
            self.assertEqual(rows[0]["inferred_role"], "unknown")
            self.assertEqual(rows[0]["alignment_backend"], "minimap2")
            self.assertEqual(rows[0]["evidence_aligner"], "miniprot")
            self.assertEqual(rows[0]["interval_scope"], "whole_locus_descriptive_fallback")
            self.assertEqual(rows[0]["homologous_dna_presence"], "unknown")
            self.assertEqual(rows[0]["anchor_interval_status"], "missing_flanking_homolog_anchor")
            self.assertEqual(rows[0]["primary_mapping_status"], "whole_locus_descriptive_candidate")
            self.assertEqual(rows[0]["predicted_exonic_role"], "unknown")
            self.assertEqual(rows[0]["supplied_annotation_role"], "intron")
            self.assertEqual(rows[0]["source_parent_occurrence_id"], "A_e1")
            self.assertEqual(rows[0]["source_parent_transcript_ids"], "tx1")
            self.assertEqual(rows[0]["target_parent_occurrence_ids"], "B_i1")
            block = json.loads(rows[0]["dna_aligned_blocks"])[0]
            self.assertEqual((block["source_start"], block["source_end"]), (1, 4))
            self.assertEqual((block["target_start"], block["target_end"]), (20, 23))
            self.assertEqual(block["source_occurrence_id"], "A_e1")
            self.assertEqual(json.loads(rows[0]["interval_candidates"])[0]["interval_scope"], "aligned_sequence")
            completed = complete_annotation(input_dir, result_dir)
            self.assertEqual(completed[0]["interval_scope"], "whole_locus_descriptive_fallback")
            self.assertEqual(completed[0]["homologous_dna_presence"], "unknown")
            self.assertEqual(completed[0]["supplied_annotation_role"], "intron")
            self.assertEqual(completed[0]["completion_call"], "ambiguous_evidence")
            align.assert_called_once_with("ACGT", "T" * 19 + "ACGT" + "T" * 17, backend="minimap2", threads=1)

    def test_short_query_uses_longer_ordered_anchor_interval_for_presence(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, bounded_target = self.write_anchored_case(Path(tmp))
            rows = generate_sequence_evidence(
                input_dir,
                result_dir,
                aligner="internal",
                short_context_max_length=4,
            )
            row = next(row for row in rows if row["homology_id"] == "Hmid")
            self.assertEqual(len(bounded_target), 20)
            self.assertEqual(row["alignment_backend"], "internal")
            self.assertEqual(row["alignment_evidence_scope"], "anchor_bounded_short_local")
            self.assertEqual(row["anchor_interval_status"], "ordered_double_flank_bounded_interval")
            self.assertEqual((row["left_anchor_id"], row["right_anchor_id"]), ("B_left", "B_right"))
            self.assertEqual(json.loads(row["search_interval"])["start"], 5)
            self.assertEqual(json.loads(row["search_interval"])["end"], 24)
            self.assertEqual(row["alignment_query_length"], 4)
            self.assertEqual(row["alignment_target_length"], 20)
            self.assertEqual(row["alignment_known_aligned_pairs"], 4)
            self.assertEqual(row["alignment_unknown_aligned_pairs"], 0)
            self.assertEqual(row["homologous_dna_presence"], "present")
            self.assertEqual((row["start"], row["end"]), (13, 16))
            self.assertEqual(row["candidate_resolution_status"], "resolved")

    def test_protein_projection_remains_predicted_candidate(self):
        row = {
            "evidence_status": "supports_hidden_segment",
            "inferred_event": "protein_cds_projection",
            "inferred_role": "predicted_CDS",
            "predicted_role": "CDS",
            "frame_status": "unknown",
        }
        self.assertEqual(completion_call(row, 0.95, 0.55), "ambiguous_evidence")
        row["anchor_interval_status"] = "ordered_double_flank_bounded_interval"
        self.assertEqual(completion_call(row, 0.95, 0.55), "predicted_exon_candidate")

    def test_truncated_candidate_enumeration_is_retained_as_search_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_minimal_case(Path(tmp))
            stat = AlignmentStats(
                identity=1.0,
                coverage=1.0,
                score=4.0,
                query_start=1,
                query_end=4,
                target_start=20,
                target_end=23,
                cigar="4M",
                backend="internal",
                query_coverage=1.0,
                aligned_pairs=4,
                aligned_blocks=[(1, 4, 20, 23)],
                enumeration_complete=False,
                incomplete_reason="max_optimal_alignments_reached",
            )
            with patch("insiphy.annotation.local_alignment_stats", return_value=stat):
                rows = generate_sequence_evidence(
                    input_dir, result_dir, aligner="internal"
                )
            row = rows[0]
            self.assertEqual(row["homologous_dna_presence"], "unknown")
            self.assertEqual(row["alignment_evidence_scope"], "whole_locus_descriptive_fallback")
            self.assertEqual(row["candidate_resolution_status"], "candidate_search_incomplete")
            self.assertEqual(row["candidate_search_complete"], "false")
            self.assertEqual(
                row["candidate_search_incomplete_reason"],
                "max_optimal_alignments_reached",
            )
            self.assertEqual(row["correspondence_status"], "unknown")

    def test_completion_table_preserves_predicted_role(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            write_tsv(
                output_dir / "sequence_synteny_evidence.tsv",
                [
                    "evidence_id",
                    "family_id",
                    "species",
                    "gene_copy_id",
                    "homology_id",
                    "annotation_status",
                    "evidence_status",
                    "inferred_role",
                    "predicted_role",
                    "anchor_interval_status",
                    "contig",
                    "start",
                    "end",
                    "strand",
                    "sequence_score",
                    "sequence_coverage",
                    "inferred_event",
                    "frame_status",
                ],
                [
                    [
                        "ev1",
                        "fam",
                        "B",
                        "B_gene",
                        "H1",
                        "protein_projection_supports_missing_cds",
                        "supports_hidden_segment",
                        "predicted_CDS",
                        "CDS",
                        "ordered_double_flank_bounded_interval",
                        "chrB",
                        "10",
                        "20",
                        "+",
                        "0.9",
                        "0.9",
                        "protein_cds_projection",
                        "unknown",
                    ]
                ],
            )
            rows = complete_annotation(input_dir, output_dir)
            self.assertEqual(rows[0]["completion_call"], "predicted_exon_candidate")
            self.assertEqual(rows[0]["inferred_role"], "predicted_CDS")
            self.assertEqual(rows[0]["predicted_role"], "CDS")
            self.assertEqual(rows[0]["interval_scope"], "unspecified")

    def test_protein_projection_cache_batches_once_per_target_locus(self):
        context1 = {"protein_id": "A|A_gene|tx1", "protein": "M" * 6, "query_start": 1, "query_end": 3}
        context2 = {"protein_id": "A|A_gene|tx2", "protein": "M" * 6, "query_start": 4, "query_end": 6}
        rows = [
            {
                "protein_id": "A|A_gene|tx1",
                "query_start": 1,
                "query_end": 3,
                "target_start": 10,
                "target_end": 18,
                "strand": "+",
                "identity": 0.95,
                "phase": "0",
                "cigar": "9M",
                "parent_id": "p1",
                "projection_status": "cds_target",
            },
            {
                "protein_id": "A|A_gene|tx2",
                "query_start": 4,
                "query_end": 6,
                "target_start": 30,
                "target_end": 38,
                "strand": "+",
                "identity": 0.95,
                "phase": "0",
                "cigar": "9M",
                "parent_id": "p2",
                "projection_status": "cds_target",
            },
        ]
        cache = {}
        protein_contexts = {"fam": {context1["protein_id"]: context1, context2["protein_id"]: context2}}
        with patch("insiphy.annotation.protein_locus_exons", return_value=rows) as project:
            hit1, reason1 = _cached_protein_projection(
                context1, "fam", ("fam", "B", "B_gene"), "B|B_gene|chrB:1-100:+", "A" * 100,
                protein_contexts, cache, 0.7, 0.6, 1
            )
            hit2, reason2 = _cached_protein_projection(
                context2, "fam", ("fam", "B", "B_gene"), "B|B_gene|chrB:1-100:+", "A" * 100,
                protein_contexts, cache, 0.7, 0.6, 1
            )
        self.assertEqual(project.call_count, 1)
        self.assertEqual(reason1, "protein_projection_supports_cds")
        self.assertEqual(reason2, "protein_projection_supports_cds")
        self.assertEqual(hit1["start"], 10)
        self.assertEqual(hit2["start"], 30)

    def test_unanchored_partial_protein_container_is_descriptive_in_tables_and_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_minimal_case(Path(tmp))
            occurrences = read_tsv(input_dir / "segment_occurrences.tsv")
            exon = occurrences[0]
            exon.update(start="10", end="15", role="CDS", coding_status="coding", cds_length="6", cds_phase="0")
            first_exon = {**exon, "occurrence_id": "A_e0", "segment_id": "s_first", "start": "1", "end": "3", "cds_length": "3"}
            write_rows(input_dir / "segment_occurrences.tsv", [first_exon, *occurrences], list(exon))
            paths = read_tsv(input_dir / "transcript_paths.tsv")
            path = paths[0]
            path.update(path_rank="2", start="10", end="15", role="CDS", coding_status="coding", cds_length="6", cds_phase="0", phase="0")
            first_path = {**path, "path_id": "A_p0", "occurrence_id": "A_e0", "path_rank": "1", "start": "1", "end": "3", "cds_length": "3"}
            write_rows(input_dir / "transcript_paths.tsv", [first_path, *paths], list(path))
            elements = read_tsv(result_dir / "element_correspondence.tsv")
            elements[0]["display_role"] = "CDS"
            write_rows(result_dir / "element_correspondence.tsv", elements, list(elements[0]))
            (input_dir / "segment_sequences.fasta").write_text(">A_e0\nATG\n>A_e1\nAAAGCC\n>B_i1\nACGT\n")
            source_locus = "ATG" + "T" * 6 + "AAAGCC"
            target_locus = "T" * 9 + "ATGAAAGCC" + "T" + "ACGT" + "T" * 17
            (input_dir / "gene_loci.fasta").write_text(f">A|A_gene|chrA:1-15:+\n{source_locus}\n>B|B_gene|chrB:1-40:+\n{target_locus}\n")
            (input_dir / "protein_sequences.fasta").write_text(">A|A_gene|tx1\nMKA\n")
            nucleotide = AlignmentStats(
                identity=1.0, coverage=1.0, score=6.0, query_start=1, query_end=6,
                target_start=13, target_end=18, backend="minimap2", cigar="6M",
                query_coverage=1.0, aligned_pairs=6, aligned_blocks=[(1, 6, 13, 18)],
            )
            projection = {
                "protein_id": "A|A_gene|tx1", "query_start": 1, "query_end": 3,
                "target_start": 10, "target_end": 18, "strand": "+", "identity": 1.0,
                "identity_scope": "cds", "phase": "0", "cigar": "3M",
                "parent_id": "p1", "projection_status": "cds_target",
            }
            with patch("insiphy.annotation.local_alignment_stats", return_value=nucleotide), patch(
                "insiphy.annotation.protein_locus_exons", return_value=[projection]
            ):
                evidence = generate_sequence_evidence(input_dir, result_dir, aligner="miniprot")
            self.assertEqual(len(evidence), 1)
            self.assertEqual((evidence[0]["start"], evidence[0]["end"]), ("NA", "NA"))
            self.assertEqual(evidence[0]["homologous_dna_presence"], "unknown")
            self.assertEqual(evidence[0]["predicted_role"], "unknown")
            completed = complete_annotation(input_dir, result_dir)
            self.assertEqual(completed[0]["completion_call"], "ambiguous_evidence")
            self.assertEqual(completed[0]["interval"], "chrB:NA-NA:+")
            expected = {
                "interval_scope": "projected_cds_container",
                "protein_cds_query_start": "1", "protein_cds_query_end": "3",
                "protein_overlap_query_start": "2", "protein_overlap_query_end": "3",
                "reference_protein_query_start": "2", "reference_protein_query_end": "3",
            }
            for filename in ("sequence_synteny_evidence.tsv", "annotation_completion_candidates.tsv"):
                with self.subTest(table=filename):
                    row = read_tsv(result_dir / filename)[0]
                    self.assertEqual(row["anchor_interval_status"], "missing_flanking_homolog_anchor")
                    candidates = json.loads(row["interval_candidates"])
                    dna, protein = candidates
                    self.assertEqual(dna["interval_scope"], "aligned_sequence")
                    self.assertEqual((dna["start"], dna["end"]), (13, 18))
                    self.assertEqual((protein["start"], protein["end"]), (10, 18))
                    for field, value in expected.items():
                        self.assertEqual(str(protein[field]), value)
            self.assertEqual(evidence[0]["predicted_role_blocks"], "NA")

    def test_projection_conflict_records_only_the_true_block_intersection(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_minimal_case(
                Path(tmp), target_role="noncoding_exon"
            )
            occurrences = read_tsv(input_dir / "segment_occurrences.tsv")
            occurrences[0].update(
                start="5", end="10", role="CDS", coding_status="coding",
                cds_length="6", cds_phase="0",
            )
            source_left = {**occurrences[0], "occurrence_id": "A_left", "segment_id": "s_left", "start": "1", "end": "4", "role": "noncoding_exon", "coding_status": "noncoding", "cds_length": "0", "cds_phase": "."}
            source_right = {**source_left, "occurrence_id": "A_right", "segment_id": "s_right", "start": "11", "end": "14"}
            target_left = {**source_left, "occurrence_id": "B_left", "species": "B", "gene_copy_id": "B_gene", "segment_id": "t_left", "contig": "chrB", "start": "10", "end": "13"}
            target_right = {**target_left, "occurrence_id": "B_right", "segment_id": "t_right", "start": "30", "end": "33"}
            write_rows(
                input_dir / "segment_occurrences.tsv",
                [source_left, occurrences[0], source_right, target_left, occurrences[1], target_right],
                list(occurrences[0]),
            )
            paths = read_tsv(input_dir / "transcript_paths.tsv")
            paths[0].update(
                path_rank="2", start="5", end="10", role="CDS", coding_status="coding",
                cds_length="6", cds_phase="0", phase="0",
            )
            source_left_path = {**paths[0], "path_id": "A_left_path", "path_rank": "1", "occurrence_id": "A_left", "start": "1", "end": "4", "role": "noncoding_exon", "coding_status": "noncoding", "cds_phase": ".", "phase": "."}
            source_right_path = {**source_left_path, "path_id": "A_right_path", "path_rank": "3", "occurrence_id": "A_right", "start": "11", "end": "14"}
            target_left_path = {**source_left_path, "path_id": "B_left_path", "species": "B", "gene_copy_id": "B_gene", "occurrence_id": "B_left", "contig": "chrB", "start": "10", "end": "13"}
            target_mid_path = {**paths[1], "path_rank": "2"}
            target_right_path = {**target_left_path, "path_id": "B_right_path", "path_rank": "3", "occurrence_id": "B_right", "start": "30", "end": "33"}
            write_rows(
                input_dir / "transcript_paths.tsv",
                [source_left_path, paths[0], source_right_path, target_left_path, target_mid_path, target_right_path],
                list(paths[0]),
            )
            elements = read_tsv(result_dir / "element_correspondence.tsv")
            elements[0]["display_role"] = "CDS"
            def anchor_element(template, element_id, homology_id, occurrence_id, species, gene_copy_id):
                return {
                    **template,
                    "element_id": element_id,
                    "homology_id": homology_id,
                    "occurrence_id": occurrence_id,
                    "species": species,
                    "gene_copy_id": gene_copy_id,
                    "display_role": "noncoding_exon",
                    "membership_score": "1",
                    "membership_call": "core_member",
                }
            anchor_rows = [
                anchor_element(elements[0], "Eleft", "Hleft", "A_left", "A", "A_gene"),
                anchor_element(elements[0], "Eleft", "Hleft", "B_left", "B", "B_gene"),
                anchor_element(elements[0], "Eright", "Hright", "A_right", "A", "A_gene"),
                anchor_element(elements[0], "Eright", "Hright", "B_right", "B", "B_gene"),
            ]
            write_rows(result_dir / "element_correspondence.tsv", [*elements, *anchor_rows], list(elements[0]))
            (input_dir / "segment_sequences.fasta").write_text(
                ">A_left\nGGGG\n>A_e1\nATGAAA\n>A_right\nCCCC\n"
                ">B_left\nGGGG\n>B_i1\nACGT\n>B_right\nCCCC\n"
            )
            (input_dir / "gene_loci.fasta").write_text(
                ">A|A_gene|chrA:1-14:+\nGGGGATGAAACCCC\n"
                ">B|B_gene|chrB:1-40:+\n" + "T" * 40 + "\n"
            )
            (input_dir / "protein_sequences.fasta").write_text(
                ">A|A_gene|tx1\nMK\n"
            )
            no_dna_hit = AlignmentStats(
                identity=0.0, coverage=0.0, score=0.0,
                query_end=6, target_end=40, backend="minimap2",
            )
            projection = {
                "protein_id": "A|A_gene|tx1", "query_start": 1, "query_end": 2,
                "target_start": 18, "target_end": 23, "strand": "+", "identity": 1.0,
                "identity_scope": "cds", "phase": "0", "cigar": "6M",
                "parent_id": "projected_1", "projection_status": "cds_target",
            }
            with patch("insiphy.annotation.local_alignment_stats", return_value=no_dna_hit), patch(
                "insiphy.annotation.protein_locus_exons", return_value=[projection]
            ):
                evidence = generate_sequence_evidence(
                    input_dir, result_dir, aligner="miniprot"
                )
            row = next(row for row in evidence if row["homology_id"] == "H1")
            self.assertEqual(row["evidence_status"], "conflicts_annotation")
            self.assertEqual(row["homologous_dna_presence"], "present")
            self.assertEqual(row["predicted_exonic_role"], "CDS")
            self.assertEqual(row["supplied_annotation_role"], "noncoding_exon")
            predicted = json.loads(row["predicted_role_blocks"])[0]
            conflict = json.loads(row["annotation_conflict_blocks"])[0]
            self.assertEqual((predicted["target_start"], predicted["target_end"]), (18, 23))
            self.assertEqual((conflict["overlap_start"], conflict["overlap_end"]), (20, 23))
            self.assertEqual(conflict["occurrence_id"], "B_i1")
            self.assertEqual(conflict["transcript_ids"], ["tx1"])
            completed = complete_annotation(input_dir, result_dir)
            completed_row = next(row for row in completed if row["homology_id"] == "H1")
            self.assertEqual(completed_row["completion_call"], "boundary_conflict_candidate")
            self.assertEqual(
                json.loads(completed_row["annotation_conflict_blocks"])[0]["overlap_start"],
                20,
            )

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
            with patch("insiphy.annotation.local_alignment_stats", return_value=no_hit):
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
        with patch("insiphy.annotation.local_alignment_stats", return_value=spanning) as align:
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
                with patch("insiphy.annotation.local_alignment_stats") as align:
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
                with patch("insiphy.annotation.local_alignment_stats") as align:
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
        with patch("insiphy.annotation.local_alignment_stats", return_value=spanning):
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

            with patch("insiphy.annotation.local_alignment_stats", side_effect=align):
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


class AlignmentEvidenceTests(unittest.TestCase):
    def test_minimap2_preserves_multiple_hits_metadata(self):
        stdout = (
            "query\t4\t0\t4\t+\ttarget\t8\t0\t4\t4\t4\t60\tcg:Z:4M\tAS:i:4\ttp:A:P\n"
            "query\t4\t0\t4\t+\ttarget\t8\t4\t8\t4\t4\t0\tcg:Z:4M\tAS:i:4\ttp:A:S\n"
        )
        with patch("insiphy.alignment.shutil.which", return_value="/usr/bin/minimap2"), patch(
            "insiphy.alignment.subprocess.run",
            return_value=SimpleNamespace(returncode=0, stdout=stdout, stderr=""),
        ):
            stats = local_alignment_stats("AAAA", "AAAAAAAA", backend="minimap2")
        self.assertEqual(stats.hit_count, 2)
        self.assertEqual(stats.ambiguous_hit_count, 1)
        self.assertEqual(stats.mapping_quality, 60)
        self.assertEqual(stats.sequence_kind, "nucleotide")
        self.assertEqual(stats.score_scheme, "minimap2_AS_tag")
        self.assertEqual(stats.known_aligned_pairs, 4)
        self.assertEqual(stats.query_covered_bases, 4)
        self.assertEqual(stats.target_covered_bases, 4)
        self.assertFalse(stats.enumeration_complete)
        self.assertEqual(
            stats.incomplete_reason,
            "backend_did_not_guarantee_complete_hit_enumeration",
        )
        self.assertEqual(len(stats.alternative_hits), 1)
        self.assertEqual(stats.alternative_hits[0]["is_secondary"], 1)
        self.assertNotEqual(stats.alternative_hits[0]["candidate_id"], stats.candidate_id)
        self.assertEqual(stats.alternative_hits[0]["nt_identity"], "1")
        self.assertEqual(stats.alternative_hits[0]["aa_identity"], "NA")


class LocalCorrespondenceContractTests(unittest.TestCase):
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
        self.assertEqual(result.ambiguous_ids, frozenset())
        self.assertEqual(result.ambiguity_status, "unique_within_reported_candidates")

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


if __name__ == "__main__":
    unittest.main()
