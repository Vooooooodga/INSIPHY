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
            with patch("intraphy.evidence.nucleotide_search.local_alignment_stats", return_value=stat) as align:
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
            with patch("intraphy.evidence.nucleotide_search.local_alignment_stats", return_value=stat):
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
        with patch("intraphy.evidence.protein.protein_locus_exons", return_value=rows) as project:
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
            with patch("intraphy.evidence.nucleotide_search.local_alignment_stats", return_value=nucleotide), patch(
                "intraphy.evidence.protein.protein_locus_exons", return_value=[projection]
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
            with patch("intraphy.evidence.nucleotide_search.local_alignment_stats", return_value=no_dna_hit), patch(
                "intraphy.evidence.protein.protein_locus_exons", return_value=[projection]
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
