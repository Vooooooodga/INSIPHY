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

    def test_singleton_annotated_exon_seeds_nucleotide_presence_in_miniprot_mode(self):
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
            self.assertEqual(rows[0]["evidence_status"], "homologous_sequence_candidate")
            self.assertEqual(rows[0]["inferred_role"], "unknown")
            self.assertEqual(rows[0]["alignment_backend"], "minimap2")
            self.assertEqual(rows[0]["evidence_aligner"], "miniprot")
            self.assertEqual(rows[0]["interval_scope"], "aligned_sequence")
            self.assertEqual(json.loads(rows[0]["interval_candidates"])[0]["interval_scope"], "aligned_sequence")
            completed = complete_annotation(input_dir, result_dir)
            self.assertEqual(completed[0]["interval_scope"], "aligned_sequence")
            align.assert_called_once_with("ACGT", "T" * 19 + "ACGT" + "T" * 17, backend="minimap2", threads=1)

    def test_protein_projection_remains_predicted_candidate(self):
        row = {
            "evidence_status": "supports_hidden_segment",
            "inferred_event": "protein_cds_projection",
            "inferred_role": "predicted_CDS",
            "predicted_role": "CDS",
            "frame_status": "unknown",
        }
        self.assertEqual(completion_call(row, 0.95, 0.55), "predicted_exon_candidate")

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

    def test_partial_protein_query_preserves_container_scope_in_tables_and_json(self):
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
            self.assertEqual((evidence[0]["start"], evidence[0]["end"]), (10, 18))
            completed = complete_annotation(input_dir, result_dir)
            self.assertEqual(completed[0]["completion_call"], "predicted_exon_candidate")
            self.assertEqual(completed[0]["interval"], "chrB:10-18:+")
            expected = {
                "interval_scope": "projected_cds_container",
                "protein_cds_query_start": "1", "protein_cds_query_end": "3",
                "protein_overlap_query_start": "2", "protein_overlap_query_end": "3",
                "reference_protein_query_start": "2", "reference_protein_query_end": "3",
            }
            for filename in ("sequence_synteny_evidence.tsv", "annotation_completion_candidates.tsv"):
                with self.subTest(table=filename):
                    row = read_tsv(result_dir / filename)[0]
                    for field, value in expected.items():
                        self.assertEqual(row[field], value)
                    candidates = json.loads(row["interval_candidates"])
                    dna, protein = candidates
                    self.assertEqual(dna["interval_scope"], "aligned_sequence")
                    self.assertEqual((dna["start"], dna["end"]), (13, 18))
                    self.assertEqual((protein["start"], protein["end"]), (10, 18))
                    for field, value in expected.items():
                        self.assertEqual(str(protein[field]), value)

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

    def test_ambiguous_primary_mapping_keeps_role_unknown(self):
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
            self.assertEqual(row["evidence_status"], "homologous_sequence_candidate")
            self.assertEqual(row["inferred_role"], "unknown")
            self.assertEqual(row["primary_mapping_status"], "ambiguous_repeated_mapping")
            self.assertEqual(row["correspondence_status"], "unknown")
            self.assertEqual((row["start"], row["end"]), ("NA", "NA"))
            self.assertEqual([candidate["start"] for candidate in json.loads(row["interval_candidates"])], [20, 30])
            completed = complete_annotation(input_dir, result_dir)
            row = next(row for row in completed if row["homology_id"] == "H1")
            self.assertEqual(row["completion_call"], "homologous_sequence_candidate")
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
        self.assertEqual(len(stats.alternative_hits), 1)
        self.assertEqual(stats.alternative_hits[0]["is_secondary"], 1)


if __name__ == "__main__":
    unittest.main()
