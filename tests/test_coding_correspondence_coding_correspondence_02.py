import json

import unittest

from unittest.mock import patch

from intraphy.alignment import AlignmentBackendError

from intraphy.coding_correspondence import CodingProjectionIndex, build_coding_transcript

from intraphy.correspondence import match_total

from intraphy.preprocess import cheap_match_evidence, cluster_segments

from intraphy.structural_sites import _mapped_reference_coordinate, _within_exon_boundary

def coding_rows(name, species, start, sequence, *, strand="+", phase="0", cds=None, transcript="tx"):
    end = start + len(sequence) - 1
    cds_left, cds_right = cds or (1, len(sequence))
    if strand == "-":
        coding_start, coding_end = end - cds_right + 1, end - cds_left + 1
    else:
        coding_start, coding_end = start + cds_left - 1, start + cds_right - 1
    occurrence = {
        "occurrence_id": name, "family_id": "fam", "species": species,
        "gene_copy_id": f"g{species}", "transcript_id": transcript,
        "role": "exon", "presence_status": "present", "contig": "chr",
        "start": str(start), "end": str(end), "strand": strand, "phase": phase,
        "transcript_order": "1", "splice_motif_score": "0.5",
    }
    path = {
        **occurrence, "path_id": f"{species}:{transcript}:{name}",
        "cds_phase": phase, "cds_length": str(coding_end - coding_start + 1),
        "cds_intervals": f"{coding_start}-{coding_end}",
        "path_status": "annotated_transcript_path", "coding_status": "coding",
    }
    return occurrence, path

def unchanged_protein_alignment(records, mode="linsi", threads=1):
    sequences = dict(records)
    if len({len(sequence) for sequence in sequences.values()}) > 1:
        raise AssertionError("This fixture requires equal-length aligned proteins")
    return sequences

def block_coordinates0(blocks):
    return tuple(
        (
            block.query.start0, block.query.end0,
            block.target.start0, block.target.end0,
        )
        for block in blocks
    )

def weak_dna_evidence(left, right, sequences, context, **kwargs):
    evidence = cheap_match_evidence(left, right, context, alignment_backend="mafft")
    evidence.update(
        alignment_score=0.4, coverage_score=1.0, sequence_score=0.58,
        total_score=evidence["total_score"] + 0.34 * 0.4 + 0.14,
        alignment_strand="+", alignment_cigar="3M", alignment_mode="overlap",
        alignment_meaning="deterministic DNA alignment fixture",
        alignment_requested_backend="mafft", query_alignment_start=1,
        query_alignment_end=3, target_alignment_start=1, target_alignment_end=3,
        projected_reference_blocks="1-3:1-3",
    )
    return evidence

from support_coding_correspondence_coding_correspondence import CodingCorrespondenceTestsSupport

class CodingCorrespondenceTests(CodingCorrespondenceTestsSupport, unittest.TestCase):
    def test_protein_blocks_enter_graph_and_junction_preserving_raw_dna(self):
        for strand in ("+", "-"):
            with self.subTest(strand=strand), patch(
                "intraphy.mapping.clustering.match_evidence", side_effect=weak_dna_evidence,
            ), patch(
                "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
                side_effect=unchanged_protein_alignment,
            ):
                occurrences, paths, sequences = self.split_fixture(strand)
                homology, matches = cluster_segments(
                    occurrences, sequences, transcript_paths=paths,
                )
                by_occurrence = {row["occurrence_id"]: row["homology_id"] for row in homology}
                self.assertEqual(len(set(by_occurrence.values())), 1)
                for match in matches:
                    self.assertEqual(match["alignment_score"], "0.4")
                    self.assertEqual(match["projected_reference_blocks"], "1-3:1-3")
                    self.assertEqual(match["alignment_cigar"], "3M")
                    self.assertEqual(match["dna_match_status"], "low_similarity")
                    self.assertEqual(match["match_status"], "mapped")
                    self.assertEqual(match["correspondence_basis"], "annotated_CDS_protein")
                    self.assertEqual(match["protein_metrics_scope"], "best_transcript_pair")
                    self.assertEqual(match["protein_best_query_transcript"], "tx")
                    self.assertEqual(match["protein_best_target_transcript"], "tx")
                    self.assertGreater(match_total(match), float(match["total_score"]))
                boundary = _within_exon_boundary(
                    "fam", "EG_1", "left", "right", {row["occurrence_id"]: row for row in occurrences},
                    {"EG_1": "ref"}, matches,
                )
                self.assertIsNotNone(boundary)
                self.assertEqual((boundary["donor_projection"], boundary["acceptor_projection"]), (6, 7))
                self.assertIsNone(_mapped_reference_coordinate(matches, "right", "ref", 12))

    def test_repeated_reference_coverage_does_not_merge_distinct_exons(self):
        first, first_path = coding_rows("first", "A", 1, "ATGGCTGCT", transcript="one")
        second, second_path = coding_rows("second", "A", 101, "ATGGCTGCT", transcript="two")
        reference, reference_path = coding_rows("ref", "B", 1, "ATGGCTGCT")
        with patch("intraphy.mapping.clustering.match_evidence", side_effect=weak_dna_evidence), patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=unchanged_protein_alignment,
        ):
            index = CodingProjectionIndex(
                [first, second, reference],
                {"first": "ATGGCTGCT", "second": "ATGGCTGCT", "ref": "ATGGCTGCT"},
                [first_path, second_path, reference_path],
            )
            evidence = index.evidence("first", "ref")
            homology, _matches = cluster_segments(
                [first, second, reference], {"first": "ATGGCTGCT", "second": "ATGGCTGCT", "ref": "ATGGCTGCT"},
                transcript_paths=[first_path, second_path, reference_path],
            )
        self.assertEqual(evidence["protein_mapping_status"], "ambiguous_mapping")
        self.assertFalse(evidence["protein_position_eligible"])
        self.assertEqual(evidence["protein_competing_occurrences"], "query:second")
        by_occurrence = {row["occurrence_id"]: row["homology_id"] for row in homology}
        self.assertNotEqual(by_occurrence["first"], by_occurrence["second"])

    def test_unavailable_metadata_preserves_dna_and_internal_backend_never_calls_protein(self):
        occurrences, paths, sequences = self.split_fixture()
        paths[0]["cds_phase"] = "."
        for backend, expected in (("mafft", "unavailable"), ("internal", "disabled_backend")):
            with self.subTest(backend=backend), patch(
                "intraphy.mapping.clustering.match_evidence", side_effect=weak_dna_evidence,
            ), patch("intraphy.coding_correspondence.alignment.protein_multiple_alignment") as aligner:
                _homology, matches = cluster_segments(occurrences, sequences, transcript_paths=paths, aligner=backend)
                aligner.assert_not_called()
                self.assertTrue(all(row["protein_status"] == expected for row in matches))
                self.assertTrue(all(row["alignment_score"] == "0.4" and row["match_status"] == "low_similarity" for row in matches))
                if backend == "mafft":
                    self.assertTrue(all("unknown_CDS_phase" in row["protein_unavailable_reason"] for row in matches))

    def test_external_aligner_error_propagates(self):
        occurrences, paths, sequences = self.split_fixture()
        with patch("intraphy.mapping.clustering.match_evidence", side_effect=weak_dna_evidence), patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=AlignmentBackendError("explicit MAFFT failure"),
        ):
            with self.assertRaisesRegex(AlignmentBackendError, "explicit MAFFT failure"):
                cluster_segments(occurrences, sequences, transcript_paths=paths)

    def test_cache_is_run_local_and_bounded_by_mapped_bases(self):
        occurrences, paths, sequences = self.split_fixture()
        index = CodingProjectionIndex(occurrences, sequences, paths)
        index.MAX_CACHED_BASE_PAIRS = 1
        with patch("intraphy.coding_correspondence.alignment.protein_multiple_alignment", side_effect=unchanged_protein_alignment):
            self.assertEqual(index.evidence("left", "ref")["protein_status"], "supported")
        self.assertFalse(index.cache)
        self.assertEqual(index.cached_bases, 0)
        self.assertFalse(CodingProjectionIndex(occurrences, sequences, paths).cache)
