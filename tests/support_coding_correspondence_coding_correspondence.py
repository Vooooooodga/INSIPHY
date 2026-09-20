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


class CodingCorrespondenceTestsSupport:
    def split_fixture(self, strand="+"):
        left_start, right_start = (101, 201) if strand == "+" else (301, 201)
        left, left_path = coding_rows("left", "A", left_start, "ATGG", strand=strand)
        right, right_path = coding_rows("right", "A", right_start, "AAGCTTAATTTT", strand=strand, phase="2", cds=(1, 8))
        reference, reference_path = coding_rows("ref", "B", 501, "CCATGGAAGCTTAAGG", cds=(3, 14))
        occurrences = [left, right, reference]
        paths = [left_path, right_path, reference_path]
        sequences = {"left": "ATGG", "right": "AAGCTTAATTTT", "ref": "CCATGGAAGCTTAAGG"}
        return occurrences, paths, sequences
