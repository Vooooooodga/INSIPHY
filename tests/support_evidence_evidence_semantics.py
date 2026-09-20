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


class EvidenceSemanticsTestsSupport:
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
