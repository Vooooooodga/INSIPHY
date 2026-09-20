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

from support_evidence_alignment_evidence import AlignmentEvidenceTestsSupport

class AlignmentEvidenceTests(AlignmentEvidenceTestsSupport, unittest.TestCase):
    def test_minimap2_preserves_multiple_hits_metadata(self):
        stdout = (
            "query\t4\t0\t4\t+\ttarget\t8\t0\t4\t4\t4\t60\tcg:Z:4M\tAS:i:4\ttp:A:P\n"
            "query\t4\t0\t4\t+\ttarget\t8\t4\t8\t4\t4\t0\tcg:Z:4M\tAS:i:4\ttp:A:S\n"
        )
        with patch("intraphy.alignment.shutil.which", return_value="/usr/bin/minimap2"), patch(
            "intraphy.alignment.subprocess.run",
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
