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


class LocalCorrespondenceContractTestsSupport:
    pass
