"""mapping / membership: explicit implementation ownership."""
from __future__ import annotations

from collections import Counter
from collections import defaultdict
from intraphy.alignment import global_alignment_stats
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0
from intraphy.coordinates import parse_legacy_blocks
from intraphy.elements import collect_element_profiles
from intraphy.elements import element_class_for_occurrence
from intraphy.elements import element_role_from_occurrences
from intraphy.elements import membership_call
from intraphy.io import norm_state
from intraphy.io import parse_fasta
from intraphy.io import read_tsv
from intraphy.io import to_float
from intraphy.io import uniq
from intraphy.io import write_tsv
from intraphy.mapping.member_intervals import _block_union_length
from intraphy.mapping.member_intervals import _blocks_overlap
from intraphy.mapping.member_intervals import _candidate_block_evidence
from intraphy.mapping.member_intervals import _edge_eligible
from intraphy.mapping.member_intervals import _genomic_block_records
from intraphy.mapping.member_intervals import _hard_membership_eligible
from intraphy.mapping.member_intervals import _hard_position_eligible
from intraphy.mapping.member_intervals import _parse_genomic_blocks
from intraphy.mapping.member_intervals import _tagged_membership_blocks
from intraphy.mapping.member_intervals import _tokens
from intraphy.mapping.member_refinement import _membership_match_details
from intraphy.mapping.membership_fields import CORRESPONDENCE_EVIDENCE_FIELDS
from intraphy.mapping.reference_coverage import _assign_reference_coverage
from intraphy.topology import SpeciesTree
import json


