"""mapping / chains: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0
from intraphy.coordinates import local_interval_to_genome
from intraphy.mapping.candidate_coordinates import _block_signature
from intraphy.mapping.candidate_serialization import _public_candidate_record
from intraphy.mapping.chain_candidates import _collect_chain_candidates
from intraphy.mapping.chain_qualification import _qualify_chain_edges
from intraphy.mapping.chain_types import ChainCandidate
from intraphy.mapping.chain_types import ChainPathMembership
from intraphy.mapping.chain_types import DEFAULT_CHAIN_CONFIGURATION
from intraphy.mapping.ordered_chains import _apply_ordered_candidate_chains
from intraphy.mapping.ordered_paths import genomic_candidate_chain
from intraphy.mapping.ordered_paths import ordered_candidate_chain
from intraphy.mapping.path_evaluation import _evaluate_candidate_paths
from intraphy.mapping.path_membership import _candidate_copy_interval
from intraphy.mapping.path_membership import _candidate_path_memberships
from intraphy.mapping.path_membership import _chain_precedes
from intraphy.mapping.path_membership import _copy_transcription_bounds
from intraphy.mapping.policies import occurrence_copy_key
from intraphy.storage.values import to_float
import json


