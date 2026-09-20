"""candidate chain: explicit implementation ownership."""
from __future__ import annotations

from dataclasses import dataclass
from intraphy.coordinates import Interval0
from intraphy.mapping.chain_graph import _ordered_before
from intraphy.mapping.chain_graph import _parent_order_compatible
from intraphy.mapping.chain_graph import _path_membership_for_context
from intraphy.mapping.chain_graph import _solve_context
from intraphy.mapping.chain_graph import _top_complete_paths
from intraphy.mapping.chain_types import CandidateChainResult
from intraphy.mapping.chain_types import ChainCandidate
from intraphy.mapping.chain_types import ChainConfiguration
from intraphy.mapping.chain_types import ChainContextSummary
from intraphy.mapping.chain_types import ChainPathMembership
from intraphy.mapping.chain_types import CoverageClassification
from intraphy.mapping.chain_types import DEFAULT_CHAIN_CONFIGURATION
from intraphy.mapping.coverage import classify_reference_coverage
from intraphy.mapping.ordered_paths import genomic_candidate_chain
from intraphy.mapping.ordered_paths import ordered_candidate_chain
from math import isfinite
from typing import FrozenSet
from typing import Iterable
from typing import Mapping
from typing import Optional
from typing import Tuple


