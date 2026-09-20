"""mapping / chain types: explicit implementation ownership."""
from __future__ import annotations

from dataclasses import dataclass
from intraphy.coordinates import Interval0
from math import isfinite
from typing import FrozenSet
from typing import Optional
from typing import Tuple


@dataclass(frozen=True)
class ChainConfiguration:
    name: str = "intraphy_nt_chain_delta_v1"
    nucleotide_match_reward: float = 2.0
    relative_delta: float = 0.05

    def score_delta(self, best_score: float, score_scheme: str) -> float:
        if str(score_scheme).startswith("nt_"):
            return max(
                2.0 * self.nucleotide_match_reward,
                self.relative_delta * abs(float(best_score)),
            )
        return 0.0

    @property
    def delta_rule(self):
        return "max(2*match_reward,relative_delta*abs(best_score))"


DEFAULT_CHAIN_CONFIGURATION = ChainConfiguration()


@dataclass(frozen=True)
class ChainPathMembership:
    """One compatible query/target transcript-path pairing for a candidate."""

    query_path_id: str
    target_path_id: str
    query_order: int
    target_order: int
    query_contig: str
    target_contig: str
    query_strand: str
    target_strand: str
    query_parent_id: str = ""
    target_parent_id: str = ""

    def __post_init__(self):
        if not self.query_path_id or not self.target_path_id:
            raise ValueError("chain path membership requires both path IDs")
        if self.query_strand not in {"+", "-"} or self.target_strand not in {"+", "-"}:
            raise ValueError("chain path membership requires explicit gene strands")

    @property
    def context(self):
        return (
            self.query_path_id,
            self.target_path_id,
            self.query_contig,
            self.target_contig,
            self.query_strand,
            self.target_strand,
        )


@dataclass(frozen=True)
class ChainCandidate:
    candidate_id: str
    query: Interval0
    target: Interval0
    score: float
    score_scheme: str
    relative_strand: str = "+"
    path_memberships: Tuple[ChainPathMembership, ...] = ()

    def __post_init__(self):
        if not self.candidate_id:
            raise ValueError("candidate_id is required")
        if not isfinite(float(self.score)):
            raise ValueError(f"candidate score must be finite: {self.candidate_id}")
        if not self.score_scheme:
            raise ValueError(f"score_scheme is required: {self.candidate_id}")
        if self.relative_strand not in {"+", "-"}:
            raise ValueError(f"invalid relative strand: {self.relative_strand}")


@dataclass(frozen=True)
class ChainContextSummary:
    context: Optional[tuple]
    best_score: float
    retained_ids: FrozenSet[str]
    best_path_member_ids: FrozenSet[str]
    retained_edges: FrozenSet[Tuple[str, str]]
    start_anchor_ids: FrozenSet[str]
    end_anchor_ids: FrozenSet[str]
    best_path_count_capped: int
    near_optimal_path_count_capped: int
    ambiguous_ids: FrozenSet[str] = frozenset()


@dataclass(frozen=True)
class CandidateChainResult:
    best_score: float
    retained_ids: FrozenSet[str]
    best_path_member_ids: FrozenSet[str]
    score_delta: float
    local_mode: bool
    retained_edges: FrozenSet[Tuple[str, str]] = frozenset()
    start_anchor_ids: FrozenSet[str] = frozenset()
    end_anchor_ids: FrozenSet[str] = frozenset()
    ambiguity_status: str = "unique_within_reported_candidates"
    configuration_name: str = DEFAULT_CHAIN_CONFIGURATION.name
    best_path_count_capped: int = 0
    near_optimal_path_count_capped: int = 0
    ambiguous_ids: FrozenSet[str] = frozenset()
    context_summaries: Tuple[ChainContextSummary, ...] = ()


@dataclass(frozen=True)
class CoverageClassification:
    relation: str
    covered_bases: int
    overlap_bases: int
    uncovered_bases: Optional[int]
