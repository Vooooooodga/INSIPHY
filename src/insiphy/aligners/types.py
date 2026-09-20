"""aligners / types: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from typing import Optional
from insiphy.coordinates import Interval0, CoordinateBlock


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


KNOWN_NT = {"A", "C", "G", "T"}


MAX_INTERNAL_DP_CELLS = 250_000


MAX_MAFFT_PAIR_CELLS = 25_000_000


MAX_OPTIMAL_ALIGNMENTS = 64


NT_BLASTN_V1 = "nt_blastn_v1"


NT_BLASTN_V1_MATCH = 2.0


NT_BLASTN_V1_MISMATCH = -3.0


NT_BLASTN_V1_GAP_OPEN = -7.0


NT_BLASTN_V1_GAP_EXTEND = -2.0


@dataclass
class AlignmentStats:
    identity: float
    coverage: float
    score: float
    query_start: int = 1
    query_end: int = 0
    target_start: int = 1
    target_end: int = 0
    cigar: str = "NA"
    backend: str = "internal"
    query_coverage: float = 0.0
    target_coverage: float = 0.0
    aligned_pairs: int = 0
    strand: str = "+"
    aligned_blocks: list[tuple[int, int, int, int]] = field(default_factory=list)
    matches: int = 0
    mismatches: int = 0
    gap_bases: int = 0
    unknown_bases: int = 0
    alignment_mode: str = "global"
    alignment_meaning: str = "nucleotide alignment"
    query_span_coverage: float = 0.0
    target_span_coverage: float = 0.0
    mapping_quality: Optional[int] = None
    is_secondary: bool = False
    hit_count: int = 0
    ambiguous_hit_count: int = 0
    alternative_hits: list[dict] = field(default_factory=list)
    score_scheme: str = "unspecified"
    enumeration_complete: bool = True
    incomplete_reason: str = ""
    sequence_kind: str = "nucleotide"
    backend_version: Optional[str] = None
    raw_score: Optional[float] = None
    nt_identity: Optional[float] = None
    aa_identity: Optional[float] = None
    known_aligned_pairs: int = 0
    unknown_aligned_pairs: int = 0
    query_covered_bases: int = 0
    target_covered_bases: int = 0
    query_length: int = 0
    target_length: int = 0
    gap_blocks: list[dict] = field(default_factory=list)
    relative_strand: str = "+"
    candidate_id: Optional[str] = None
    alternative_candidate_ids: list[str] = field(default_factory=list)
    query_occurrence_id: Optional[str] = None
    target_occurrence_id: Optional[str] = None
    query_transcript_id: Optional[str] = None
    target_transcript_id: Optional[str] = None
    left_anchor_id: Optional[str] = None
    right_anchor_id: Optional[str] = None
    search_interval: Optional[dict] = None


@dataclass
class AlignmentGap:
    """One gap placement; exactly one half-open interval is empty."""

    query: Interval0
    target: Interval0

    def __post_init__(self):
        if (self.query.length == 0) == (self.target.length == 0):
            raise ValueError("an alignment gap must consume exactly one sequence")

    @property
    def gap_in(self) -> str:
        return "query" if self.query.length == 0 else "target"

    @property
    def query_cut0(self) -> int:
        return self.query.start0

    @property
    def target_cut0(self) -> int:
        return self.target.start0


@dataclass
class AlignmentCandidate:
    """One optimal alignment in 0-based, half-open internal coordinates."""

    query_interval: Interval0
    target_interval: Interval0
    aligned_blocks: tuple[CoordinateBlock, ...]
    gap_blocks: tuple[AlignmentGap, ...]
    identity: float
    coverage: float
    score: float
    cigar: str
    query_coverage: float
    target_coverage: float
    aligned_pairs: int
    matches: int
    mismatches: int
    gap_bases: int
    unknown_bases: int
    query_span_coverage: float
    target_span_coverage: float
    strand: str = "+"
    backend: str = "internal"
    score_scheme: str = NT_BLASTN_V1
    alignment_mode: str = "global"
    candidate_id: Optional[str] = None
    query_occurrence_id: Optional[str] = None
    target_occurrence_id: Optional[str] = None
    query_transcript_id: Optional[str] = None
    target_transcript_id: Optional[str] = None
    sequence_kind: str = "nucleotide"
    backend_version: Optional[str] = None
    raw_score: Optional[float] = None
    nt_identity: Optional[float] = None
    aa_identity: Optional[float] = None
    known_aligned_pairs: int = 0
    unknown_aligned_pairs: int = 0
    query_covered_bases: int = 0
    target_covered_bases: int = 0
    query_length: int = 0
    target_length: int = 0
    relative_strand: str = "+"
    mapping_quality: Optional[int] = None
    is_secondary: bool = False
    hit_count: int = 0
    alternative_candidate_ids: tuple[str, ...] = ()
    left_anchor_id: Optional[str] = None
    right_anchor_id: Optional[str] = None
    search_interval: Optional[dict] = None
    enumeration_complete: bool = True
    incomplete_reason: str = ""


@dataclass
class AlignmentCandidateSet:
    """Distinct optimal alignments reported for one bounded sequence pair."""

    candidates: list[AlignmentCandidate] = field(default_factory=list)
    enumeration_complete: bool = True
    incomplete_reason: str = ""
    score_scheme: str = NT_BLASTN_V1
    alignment_mode: str = "global"
    sequence_kind: str = "nucleotide"
    backend: str = "internal"
    backend_version: Optional[str] = None
    query_length: int = 0
    target_length: int = 0
    query_occurrence_id: Optional[str] = None
    target_occurrence_id: Optional[str] = None
    query_transcript_id: Optional[str] = None
    target_transcript_id: Optional[str] = None
    left_anchor_id: Optional[str] = None
    right_anchor_id: Optional[str] = None
    search_interval: Optional[dict] = None

    @property
    def primary(self) -> Optional[AlignmentCandidate]:
        return self.candidates[0] if self.candidates else None

    @property
    def candidate_ids(self) -> tuple[str, ...]:
        return tuple(candidate.candidate_id for candidate in self.candidates if candidate.candidate_id)

    @property
    def alternative_candidate_ids(self) -> tuple[str, ...]:
        return self.candidate_ids[1:]

    @property
    def hit_count(self) -> int:
        return len(self.candidates)


class AlignmentBackendError(RuntimeError):
    """Raised when an explicitly requested external alignment backend fails."""


def _gap_as_dict(gap: AlignmentGap) -> dict:
    if gap.gap_in == "query":
        return {
            "gap_in": "query",
            "query_cut0": gap.query_cut0,
            "target_start0": gap.target.start0,
            "target_end0": gap.target.end0,
        }
    return {
        "gap_in": "target",
        "query_start0": gap.query.start0,
        "query_end0": gap.query.end0,
        "target_cut0": gap.target_cut0,
    }
