"""aligners / pairwise: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.aligners.columns import _empty_stats
from intraphy.aligners.columns import _stats_from_alignment_columns
from intraphy.aligners.external import _external_lastz_stats
from intraphy.aligners.external import _external_mafft_stats
from intraphy.aligners.external import _external_minimap2_stats
from intraphy.aligners.external import _external_miniprot_stats
from intraphy.aligners.external import _mafft_pair_alignment
from intraphy.aligners.legacy import _candidate_as_legacy_hit
from intraphy.aligners.legacy import _legacy_stats_from_candidate
from intraphy.aligners.short import anchored_short_alignment
from intraphy.aligners.types import AlignmentBackendError
from intraphy.aligners.types import AlignmentStats
from intraphy.aligners.types import MAX_INTERNAL_DP_CELLS
from intraphy.aligners.types import NT_BLASTN_V1
from typing import Optional
import shutil


def _trim_terminal_overhangs(query_aligned: str, target_aligned: str) -> tuple[str, str, int, int]:
    paired = [idx for idx, (q_char, t_char) in enumerate(zip(query_aligned, target_aligned)) if q_char != "-" and t_char != "-"]
    if not paired:
        return "", "", 0, 0
    first = paired[0]
    last = paired[-1]
    query_offset = sum(1 for char in query_aligned[:first] if char != "-")
    target_offset = sum(1 for char in target_aligned[:first] if char != "-")
    return query_aligned[first : last + 1], target_aligned[first : last + 1], query_offset, target_offset


def overlap_alignment_stats(query: str, target: str, backend: str = "mafft", threads: int = 1) -> AlignmentStats:
    """Project a pairwise overlap by trimming terminal overhangs from an external global alignment."""
    backend = (backend or "mafft").lower()
    if backend != "mafft":
        raise AlignmentBackendError(f"unsupported overlap alignment backend: {backend}")
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return _empty_stats(
            len(query),
            len(target),
            "mafft_overlap",
            alignment_mode="overlap_projection",
            alignment_meaning="terminal-overhang-trimmed MAFFT overlap projection",
        )
    left, right = _mafft_pair_alignment(query, target, threads)
    trimmed_query, trimmed_target, query_offset, target_offset = _trim_terminal_overhangs(left, right)
    return _stats_from_alignment_columns(
        trimmed_query,
        trimmed_target,
        len(query),
        len(target),
        query_start0=query_offset,
        target_start0=target_offset,
        score=0.0,
        backend="mafft_overlap",
        mode="overlap",
        alignment_mode="overlap_projection",
        alignment_meaning="terminal-overhang-trimmed MAFFT overlap projection; not a true local algorithm",
    )


def _pairwise_aligner_stats(seq_a: str, seq_b: str, mode: str, **candidate_context) -> AlignmentStats:
    result = anchored_short_alignment(seq_a, seq_b, mode=mode, **candidate_context)
    if result.primary is None:
        stats = _empty_stats(len(seq_a or ""), len(seq_b or ""), "internal")
        stats.score_scheme = result.score_scheme
        stats.enumeration_complete = result.enumeration_complete
        stats.incomplete_reason = result.incomplete_reason
        return stats
    stats = _legacy_stats_from_candidate(result.primary)
    stats.hit_count = len(result.candidates)
    stats.ambiguous_hit_count = max(0, len(result.candidates) - 1)
    stats.alternative_hits = [
        _candidate_as_legacy_hit(candidate, rank)
        for rank, candidate in enumerate(result.candidates[1:], start=2)
    ]
    stats.enumeration_complete = result.enumeration_complete
    stats.incomplete_reason = result.incomplete_reason
    return stats


def global_alignment_stats(seq_a: str, seq_b: str, backend: str = "internal", threads: int = 1) -> AlignmentStats:
    backend = (backend or "internal").lower()
    if backend == "auto":
        return global_alignment_stats(seq_a, seq_b, "mafft", threads)
    if backend == "minimap2":
        return _external_minimap2_stats(seq_a, seq_b, "global", threads)
    if backend == "miniprot":
        return _external_miniprot_stats(seq_a, seq_b, threads)
    if backend == "mafft":
        return _external_mafft_stats(seq_a, seq_b, threads)
    if backend == "lastz":
        raise AlignmentBackendError("LASTZ is a local nucleotide backend; use local_alignment_stats")
    if backend != "internal":
        raise AlignmentBackendError(f"unsupported alignment backend: {backend}")
    seq_a = (seq_a or "").upper()
    seq_b = (seq_b or "").upper()
    if not seq_a or not seq_b:
        stats = _empty_stats(len(seq_a), len(seq_b), "internal")
        stats.score_scheme = NT_BLASTN_V1
        return stats
    if len(seq_a) * len(seq_b) > MAX_INTERNAL_DP_CELLS:
        raise AlignmentBackendError(
            f"internal global alignment requires {len(seq_a) * len(seq_b)} DP cells; "
            "select MAFFT or minimap2"
        )
    return _pairwise_aligner_stats(seq_a, seq_b, "global")


def local_alignment_stats(
    query: str,
    target: str,
    backend: str = "internal",
    threads: int = 1,
    *,
    query_occurrence_id: Optional[str] = None,
    target_occurrence_id: Optional[str] = None,
    query_transcript_id: Optional[str] = None,
    target_transcript_id: Optional[str] = None,
    left_anchor_id: Optional[str] = None,
    right_anchor_id: Optional[str] = None,
    search_interval: Optional[dict] = None,
) -> AlignmentStats:
    backend = (backend or "internal").lower()
    if backend == "auto":
        if shutil.which("lastz"):
            return _external_lastz_stats(query, target, "local", threads)
        if shutil.which("minimap2"):
            return _external_minimap2_stats(query, target, "local", threads)
        raise AlignmentBackendError("local auto requires LASTZ or minimap2 on PATH; request internal explicitly for small Biopython alignment")
    if backend == "minimap2":
        return _external_minimap2_stats(query, target, "local", threads)
    if backend == "miniprot":
        return _external_miniprot_stats(query, target, threads)
    if backend == "mafft":
        return overlap_alignment_stats(query, target, backend="mafft", threads=threads)
    if backend == "lastz":
        return _external_lastz_stats(query, target, "local", threads)
    if backend != "internal":
        raise AlignmentBackendError(f"unsupported alignment backend: {backend}")
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        stats = _empty_stats(len(query), len(target), "internal")
        stats.score_scheme = NT_BLASTN_V1
        return stats
    if len(query) * len(target) > MAX_INTERNAL_DP_CELLS:
        raise AlignmentBackendError(
            f"internal local alignment requires {len(query) * len(target)} DP cells; select LASTZ or minimap2"
        )
    return _pairwise_aligner_stats(
        query,
        target,
        "local",
        query_occurrence_id=query_occurrence_id,
        target_occurrence_id=target_occurrence_id,
        query_transcript_id=query_transcript_id,
        target_transcript_id=target_transcript_id,
        left_anchor_id=left_anchor_id,
        right_anchor_id=right_anchor_id,
        search_interval=search_interval,
    )


def best_ungapped_hit(query: str, target: str) -> AlignmentStats:
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return AlignmentStats(0.0, 0.0, 0.0)
    if len(target) < len(query):
        span = len(target)
        matches = sum(1 for a, b in zip(query[:span], target) if a == b and a != "N" and b != "N")
        return AlignmentStats(matches / max(1, span), span / max(1, len(query)), float(matches), 1, span, 1, span, f"{span}M")
    best = AlignmentStats(-1.0, 1.0, 0.0, 1, len(query), 1, len(query), f"{len(query)}M")
    qlen = len(query)
    for offset in range(0, len(target) - qlen + 1):
        window = target[offset : offset + qlen]
        matches = sum(1 for a, b in zip(query, window) if a == b and a != "N" and b != "N")
        identity = matches / qlen
        if identity > best.identity:
            best = AlignmentStats(identity, 1.0, float(matches), 1, qlen, offset + 1, offset + qlen, f"{qlen}M")
    return best
