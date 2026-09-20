"""aligners / legacy: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.aligners.types import AlignmentCandidate
from intraphy.aligners.types import AlignmentStats
from intraphy.aligners.types import _gap_as_dict
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import CoordinateBlock


def _candidate_from_legacy_stats(stats: AlignmentStats, gap_blocks=()) -> AlignmentCandidate:
    blocks = tuple(
        CoordinateBlock(
            ClosedInterval1(int(q_start), int(q_end)).to_interval0(),
            ClosedInterval1(int(t_start), int(t_end)).to_interval0(),
        )
        for q_start, q_end, t_start, t_end in stats.aligned_blocks
    )
    return AlignmentCandidate(
        query_interval=ClosedInterval1(
            int(stats.query_start), int(stats.query_end),
        ).to_interval0(),
        target_interval=ClosedInterval1(
            int(stats.target_start), int(stats.target_end),
        ).to_interval0(),
        aligned_blocks=blocks,
        gap_blocks=tuple(gap_blocks),
        identity=stats.identity,
        coverage=stats.coverage,
        score=stats.score,
        cigar=stats.cigar,
        query_coverage=stats.query_coverage,
        target_coverage=stats.target_coverage,
        aligned_pairs=stats.aligned_pairs,
        matches=stats.matches,
        mismatches=stats.mismatches,
        gap_bases=stats.gap_bases,
        unknown_bases=stats.unknown_bases,
        query_span_coverage=stats.query_span_coverage,
        target_span_coverage=stats.target_span_coverage,
        strand=stats.strand,
        backend=stats.backend,
        score_scheme=stats.score_scheme,
        alignment_mode=stats.alignment_mode,
        candidate_id=stats.candidate_id,
        query_occurrence_id=stats.query_occurrence_id,
        target_occurrence_id=stats.target_occurrence_id,
        query_transcript_id=stats.query_transcript_id,
        target_transcript_id=stats.target_transcript_id,
        sequence_kind=stats.sequence_kind,
        backend_version=stats.backend_version,
        raw_score=stats.raw_score if stats.raw_score is not None else stats.score,
        nt_identity=stats.nt_identity if stats.nt_identity is not None else stats.identity,
        aa_identity=stats.aa_identity,
        known_aligned_pairs=stats.known_aligned_pairs or stats.aligned_pairs,
        unknown_aligned_pairs=stats.unknown_aligned_pairs,
        query_covered_bases=stats.query_covered_bases or stats.aligned_pairs,
        target_covered_bases=stats.target_covered_bases or stats.aligned_pairs,
        query_length=stats.query_length,
        target_length=stats.target_length,
        relative_strand=stats.relative_strand or stats.strand,
        mapping_quality=stats.mapping_quality,
        is_secondary=stats.is_secondary,
        hit_count=stats.hit_count,
        alternative_candidate_ids=tuple(stats.alternative_candidate_ids),
        left_anchor_id=stats.left_anchor_id,
        right_anchor_id=stats.right_anchor_id,
        search_interval=stats.search_interval,
        enumeration_complete=stats.enumeration_complete,
        incomplete_reason=stats.incomplete_reason,
    )


def _legacy_coordinate_tuple(block: CoordinateBlock) -> tuple[int, int, int, int]:
    query = ClosedInterval1.from_interval0(block.query)
    target = ClosedInterval1.from_interval0(block.target)
    return query.start, query.end, target.start, target.end


def _legacy_stats_from_candidate(candidate: AlignmentCandidate) -> AlignmentStats:
    query_interval = ClosedInterval1.from_interval0(candidate.query_interval)
    target_interval = ClosedInterval1.from_interval0(candidate.target_interval)
    stats = AlignmentStats(
        identity=candidate.identity,
        coverage=candidate.coverage,
        score=candidate.score,
        query_start=query_interval.start,
        query_end=query_interval.end,
        target_start=target_interval.start,
        target_end=target_interval.end,
        cigar=candidate.cigar,
        backend=candidate.backend,
        query_coverage=candidate.query_coverage,
        target_coverage=candidate.target_coverage,
        aligned_pairs=candidate.aligned_pairs,
        strand=candidate.strand,
        aligned_blocks=[
            _legacy_coordinate_tuple(block) for block in candidate.aligned_blocks
        ],
        matches=candidate.matches,
        mismatches=candidate.mismatches,
        gap_bases=candidate.gap_bases,
        unknown_bases=candidate.unknown_bases,
        alignment_mode=candidate.alignment_mode,
        alignment_meaning=f"bounded short nucleotide {candidate.alignment_mode} alignment",
        query_span_coverage=candidate.query_span_coverage,
        target_span_coverage=candidate.target_span_coverage,
        mapping_quality=candidate.mapping_quality,
        is_secondary=candidate.is_secondary,
        hit_count=candidate.hit_count,
        score_scheme=candidate.score_scheme,
        enumeration_complete=candidate.enumeration_complete,
        incomplete_reason=candidate.incomplete_reason,
        sequence_kind=candidate.sequence_kind,
        backend_version=candidate.backend_version,
        raw_score=candidate.raw_score,
        nt_identity=candidate.nt_identity,
        aa_identity=candidate.aa_identity,
        known_aligned_pairs=candidate.known_aligned_pairs,
        unknown_aligned_pairs=candidate.unknown_aligned_pairs,
        query_covered_bases=candidate.query_covered_bases,
        target_covered_bases=candidate.target_covered_bases,
        query_length=candidate.query_length,
        target_length=candidate.target_length,
        gap_blocks=[_gap_as_dict(gap) for gap in candidate.gap_blocks],
        relative_strand=candidate.relative_strand,
        candidate_id=candidate.candidate_id,
        alternative_candidate_ids=list(candidate.alternative_candidate_ids),
        query_occurrence_id=candidate.query_occurrence_id,
        target_occurrence_id=candidate.target_occurrence_id,
        query_transcript_id=candidate.query_transcript_id,
        target_transcript_id=candidate.target_transcript_id,
        left_anchor_id=candidate.left_anchor_id,
        right_anchor_id=candidate.right_anchor_id,
        search_interval=candidate.search_interval,
    )
    return stats


def _candidate_as_legacy_hit(candidate: AlignmentCandidate, rank: int) -> dict:
    stats = _legacy_stats_from_candidate(candidate)
    return {
        "rank": rank,
        "identity": f"{stats.identity:.6g}",
        "coverage": f"{stats.coverage:.6g}",
        "query_start": stats.query_start,
        "query_end": stats.query_end,
        "target_start": stats.target_start,
        "target_end": stats.target_end,
        "strand": stats.strand,
        "mapping_quality": stats.mapping_quality if stats.mapping_quality is not None else "NA",
        "is_secondary": int(stats.is_secondary or rank > 1),
        "score": f"{stats.score:.6g}",
        "raw_score": f"{stats.raw_score:.6g}" if stats.raw_score is not None else "NA",
        "cigar": stats.cigar,
        "aligned_blocks": list(stats.aligned_blocks),
        "gap_blocks": list(stats.gap_blocks),
        "score_scheme": stats.score_scheme,
        "candidate_id": stats.candidate_id or "NA",
        "query_occurrence_id": stats.query_occurrence_id or "NA",
        "target_occurrence_id": stats.target_occurrence_id or "NA",
        "query_transcript_id": stats.query_transcript_id or "NA",
        "target_transcript_id": stats.target_transcript_id or "NA",
        "sequence_kind": stats.sequence_kind,
        "backend": stats.backend,
        "backend_version": stats.backend_version or "NA",
        "nt_identity": f"{stats.nt_identity:.6g}" if stats.nt_identity is not None else "NA",
        "aa_identity": f"{stats.aa_identity:.6g}" if stats.aa_identity is not None else "NA",
        "known_aligned_pairs": stats.known_aligned_pairs,
        "unknown_aligned_pairs": stats.unknown_aligned_pairs,
        "query_covered_bases": stats.query_covered_bases,
        "target_covered_bases": stats.target_covered_bases,
        "query_length": stats.query_length,
        "target_length": stats.target_length,
        "relative_strand": stats.relative_strand,
        "hit_count": stats.hit_count,
        "alternative_candidate_ids": list(stats.alternative_candidate_ids),
        "left_anchor_id": stats.left_anchor_id or "NA",
        "right_anchor_id": stats.right_anchor_id or "NA",
        "search_interval": stats.search_interval or "NA",
        "enumeration_complete": stats.enumeration_complete,
        "incomplete_reason": stats.incomplete_reason or "NA",
    }
