"""mapping / alternative: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from insiphy.mapping.fields import EXON_LIKE_ROLES
from insiphy.mapping.fields import _transcript_id_set
from insiphy.mapping.matches import context_score
from insiphy.mapping.matches import phase_score
from insiphy.mapping.matches import role_boundary_score
from insiphy.mapping.matches import segment_length
from insiphy.mapping.policies import occurrence_copy_key
from insiphy.storage.values import to_float


def _relative_overlap_interval(row, overlap_start, overlap_end):
    start = int(row["start"])
    end = int(row["end"])
    if row.get("strand") == "-":
        rel_start = end - int(overlap_end) + 1
        rel_end = end - int(overlap_start) + 1
    else:
        rel_start = int(overlap_start) - start + 1
        rel_end = int(overlap_end) - start + 1
    return min(rel_start, rel_end), max(rel_start, rel_end)


def alternative_overlap_evidence(left, right, context):
    if occurrence_copy_key(left) != occurrence_copy_key(right):
        return None
    if left.get("role") not in EXON_LIKE_ROLES or right.get("role") not in EXON_LIKE_ROLES:
        return None
    if left.get("contig") != right.get("contig") or left.get("strand") != right.get("strand"):
        return None
    left_tx = _transcript_id_set(left)
    right_tx = _transcript_id_set(right)
    if not left_tx or not right_tx or left_tx & right_tx:
        return None
    try:
        overlap_start = max(int(left["start"]), int(right["start"]))
        overlap_end = min(int(left["end"]), int(right["end"]))
    except (KeyError, TypeError, ValueError):
        return None
    if overlap_start > overlap_end:
        return None
    query_start, query_end = _relative_overlap_interval(left, overlap_start, overlap_end)
    target_start, target_end = _relative_overlap_interval(right, overlap_start, overlap_end)
    overlap_len = overlap_end - overlap_start + 1
    left_len = segment_length(left)
    right_len = segment_length(right)
    coverage = overlap_len / max(1, min(left_len, right_len))
    left_ctx = context.get(left["occurrence_id"], {})
    right_ctx = context.get(right["occurrence_id"], {})
    order = 1.0 - abs(to_float(left_ctx.get("scaled_index"), 0.5) - to_float(right_ctx.get("scaled_index"), 0.5))
    left_context = context_score(left_ctx, right_ctx, "left")
    right_context = context_score(left_ctx, right_ctx, "right")
    boundary = role_boundary_score(left, right)
    phase = phase_score(left, right)
    splice = 1.0 - abs(to_float(left.get("splice_motif_score"), 0.5) - to_float(right.get("splice_motif_score"), 0.5))
    total = (
        0.34 * 1.0
        + 0.14 * coverage
        + 0.10 * left_context
        + 0.10 * right_context
        + 0.10 * boundary
        + 0.08 * phase
        + 0.06 * order
        + 0.04 * 1.0
        + 0.04 * splice
    )
    return {
        "alignment_score": 1.0,
        "coverage_score": coverage,
        "sequence_score": 0.70 * 1.0 + 0.30 * coverage,
        "structural_context_score": 0.5 * left_context + 0.5 * right_context,
        "query_coverage": overlap_len / max(1, left_len),
        "target_coverage": overlap_len / max(1, right_len),
        "aligned_pairs": overlap_len,
        "query_alignment_start": query_start,
        "query_alignment_end": query_end,
        "target_alignment_start": target_start,
        "target_alignment_end": target_end,
        "query_mapped_contig": left.get("contig", "NA"),
        "query_mapped_start": overlap_start,
        "query_mapped_end": overlap_end,
        "query_mapped_strand": left.get("strand", "NA"),
        "query_mapped_length": overlap_len,
        "subject_mapped_contig": right.get("contig", "NA"),
        "subject_mapped_start": overlap_start,
        "subject_mapped_end": overlap_end,
        "subject_mapped_strand": right.get("strand", "NA"),
        "subject_mapped_length": overlap_len,
        "alignment_strand": "+",
        "projected_reference_occurrence_id": right["occurrence_id"],
        "projected_reference_start": target_start,
        "projected_reference_end": target_end,
        "projected_reference_blocks": f"{query_start}-{query_end}:{target_start}-{target_end}",
        "matched_blocks": f"{query_start}-{query_end}:{target_start}-{target_end}",
        "query_genomic_matched_blocks": f"{left.get('contig', 'NA')}:{overlap_start}-{overlap_end}:{left.get('strand', 'NA')}",
        "subject_genomic_matched_blocks": f"{right.get('contig', 'NA')}:{overlap_start}-{overlap_end}:{right.get('strand', 'NA')}",
        "query_parent_feature_ids": left.get("source_feature_id", "NA"),
        "subject_parent_feature_ids": right.get("source_feature_id", "NA"),
        "query_transcript_ids": left.get("transcript_id", "NA"),
        "subject_transcript_ids": right.get("transcript_id", "NA"),
        "left_context_score": left_context,
        "right_context_score": right_context,
        "boundary_score": boundary,
        "phase_score": phase,
        "order_score": order,
        "strand_score": 1.0,
        "splice_score": splice,
        "size_ratio": min(left_len, right_len) / max(left_len, right_len),
        "alignment_cigar": f"{overlap_len}=",
        "alignment_backend": "genomic_overlap",
        "alignment_mode": "coordinate_overlap",
        "alignment_meaning": "same-copy alternative isoform shared genomic interval",
        "alignment_requested_backend": "genomic_overlap",
        "mapping_quality": "NA",
        "hit_count": 1,
        "ambiguous_hit_count": 0,
        "alternative_hits": [],
        "score_scheme": "genomic_coordinate_overlap",
        "raw_alignment_score": overlap_len,
        "enumeration_complete": True,
        "candidate_enumeration_status": "complete",
        "incomplete_reason": "NA",
        "short_context_route": "not_used",
        "local_boundary_range": "available",
        "flanking_anchor_status": "same_locus_annotation_overlap",
        "true_absence_eligible": 0,
        "true_absence_evidence_status": "not_applicable",
        "true_absence_reason": "same_locus_annotation_overlap_is_not_deletion_evidence",
        "total_score": total,
    }
