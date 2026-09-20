"""evidence / alignment provenance: explicit implementation ownership."""
from __future__ import annotations

from intraphy.evidence.projection import _project_locus_interval
import json


def _offset_alignment_target(alignment, target_offset0):
    if alignment is None or not target_offset0:
        return alignment
    if alignment.target_start > 0:
        alignment.target_start += target_offset0
    if alignment.target_end > 0:
        alignment.target_end += target_offset0
    alignment.aligned_blocks = [
        (query_start, query_end, target_start + target_offset0, target_end + target_offset0)
        for query_start, query_end, target_start, target_end in alignment.aligned_blocks
    ]
    for gap in alignment.gap_blocks:
        if gap.get("target_start0") is not None:
            gap["target_start0"] = int(gap["target_start0"]) + target_offset0
        if gap.get("target_end0") is not None:
            gap["target_end0"] = int(gap["target_end0"]) + target_offset0
    for hit in alignment.alternative_hits:
        if hit.get("target_start") not in {None, "", "NA"}:
            hit["target_start"] = int(hit["target_start"]) + target_offset0
        if hit.get("target_end") not in {None, "", "NA"}:
            hit["target_end"] = int(hit["target_end"]) + target_offset0
        if hit.get("aligned_blocks"):
            hit["aligned_blocks"] = [
                (int(block[0]), int(block[1]), int(block[2]) + target_offset0, int(block[3]) + target_offset0)
                for block in hit["aligned_blocks"]
            ]
        for gap in hit.get("gap_blocks", ()):
            if gap.get("target_start0") is not None:
                gap["target_start0"] = int(gap["target_start0"]) + target_offset0
            if gap.get("target_end0") is not None:
                gap["target_end0"] = int(gap["target_end0"]) + target_offset0
    return alignment


def _set_alignment_context(
    alignment,
    query_occurrence_id,
    target_occurrence_id,
    query_transcript_id,
    left_anchor_id,
    right_anchor_id,
    search_interval,
):
    if alignment is None:
        return None
    alignment.query_occurrence_id = query_occurrence_id
    alignment.target_occurrence_id = target_occurrence_id
    alignment.query_transcript_id = query_transcript_id
    alignment.left_anchor_id = left_anchor_id
    alignment.right_anchor_id = right_anchor_id
    alignment.search_interval = search_interval
    query_candidate_scope = query_occurrence_id or "query"
    target_candidate_scope = target_occurrence_id or "target_interval"
    candidate_ids = [
        alignment.candidate_id or f"{query_candidate_scope}->{target_candidate_scope}.candidate_001"
    ]
    for index, hit in enumerate(alignment.alternative_hits, start=2):
        candidate_ids.append(
            hit.get("candidate_id") or f"{query_candidate_scope}->{target_candidate_scope}.candidate_{index:03d}"
        )
    alignment.candidate_id = candidate_ids[0]
    alignment.hit_count = max(int(alignment.hit_count or 0), len(candidate_ids))
    alignment.alternative_candidate_ids = candidate_ids[1:]
    for index, hit in enumerate(alignment.alternative_hits, start=1):
        hit["candidate_id"] = candidate_ids[index]
        hit.setdefault("query_occurrence_id", query_occurrence_id)
        hit.setdefault("target_occurrence_id", target_occurrence_id or "NA")
        hit.setdefault("query_transcript_id", query_transcript_id or "NA")
        hit.setdefault("target_transcript_id", "NA")
        hit.setdefault("left_anchor_id", left_anchor_id)
        hit.setdefault("right_anchor_id", right_anchor_id)
        hit.setdefault("search_interval", search_interval)
        hit.setdefault("alternative_candidate_ids", [value for value in candidate_ids if value != hit["candidate_id"]])
    return alignment


def _projection_within_anchor_interval(projection, anchor_context):
    interval = anchor_context.get("search_interval") or {}
    if not projection or not anchor_context.get("valid_double_flank"):
        return False
    try:
        return (
            projection.get("contig") == interval["contig"]
            and projection.get("strand") == interval["strand"]
            and int(projection["start"]) >= int(interval["start"])
            and int(projection["end"]) <= int(interval["end"])
        )
    except (KeyError, TypeError, ValueError):
        return False


def _format_alternative_hits(alignment):
    if alignment is None or not alignment.alternative_hits:
        return "NA"
    tokens = []
    for hit in alignment.alternative_hits:
        tokens.append(
            f"rank{hit.get('rank', 'NA')}:{hit.get('target_start', 'NA')}-{hit.get('target_end', 'NA')}:"
            f"{hit.get('strand', 'NA')}:id{hit.get('identity', 'NA')}:cov{hit.get('coverage', 'NA')}:"
            f"mapq{hit.get('mapping_quality', 'NA')}:secondary{hit.get('is_secondary', 'NA')}:"
            f"candidate{hit.get('candidate_id', 'NA')}"
        )
    return ";".join(tokens)


def _alignment_evidence_provenance(alignment, query_length=0, target_length=0):
    if alignment is None:
        return {
            "alignment_sequence_kind": "unknown",
            "alignment_backend_version": "NA",
            "alignment_score_scheme": "NA",
            "alignment_nt_identity": "NA",
            "alignment_aa_identity": "NA",
            "alignment_known_aligned_pairs": "NA",
            "alignment_unknown_aligned_pairs": "NA",
            "alignment_query_covered_bases": "NA",
            "alignment_target_covered_bases": "NA",
            "alignment_query_length": query_length or "NA",
            "alignment_target_length": target_length or "NA",
            "alignment_gap_blocks": "NA",
            "alignment_relative_strand": "unknown",
            "alignment_candidate_id": "NA",
            "alignment_alternative_candidate_ids": "NA",
            "alignment_enumeration_complete": "unknown",
            "alignment_incomplete_reason": "NA",
        }
    known_pairs = int(getattr(alignment, "known_aligned_pairs", 0) or alignment.aligned_pairs or 0)
    query_covered = int(getattr(alignment, "query_covered_bases", 0) or known_pairs)
    target_covered = int(getattr(alignment, "target_covered_bases", 0) or known_pairs)
    candidate_ids = [
        hit.get("candidate_id")
        for hit in alignment.alternative_hits
        if hit.get("candidate_id") not in {None, "", "NA"}
    ]
    candidate_ids = list(getattr(alignment, "alternative_candidate_ids", ()) or candidate_ids)
    return {
        "alignment_sequence_kind": getattr(alignment, "sequence_kind", "nucleotide"),
        "alignment_backend_version": getattr(alignment, "backend_version", None) or "NA",
        "alignment_score_scheme": getattr(alignment, "score_scheme", None) or "NA",
        "alignment_nt_identity": (
            f"{alignment.nt_identity:.6g}"
            if getattr(alignment, "nt_identity", None) is not None
            else f"{alignment.identity:.6g}"
        ),
        "alignment_aa_identity": (
            f"{alignment.aa_identity:.6g}"
            if getattr(alignment, "aa_identity", None) is not None
            else "NA"
        ),
        "alignment_known_aligned_pairs": known_pairs,
        "alignment_unknown_aligned_pairs": int(getattr(alignment, "unknown_aligned_pairs", 0) or 0),
        "alignment_query_covered_bases": query_covered,
        "alignment_target_covered_bases": target_covered,
        "alignment_query_length": int(getattr(alignment, "query_length", 0) or query_length),
        "alignment_target_length": int(getattr(alignment, "target_length", 0) or target_length),
        "alignment_gap_blocks": json.dumps(
            getattr(alignment, "gap_blocks", ()) or (), separators=(",", ":")
        ),
        "alignment_relative_strand": getattr(alignment, "relative_strand", None) or alignment.strand,
        "alignment_candidate_id": getattr(alignment, "candidate_id", None) or "NA",
        "alignment_alternative_candidate_ids": ";".join(candidate_ids) or "NA",
        "alignment_enumeration_complete": str(bool(alignment.enumeration_complete)).lower(),
        "alignment_incomplete_reason": alignment.incomplete_reason or "NA",
    }


def _ambiguous_repeated_mapping(alignment):
    if alignment is None:
        return False
    if int(alignment.ambiguous_hit_count or 0) > 0:
        return True
    return int(alignment.hit_count or 0) > 1 and (
        alignment.is_secondary
        or (alignment.mapping_quality is not None and int(alignment.mapping_quality) == 0)
    )


def _nucleotide_interval_candidates(alignment, locus_header):
    if (
        alignment is None
        or alignment.target_start <= 0
        or alignment.target_end <= 0
        or (not alignment.aligned_blocks and int(alignment.aligned_pairs or 0) <= 0)
    ):
        return []
    hits = [{
        "target_start": alignment.target_start,
        "target_end": alignment.target_end,
        "strand": alignment.strand,
        "identity": alignment.identity,
        "coverage": alignment.query_coverage or alignment.coverage,
        "mapping_quality": alignment.mapping_quality if alignment.mapping_quality is not None else "NA",
        "is_secondary": alignment.is_secondary,
        "candidate_id": alignment.candidate_id or "NA",
        "left_anchor_id": alignment.left_anchor_id or "NA",
        "right_anchor_id": alignment.right_anchor_id or "NA",
        "search_interval": alignment.search_interval or "NA",
        "enumeration_complete": alignment.enumeration_complete,
        "incomplete_reason": alignment.incomplete_reason or "NA",
    }, *alignment.alternative_hits]
    candidates = []
    for hit in hits:
        contig, start, end, strand = _project_locus_interval(hit["target_start"], hit["target_end"], locus_header)
        candidates.append({
            "contig": contig, "start": start, "end": end,
            "strand": strand if hit["strand"] == "+" else ("-" if strand == "+" else "+"),
            "backend": alignment.backend,
            "interval_scope": "aligned_sequence",
            "identity": hit["identity"], "coverage": hit["coverage"],
            "mapping_quality": hit.get("mapping_quality", "NA"), "is_secondary": hit.get("is_secondary", "NA"),
            "candidate_id": hit.get("candidate_id", "NA"),
            "left_anchor_id": hit.get("left_anchor_id", "NA"),
            "right_anchor_id": hit.get("right_anchor_id", "NA"),
            "search_interval": hit.get("search_interval", "NA"),
            "enumeration_complete": hit.get("enumeration_complete", alignment.enumeration_complete),
            "incomplete_reason": hit.get("incomplete_reason", alignment.incomplete_reason or "NA"),
        })
    return candidates
