"""evidence / deletion support: explicit implementation ownership."""
from __future__ import annotations

from intraphy.aligners.pairwise import local_alignment_stats
from intraphy.aligners.types import AlignmentBackendError
from intraphy.evidence.alignment_provenance import _ambiguous_repeated_mapping
from intraphy.evidence.projection import _oriented_locus_slice
from intraphy.evidence.projection import _relative_interval_within_span
from intraphy.evidence.search_context import _occurrence_for_element
from intraphy.evidence.search_context import _occurrences_for_element
from intraphy.evidence.search_context import _ordered_occurrence_pair
from intraphy.evidence.search_context import _ordered_occurrence_triplet
from intraphy.storage.values import to_float


def _cigar_ops(cigar):
    import re

    if not cigar or cigar == "NA":
        return []
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def _supported_alignment_blocks(alignment):
    return [
        {
            "query_start": int(query_start),
            "query_end": int(query_end),
            "target_start": int(target_start),
            "target_end": int(target_end),
        }
        for query_start, query_end, target_start, target_end in alignment.aligned_blocks
    ]


def _query_only_gaps(cigar, query_start=1, target_start=1):
    query_pos = int(query_start or 1)
    target_pos = int(target_start or 1)
    query_only_gaps = []
    for length, op in _cigar_ops(cigar):
        if op in {"M", "=", "X"}:
            query_pos += length
            target_pos += length
        elif op == "I":
            query_only_gaps.append(
                {
                    "query_start": query_pos,
                    "query_end": query_pos + length - 1,
                    "previous_query_end": query_pos - 1,
                    "next_query_start": query_pos + length,
                    "previous_target_end": target_pos - 1,
                    "next_target_start": target_pos,
                }
            )
            query_pos += length
        elif op in {"D", "N"}:
            target_pos += length
        elif op in {"S", "H", "P"}:
            if op == "S":
                query_pos += length
    return query_only_gaps


def _alternative_supports_expected_deletion(hit, expected_interval):
    """Return True/False for verifiable CIGAR evidence, or None if unavailable."""
    cigar = hit.get("cigar")
    raw_blocks = hit.get("aligned_blocks")
    if not cigar or cigar == "NA" or not raw_blocks or hit.get("strand") not in {"+", "-"}:
        return None
    if hit["strand"] != "+":
        return False
    try:
        query_start = int(hit["query_start"])
        target_start = int(hit["target_start"])
        blocks = [
            {
                "query_start": int(block[0]),
                "query_end": int(block[1]),
                "target_start": int(block[2]),
                "target_end": int(block[3]),
            }
            for block in raw_blocks
        ]
    except (KeyError, TypeError, ValueError, IndexError):
        return None
    for deletion in _query_only_gaps(cigar, query_start, target_start):
        if not (
            deletion["query_start"] <= expected_interval[0]
            and deletion["query_end"] >= expected_interval[1]
        ):
            continue
        previous = next(
            (
                block for block in blocks
                if block["query_end"] == deletion["previous_query_end"]
                and block["target_end"] == deletion["previous_target_end"]
            ),
            None,
        )
        following = next(
            (
                block for block in blocks
                if block["query_start"] == deletion["next_query_start"]
                and block["target_start"] == deletion["next_target_start"]
            ),
            None,
        )
        if (
            previous is not None
            and following is not None
            and following["target_start"] == previous["target_end"] + 1
        ):
            return True
    return False


def _projected_interval_coverage(blocks, query_interval, target_interval):
    q_start, q_end = query_interval
    t_start, t_end = target_interval
    covered = 0
    for block in blocks:
        left = max(q_start, block["query_start"])
        right = min(q_end, block["query_end"])
        if left > right:
            continue
        offset = left - block["query_start"]
        projected_left = block["target_start"] + offset
        projected_right = projected_left + (right - left)
        overlap_left = max(min(projected_left, projected_right), t_start)
        overlap_right = min(max(projected_left, projected_right), t_end)
        if overlap_left <= overlap_right:
            covered += overlap_right - overlap_left + 1
    return covered / max(1, q_end - q_start + 1)


def _flank_support_from_spanning_blocks(source_span, target_span, blocks, source_interval, target_interval):
    source_start, source_end = int(source_interval[0]), int(source_interval[1])
    target_start, target_end = int(target_interval[0]), int(target_interval[1])
    source_low, source_high = min(source_start, source_end), max(source_start, source_end)
    target_low, target_high = min(target_start, target_end), max(target_start, target_end)
    matches = mismatches = paired = 0
    for block in blocks:
        q_left = max(source_low, int(block["query_start"]))
        q_right = min(source_high, int(block["query_end"]))
        if q_left > q_right:
            continue
        offset = q_left - int(block["query_start"])
        t_left = int(block["target_start"]) + offset
        for q_pos in range(q_left, q_right + 1):
            t_pos = t_left + (q_pos - q_left)
            if t_pos < target_low or t_pos > target_high:
                continue
            q_base = source_span[q_pos - 1].upper()
            t_base = target_span[t_pos - 1].upper()
            if q_base not in {"A", "C", "G", "T"} or t_base not in {"A", "C", "G", "T"}:
                continue
            paired += 1
            if q_base == t_base:
                matches += 1
            else:
                mismatches += 1
    identity = matches / max(1, matches + mismatches)
    coverage = paired / max(1, source_high - source_low + 1)
    return {
        "identity": identity,
        "coverage": coverage,
        "matches": matches,
        "mismatches": mismatches,
        "paired_bases": paired,
    }


def _ordered_anchor_deletion_support(
    source_occurrences,
    target_occurrences,
    element_by_occurrence,
    expected_occurrence,
    left_element,
    right_element,
    source_locus_header,
    source_locus_sequence,
    target_locus_header,
    target_locus_sequence,
    min_identity,
    min_coverage,
    threads=1,
):
    provenance = {
        "backend": "NA",
        "cigar": "NA",
        "left_flank_backend": "NA",
        "right_flank_backend": "NA",
    }
    if not left_element or not right_element:
        return False, "missing_flanking_homolog_anchor", provenance
    if not source_locus_header or not source_locus_sequence:
        return False, "source_locus_sequence_unavailable", provenance
    if not target_locus_header or not target_locus_sequence:
        return False, "target_locus_sequence_unavailable", provenance
    source_left = _occurrence_for_element(source_occurrences, element_by_occurrence, left_element)
    source_right = _occurrence_for_element(source_occurrences, element_by_occurrence, right_element)
    target_left_matches = _occurrences_for_element(target_occurrences, element_by_occurrence, left_element)
    target_right_matches = _occurrences_for_element(target_occurrences, element_by_occurrence, right_element)
    if not source_left or not source_right or not target_left_matches or not target_right_matches:
        return False, "flanking_homolog_anchor_not_observed", provenance
    if len(target_left_matches) != 1 or len(target_right_matches) != 1:
        return False, "flanking_homolog_anchor_mapping_ambiguous", provenance
    target_left = target_left_matches[0]
    target_right = target_right_matches[0]
    if not _ordered_occurrence_triplet(source_left, expected_occurrence, source_right):
        return False, "source_expected_exon_not_between_ordered_flanks", provenance
    if not _ordered_occurrence_pair(target_left, target_right):
        return False, "target_flanking_homolog_anchors_not_ordered", provenance
    try:
        source_span_start = min(int(source_left["start"]), int(source_right["start"]))
        source_span_end = max(int(source_left["end"]), int(source_right["end"]))
        target_span_start = min(int(target_left["start"]), int(target_right["start"]))
        target_span_end = max(int(target_left["end"]), int(target_right["end"]))
    except (KeyError, TypeError, ValueError):
        return False, "flanking_homolog_anchor_coordinates_unusable", provenance
    source_span, source_status = _oriented_locus_slice(
        source_locus_header, source_locus_sequence, source_span_start, source_span_end
    )
    target_span, target_status = _oriented_locus_slice(
        target_locus_header, target_locus_sequence, target_span_start, target_span_end
    )
    if source_status != "ok":
        return False, f"source_{source_status}", provenance
    if target_status != "ok":
        return False, f"target_{target_status}", provenance
    expected_interval = _relative_interval_within_span(
        source_locus_header,
        source_span_start,
        source_span_end,
        expected_occurrence["start"],
        expected_occurrence["end"],
    )
    source_left_interval = _relative_interval_within_span(
        source_locus_header, source_span_start, source_span_end, source_left["start"], source_left["end"]
    )
    source_right_interval = _relative_interval_within_span(
        source_locus_header, source_span_start, source_span_end, source_right["start"], source_right["end"]
    )
    target_left_interval = _relative_interval_within_span(
        target_locus_header, target_span_start, target_span_end, target_left["start"], target_left["end"]
    )
    target_right_interval = _relative_interval_within_span(
        target_locus_header, target_span_start, target_span_end, target_right["start"], target_right["end"]
    )
    local_target_start = max(1, min(target_left_interval[1], target_right_interval[0]) - 10)
    local_target_end = min(len(target_span), max(target_left_interval[1], target_right_interval[0]) + 10)
    expected_sequence = source_span[expected_interval[0] - 1 : expected_interval[1]]
    local_target_sequence = target_span[local_target_start - 1 : local_target_end]
    target_gap_left = min(target_left_interval[1], target_right_interval[1]) + 1
    target_gap_right = max(target_left_interval[0], target_right_interval[0]) - 1
    if any(base not in "ACGT" for base in expected_sequence.upper()):
        return False, "deletion_interval_contains_non_acgt_bases", provenance
    if target_gap_left <= target_gap_right and any(
        base not in "ACGT"
        for base in target_span[target_gap_left - 1 : target_gap_right].upper()
    ):
        return False, "target_sequence_between_flanking_homologs_contains_non_acgt_bases", provenance
    if any(base not in "ACGT" for base in local_target_sequence.upper()):
        return False, "local_target_sequence_contains_non_acgt_bases", provenance
    try:
        alignment = local_alignment_stats(source_span, target_span, backend="minimap2", threads=threads)
    except AlignmentBackendError as exc:
        return False, f"deletion_spanning_alignment_unresolved:{exc}", provenance
    provenance.update({"backend": alignment.backend, "cigar": alignment.cigar})
    if alignment.strand != "+":
        return False, "deletion_spanning_alignment_not_relative_plus_strand", provenance
    acceptable_alternatives = []
    for hit in alignment.alternative_hits:
        identity = to_float(hit.get("identity"), 0.0)
        coverage = to_float(hit.get("coverage"), 0.0)
        if identity >= min_identity and coverage >= min_coverage:
            acceptable_alternatives.append(hit)
    alternative_support = [
        _alternative_supports_expected_deletion(hit, expected_interval)
        for hit in acceptable_alternatives
    ]
    if any(support is None for support in alternative_support):
        return False, "alternative_alignment_concordance_unverifiable", provenance
    if any(not support for support in alternative_support):
        return False, "alternative_alignment_does_not_support_same_deletion", provenance
    if _ambiguous_repeated_mapping(alignment) and not acceptable_alternatives:
        return False, "deletion_spanning_alignment_ambiguous_repeated_mapping", provenance
    blocks = _supported_alignment_blocks(alignment)
    left_support = _flank_support_from_spanning_blocks(
        source_span, target_span, blocks, source_left_interval, target_left_interval
    )
    right_support = _flank_support_from_spanning_blocks(
        source_span, target_span, blocks, source_right_interval, target_right_interval
    )
    provenance.update(
        {
            "left_flank_backend": alignment.backend,
            "left_flank_cigar": alignment.cigar,
            "left_flank_identity": f"{left_support['identity']:.6g}",
            "left_flank_coverage": f"{left_support['coverage']:.6g}",
            "left_flank_paired_bases": left_support["paired_bases"],
            "left_flank_matches": left_support["matches"],
            "left_flank_mismatches": left_support["mismatches"],
            "right_flank_backend": alignment.backend,
            "right_flank_cigar": alignment.cigar,
            "right_flank_identity": f"{right_support['identity']:.6g}",
            "right_flank_coverage": f"{right_support['coverage']:.6g}",
            "right_flank_paired_bases": right_support["paired_bases"],
            "right_flank_matches": right_support["matches"],
            "right_flank_mismatches": right_support["mismatches"],
        }
    )
    if left_support["identity"] < min_identity:
        return False, "left_flanking_homolog_identity_insufficient_in_spanning_alignment", provenance
    if left_support["coverage"] < min_coverage:
        return False, "left_flanking_homolog_coverage_insufficient_in_spanning_alignment", provenance
    if right_support["identity"] < min_identity:
        return False, "right_flanking_homolog_identity_insufficient_in_spanning_alignment", provenance
    if right_support["coverage"] < min_coverage:
        return False, "right_flanking_homolog_coverage_insufficient_in_spanning_alignment", provenance
    deletions = _query_only_gaps(
        alignment.cigar,
        query_start=alignment.query_start,
        target_start=alignment.target_start,
    )
    left_coverage = _projected_interval_coverage(blocks, source_left_interval, target_left_interval)
    right_coverage = _projected_interval_coverage(blocks, source_right_interval, target_right_interval)
    if left_coverage < min_coverage or right_coverage < min_coverage:
        return False, "both_flanking_homologs_not_projected_in_spanning_alignment", provenance
    expected_start, expected_end = expected_interval
    for deletion in deletions:
        covers_expected = deletion["query_start"] <= expected_start and deletion["query_end"] >= expected_end
        previous_block = next(
            (
                block for block in blocks
                if block["query_end"] == deletion["previous_query_end"]
                and block["target_end"] == deletion["previous_target_end"]
            ),
            None,
        )
        next_block = next(
            (
                block for block in blocks
                if block["query_start"] == deletion["next_query_start"]
                and block["target_start"] == deletion["next_target_start"]
            ),
            None,
        )
        continuous_target = (
            previous_block is not None
            and next_block is not None
            and next_block["target_start"] == previous_block["target_end"] + 1
        )
        if covers_expected and continuous_target:
            return True, "source_expected_exon_deleted_in_target_spanning_alignment", provenance
    return False, "no_query_only_gap_covers_expected_exon", provenance
