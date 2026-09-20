"""evidence / search_context: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.aligners.columns import revcomp
from insiphy.elements import EXON_LIKE_ROLES
from insiphy.evidence.projection import _locus_geometry
from insiphy.evidence.projection import _oriented_locus_slice
from insiphy.evidence.projection import _relative_interval_within_span
from insiphy.evidence.projection import _split_transcript_ids
from insiphy.storage.fasta import fasta_record_length
from insiphy.storage.fasta import read_fasta_interval
from insiphy.storage.tabular import read_tsv
from pathlib import Path


def _read_locus_metadata(input_dir):
    rows = read_tsv(Path(input_dir) / "gene_loci.tsv", optional=True)
    return {(row.get("species"), row.get("gene_copy_id")): row for row in rows}


def _search_provenance(metadata, locus_header, terminal_need, candidates):
    try:
        _contig, bound_start, bound_end, _strand = _locus_geometry(locus_header)
    except (IndexError, TypeError, ValueError):
        bound_start = bound_end = None
    annotation_start = _safe_int(metadata.get("annotation_start"), 0)
    annotation_end = _safe_int(metadata.get("annotation_end"), 0)
    linked_start = _safe_int(metadata.get("linked_start"), 0)
    linked_end = _safe_int(metadata.get("linked_end"), 0)
    contig_length = _safe_int(metadata.get("contig_length"), 0)
    max_extension = _safe_int(metadata.get("max_extension"), 0)

    left_limit = "unassessed"
    right_limit = "unassessed"
    if bound_start is not None:
        if bound_start == 1:
            left_limit = "contig_boundary_reached"
        elif annotation_start and max_extension and bound_start <= annotation_start - max_extension:
            left_limit = "max_extension_reached"
        else:
            left_limit = "within_declared_search_limit"
        if contig_length and bound_end == contig_length:
            right_limit = "contig_boundary_reached"
        elif annotation_end and max_extension and bound_end >= annotation_end + max_extension:
            right_limit = "max_extension_reached"
        else:
            right_limit = "within_declared_search_limit"

    touches_left = False
    touches_right = False
    for candidate in candidates:
        try:
            touches_left = touches_left or int(candidate["start"]) <= int(bound_start)
            touches_right = touches_right or int(candidate["end"]) >= int(bound_end)
        except (KeyError, TypeError, ValueError):
            continue
    if touches_left and touches_right:
        hit_limit = "candidate_touches_both_search_bounds"
    elif touches_left:
        hit_limit = "candidate_touches_left_search_bound"
    elif touches_right:
        hit_limit = "candidate_touches_right_search_bound"
    elif candidates:
        hit_limit = "candidate_inside_search_bounds"
    elif terminal_need == "left" and left_limit in {
        "contig_boundary_reached", "max_extension_reached"
    }:
        hit_limit = "no_hit_with_left_search_limit_reached"
    elif terminal_need == "right" and right_limit in {
        "contig_boundary_reached", "max_extension_reached"
    }:
        hit_limit = "no_hit_with_right_search_limit_reached"
    else:
        hit_limit = "no_candidate_reported"
    return {
        "original_annotation_start": annotation_start or "NA",
        "original_annotation_end": annotation_end or "NA",
        "linked_feature_start": linked_start or "NA",
        "linked_feature_end": linked_end or "NA",
        "search_bound_start": bound_start if bound_start is not None else "NA",
        "search_bound_end": bound_end if bound_end is not None else "NA",
        "search_left_limit_status": left_limit,
        "search_right_limit_status": right_limit,
        "hit_search_limit_status": hit_limit,
    }


def _extended_locus(input_dir, locus_header, locus, metadata, need_side):
    if need_side not in {"left", "right"} or not metadata:
        return locus_header, locus, {"search_expansion_status": "not_requested"}
    genome = metadata.get("genome_fasta")
    contig = metadata.get("contig")
    strand = metadata.get("strand", "+")
    try:
        max_extension = int(metadata.get("max_extension", 0) or 0)
        annotation_start = int(metadata.get("annotation_start"))
        annotation_end = int(metadata.get("annotation_end"))
        search_start = int(metadata.get("search_start") or metadata.get("annotation_start"))
        search_end = int(metadata.get("search_end") or metadata.get("annotation_end"))
        contig_length = int(metadata.get("contig_length", 0) or 0)
    except (TypeError, ValueError):
        return locus_header, locus, {"search_expansion_status": "metadata_unusable"}
    if not genome or not contig or max_extension <= 0:
        return locus_header, locus, {"search_expansion_status": "metadata_unusable"}
    if contig_length <= 0:
        contig_length = int(fasta_record_length(genome, contig))
    if contig_length <= 0:
        return locus_header, locus, {"search_expansion_status": "contig_length_unknown"}

    new_start, new_end = search_start, search_end
    if (need_side == "left" and strand != "-") or (need_side == "right" and strand == "-"):
        new_start = min(search_start, max(1, annotation_start - max_extension))
    else:
        new_end = max(search_end, min(contig_length, annotation_end + max_extension))
    if (new_start, new_end) == (search_start, search_end):
        boundary = (
            new_start == 1
            if (need_side == "left" and strand != "-") or (need_side == "right" and strand == "-")
            else new_end == contig_length
        )
        return locus_header, locus, {
            "search_expansion_status": (
                "at_contig_boundary" if boundary else "at_declared_search_limit"
            )
        }
    sequence = read_fasta_interval(genome, contig, new_start, new_end)
    if strand == "-":
        sequence = revcomp(sequence)
    parts = locus_header.split("|", 2)
    header = f"{parts[0]}|{parts[1]}|{contig}:{new_start}-{new_end}:{strand}" if len(parts) >= 2 else locus_header
    return header, sequence, {
        "search_expansion_status": "expanded",
        "search_original_start": search_start,
        "search_original_end": search_end,
        "search_expanded_start": new_start,
        "search_expanded_end": new_end,
        "range_status": metadata.get("range_status", "NA"),
    }


def _occurrence_for_element(copy_occurrences, element_by_occurrence, element):
    for occurrence in copy_occurrences:
        if (
            element_by_occurrence.get(occurrence.get("occurrence_id")) == element
            and occurrence.get("presence_status") == "present"
        ):
            return occurrence
    return None


def _occurrences_for_element(copy_occurrences, element_by_occurrence, element):
    return [
        occurrence for occurrence in copy_occurrences
        if element_by_occurrence.get(occurrence.get("occurrence_id")) == element
        and occurrence.get("presence_status") == "present"
    ]


def _safe_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _path_occurrences(path_rows, occ_by_id, element_by_occurrence):
    ordered = []
    seen = set()
    for row in sorted(path_rows, key=lambda item: (_safe_int(item.get("path_rank")), item.get("occurrence_id", ""))):
        occurrence_id = row.get("occurrence_id")
        occurrence = occ_by_id.get(occurrence_id)
        if (
            not occurrence
            or occurrence_id in seen
            or occurrence.get("presence_status") != "present"
            or occurrence.get("role") not in EXON_LIKE_ROLES
            or occurrence_id not in element_by_occurrence
        ):
            continue
        ordered.append(occurrence)
        seen.add(occurrence_id)
    return ordered


def _disjoint_genomic_occurrences(occurrences, core_exonic_occurrences):
    candidates = [
        row for row in occurrences
        if row.get("occurrence_id") in core_exonic_occurrences
        and row.get("presence_status") == "present"
        and row.get("role") in EXON_LIKE_ROLES
    ]
    if not candidates:
        return [], "source_path_unknown"
    contigs = {row.get("contig") for row in candidates}
    strands = {row.get("strand") for row in candidates}
    if len(contigs) != 1 or len(strands) != 1:
        return [], "source_path_ambiguous_genomic_order"
    intervals = []
    for row in candidates:
        try:
            intervals.append((int(row["start"]), int(row["end"]), row))
        except (KeyError, TypeError, ValueError):
            return [], "source_path_ambiguous_genomic_order"
    intervals.sort(key=lambda item: (item[0], item[1], item[2].get("occurrence_id", "")))
    for (_left_start, left_end, _left), (right_start, _right_end, _right) in zip(intervals, intervals[1:]):
        if left_end >= right_start:
            return [], "source_path_ambiguous_genomic_order"
    if next(iter(strands)) == "-":
        intervals.reverse()
    return [row for _start, _end, row in intervals], "source_path_genomic_order"


def _source_occurrence_path(
    family,
    representative,
    transcript_paths,
    occurrences_by_copy,
    occ_by_id,
    element_by_occurrence,
    core_exonic_occurrences,
):
    copy_key = (family, representative.get("species"), representative.get("gene_copy_id"))
    representative_id = representative.get("occurrence_id")
    copy_path_rows = [
        row for row in transcript_paths
        if row.get("family_id") == family
        and row.get("species") == representative.get("species")
        and row.get("gene_copy_id") == representative.get("gene_copy_id")
    ]
    rows_by_transcript = defaultdict(list)
    for row in copy_path_rows:
        rows_by_transcript[row.get("transcript_id", "")].append(row)
    transcript_candidates = []
    representative_transcripts = set(_split_transcript_ids(representative.get("transcript_id")))
    for transcript_id, rows in rows_by_transcript.items():
        if any(row.get("occurrence_id") == representative_id for row in rows):
            preferred = 0 if transcript_id in representative_transcripts else 1
            rank = min(_safe_int(row.get("path_rank")) for row in rows)
            transcript_candidates.append((preferred, transcript_id, rank, rows))
    if transcript_candidates:
        _preferred, transcript_id, _rank, rows = min(transcript_candidates)
        ordered = _path_occurrences(rows, occ_by_id, element_by_occurrence)
        if any(row.get("occurrence_id") == representative_id for row in ordered):
            return ordered, f"source_path_transcript:{transcript_id}"
    ordered, status = _disjoint_genomic_occurrences(occurrences_by_copy[copy_key], core_exonic_occurrences)
    if any(row.get("occurrence_id") == representative_id for row in ordered):
        return ordered, status
    return [], status


def _ordered_occurrence_triplet(left, middle, right):
    if not left or not middle or not right:
        return False
    if len({left.get("contig"), middle.get("contig"), right.get("contig")}) != 1:
        return False
    if len({left.get("strand"), middle.get("strand"), right.get("strand")}) != 1:
        return False
    try:
        left_start, left_end = int(left["start"]), int(left["end"])
        middle_start, middle_end = int(middle["start"]), int(middle["end"])
        right_start, right_end = int(right["start"]), int(right["end"])
    except (KeyError, TypeError, ValueError):
        return False
    if left.get("strand") == "-":
        return right_end < middle_start and middle_end < left_start
    return left_end < middle_start and middle_end < right_start


def _ordered_occurrence_pair(left, right):
    if not left or not right:
        return False
    if left.get("contig") != right.get("contig") or left.get("strand") != right.get("strand"):
        return False
    try:
        left_start, left_end = int(left["start"]), int(left["end"])
        right_start, right_end = int(right["start"]), int(right["end"])
    except (KeyError, TypeError, ValueError):
        return False
    if left.get("strand") == "-":
        return right_end < left_start
    return left_end < right_start


def _ordered_anchor_search_context(
    source_occurrences,
    target_occurrences,
    element_by_occurrence,
    expected_occurrence,
    left_element,
    right_element,
    target_locus_header,
    target_locus_sequence,
):
    context = {
        "valid_double_flank": False,
        "status": "missing_flanking_homolog_anchor",
        "left_anchor_id": "NA",
        "right_anchor_id": "NA",
        "search_interval": None,
        "sequence": "",
        "target_offset0": 0,
    }
    if not left_element or not right_element:
        return context
    source_left_matches = _occurrences_for_element(
        source_occurrences, element_by_occurrence, left_element
    )
    source_right_matches = _occurrences_for_element(
        source_occurrences, element_by_occurrence, right_element
    )
    target_left_matches = _occurrences_for_element(target_occurrences, element_by_occurrence, left_element)
    target_right_matches = _occurrences_for_element(target_occurrences, element_by_occurrence, right_element)
    if not source_left_matches or not source_right_matches or not target_left_matches or not target_right_matches:
        context["status"] = "flanking_homolog_anchor_not_observed"
        return context
    if len(source_left_matches) != 1 or len(source_right_matches) != 1:
        context["status"] = "source_flanking_homolog_anchor_mapping_ambiguous"
        return context
    if len(target_left_matches) != 1 or len(target_right_matches) != 1:
        context["status"] = "flanking_homolog_anchor_mapping_ambiguous"
        return context
    source_left = source_left_matches[0]
    source_right = source_right_matches[0]
    target_left = target_left_matches[0]
    target_right = target_right_matches[0]
    context["left_anchor_id"] = target_left.get("occurrence_id", "NA")
    context["right_anchor_id"] = target_right.get("occurrence_id", "NA")
    if not _ordered_occurrence_triplet(source_left, expected_occurrence, source_right):
        context["status"] = "source_expected_exon_not_between_ordered_flanks"
        return context
    if not _ordered_occurrence_pair(target_left, target_right):
        context["status"] = "target_flanking_homolog_anchors_not_ordered"
        return context
    try:
        locus_contig, locus_start, locus_end, locus_strand = _locus_geometry(target_locus_header)
        if target_left.get("contig") != locus_contig or target_left.get("strand") != locus_strand:
            context["status"] = "target_flanking_homolog_anchors_incompatible_with_locus"
            return context
        if locus_strand == "-":
            interval_start = int(target_right["end"]) + 1
            interval_end = int(target_left["start"]) - 1
        else:
            interval_start = int(target_left["end"]) + 1
            interval_end = int(target_right["start"]) - 1
    except (KeyError, TypeError, ValueError):
        context["status"] = "flanking_homolog_anchor_coordinates_unusable"
        return context

    search_interval = {
        "contig": locus_contig,
        "start": interval_start,
        "end": interval_end,
        "strand": locus_strand,
        "coordinate_system": "1-based-closed",
    }
    context["search_interval"] = search_interval
    context["valid_double_flank"] = True
    if interval_start > interval_end:
        context["status"] = "ordered_double_flank_empty_interval"
        return context
    sequence, interval_status = _oriented_locus_slice(
        target_locus_header,
        target_locus_sequence,
        interval_start,
        interval_end,
    )
    if interval_status != "ok":
        context["valid_double_flank"] = False
        context["status"] = f"target_{interval_status}"
        return context
    local_interval = _relative_interval_within_span(
        target_locus_header,
        locus_start,
        locus_end,
        interval_start,
        interval_end,
    )
    context.update(
        {
            "status": "ordered_double_flank_bounded_interval",
            "sequence": sequence,
            "target_offset0": local_interval[0] - 1,
        }
    )
    return context
