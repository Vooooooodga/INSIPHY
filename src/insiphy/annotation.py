"""Sequence-supported annotation completion."""

import json
from collections import Counter, defaultdict
from pathlib import Path

from .alignment import (
    AlignmentBackendError,
    local_alignment_stats,
    protein_locus_exons,
    revcomp,
)
from .coordinates import ClosedInterval1, local_interval_to_genome
from .elements import EXON_LIKE_ROLES
from .io import fasta_record_length, parse_fasta, read_fasta_interval, read_tsv, to_float, write_tsv


_PROTEIN_QUERY_INTERVAL_FIELDS = (
    "protein_cds_query_start", "protein_cds_query_end",
    "protein_overlap_query_start", "protein_overlap_query_end",
)

_EVIDENCE_PROVENANCE_FIELDS = (
    "homologous_dna_presence",
    "homologous_dna_evidence",
    "predicted_exonic_role",
    "supplied_annotation_role",
    "supplied_annotation_roles",
    "source_parent_occurrence_id",
    "source_parent_transcript_ids",
    "target_parent_occurrence_ids",
    "target_parent_transcript_ids",
    "dna_aligned_blocks",
    "predicted_role_blocks",
    "supplied_annotation_overlaps",
    "annotation_conflict_blocks",
    "candidate_resolution_status",
    "candidate_search_complete",
    "candidate_search_incomplete_reason",
    "left_anchor_id",
    "right_anchor_id",
    "anchor_interval_status",
    "search_interval",
    "alignment_evidence_scope",
    "alignment_sequence_kind",
    "alignment_backend_version",
    "alignment_score_scheme",
    "alignment_nt_identity",
    "alignment_aa_identity",
    "alignment_known_aligned_pairs",
    "alignment_unknown_aligned_pairs",
    "alignment_query_covered_bases",
    "alignment_target_covered_bases",
    "alignment_query_length",
    "alignment_target_length",
    "alignment_gap_blocks",
    "alignment_relative_strand",
    "alignment_candidate_id",
    "alignment_alternative_candidate_ids",
    "alignment_enumeration_complete",
    "alignment_incomplete_reason",
    "original_annotation_start",
    "original_annotation_end",
    "linked_feature_start",
    "linked_feature_end",
    "search_bound_start",
    "search_bound_end",
    "search_left_limit_status",
    "search_right_limit_status",
    "hit_search_limit_status",
)

_SEARCH_COMPAT_FIELDS = (
    "search_expansion_status",
    "search_original_start",
    "search_original_end",
    "search_expanded_start",
    "search_expanded_end",
    "range_status",
)


def support_score(row):
    if row.get("evidence_status") == "supports_absence":
        return None
    identity = to_float(row.get("sequence_score"))
    coverage = to_float(row.get("sequence_coverage"), 1.0)
    return min(identity, coverage)


def completion_call(row, score, threshold):
    status = row.get("evidence_status", "ambiguous")
    anchor_status = row.get("anchor_interval_status", "unknown")
    if anchor_status not in {
        "ordered_double_flank_bounded_interval",
        "ordered_double_flank_empty_interval",
    }:
        return "ambiguous_evidence"
    event = row.get("inferred_event", "")
    frame = row.get("frame_status", "")
    inferred_role = row.get("inferred_role", "unknown")
    if (
        status == "supports_hidden_segment"
        and event == "protein_cds_projection"
        and row.get("predicted_role") in EXON_LIKE_ROLES
        and score >= threshold
    ):
        return "predicted_exon_candidate"
    if (
        status == "conflicts_annotation"
        and row.get("annotation_status") == "protein_projection_boundary_conflict"
        and row.get("predicted_role") in EXON_LIKE_ROLES
        and score >= threshold
    ):
        return "boundary_conflict_candidate"
    if status == "supports_hidden_segment" and inferred_role in EXON_LIKE_ROLES and score >= threshold:
        if event == "shifted_splice_site":
            return "shifted_splice_site_candidate"
        if event == "intron_deletion_joined_exon":
            return "joined_exon_candidate"
        if frame == "frameshift_or_stop_risk":
            return "hidden_segment_with_frame_disruption"
        return "hidden_segment_candidate"
    if status == "conflicts_annotation" and inferred_role in EXON_LIKE_ROLES and score >= threshold:
        if row.get("annotation_status") == "protein_projection_boundary_conflict":
            return "boundary_conflict_candidate"
        return "annotation_conflict_candidate"
    if status == "supports_annotation":
        return "supports_annotation"
    if status == "supports_absence":
        return "supports_true_absence"
    if status == "homologous_sequence_candidate":
        return "homologous_sequence_candidate"
    return "ambiguous_evidence"


def _locus_key(header):
    parts = header.split("|", 2)
    return tuple(parts[:2]) if len(parts) >= 2 else ("NA", header)


def _locus_geometry(locus_header):
    locus = locus_header.split("|", 2)[2]
    contig, interval, strand = locus.rsplit(":", 2)
    lower, upper = (int(value) for value in interval.split("-", 1))
    return contig, lower, upper, strand


def _project_locus_interval(start, end, locus_header):
    contig, lower, upper, strand = _locus_geometry(locus_header)
    local = ClosedInterval1(int(start), int(end)).to_interval0()
    locus = ClosedInterval1(lower, upper).to_interval0()
    projected = ClosedInterval1.from_interval0(local_interval_to_genome(local, locus, strand))
    return contig, projected.start, projected.end, strand


def _project_occurrence_interval(start, end, occurrence):
    local = ClosedInterval1(int(start), int(end)).to_interval0()
    parent = ClosedInterval1(int(occurrence["start"]), int(occurrence["end"])).to_interval0()
    projected = ClosedInterval1.from_interval0(
        local_interval_to_genome(local, parent, occurrence.get("strand", "+"))
    )
    return projected.start, projected.end


def _json_records(records):
    return json.dumps(records, separators=(",", ":"), sort_keys=True) if records else "NA"


def _transcripts_for_occurrence(occurrence, transcript_ids_by_occurrence):
    transcript_ids = transcript_ids_by_occurrence.get(occurrence.get("occurrence_id"), set())
    if transcript_ids:
        return sorted(transcript_ids)
    return sorted(set(_split_transcript_ids(occurrence.get("transcript_id"))))


def _alignment_genome_blocks(alignment, locus_header, source_occurrence, source_transcript_ids):
    if alignment is None:
        return []
    blocks = []
    for rank, (query_start, query_end, target_start, target_end) in enumerate(
        alignment.aligned_blocks, start=1
    ):
        source_start, source_end = _project_occurrence_interval(
            query_start, query_end, source_occurrence
        )
        target_contig, target_genome_start, target_genome_end, target_locus_strand = (
            _project_locus_interval(target_start, target_end, locus_header)
        )
        target_strand = target_locus_strand
        if alignment.strand in {"+", "-"}:
            target_strand = (
                target_locus_strand
                if alignment.strand == "+"
                else ("-" if target_locus_strand == "+" else "+")
            )
        blocks.append(
            {
                "block_id": f"dna_block_{rank}",
                "block_resolution": "aligned_block",
                "query_start": int(query_start),
                "query_end": int(query_end),
                "target_locus_start": int(target_start),
                "target_locus_end": int(target_end),
                "source_contig": source_occurrence.get("contig", "NA"),
                "source_start": source_start,
                "source_end": source_end,
                "source_strand": source_occurrence.get("strand", "NA"),
                "source_occurrence_id": source_occurrence.get("occurrence_id", "NA"),
                "source_transcript_ids": source_transcript_ids,
                "target_contig": target_contig,
                "target_start": target_genome_start,
                "target_end": target_genome_end,
                "target_strand": target_strand,
                "relative_strand": alignment.strand,
            }
        )
    return blocks


def _protein_projection_blocks(candidates, source_occurrence, source_transcript_ids):
    blocks = []
    for rank, candidate in enumerate(candidates, start=1):
        try:
            start = int(candidate["start"])
            end = int(candidate["end"])
        except (KeyError, TypeError, ValueError):
            continue
        blocks.append(
            {
                "block_id": f"protein_projection_{rank}",
                "block_resolution": "projected_cds_block",
                "source_occurrence_id": source_occurrence.get("occurrence_id", "NA"),
                "source_transcript_ids": (
                    [candidate["reference_transcript_id"]]
                    if candidate.get("reference_transcript_id") not in {None, "", "NA"}
                    else source_transcript_ids
                ),
                "source_protein_id": candidate.get("reference_protein_id", "NA"),
                "source_protein_start": candidate.get("protein_overlap_query_start", "NA"),
                "source_protein_end": candidate.get("protein_overlap_query_end", "NA"),
                "projection_parent_id": candidate.get("parent_id", "NA"),
                "target_contig": candidate.get("contig", "NA"),
                "target_start": start,
                "target_end": end,
                "target_strand": candidate.get("strand", "NA"),
            }
        )
    return blocks


def _candidate_span_blocks(candidates, source_occurrence, source_transcript_ids, backend=None):
    blocks = []
    for rank, candidate in enumerate(candidates, start=1):
        if backend and candidate.get("backend") != backend:
            continue
        try:
            start = int(candidate["start"])
            end = int(candidate["end"])
        except (KeyError, TypeError, ValueError):
            continue
        blocks.append(
            {
                "block_id": f"candidate_span_{rank}",
                "block_resolution": "candidate_span",
                "source_occurrence_id": source_occurrence.get("occurrence_id", "NA"),
                "source_transcript_ids": source_transcript_ids,
                "target_contig": candidate.get("contig", "NA"),
                "target_start": start,
                "target_end": end,
                "target_strand": candidate.get("strand", "NA"),
            }
        )
    return blocks


def _annotation_overlaps_for_blocks(
    blocks,
    occurrences,
    transcript_ids_by_occurrence,
):
    records = []
    seen = set()
    for block in blocks:
        try:
            block_interval = ClosedInterval1(
                int(block["target_start"]), int(block["target_end"])
            ).to_interval0()
        except (KeyError, TypeError, ValueError):
            continue
        for occurrence in occurrences:
            if occurrence.get("presence_status") != "present":
                continue
            if occurrence.get("contig") != block.get("target_contig"):
                continue
            try:
                occurrence_interval = ClosedInterval1(
                    int(occurrence["start"]), int(occurrence["end"])
                ).to_interval0()
            except (KeyError, TypeError, ValueError):
                continue
            intersection = block_interval.intersection(occurrence_interval)
            if intersection is None:
                continue
            closed = ClosedInterval1.from_interval0(intersection)
            transcript_ids = _transcripts_for_occurrence(
                occurrence, transcript_ids_by_occurrence
            )
            key = (
                block.get("block_id", "NA"),
                occurrence.get("occurrence_id", "NA"),
                closed.start,
                closed.end,
            )
            if key in seen:
                continue
            seen.add(key)
            records.append(
                {
                    "block_id": block.get("block_id", "NA"),
                    "occurrence_id": occurrence.get("occurrence_id", "NA"),
                    "transcript_ids": transcript_ids,
                    "source_feature_ids": _split_transcript_ids(
                        occurrence.get("source_feature_id")
                    ),
                    "source_parent_ids": _split_transcript_ids(
                        occurrence.get("source_parents")
                        or occurrence.get("source_parent")
                    ),
                    "role": occurrence.get("role", "unknown"),
                    "coding_status": occurrence.get("coding_status", "unknown"),
                    "strand_relation": (
                        "sense"
                        if occurrence.get("strand") == block.get("target_strand")
                        else "antisense"
                    ),
                    "overlap_start": closed.start,
                    "overlap_end": closed.end,
                    "contains_block": (
                        occurrence_interval.start0 <= block_interval.start0
                        and block_interval.end0 <= occurrence_interval.end0
                    ),
                }
            )
    return records


def _supplied_annotation_role(records):
    if not records:
        return "unknown"
    sense_roles = sorted(
        {record["role"] for record in records if record["strand_relation"] == "sense"}
    )
    if len(sense_roles) == 1:
        return sense_roles[0]
    if len(sense_roles) > 1:
        return "multiple_supplied_roles"
    return "antisense_annotation_only"


def _overlapping_annotation_role(start, end, locus_header, occurrences, hit_strand):
    contig, hit_start, hit_end, _strand = _project_locus_interval(start, end, locus_header)
    roles = []
    for row in occurrences:
        if row.get("contig") != contig:
            continue
        try:
            row_start = int(row["start"])
            row_end = int(row["end"])
        except (KeyError, TypeError, ValueError):
            continue
        if row_start <= hit_end and row_end >= hit_start:
            roles.append((row.get("role", "unknown"), row.get("strand", "NA")))
    role_values = {role for role, _strand in roles}
    if role_values & EXON_LIKE_ROLES:
        sense_roles = {role for role, strand in roles if role in EXON_LIKE_ROLES and strand == hit_strand}
        if "CDS" in sense_roles:
            return "CDS", True
        if sense_roles:
            return next(iter(sorted(sense_roles))), True
        if "CDS" in role_values:
            return "CDS", False
        return next(iter(sorted(role_values & EXON_LIKE_ROLES))), False
    if roles:
        return next(iter(sorted(role_values))), True
    return "unknown", True


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


def _split_transcript_ids(value):
    if not value:
        return []
    return [part for part in str(value).replace(",", ";").split(";") if part]


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


def _oriented_locus_slice(locus_header, locus_sequence, start, end):
    contig, locus_start, locus_end, strand = _locus_geometry(locus_header)
    start = int(start)
    end = int(end)
    if start < locus_start or end > locus_end:
        return None, "interval_outside_locus_sequence"
    if strand == "-":
        rel_start = locus_end - end + 1
        rel_end = locus_end - start + 1
    else:
        rel_start = start - locus_start + 1
        rel_end = end - locus_start + 1
    return locus_sequence[rel_start - 1 : rel_end], "ok"


def _relative_interval_within_span(locus_header, span_start, span_end, feature_start, feature_end):
    _contig, _locus_start, _locus_end, strand = _locus_geometry(locus_header)
    if strand == "-":
        return span_end - int(feature_end) + 1, span_end - int(feature_start) + 1
    return int(feature_start) - span_start + 1, int(feature_end) - span_start + 1


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


def _protein_context_provenance(context):
    if not context:
        return {
            "reference_protein_id": "NA",
            "reference_transcript_id": "NA",
            "reference_protein_query_start": "NA",
            "reference_protein_query_end": "NA",
            "reference_cds_length": "NA",
            "reference_cds_phase": "NA",
        }
    return {
        "reference_protein_id": context.get("protein_id", "NA"),
        "reference_transcript_id": context.get("transcript_id", "NA"),
        "reference_protein_query_start": context.get("query_start", "NA"),
        "reference_protein_query_end": context.get("query_end", "NA"),
        "reference_cds_length": context.get("cds_length", "NA"),
        "reference_cds_phase": context.get("cds_phase", "NA"),
    }


def _coding_length_and_phase(path_row, occurrence):
    length = 0
    for record in (path_row, occurrence):
        if record.get("coding_status") in {"noncoding", "not_applicable"}:
            break
        value = next((record[key] for key in ("cds_length", "coding_length")
                      if record.get(key) not in {None, "", "NA", "."}), None)
        if value is not None:
            length = int(value)
            break
        intervals = record.get("cds_intervals", "NA")
        if intervals not in {None, "", "NA", "."}:
            previous_end = 0
            for start, end in sorted(tuple(map(int, interval.split("-"))) for interval in intervals.split(";")):
                length += max(0, end - max(start, previous_end + 1) + 1)
                previous_end = max(previous_end, end)
            break
    phase = next((record[key] for record in (path_row, occurrence) for key in ("cds_phase", "phase")
                  if str(record.get(key)) in {"0", "1", "2"}), ".")
    return length, phase


def _reference_protein_context(representative, transcript_paths, proteins, occ_by_id):
    occurrence_id = representative.get("occurrence_id")
    path_rows_by_transcript = defaultdict(list)
    containing_transcripts = []
    for row in transcript_paths:
        if (
            row.get("species") == representative.get("species")
            and row.get("gene_copy_id") == representative.get("gene_copy_id")
            and row.get("transcript_id")
        ):
            path_rows_by_transcript[row.get("transcript_id")].append(row)
            if row.get("occurrence_id") == occurrence_id:
                containing_transcripts.append(row.get("transcript_id"))

    if not containing_transcripts:
        for transcript_id in _split_transcript_ids(representative.get("transcript_id")):
            if transcript_id:
                containing_transcripts.append(transcript_id)
    if not containing_transcripts:
        return None

    seen = set()
    contexts = []
    for transcript_id in containing_transcripts:
        if transcript_id in seen:
            continue
        seen.add(transcript_id)
        protein_id = f"{representative['species']}|{representative['gene_copy_id']}|{transcript_id}"
        protein = proteins.get(protein_id, "")
        if not protein:
            continue
        rows = path_rows_by_transcript.get(transcript_id, [])
        if not rows:
            continue
        ordered = sorted(rows, key=lambda item: _safe_int(item.get("path_rank"), 0))
        coding_rows = []
        initial_phase = None
        for row in ordered:
            occurrence = occ_by_id.get(row.get("occurrence_id"), {})
            if row.get("coding_status", occurrence.get("coding_status")) != "coding":
                continue
            cds_length, phase = _coding_length_and_phase(row, occurrence)
            if cds_length <= 0:
                continue
            if initial_phase is None and str(phase) in {"0", "1", "2"}:
                initial_phase = int(phase)
            coding_rows.append((row, cds_length, phase))
        if initial_phase is None:
            initial_phase = 0
        cds_offset = 0
        for row, cds_length, phase in coding_rows:
            if row.get("occurrence_id") == occurrence_id:
                translated_start = max(0, cds_offset - initial_phase)
                translated_end = min(len(protein) * 3, max(0, cds_offset + cds_length - initial_phase)) - 1
                if translated_end >= translated_start:
                    contexts.append(
                        {
                            "protein_id": protein_id,
                            "protein": protein,
                            "query_start": translated_start // 3 + 1,
                            "query_end": translated_end // 3 + 1,
                            "transcript_id": transcript_id,
                            "cds_length": cds_length,
                            "cds_phase": phase,
                        }
                    )
                break
            cds_offset += cds_length
    if contexts:
        return max(contexts, key=lambda item: (item["cds_length"], item["query_end"] - item["query_start"], item["transcript_id"]))
    return None


def _protein_projection(
    representative,
    locus,
    locus_header,
    transcript_paths,
    proteins,
    occ_by_id,
    min_identity,
    min_coverage,
    threads,
):
    context = _reference_protein_context(representative, transcript_paths, proteins, occ_by_id)
    if context is None:
        return None, "protein_reference_context_unavailable"
    try:
        projections = protein_locus_exons({context["protein_id"]: context["protein"]}, locus, threads=threads)
    except AlignmentBackendError as exc:
        return None, f"protein_projection_unresolved:{exc}"
    return _protein_projection_from_rows(context, projections, locus_header, min_identity, min_coverage)


def _protein_projection_from_rows(context, projections, locus_header, min_identity, min_coverage, candidates=None):
    matches = []
    exon_aa_length = max(1, context["query_end"] - context["query_start"] + 1)
    for row in projections:
        if row.get("protein_id") != context["protein_id"]:
            continue
        if row.get("projection_status") == "missing_cds_target":
            continue
        try:
            query_start = int(row.get("query_start", row.get("query_aa_start", 0)) or 0)
            query_end = int(row.get("query_end", row.get("query_aa_end", 0)) or 0)
        except (TypeError, ValueError):
            continue
        overlap_start = max(query_start, context["query_start"])
        overlap_end = min(query_end, context["query_end"])
        per_exon_coverage = max(0, overlap_end - overlap_start + 1) / exon_aa_length
        if (
            overlap_start <= overlap_end
            and to_float(row.get("identity")) >= min_identity
            and per_exon_coverage >= min_coverage
        ):
            matches.append({
                **row,
                "per_exon_coverage": per_exon_coverage,
                "protein_cds_query_start": query_start,
                "protein_cds_query_end": query_end,
                "protein_overlap_query_start": overlap_start,
                "protein_overlap_query_end": overlap_end,
            })
    if not matches:
        return None, "protein_projection_did_not_unambiguously_cover_reference_cds_interval"
    if candidates is not None:
        for row in matches:
            contig, start, end, strand = _project_locus_interval(row["target_start"], row["target_end"], locus_header)
            candidates.append({
                "contig": contig, "start": start, "end": end,
                "strand": strand if row.get("strand", "+") == "+" else ("-" if strand == "+" else "+"),
                "backend": "miniprot", "identity": to_float(row.get("identity")),
                "interval_scope": "projected_cds_container",
                **{field: row[field] for field in _PROTEIN_QUERY_INTERVAL_FIELDS},
                "reference_protein_id": context["protein_id"],
                "reference_transcript_id": context.get("transcript_id", "NA"),
                "reference_protein_query_start": context["query_start"],
                "reference_protein_query_end": context["query_end"],
                "coverage": row["per_exon_coverage"], "parent_id": row.get("parent_id", "NA"),
                "identity_scope": row.get("identity_scope", "unknown"),
                "phase": row.get("phase", "NA"), "cigar": row.get("cigar", "NA"),
            })
    signatures = {
        (
            row.get("parent_id", "NA"),
            int(row.get("target_start", 0) or 0),
            int(row.get("target_end", 0) or 0),
            row.get("strand", "NA"),
        )
        for row in matches
    }
    if len(signatures) != 1:
        return None, "protein_projection_ambiguous_multiple_possible_mappings"
    best = max(matches, key=lambda row: (to_float(row.get("identity")), row.get("per_exon_coverage", 0.0)))
    contig, start, end, strand = _project_locus_interval(best["target_start"], best["target_end"], locus_header)
    projected_strand = best.get("strand", "+")
    if projected_strand in {"+", "-"} and strand in {"+", "-"}:
        strand = "+" if projected_strand == strand else "-"
    # Query overlap identifies a CDS container, without localizing exact segment bounds inside it.
    return {
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "identity": to_float(best.get("identity")),
        "coverage": best.get("per_exon_coverage", 0.0),
        "phase": best.get("phase", "NA"),
        "cigar": best.get("cigar", "NA"),
        "parent_id": best.get("parent_id", "NA"),
        "identity_scope": best.get("identity_scope", "unknown"),
        "interval_scope": "projected_cds_container",
        **{field: best[field] for field in _PROTEIN_QUERY_INTERVAL_FIELDS},
    }, "protein_projection_supports_cds"


def _cached_protein_projection(
    context,
    family,
    copy_key,
    locus_header,
    locus,
    protein_contexts_by_family,
    projection_cache,
    min_identity,
    min_coverage,
    threads,
    candidates=None,
):
    if context is None:
        return None, "protein_reference_context_unavailable"
    cache_key = (family, copy_key[1], copy_key[2], locus_header)
    if cache_key not in projection_cache:
        proteins = {
            protein_id: protein_context["protein"]
            for protein_id, protein_context in protein_contexts_by_family.get(family, {}).items()
        }
        try:
            projection_cache[cache_key] = (protein_locus_exons(proteins, locus, threads=threads), "")
        except AlignmentBackendError as exc:
            projection_cache[cache_key] = ([], f"protein_projection_unresolved:{exc}")
    projections, error = projection_cache[cache_key]
    if error:
        return None, error
    return _protein_projection_from_rows(context, projections, locus_header, min_identity, min_coverage, candidates)


def generate_sequence_evidence(
    input_dir,
    result_dir,
    min_identity=0.70,
    min_coverage=0.60,
    threads=1,
    aligner="internal",
    short_context_max_length=300,
):
    """Search missing homologous exon sequences inside supplied homologous gene loci."""
    input_dir = Path(input_dir)
    result_dir = Path(result_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv")
    elements = read_tsv(result_dir / "element_correspondence.tsv", optional=True)
    sequences = parse_fasta(input_dir / "segment_sequences.fasta")
    locus_records = parse_fasta(input_dir / "gene_loci.fasta")
    locus_metadata = _read_locus_metadata(input_dir)
    transcript_paths = read_tsv(input_dir / "transcript_paths.tsv", optional=True)
    transcript_ids_by_occurrence = defaultdict(set)
    for path_row in transcript_paths:
        if path_row.get("occurrence_id") and path_row.get("transcript_id"):
            transcript_ids_by_occurrence[path_row["occurrence_id"]].add(path_row["transcript_id"])
    proteins = parse_fasta(input_dir / "protein_sequences.fasta") if aligner == "miniprot" else {}
    loci = {_locus_key(name): (name, sequence) for name, sequence in locus_records.items()}
    if not elements or not loci:
        return []

    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    element_by_occ = {row["occurrence_id"]: row["element_id"] for row in elements}
    rows_by_element = defaultdict(list)
    copies_by_family_species = defaultdict(set)
    occurrences_by_copy = defaultdict(list)
    present_elements = defaultdict(set)
    core_exonic_occurrences = set()
    annotated_exon_counts = Counter()
    for occurrence in occurrences:
        key = (occurrence["family_id"], occurrence["species"])
        copies_by_family_species[key].add(occurrence["gene_copy_id"])
        occurrences_by_copy[(occurrence["family_id"], occurrence["species"], occurrence["gene_copy_id"])].append(occurrence)
        element = element_by_occ.get(occurrence["occurrence_id"])
        if element and occurrence.get("role") in EXON_LIKE_ROLES:
            present_elements[(occurrence["family_id"], occurrence["species"], occurrence["gene_copy_id"])].add(element)
    for row in elements:
        occurrence = occ_by_id.get(row.get("occurrence_id", ""))
        if (
            occurrence
            and row.get("element_class") == "exon_like"
            and occurrence.get("role") in EXON_LIKE_ROLES
        ):
            annotated_exon_counts[(occurrence["family_id"], row["element_id"])] += 1
    for row in elements:
        occurrence = occ_by_id.get(row.get("occurrence_id", ""))
        if (
            occurrence
            and row.get("element_class") == "exon_like"
            and occurrence.get("role") in EXON_LIKE_ROLES
        ):
            key = (occurrence["family_id"], row["element_id"])
            if row.get("membership_call", "core_member") == "core_member" or annotated_exon_counts[key] == 1:
                rows_by_element[key].append((row, occurrence))
            if row.get("membership_call", "core_member") == "core_member":
                core_exonic_occurrences.add(row.get("occurrence_id"))

    representative_infos = []
    for (family, element), members in sorted(rows_by_element.items()):
        representative_row, representative = max(
            members,
            key=lambda pair: len(sequences.get(pair[1]["occurrence_id"], "")),
        )
        if aligner == "miniprot":
            protein_members = [
                pair for pair in members
                if _reference_protein_context(pair[1], transcript_paths, proteins, occ_by_id) is not None
            ]
            if protein_members:
                representative_row, representative = min(
                    protein_members,
                    key=lambda pair: (-_safe_int(pair[1].get("cds_length"), 0), pair[1]["occurrence_id"]),
                )
        query = sequences.get(representative["occurrence_id"], "")
        if not query:
            continue
        representative_infos.append(
            {
                "family": family,
                "element": element,
                "members": members,
                "representative_row": representative_row,
                "representative": representative,
                "query": query,
                "protein_context": _reference_protein_context(representative, transcript_paths, proteins, occ_by_id)
                if aligner == "miniprot"
                else None,
            }
        )

    protein_contexts_by_family = defaultdict(dict)
    if aligner == "miniprot":
        for info in representative_infos:
            context = info.get("protein_context")
            if context is not None:
                protein_contexts_by_family[info["family"]][context["protein_id"]] = context

    protein_projection_cache = {}
    evidence = []
    for info in representative_infos:
        family = info["family"]
        element = info["element"]
        representative_row = info["representative_row"]
        representative = info["representative"]
        query = info["query"]
        protein_context = info.get("protein_context")
        protein_context_provenance = _protein_context_provenance(protein_context)
        source_transcript_ids = _transcripts_for_occurrence(
            representative, transcript_ids_by_occurrence
        )
        source_occurrences, source_path_status = _source_occurrence_path(
            family,
            representative,
            transcript_paths,
            occurrences_by_copy,
            occ_by_id,
            element_by_occ,
            core_exonic_occurrences,
        )
        source_elements = [element_by_occ.get(row.get("occurrence_id")) for row in source_occurrences]
        source_index = next(
            (
                index for index, row in enumerate(source_occurrences)
                if row.get("occurrence_id") == representative.get("occurrence_id")
            ),
            -1,
        )
        left_element = source_elements[source_index - 1] if source_index > 0 else None
        right_element = source_elements[source_index + 1] if 0 <= source_index < len(source_elements) - 1 else None
        terminal_need = None
        if source_elements and source_index == 0:
            terminal_need = "left"
        elif source_index >= 0 and source_index == len(source_elements) - 1:
            terminal_need = "right"
        for (candidate_family, species), copies in sorted(copies_by_family_species.items()):
            if candidate_family != family or len(copies) != 1:
                continue
            gene_copy = next(iter(copies))
            copy_key = (family, species, gene_copy)
            if element in present_elements[copy_key]:
                continue
            locus_header, locus = loci.get((species, gene_copy), ("", ""))
            if not locus:
                continue
            metadata = locus_metadata.get((species, gene_copy), {})
            locus_header, locus, expansion = _extended_locus(
                input_dir, locus_header, locus, metadata, terminal_need
            )
            anchor_context = _ordered_anchor_search_context(
                source_occurrences,
                occurrences_by_copy[copy_key],
                element_by_occ,
                representative,
                left_element,
                right_element,
                locus_header,
                locus,
            )
            alignment = None
            alignment_error = ""
            alignment_evidence_scope = "unavailable"
            identity = coverage = 0.0
            dna_identity = dna_coverage = 0.0
            nucleotide_aligner = "minimap2" if aligner == "miniprot" else aligner
            bounded_target = anchor_context.get("sequence", "")
            bounded_backend = (
                "internal"
                if len(query) <= int(short_context_max_length)
                else nucleotide_aligner
            )
            if anchor_context["valid_double_flank"] and bounded_target:
                try:
                    alignment = local_alignment_stats(
                        query,
                        bounded_target,
                        backend=bounded_backend,
                        threads=threads,
                        query_occurrence_id=representative.get("occurrence_id"),
                        target_occurrence_id=None,
                        query_transcript_id=(source_transcript_ids[0] if len(source_transcript_ids) == 1 else None),
                        left_anchor_id=anchor_context["left_anchor_id"],
                        right_anchor_id=anchor_context["right_anchor_id"],
                        search_interval=anchor_context["search_interval"],
                    )
                    alignment_evidence_scope = (
                        "anchor_bounded_short_local"
                        if bounded_backend == "internal"
                        else "anchor_bounded_external_local"
                    )
                except AlignmentBackendError as exc:
                    if bounded_backend == "internal" and nucleotide_aligner != "internal":
                        try:
                            alignment = local_alignment_stats(
                                query,
                                bounded_target,
                                backend=nucleotide_aligner,
                                threads=threads,
                            )
                            alignment_evidence_scope = "anchor_bounded_external_local"
                        except AlignmentBackendError as fallback_exc:
                            alignment_error = f"bounded_alignment_unavailable:{fallback_exc}"
                    else:
                        alignment_error = f"bounded_alignment_unavailable:{exc}"
                if alignment is not None:
                    _set_alignment_context(
                        alignment,
                        representative.get("occurrence_id", "NA"),
                        None,
                        source_transcript_ids[0] if len(source_transcript_ids) == 1 else None,
                        anchor_context["left_anchor_id"],
                        anchor_context["right_anchor_id"],
                        anchor_context["search_interval"],
                    )
                    _offset_alignment_target(alignment, anchor_context["target_offset0"])
            if alignment is None and not anchor_context["valid_double_flank"]:
                try:
                    alignment = local_alignment_stats(
                        query, locus, backend=nucleotide_aligner, threads=threads
                    )
                    alignment_evidence_scope = "whole_locus_descriptive_fallback"
                except AlignmentBackendError as exc:
                    alignment_error = str(exc)
            if alignment is not None:
                dna_identity = alignment.identity
                dna_coverage = alignment.query_coverage or alignment.coverage
                identity = dna_identity
                coverage = dna_coverage
            supported = (
                alignment_evidence_scope.startswith("anchor_bounded_")
                and alignment is not None
                and (bool(alignment.aligned_blocks) or int(alignment.aligned_pairs or 0) > 0)
                and identity >= min_identity
                and coverage >= min_coverage
            )
            projected = None
            projected_candidate = None
            protein_candidates = []
            projection_reason = "not_requested"
            protein_identity = protein_coverage = 0.0
            if aligner == "miniprot":
                projected_candidate, projection_reason = _cached_protein_projection(
                    protein_context,
                    family,
                    copy_key,
                    locus_header,
                    locus,
                    protein_contexts_by_family,
                    protein_projection_cache,
                    min_identity,
                    min_coverage,
                    threads,
                    candidates=protein_candidates,
                )
                projected = (
                    projected_candidate
                    if _projection_within_anchor_interval(projected_candidate, anchor_context)
                    else None
                )
                if projected:
                    protein_identity = projected["identity"]
                    protein_coverage = projected["coverage"]
                    identity = protein_identity
                    coverage = protein_coverage
            if not source_occurrences:
                deletion_supported = False
                deletion_reason = source_path_status
                deletion_provenance = {"backend": "NA", "cigar": "NA"}
            else:
                deletion_supported, deletion_reason, deletion_provenance = _ordered_anchor_deletion_support(
                    source_occurrences,
                    occurrences_by_copy[copy_key],
                    element_by_occ,
                    representative,
                    left_element,
                    right_element,
                    loci.get((representative["species"], representative["gene_copy_id"]), ("", ""))[0],
                    loci.get((representative["species"], representative["gene_copy_id"]), ("", ""))[1],
                    locus_header,
                    locus,
                    min_identity,
                    min_coverage,
                    threads=threads,
                )
            bounded_protein_candidates = [
                candidate for candidate in protein_candidates
                if _projection_within_anchor_interval(candidate, anchor_context)
            ]
            if projected is None and len(bounded_protein_candidates) == 1:
                projected = bounded_protein_candidates[0]
                projection_reason = "protein_projection_resolved_by_ordered_anchor_interval"
                protein_identity = projected["identity"]
                protein_coverage = projected["coverage"]
                identity = protein_identity
                coverage = protein_coverage
            elif projected_candidate is not None and projected is None:
                projection_reason = "protein_projection_outside_ordered_anchor_interval"
            ambiguous_dna = supported and _ambiguous_repeated_mapping(alignment)
            ambiguous_protein = len(bounded_protein_candidates) > 1
            interval_candidates = _nucleotide_interval_candidates(alignment, locus_header) + protein_candidates
            dna_blocks = _alignment_genome_blocks(
                alignment, locus_header, representative, source_transcript_ids
            )
            dna_annotation_overlaps = _annotation_overlaps_for_blocks(
                dna_blocks,
                occurrences_by_copy[copy_key],
                transcript_ids_by_occurrence,
            )
            protein_blocks = _protein_projection_blocks(
                bounded_protein_candidates, representative, source_transcript_ids
            )
            projection_annotation_overlaps = _annotation_overlaps_for_blocks(
                protein_blocks,
                occurrences_by_copy[copy_key],
                transcript_ids_by_occurrence,
            )
            conflict_blocks = []
            correspondence_status = "resolved"
            if not anchor_context["valid_double_flank"]:
                has_descriptive_candidate = (
                    (
                        alignment is not None
                        and (bool(alignment.aligned_blocks) or int(alignment.aligned_pairs or 0) > 0)
                        and dna_identity >= min_identity
                        and dna_coverage >= min_coverage
                    )
                    or bool(protein_candidates)
                )
                primary_mapping_status = (
                    "whole_locus_descriptive_candidate" if has_descriptive_candidate else "unresolved"
                )
                interval_scope = alignment_evidence_scope
                correspondence_status = "unknown"
                primary_backend = (
                    "miniprot"
                    if protein_candidates
                    else alignment.backend if alignment is not None else nucleotide_aligner
                )
                primary_cigar = alignment.cigar if alignment is not None else "NA"
                primary_score = dna_identity if alignment is not None else 0.0
                primary_coverage = dna_coverage if alignment is not None else 0.0
                predicted_role = "unknown"
                inferred_role = "unknown"
                status = "ambiguous"
                annotation_status = anchor_context["status"]
                conclusion = "unknown"
                contig = metadata.get("contig") or _locus_geometry(locus_header)[0]
                start = end = "NA"
                strand = metadata.get("strand") or _locus_geometry(locus_header)[3]
                confidence = "low"
            elif ambiguous_dna or ambiguous_protein:
                primary_mapping_status = "ambiguous_repeated_mapping"
                interval_scope = "ambiguous_candidates"
                correspondence_status = "unknown"
                primary_backend = "miniprot" if bounded_protein_candidates else alignment.backend
                primary_cigar = "NA"
                if bounded_protein_candidates:
                    primary_candidate = max(bounded_protein_candidates, key=lambda row: (row["coverage"], row["identity"]))
                    primary_score = primary_candidate["identity"]
                    primary_coverage = primary_candidate["coverage"]
                else:
                    primary_score, primary_coverage = dna_identity, dna_coverage
                predicted_role = "CDS" if bounded_protein_candidates else "unknown"
                inferred_role = "unknown"
                status = "homologous_sequence_candidate"
                annotation_status = "alignment_ambiguous_repeated_mapping"
                conclusion = "sequence_present_repeated_mapping_ambiguous"
                contig = _locus_geometry(locus_header)[0]
                start = end = strand = "NA"
                confidence = "low"
            elif projected:
                primary_mapping_status = "protein_projection"
                interval_scope = "projected_cds_container"
                primary_backend = "miniprot"
                primary_cigar = projected.get("cigar", "NA")
                primary_score = protein_identity
                primary_coverage = protein_coverage
                predicted_role = "CDS"
                inferred_role = "predicted_CDS"
                contig = projected["contig"]
                start = projected["start"]
                end = projected["end"]
                strand = projected["strand"]
                overlapping_exons = [
                    row for row in projection_annotation_overlaps
                    if row.get("role") in EXON_LIKE_ROLES
                    and row.get("strand_relation") == "sense"
                ]
                contained_roles = {
                    row["role"] for row in overlapping_exons if row.get("contains_block")
                }
                if strand != _locus_geometry(locus_header)[3]:
                    status = "homologous_sequence_candidate"
                    annotation_status = "protein_projection_antisense_to_gene_locus"
                    inferred_role = "unknown"
                    predicted_role = "CDS"
                    conclusion = "sequence_present_role_unknown"
                elif contained_roles:
                    status = "supports_annotation"
                    annotation_status = "protein_projection_contained_in_annotated_exon"
                    inferred_role = "CDS" if "CDS" in contained_roles else sorted(contained_roles)[0]
                    conclusion = "annotated_exon_sequence_present"
                elif overlapping_exons:
                    status = "conflicts_annotation"
                    annotation_status = "protein_projection_boundary_conflict"
                    inferred_role = "predicted_CDS"
                    conclusion = "predicted_cds_boundary_conflict"
                    conflict_blocks = overlapping_exons
                else:
                    status = "supports_hidden_segment"
                    annotation_status = "protein_projection_supports_missing_cds"
                    conclusion = "predicted_exon_candidate"
                confidence = "high" if identity >= min_identity and coverage >= min_coverage else "low"
            elif supported:
                interval_scope = "aligned_sequence"
                primary_backend = alignment.backend
                primary_cigar = alignment.cigar
                primary_score = dna_identity
                primary_coverage = dna_coverage
                predicted_role = "unknown"
                contig, start, end, strand = _project_locus_interval(
                    alignment.target_start, alignment.target_end, locus_header
                )
                if alignment.strand in {"+", "-"} and strand in {"+", "-"}:
                    strand = "+" if alignment.strand == strand else "-"
                primary_mapping_status = "unique_nucleotide_mapping"
                sense_roles = {
                    row["role"] for row in dna_annotation_overlaps
                    if row.get("strand_relation") == "sense"
                }
                antisense_exonic = any(
                    row.get("role") in EXON_LIKE_ROLES
                    and row.get("strand_relation") == "antisense"
                    for row in dna_annotation_overlaps
                )
                exonic_roles = sense_roles & EXON_LIKE_ROLES
                if exonic_roles:
                    status = "supports_annotation"
                    annotation_status = "alignment_overlaps_annotated_exon"
                    inferred_role = "CDS" if "CDS" in exonic_roles else sorted(exonic_roles)[0]
                    conclusion = "annotated_exon_sequence_present"
                else:
                    status = "homologous_sequence_candidate"
                    if antisense_exonic:
                        annotation_status = "homologous_sequence_overlaps_antisense_exon_annotation"
                    elif sense_roles:
                        annotation_status = "homologous_sequence_overlaps_non_exonic_annotation"
                    else:
                        annotation_status = "unannotated_homologous_sequence_candidate"
                    inferred_role = "unknown"
                    conclusion = "sequence_present_role_unknown"
                confidence = "medium"
            elif deletion_supported:
                primary_mapping_status = "ordered_flank_deletion"
                interval_scope = "not_applicable"
                primary_backend = deletion_provenance.get("backend", "NA")
                primary_cigar = deletion_provenance.get("cigar", "NA")
                primary_score = None
                primary_coverage = None
                predicted_role = "unknown"
                status = "supports_absence"
                annotation_status = "sequence_absence_between_ordered_flanking_homologs"
                inferred_role = "unknown"
                contig = metadata.get("contig") or "supplied_gene_locus"
                start = end = "NA"
                strand = metadata.get("strand") or "+"
                conclusion = "targeted_deletion_supported"
                confidence = "medium"
            else:
                primary_mapping_status = "unresolved"
                interval_scope = "unresolved"
                correspondence_status = "unknown"
                primary_backend = nucleotide_aligner
                primary_cigar = "NA"
                primary_score = 0.0
                primary_coverage = 0.0
                predicted_role = "unknown"
                status = "ambiguous"
                annotation_status = (
                    "alignment_backend_error" if alignment_error else
                    projection_reason if aligner == "miniprot" and projection_reason != "not_requested" else
                    deletion_reason
                )
                inferred_role = "unknown"
                contig = metadata.get("contig") or "supplied_gene_locus"
                start = end = "NA"
                strand = metadata.get("strand") or "+"
                conclusion = "unknown"
                confidence = "low"

            if protein_blocks and (projected or ambiguous_protein):
                observation_blocks = protein_blocks
                annotation_overlaps = projection_annotation_overlaps
            elif ambiguous_dna:
                observation_blocks = _candidate_span_blocks(
                    interval_candidates,
                    representative,
                    source_transcript_ids,
                    backend=alignment.backend,
                )
                annotation_overlaps = _annotation_overlaps_for_blocks(
                    observation_blocks,
                    occurrences_by_copy[copy_key],
                    transcript_ids_by_occurrence,
                )
            elif dna_blocks:
                observation_blocks = dna_blocks
                annotation_overlaps = dna_annotation_overlaps
            else:
                observation_blocks = []
                annotation_overlaps = []

            supplied_role = (
                "ambiguous_candidate_roles"
                if ambiguous_dna or ambiguous_protein
                else _supplied_annotation_role(annotation_overlaps)
            )
            supplied_roles = sorted({row["role"] for row in annotation_overlaps})
            target_parent_occurrence_ids = sorted(
                {row["occurrence_id"] for row in annotation_overlaps}
            )
            target_parent_transcript_ids = sorted(
                {
                    transcript_id
                    for row in annotation_overlaps
                    for transcript_id in row.get("transcript_ids", [])
                }
            )

            if status == "supports_absence":
                homologous_dna_presence = "absent"
                homologous_dna_evidence = "ordered_flank_deletion"
            elif supported or projected or bounded_protein_candidates:
                homologous_dna_presence = "present"
                dna_evidence_types = []
                if supported:
                    dna_evidence_types.append("anchor_bounded_nucleotide_alignment")
                if projected or bounded_protein_candidates:
                    dna_evidence_types.append("anchor_bounded_protein_coding_projection")
                homologous_dna_evidence = ";".join(dna_evidence_types)
            else:
                homologous_dna_presence = "unknown"
                homologous_dna_evidence = "no_resolved_sequence_evidence"

            search_provenance = _search_provenance(
                metadata, locus_header, terminal_need, interval_candidates
            )
            incomplete_reasons = []
            if alignment is not None and not alignment.enumeration_complete:
                incomplete_reasons.append(
                    alignment.incomplete_reason or "nucleotide_candidate_enumeration_truncated"
                )
            if search_provenance["hit_search_limit_status"].startswith("candidate_touches_"):
                incomplete_reasons.append(
                    search_provenance["hit_search_limit_status"]
                )
            if search_provenance["hit_search_limit_status"].startswith("no_hit_with_"):
                incomplete_reasons.append(
                    search_provenance["hit_search_limit_status"]
                )
            if incomplete_reasons:
                candidate_resolution_status = "candidate_search_incomplete"
                candidate_search_complete = "false"
                correspondence_status = "unknown"
            elif ambiguous_dna or ambiguous_protein:
                candidate_resolution_status = "ambiguous"
                candidate_search_complete = (
                    "true"
                    if not bounded_protein_candidates
                    and alignment is not None
                    and alignment.enumeration_complete
                    else "unknown"
                )
            elif correspondence_status == "resolved" and homologous_dna_presence in {"present", "absent"}:
                candidate_resolution_status = "resolved"
                candidate_search_complete = (
                    "true"
                    if supported
                    and not bounded_protein_candidates
                    and alignment is not None
                    and alignment.enumeration_complete
                    else "unknown"
                )
            else:
                candidate_resolution_status = "unresolved"
                candidate_search_complete = "unknown"

            evidence.append(
                {
                    "evidence_id": f"completion_{family}_{element}_{species}",
                    "family_id": family,
                    "species": species,
                    "gene_copy_id": gene_copy,
                    "homology_id": representative_row.get("homology_id", "NA"),
                    "annotation_status": annotation_status,
                    "evidence_status": status,
                    "inferred_role": inferred_role,
                    "predicted_role": predicted_role,
                    "homologous_dna_presence": homologous_dna_presence,
                    "homologous_dna_evidence": homologous_dna_evidence,
                    "predicted_exonic_role": predicted_role,
                    "supplied_annotation_role": supplied_role,
                    "supplied_annotation_roles": ";".join(supplied_roles) or "NA",
                    "source_parent_occurrence_id": representative.get("occurrence_id", "NA"),
                    "source_parent_transcript_ids": ";".join(source_transcript_ids) or "NA",
                    "target_parent_occurrence_ids": ";".join(target_parent_occurrence_ids) or "NA",
                    "target_parent_transcript_ids": ";".join(target_parent_transcript_ids) or "NA",
                    "dna_aligned_blocks": _json_records(dna_blocks),
                    "predicted_role_blocks": _json_records(protein_blocks if predicted_role in EXON_LIKE_ROLES else []),
                    "supplied_annotation_overlaps": _json_records(annotation_overlaps),
                    "annotation_conflict_blocks": _json_records(conflict_blocks),
                    "candidate_resolution_status": candidate_resolution_status,
                    "candidate_search_complete": candidate_search_complete,
                    "candidate_search_incomplete_reason": ";".join(sorted(set(incomplete_reasons))) or "NA",
                    "left_anchor_id": anchor_context["left_anchor_id"],
                    "right_anchor_id": anchor_context["right_anchor_id"],
                    "anchor_interval_status": anchor_context["status"],
                    "search_interval": (
                        json.dumps(anchor_context["search_interval"], separators=(",", ":"))
                        if anchor_context["search_interval"] is not None
                        else "NA"
                    ),
                    "alignment_evidence_scope": alignment_evidence_scope,
                    **_alignment_evidence_provenance(
                        alignment,
                        query_length=len(query),
                        target_length=(
                            len(bounded_target)
                            if alignment_evidence_scope.startswith("anchor_bounded_")
                            else len(locus)
                        ),
                    ),
                    **search_provenance,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "sequence_score": "NA" if primary_score is None else f"{primary_score:.6g}",
                    "sequence_coverage": "NA" if primary_coverage is None else f"{primary_coverage:.6g}",
                    "left_synteny_score": "1" if left_element in present_elements[copy_key] else "0",
                    "right_synteny_score": "1" if right_element in present_elements[copy_key] else "0",
                    "splice_motif_score": "NA",
                    "phase_compatibility": "unknown",
                    "inferred_event": "protein_cds_projection" if projected else "homologous_exon_sequence_search",
                    "frame_status": "unknown",
                    "evidence_conclusion": conclusion,
                    "confidence_flag": confidence,
                    "absence_evidence": deletion_reason,
                    "alignment_backend": primary_backend,
                    "primary_mapping_status": primary_mapping_status,
                    "correspondence_status": correspondence_status,
                    "interval_scope": interval_scope,
                    "interval_candidates": json.dumps(interval_candidates, separators=(",", ":")) if interval_candidates else "NA",
                    "evidence_aligner": aligner,
                    "alignment_error": alignment_error or "NA",
                    "alignment_cigar": primary_cigar,
                    "dna_sequence_score": f"{dna_identity:.6g}" if alignment is not None else "NA",
                    "dna_sequence_coverage": f"{dna_coverage:.6g}" if alignment is not None else "NA",
                    "dna_alignment_backend": alignment.backend if alignment is not None else nucleotide_aligner,
                    "dna_alignment_cigar": alignment.cigar if alignment is not None else "NA",
                    "protein_sequence_score": f"{protein_identity:.6g}" if projected else "NA",
                    "protein_sequence_coverage": f"{protein_coverage:.6g}" if projected else "NA",
                    "protein_alignment_backend": "miniprot" if projected else "NA",
                    "protein_alignment_cigar": projected.get("cigar", "NA") if projected else "NA",
                    "protein_projection_parent_id": projected.get("parent_id", "NA") if projected else "NA",
                    "protein_projection_reason": projection_reason,
                    **{field: projected.get(field, "NA") if projected else "NA" for field in _PROTEIN_QUERY_INTERVAL_FIELDS},
                    "protein_identity_scope": projected.get("identity_scope", "unknown") if projected else "NA",
                    "protein_cigar_scope": "parent_alignment" if projected and projected.get("cigar", "NA") != "NA" else "NA",
                    "protein_cds_phase": projected.get("phase", "NA") if projected else "NA",
                    **protein_context_provenance,
                    "alignment_hit_count": len(protein_candidates) if primary_backend == "miniprot" else alignment.hit_count if alignment is not None and primary_mapping_status != "ordered_flank_deletion" else "NA",
                    "alignment_ambiguous_hit_count": max(0, len(protein_candidates) - 1) if primary_backend == "miniprot" else alignment.ambiguous_hit_count if alignment is not None and primary_mapping_status != "ordered_flank_deletion" else "NA",
                    "alignment_mapping_quality": alignment.mapping_quality if alignment is not None and alignment.mapping_quality is not None and primary_backend != "miniprot" and primary_mapping_status != "ordered_flank_deletion" else "NA",
                    "alignment_alternative_hits": _format_alternative_hits(alignment) if primary_backend != "miniprot" and primary_mapping_status != "ordered_flank_deletion" else "NA",
                    "dna_alignment_hit_count": alignment.hit_count if alignment is not None else 0,
                    "dna_alignment_ambiguous_hit_count": alignment.ambiguous_hit_count if alignment is not None else 0,
                    "dna_alignment_mapping_quality": alignment.mapping_quality if alignment is not None and alignment.mapping_quality is not None else "NA",
                    "dna_alignment_alternative_hits": _format_alternative_hits(alignment),
                    "deletion_alignment_backend": deletion_provenance.get("backend", "NA"),
                    "deletion_alignment_cigar": deletion_provenance.get("cigar", "NA"),
                    "left_flank_identity": deletion_provenance.get("left_flank_identity", "NA"),
                    "left_flank_coverage": deletion_provenance.get("left_flank_coverage", "NA"),
                    "left_flank_paired_bases": deletion_provenance.get("left_flank_paired_bases", "NA"),
                    "left_flank_matches": deletion_provenance.get("left_flank_matches", "NA"),
                    "left_flank_mismatches": deletion_provenance.get("left_flank_mismatches", "NA"),
                    "right_flank_identity": deletion_provenance.get("right_flank_identity", "NA"),
                    "right_flank_coverage": deletion_provenance.get("right_flank_coverage", "NA"),
                    "right_flank_paired_bases": deletion_provenance.get("right_flank_paired_bases", "NA"),
                    "right_flank_matches": deletion_provenance.get("right_flank_matches", "NA"),
                    "right_flank_mismatches": deletion_provenance.get("right_flank_mismatches", "NA"),
                    "search_expansion_status": expansion.get("search_expansion_status", "not_requested"),
                    "search_original_start": expansion.get("search_original_start", metadata.get("search_start", "NA")),
                    "search_original_end": expansion.get("search_original_end", metadata.get("search_end", "NA")),
                    "search_expanded_start": expansion.get("search_expanded_start", "NA"),
                    "search_expanded_end": expansion.get("search_expanded_end", "NA"),
                    "range_status": expansion.get("range_status", metadata.get("range_status", "NA")),
                }
            )
    fields = [
        "evidence_id", "family_id", "species", "gene_copy_id", "homology_id",
        "annotation_status", "evidence_status", "inferred_role", "predicted_role", "contig", "start", "end",
        *_EVIDENCE_PROVENANCE_FIELDS,
        "strand", "sequence_score", "sequence_coverage", "left_synteny_score",
        "right_synteny_score", "splice_motif_score", "phase_compatibility",
        "inferred_event", "frame_status", "evidence_conclusion", "confidence_flag",
        "absence_evidence", "alignment_backend", "primary_mapping_status", "evidence_aligner", "alignment_error", "alignment_cigar",
        "correspondence_status", "interval_scope", "interval_candidates",
        "dna_sequence_score", "dna_sequence_coverage", "dna_alignment_backend", "dna_alignment_cigar",
        "protein_sequence_score", "protein_sequence_coverage", "protein_alignment_backend", "protein_alignment_cigar",
        "protein_projection_parent_id", "protein_projection_reason",
        *_PROTEIN_QUERY_INTERVAL_FIELDS,
        "protein_identity_scope", "protein_cigar_scope", "protein_cds_phase",
        "reference_protein_id", "reference_transcript_id", "reference_protein_query_start",
        "reference_protein_query_end", "reference_cds_length", "reference_cds_phase",
        "alignment_hit_count", "alignment_ambiguous_hit_count", "alignment_mapping_quality", "alignment_alternative_hits",
        "dna_alignment_hit_count", "dna_alignment_ambiguous_hit_count", "dna_alignment_mapping_quality", "dna_alignment_alternative_hits",
        "deletion_alignment_backend", "deletion_alignment_cigar",
        "left_flank_identity", "left_flank_coverage", "left_flank_paired_bases",
        "left_flank_matches", "left_flank_mismatches",
        "right_flank_identity", "right_flank_coverage", "right_flank_paired_bases",
        "right_flank_matches", "right_flank_mismatches",
        "search_expansion_status", "search_original_start", "search_original_end",
        "search_expanded_start", "search_expanded_end", "range_status",
    ]
    write_tsv(result_dir / "sequence_synteny_evidence.tsv", evidence, fields)
    return evidence


def complete_annotation(input_dir, output_dir, threshold=0.55):
    generated = Path(output_dir) / "sequence_synteny_evidence.tsv"
    evidence_path = generated if generated.exists() else Path(input_dir) / "sequence_synteny_evidence.tsv"
    evidence = read_tsv(
        evidence_path,
        ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "annotation_status", "evidence_status"],
        optional=True,
    )
    rows = []
    summary = Counter()
    for row in evidence:
        score = support_score(row)
        status = row.get("evidence_status", "ambiguous")
        call = completion_call(row, score, threshold)
        summary[call] += 1
        rows.append(
            {
                "evidence_id": row["evidence_id"],
                "family_id": row["family_id"],
                "species": row["species"],
                "gene_copy_id": row["gene_copy_id"],
                "homology_id": row["homology_id"],
                "interval": f"{row.get('contig', 'NA')}:{row.get('start', 'NA')}-{row.get('end', 'NA')}:{row.get('strand', 'NA')}",
                "annotation_status": row.get("annotation_status", "unknown"),
                "inferred_role": row.get("inferred_role", "unknown"),
                "predicted_role": row.get("predicted_role", "unknown"),
                **{
                    field: row.get(
                        field,
                        "unknown"
                        if field in {
                            "homologous_dna_presence",
                            "homologous_dna_evidence",
                            "predicted_exonic_role",
                            "supplied_annotation_role",
                            "candidate_resolution_status",
                            "candidate_search_complete",
                        }
                        else "NA",
                    )
                    for field in _EVIDENCE_PROVENANCE_FIELDS
                },
                **{field: row.get(field, "NA") for field in _SEARCH_COMPAT_FIELDS},
                "primary_mapping_status": row.get("primary_mapping_status", "unknown"),
                # Legacy evidence lacks this field; only explicit unknown marks unresolved correspondence.
                "correspondence_status": row.get("correspondence_status", "unspecified"),
                "interval_scope": row.get("interval_scope", "unspecified"),
                **{field: row.get(field, "NA") for field in _PROTEIN_QUERY_INTERVAL_FIELDS},
                "reference_protein_id": row.get("reference_protein_id", "NA"),
                "reference_transcript_id": row.get("reference_transcript_id", "NA"),
                "reference_protein_query_start": row.get("reference_protein_query_start", "NA"),
                "reference_protein_query_end": row.get("reference_protein_query_end", "NA"),
                "interval_candidates": row.get("interval_candidates", "NA"),
                "support_score": "NA" if score is None else f"{score:.6g}",
                "completion_call": call,
                "evidence_status": status,
                "inferred_event": row.get("inferred_event", "NA"),
                "frame_status": row.get("frame_status", "NA"),
                "evidence_conclusion": row.get("evidence_conclusion", "unknown"),
                "confidence_flag": row.get("confidence_flag", "low"),
                "absence_evidence": row.get("absence_evidence", "NA"),
            }
        )
    write_tsv(
        f"{output_dir}/annotation_completion_candidates.tsv",
        rows,
        [
            "evidence_id", "family_id", "species", "gene_copy_id", "homology_id",
            "interval", "annotation_status", "inferred_role", "predicted_role", "support_score",
            *_EVIDENCE_PROVENANCE_FIELDS,
            *_SEARCH_COMPAT_FIELDS,
            "primary_mapping_status", "correspondence_status", "interval_scope", "interval_candidates",
            *_PROTEIN_QUERY_INTERVAL_FIELDS,
            "reference_protein_id", "reference_transcript_id",
            "reference_protein_query_start", "reference_protein_query_end",
            "completion_call", "evidence_status", "inferred_event", "frame_status",
            "evidence_conclusion", "confidence_flag", "absence_evidence",
        ],
    )
    write_tsv(
        f"{output_dir}/annotation_completion_summary.tsv",
        [{"completion_call": key, "count": value} for key, value in sorted(summary.items())],
        ["completion_call", "count"],
    )
    return rows
