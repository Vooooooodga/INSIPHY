"""Sequence-supported annotation completion."""

from collections import Counter, defaultdict
from pathlib import Path

from .alignment import AlignmentBackendError, local_alignment_stats, revcomp
from .elements import EXON_LIKE_ROLES
from .io import fasta_record_length, parse_fasta, read_fasta_interval, read_tsv, to_float, write_tsv


def support_score(row):
    if row.get("evidence_status") == "supports_absence":
        return None
    identity = to_float(row.get("sequence_score"))
    coverage = to_float(row.get("sequence_coverage"), 1.0)
    return min(identity, coverage)


def completion_call(row, score, threshold):
    status = row.get("evidence_status", "ambiguous")
    event = row.get("inferred_event", "")
    frame = row.get("frame_status", "")
    inferred_role = row.get("inferred_role", "unknown")
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
    start = int(start)
    end = int(end)
    if strand == "-":
        hit_start, hit_end = upper - end + 1, upper - start + 1
    else:
        hit_start, hit_end = lower + start - 1, lower + end - 1
    return contig, min(hit_start, hit_end), max(hit_start, hit_end), strand


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
        return locus_header, locus, {"search_expansion_status": "at_contig_boundary"}
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
    provenance = {"backend": "NA", "cigar": "NA"}
    if not left_element or not right_element:
        return False, "missing_flanking_homolog_anchor", provenance
    if not source_locus_header or not source_locus_sequence:
        return False, "source_locus_sequence_unavailable", provenance
    if not target_locus_header or not target_locus_sequence:
        return False, "target_locus_sequence_unavailable", provenance
    source_left = _occurrence_for_element(source_occurrences, element_by_occurrence, left_element)
    source_right = _occurrence_for_element(source_occurrences, element_by_occurrence, right_element)
    target_left = _occurrence_for_element(target_occurrences, element_by_occurrence, left_element)
    target_right = _occurrence_for_element(target_occurrences, element_by_occurrence, right_element)
    if not source_left or not source_right or not target_left or not target_right:
        return False, "flanking_homolog_anchor_not_observed", provenance
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
    if "N" in expected_sequence.upper() or "N" in local_target_sequence.upper():
        return False, "deletion_interval_or_local_target_sequence_contains_N", provenance
    try:
        alignment = local_alignment_stats(source_span, target_span, backend="minimap2", threads=threads)
    except AlignmentBackendError as exc:
        return False, f"deletion_spanning_alignment_unresolved:{exc}", provenance
    provenance = {"backend": alignment.backend, "cigar": alignment.cigar}
    if alignment.strand != "+":
        return False, "deletion_spanning_alignment_not_relative_plus_strand", provenance
    if alignment.identity < min_identity:
        return False, "deletion_spanning_alignment_identity_insufficient", provenance
    blocks = _supported_alignment_blocks(alignment)
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


def _reference_protein_context(representative, transcript_paths, proteins, occ_by_id):
    transcript_id = representative.get("transcript_id")
    if not transcript_id:
        for row in transcript_paths:
            if row.get("occurrence_id") == representative.get("occurrence_id"):
                transcript_id = row.get("transcript_id")
                break
    if not transcript_id:
        return None
    protein_id = f"{representative['species']}|{representative['gene_copy_id']}|{transcript_id}"
    protein = proteins.get(protein_id, "")
    if not protein:
        return None
    rows = [
        row for row in transcript_paths
        if row.get("species") == representative.get("species")
        and row.get("gene_copy_id") == representative.get("gene_copy_id")
        and row.get("transcript_id") == transcript_id
    ]
    ordered = sorted(rows, key=lambda item: int(item.get("path_rank", 0) or 0))
    coding_rows = []
    initial_phase = None
    for row in ordered:
        occurrence = occ_by_id.get(row.get("occurrence_id"), {})
        if occurrence.get("coding_status") != "coding" and row.get("coding_status") != "coding":
            continue
        cds_length = int(occurrence.get("cds_length", 0) or 0)
        if cds_length <= 0:
            continue
        phase = occurrence.get("cds_phase") or row.get("cds_phase") or "."
        if initial_phase is None and phase in {"0", "1", "2"}:
            initial_phase = int(phase)
        coding_rows.append((row, cds_length))
    if initial_phase is None:
        initial_phase = 0
    cds_offset = 0
    for row, cds_length in coding_rows:
        if row.get("occurrence_id") == representative.get("occurrence_id"):
            translated_start = max(0, cds_offset - initial_phase)
            translated_end = max(0, cds_offset + cds_length - initial_phase) - 1
            if translated_end < translated_start:
                return None
            return {
                "protein_id": protein_id,
                "protein": protein,
                "query_start": translated_start // 3 + 1,
                "query_end": translated_end // 3 + 1,
            }
        cds_offset += cds_length
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
    from .alignment import protein_locus_exons
    try:
        projections = protein_locus_exons({context["protein_id"]: context["protein"]}, locus, threads=threads)
    except AlignmentBackendError as exc:
        return None, f"protein_projection_unresolved:{exc}"
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
            matches.append({**row, "per_exon_coverage": per_exon_coverage})
    if not matches:
        return None, "protein_projection_did_not_unambiguously_cover_reference_cds_interval"
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
    }, "protein_projection_supports_cds"


def generate_sequence_evidence(input_dir, result_dir, min_identity=0.70, min_coverage=0.60, threads=1, aligner="internal"):
    """Search missing homologous exon sequences inside supplied homologous gene loci."""
    input_dir = Path(input_dir)
    result_dir = Path(result_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv")
    elements = read_tsv(result_dir / "element_correspondence.tsv", optional=True)
    sequences = parse_fasta(input_dir / "segment_sequences.fasta")
    locus_records = parse_fasta(input_dir / "gene_loci.fasta")
    locus_metadata = _read_locus_metadata(input_dir)
    transcript_paths = read_tsv(input_dir / "transcript_paths.tsv", optional=True)
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
            and row.get("membership_call", "core_member") == "core_member"
            and occurrence.get("role") in EXON_LIKE_ROLES
        ):
            rows_by_element[(occurrence["family_id"], row["element_id"])].append((row, occurrence))
            core_exonic_occurrences.add(row.get("occurrence_id"))

    evidence = []
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
                    key=lambda pair: (-int(pair[1]["cds_length"]), pair[1]["occurrence_id"]),
                )
        query = sequences.get(representative["occurrence_id"], "")
        if not query:
            continue
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
            alignment = None
            alignment_error = ""
            identity = coverage = 0.0
            if aligner != "miniprot":
                try:
                    alignment = local_alignment_stats(query, locus, backend=aligner, threads=threads)
                    identity = alignment.identity
                    coverage = alignment.query_coverage or alignment.coverage
                except AlignmentBackendError as exc:
                    alignment_error = str(exc)
            supported = alignment is not None and identity >= min_identity and coverage >= min_coverage
            projected = None
            projection_reason = "not_requested"
            if aligner == "miniprot":
                projected, projection_reason = _protein_projection(
                    representative, locus, locus_header, transcript_paths, proteins,
                    occ_by_id, min_identity, min_coverage, threads,
                )
                if projected:
                    identity = projected["identity"]
                    coverage = projected["coverage"]
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
            if projected:
                inferred_role = "CDS"
                contig = projected["contig"]
                start = projected["start"]
                end = projected["end"]
                strand = projected["strand"]
                overlapping_exons = [
                    row for row in occurrences_by_copy[copy_key]
                    if row.get("presence_status") == "present"
                    and row.get("role") in EXON_LIKE_ROLES
                    and row.get("contig") == contig
                    and row.get("strand") == strand
                    and int(row["start"]) <= end
                    and int(row["end"]) >= start
                ]
                contained = any(
                    int(row["start"]) <= start and end <= int(row["end"])
                    for row in overlapping_exons
                )
                if strand != _locus_geometry(locus_header)[3]:
                    status = "homologous_sequence_candidate"
                    annotation_status = "protein_projection_antisense_to_gene_locus"
                    inferred_role = "unknown"
                    conclusion = "sequence_present_role_unknown"
                elif contained:
                    status = "supports_annotation"
                    annotation_status = "protein_projection_contained_in_annotated_exon"
                    conclusion = "annotated_exon_sequence_present"
                elif overlapping_exons:
                    status = "conflicts_annotation"
                    annotation_status = "protein_projection_boundary_conflict"
                    conclusion = "predicted_cds_boundary_conflict"
                else:
                    status = "supports_hidden_segment"
                    annotation_status = "protein_projection_supports_missing_cds"
                    conclusion = "predicted_exon"
                confidence = "high" if identity >= min_identity and coverage >= min_coverage else "low"
            elif supported:
                contig, start, end, strand = _project_locus_interval(
                    alignment.target_start, alignment.target_end, locus_header
                )
                if alignment.strand in {"+", "-"} and strand in {"+", "-"}:
                    strand = "+" if alignment.strand == strand else "-"
                overlap_role, sense_overlap = _overlapping_annotation_role(
                    alignment.target_start, alignment.target_end,
                    locus_header, occurrences_by_copy[copy_key], strand,
                )
                if overlap_role in EXON_LIKE_ROLES and sense_overlap:
                    status = "supports_annotation"
                    annotation_status = "alignment_overlaps_annotated_exon"
                    inferred_role = overlap_role
                    conclusion = "annotated_exon_sequence_present"
                else:
                    status = "homologous_sequence_candidate"
                    if overlap_role in EXON_LIKE_ROLES and not sense_overlap:
                        annotation_status = "homologous_sequence_overlaps_antisense_exon_annotation"
                    elif overlap_role != "unknown":
                        annotation_status = "homologous_sequence_overlaps_non_exonic_annotation"
                    else:
                        annotation_status = "unannotated_homologous_sequence_candidate"
                    inferred_role = "unknown"
                    conclusion = "sequence_present_role_unknown"
                confidence = "medium"
            elif deletion_supported:
                status = "supports_absence"
                annotation_status = "sequence_absence_between_ordered_flanking_homologs"
                inferred_role = "unknown"
                contig = metadata.get("contig") or "supplied_gene_locus"
                start = end = "NA"
                strand = metadata.get("strand") or "+"
                conclusion = "targeted_deletion_supported"
                confidence = "medium"
            else:
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
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "sequence_score": "NA" if status == "supports_absence" else f"{identity:.6g}",
                    "sequence_coverage": "NA" if status == "supports_absence" else f"{coverage:.6g}",
                    "left_synteny_score": "1" if left_element in present_elements[copy_key] else "0",
                    "right_synteny_score": "1" if right_element in present_elements[copy_key] else "0",
                    "splice_motif_score": "NA",
                    "phase_compatibility": projected.get("phase", "unknown") if projected else "unknown",
                    "inferred_event": "protein_cds_projection" if projected else "homologous_exon_sequence_search",
                    "frame_status": "unknown",
                    "evidence_conclusion": conclusion,
                    "confidence_flag": confidence,
                    "absence_evidence": deletion_reason,
                    "alignment_backend": aligner,
                    "alignment_error": alignment_error or "NA",
                    "alignment_cigar": (
                        projected.get("cigar", "NA") if projected else
                        alignment.cigar if alignment is not None else "NA"
                    ),
                    "deletion_alignment_backend": deletion_provenance.get("backend", "NA"),
                    "deletion_alignment_cigar": deletion_provenance.get("cigar", "NA"),
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
        "annotation_status", "evidence_status", "inferred_role", "contig", "start", "end",
        "strand", "sequence_score", "sequence_coverage", "left_synteny_score",
        "right_synteny_score", "splice_motif_score", "phase_compatibility",
        "inferred_event", "frame_status", "evidence_conclusion", "confidence_flag",
        "absence_evidence", "alignment_backend", "alignment_error", "alignment_cigar",
        "deletion_alignment_backend", "deletion_alignment_cigar",
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
            "interval", "annotation_status", "inferred_role", "support_score",
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
