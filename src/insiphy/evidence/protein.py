"""evidence / protein: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.aligners.projection import protein_locus_exons
from insiphy.aligners.types import AlignmentBackendError
from insiphy.evidence.fields import _PROTEIN_QUERY_INTERVAL_FIELDS
from insiphy.evidence.projection import _project_locus_interval
from insiphy.evidence.projection import _split_transcript_ids
from insiphy.evidence.search_context import _safe_int
from insiphy.storage.values import to_float


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
