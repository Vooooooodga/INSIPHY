"""evidence / projection: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from insiphy.coordinates import ClosedInterval1
from insiphy.coordinates import local_interval_to_genome
from insiphy.elements import EXON_LIKE_ROLES
import json


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


def _split_transcript_ids(value):
    if not value:
        return []
    return [part for part in str(value).replace(",", ";").split(";") if part]


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
