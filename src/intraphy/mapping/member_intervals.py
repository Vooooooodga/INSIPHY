"""mapping / member intervals: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import parse_legacy_blocks
import json


def _tokens(value):
    return {
        token
        for token in str(value or "").replace(",", ";").split(";")
        if token and token != "NA"
    }


def _edge_eligible(match, field, legacy_default=False):
    value = match.get(field)
    if value in {1, "1", True}:
        return True
    if value in {0, "0", False}:
        return False
    return bool(legacy_default)


def _hard_membership_eligible(match):
    legacy_resolved = (
        match.get("membership_edge_eligible") in {None, "", "NA"}
        and match.get("match_status") == "mapped"
        and match.get("candidate_resolution", "resolved") == "resolved"
    )
    if not _edge_eligible(match, "membership_edge_eligible", legacy_resolved):
        return False
    if "annotated_CDS_protein" not in str(match.get("correspondence_basis", "")):
        return True
    return _edge_eligible(match, "protein_hard_observation_eligible", False)


def _hard_position_eligible(match):
    legacy_resolved = (
        match.get("position_edge_eligible") in {None, "", "NA"}
        and match.get("match_status") == "mapped"
        and match.get("candidate_resolution", "resolved") == "resolved"
    )
    return (
        _hard_membership_eligible(match)
        and _edge_eligible(match, "position_edge_eligible", legacy_resolved)
    )


def _genomic_block_records(blocks):
    records = []
    for index, token in enumerate(sorted(blocks), start=1):
        try:
            contig, bounds, strand = token.rsplit(":", 2)
            start, end = (int(value) for value in bounds.split("-", 1))
        except (TypeError, ValueError):
            continue
        records.append(
            {
                "block_id": f"membership_block_{index:03d}",
                "target_contig": contig,
                "target_start": min(start, end),
                "target_end": max(start, end),
                "target_strand": strand,
            }
        )
    return records


def _parse_genomic_blocks(value):
    records = []
    for token in str(value or "").split(";"):
        if not token or token == "NA":
            continue
        try:
            contig, bounds, strand = token.rsplit(":", 2)
            start, end = (int(item) for item in bounds.split("-", 1))
            interval = ClosedInterval1(
                min(start, end), max(start, end),
            ).to_interval0()
        except (TypeError, ValueError):
            continue
        records.append((contig, strand, interval))
    return tuple(records)


def _blocks_overlap(left, right):
    return any(
        left_contig == right_contig
        and left_strand == right_strand
        and left_interval.overlaps(right_interval)
        for left_contig, left_strand, left_interval in left
        for right_contig, right_strand, right_interval in right
    )


def _block_union_length(blocks):
    grouped = defaultdict(list)
    for contig, strand, interval in blocks:
        grouped[(contig, strand)].append(interval)
    total = 0
    for intervals in grouped.values():
        intervals = sorted(intervals)
        start0, end0 = intervals[0].start0, intervals[0].end0
        for interval in intervals[1:]:
            if interval.start0 <= end0:
                end0 = max(end0, interval.end0)
            else:
                total += end0 - start0
                start0, end0 = interval.start0, interval.end0
        total += end0 - start0
    return total


def _tagged_membership_blocks(nodes, node_info):
    """Preserve each alignment block and the occurrence side it supports."""

    blocks_by_match = defaultdict(set)
    for node in nodes:
        match = node_info[node]["match"]
        side = node_info[node]["side"]
        try:
            blocks = parse_legacy_blocks(match.get("matched_blocks"))
        except (TypeError, ValueError):
            continue
        for block in blocks:
            query = ClosedInterval1.from_interval0(block.query)
            target = ClosedInterval1.from_interval0(block.target)
            blocks_by_match[(match["match_id"], side)].add(
                (query.start, query.end, target.start, target.end)
            )
    tokens = []
    for (match_id, side), blocks in sorted(blocks_by_match.items()):
        for index, (query_start, query_end, target_start, target_end) in enumerate(
            sorted(blocks)
        ):
            prefix = f"{match_id}:{side}:" if index == 0 else ""
            tokens.append(
                f"{prefix}{query_start}-{query_end}:"
                f"{target_start}-{target_end}"
            )
    return ";".join(tokens) or "NA"


def _candidate_block_evidence(match, side):
    field = (
        "query_genomic_matched_blocks"
        if side == "query"
        else "subject_genomic_matched_blocks"
    )
    evidence = []
    retained_ids = _tokens(match.get("retained_candidate_ids"))
    for encoded_field in ("candidate_assessments", "dna_candidate_assessments"):
        encoded = match.get(encoded_field)
        if encoded is None or encoded == "" or encoded == "NA":
            continue
        try:
            records = json.loads(encoded) if isinstance(encoded, str) else encoded
        except (TypeError, ValueError):
            continue
        if not isinstance(records, list):
            continue
        for record in records:
            if not isinstance(record, dict):
                continue
            if record.get("accepted", 1) not in {1, "1", True}:
                continue
            record_ids = _tokens(record.get("candidate_id"))
            if retained_ids and record_ids.isdisjoint(retained_ids):
                continue
            blocks = _parse_genomic_blocks(record.get(field))
            if not blocks:
                continue
            evidence.append((blocks, record_ids))
    if not evidence:
        top_blocks = _parse_genomic_blocks(match.get(field))
        if top_blocks:
            evidence.append((top_blocks, _tokens(match.get("candidate_ids"))))
    distinct = {}
    for blocks, candidate_ids in evidence:
        signature = tuple(
            (contig, strand, interval.start0, interval.end0)
            for contig, strand, interval in blocks
        )
        distinct.setdefault(signature, [blocks, set()])[1].update(candidate_ids)
    return tuple(
        (blocks, candidate_ids)
        for blocks, candidate_ids in distinct.values()
    )
