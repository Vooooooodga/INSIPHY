"""mapping / flank context: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0
from intraphy.coordinates import genome_interval_to_local
from intraphy.coordinates import local_interval_to_genome
from intraphy.mapping.path_membership import _candidate_path_memberships
from intraphy.mapping.policies import occurrence_copy_key
from intraphy.storage.fasta import parse_fasta
from intraphy.storage.tabular import read_tsv
from pathlib import Path


def _gene_locus_records(input_dir):
    input_dir = Path(input_dir)
    sequences = parse_fasta(input_dir / "gene_loci.fasta")
    metadata = {
        (row.get("species"), row.get("gene_copy_id")): row
        for row in read_tsv(input_dir / "gene_loci.tsv", optional=True)
    }
    records = {}
    for header, sequence in sequences.items():
        try:
            species, gene_copy_id, geometry = header.split("|", 2)
            contig, bounds, strand = geometry.rsplit(":", 2)
            start, end = (int(value) for value in bounds.split("-", 1))
            interval = ClosedInterval1(start, end).to_interval0()
        except (TypeError, ValueError):
            continue
        if strand not in {"+", "-"} or interval.length != len(sequence):
            continue
        row = metadata.get((species, gene_copy_id), {})
        records[(species, gene_copy_id)] = {
            "header": header,
            "sequence": sequence,
            "contig": contig,
            "strand": strand,
            "interval": interval,
            "range_status": row.get("range_status", "NA"),
        }
    return records


def _value_tokens(value):
    return {
        token
        for token in str(value or "").replace(",", ";").split(";")
        if token and token != "NA"
    }


def _candidate_blocks_for_copy(record, owner, copy_key, occurrence_by_id):
    query = occurrence_by_id.get(owner.get("query_occurrence_id"), {})
    subject = occurrence_by_id.get(owner.get("subject_occurrence_id"), {})
    if occurrence_copy_key(query) == copy_key:
        return tuple(record.get("query_genomic_blocks0", ()) or ())
    if occurrence_copy_key(subject) == copy_key:
        return tuple(record.get("target_genomic_blocks0", ()) or ())
    return tuple()


def _owner_occurrence_for_copy(owner, copy_key, occurrence_by_id):
    query = occurrence_by_id.get(owner.get("query_occurrence_id"), {})
    subject = occurrence_by_id.get(owner.get("subject_occurrence_id"), {})
    if occurrence_copy_key(query) == copy_key:
        return query
    if occurrence_copy_key(subject) == copy_key:
        return subject
    return None


def _interval_between_transcript_flanks(left_blocks, right_blocks, strand):
    if not left_blocks or not right_blocks:
        return None
    left = Interval0(
        min(block.start0 for block in left_blocks),
        max(block.end0 for block in left_blocks),
    )
    right = Interval0(
        min(block.start0 for block in right_blocks),
        max(block.end0 for block in right_blocks),
    )
    if strand == "+" and left.end0 <= right.start0:
        return Interval0(left.end0, right.start0)
    if strand == "-" and right.end0 <= left.start0:
        return Interval0(right.end0, left.start0)
    return None


def _cut0_between_loci(cut0, inner_locus, outer_locus, strand):
    cut0 = int(cut0)
    if cut0 < 0 or cut0 > inner_locus.length:
        raise ValueError("cut lies outside the inner locus")
    genome_cut0 = (
        inner_locus.start0 + cut0
        if strand == "+"
        else inner_locus.end0 - cut0
    )
    if genome_cut0 < outer_locus.start0 or genome_cut0 > outer_locus.end0:
        raise ValueError("cut lies outside the parent feature")
    return (
        genome_cut0 - outer_locus.start0
        if strand == "+"
        else outer_locus.end0 - genome_cut0
    )


def _anchor_bounded_gap_blocks(gaps, search_interval, parent_interval, strand):
    projected = []
    for gap in gaps or ():
        if gap.get("gap_in") == "query":
            target = genome_interval_to_local(
                local_interval_to_genome(
                    Interval0(
                        int(gap["target_start0"]), int(gap["target_end0"]),
                    ),
                    search_interval,
                    strand,
                ),
                parent_interval,
                strand,
            )
            projected.append({
                "gap_in": "query",
                "query_cut0": int(gap["query_cut0"]),
                "target_start0": target.start0,
                "target_end0": target.end0,
            })
        elif gap.get("gap_in") == "target":
            projected.append({
                "gap_in": "target",
                "query_start0": int(gap["query_start0"]),
                "query_end0": int(gap["query_end0"]),
                "target_cut0": _cut0_between_loci(
                    gap["target_cut0"], search_interval, parent_interval, strand,
                ),
            })
    return projected


def _context_flank_pairs(
    focus_row,
    source,
    bounded,
    candidate_owner,
    candidate_record,
    occurrence_by_id,
    transcript_paths,
):
    if not transcript_paths:
        return set()
    if (
        source.get("contig") in {None, "", "NA"}
        or bounded.get("contig") in {None, "", "NA"}
        or source.get("strand") not in {"+", "-"}
        or bounded.get("strand") not in {"+", "-"}
    ):
        return set()
    focal_memberships = _candidate_path_memberships(
        source, bounded, transcript_paths,
    )
    if not focal_memberships:
        return set()
    source_copy = occurrence_copy_key(source)
    bounded_copy = occurrence_copy_key(bounded)

    def distinct_placement_ids(candidate_ids):
        by_placement = defaultdict(list)
        for candidate_id in candidate_ids:
            owner = candidate_owner[candidate_id]
            record = candidate_record[candidate_id]
            source_blocks = _candidate_blocks_for_copy(
                record, owner, source_copy, occurrence_by_id,
            )
            bounded_blocks = _candidate_blocks_for_copy(
                record, owner, bounded_copy, occurrence_by_id,
            )
            signature = (
                tuple((block.start0, block.end0) for block in source_blocks),
                tuple((block.start0, block.end0) for block in bounded_blocks),
            )
            by_placement[signature].append(candidate_id)
        return tuple(
            min(candidate_ids)
            for _signature, candidate_ids in sorted(by_placement.items())
        )

    anchors_by_context = defaultdict(list)
    for candidate_id, record in candidate_record.items():
        owner = candidate_owner[candidate_id]
        if owner is focus_row or not owner.get("_position_edge_eligible"):
            continue
        if candidate_id not in _value_tokens(owner.get("retained_candidate_ids")):
            continue
        query = occurrence_by_id.get(owner.get("query_occurrence_id"), {})
        subject = occurrence_by_id.get(owner.get("subject_occurrence_id"), {})
        query_copy = occurrence_copy_key(query)
        subject_copy = occurrence_copy_key(subject)
        if query_copy == source_copy and subject_copy == bounded_copy:
            anchor_source, anchor_bounded = query, subject
            source_transcript = record.get("query_transcript_id")
            bounded_transcript = record.get("target_transcript_id")
        elif subject_copy == source_copy and query_copy == bounded_copy:
            anchor_source, anchor_bounded = subject, query
            source_transcript = record.get("target_transcript_id")
            bounded_transcript = record.get("query_transcript_id")
        else:
            continue
        if (
            anchor_source.get("occurrence_id") == source.get("occurrence_id")
            or anchor_bounded.get("occurrence_id") == bounded.get("occurrence_id")
        ):
            continue
        if (
            anchor_source.get("contig") != source.get("contig")
            or anchor_source.get("strand") != source.get("strand")
            or anchor_bounded.get("contig") != bounded.get("contig")
            or anchor_bounded.get("strand") != bounded.get("strand")
        ):
            continue
        for membership in _candidate_path_memberships(
            anchor_source,
            anchor_bounded,
            transcript_paths,
            query_transcript_ids=source_transcript,
            subject_transcript_ids=bounded_transcript,
        ):
            anchors_by_context[membership.context].append((candidate_id, membership))

    flank_pairs = set()
    for focal in focal_memberships:
        contextual = anchors_by_context.get(focal.context, ())
        left = [
            (candidate_id, membership)
            for candidate_id, membership in contextual
            if membership.query_order < focal.query_order
            and membership.target_order < focal.target_order
        ]
        right = [
            (candidate_id, membership)
            for candidate_id, membership in contextual
            if membership.query_order > focal.query_order
            and membership.target_order > focal.target_order
        ]
        nearest_left = {
            candidate_id
            for candidate_id, membership in left
            if not any(
                (other.query_order >= membership.query_order)
                and (other.target_order >= membership.target_order)
                and (
                    other.query_order > membership.query_order
                    or other.target_order > membership.target_order
                )
                for _other_id, other in left
            )
        }
        nearest_right = {
            candidate_id
            for candidate_id, membership in right
            if not any(
                (other.query_order <= membership.query_order)
                and (other.target_order <= membership.target_order)
                and (
                    other.query_order < membership.query_order
                    or other.target_order < membership.target_order
                )
                for _other_id, other in right
            )
        }
        if nearest_left and nearest_right:
            flank_pairs.add((
                distinct_placement_ids(nearest_left),
                distinct_placement_ids(nearest_right),
            ))
    return flank_pairs
