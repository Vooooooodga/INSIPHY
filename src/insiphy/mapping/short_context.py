"""mapping / short_context: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.aligners.short import anchored_short_alignment
from insiphy.aligners.types import AlignmentBackendError
from insiphy.coordinates import ClosedInterval1
from insiphy.coordinates import CoordinateBlock
from insiphy.coordinates import Interval0
from insiphy.coordinates import genome_interval_to_local
from insiphy.coordinates import local_interval_to_genome
from insiphy.mapping.candidate_codec import _candidate_record
from insiphy.mapping.candidate_codec import _coordinate_block0
from insiphy.mapping.candidate_codec import _format_alignment_blocks
from insiphy.mapping.candidate_codec import _format_genomic_blocks
from insiphy.mapping.candidate_codec import _genomic_blocks0
from insiphy.mapping.candidate_codec import _public_candidate_record
from insiphy.mapping.candidate_codec import _public_gap_blocks
from insiphy.mapping.candidate_codec import _public_interval
from insiphy.mapping.candidate_codec import _transpose_candidate_record
from insiphy.mapping.candidate_codec import _transpose_gap_blocks
from insiphy.mapping.chains import _candidate_path_memberships
from insiphy.mapping.match_records import _format_contract_value
from insiphy.mapping.match_records import _format_optional_number
from insiphy.mapping.policies import _candidate_sequence_accepted
from insiphy.mapping.policies import occurrence_copy_key
from insiphy.storage.fasta import parse_fasta
from insiphy.storage.tabular import read_tsv
from insiphy.storage.values import to_float
from pathlib import Path
import json


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


def _rerun_anchor_bounded_short_candidates(
    rows,
    occurrence_by_id,
    seqs,
    gene_loci,
    transcript_paths=None,
):
    candidate_owner = {}
    candidate_record = {}
    for owner in rows:
        for record in owner.get("_candidate_records", ()):
            candidate_id = record.get("candidate_id")
            if candidate_id not in {None, "", "NA"}:
                candidate_owner[candidate_id] = owner
                candidate_record[candidate_id] = record

    changed = False
    for row in rows:
        if row.get("short_context_route") != "feature_bounded_candidate":
            continue
        row["membership_edge_eligible"] = 0
        row["position_edge_eligible"] = 0
        row["_membership_edge_eligible"] = False
        row["_position_edge_eligible"] = False
        input_transposed = row.get("alignment_input_transposed") in {1, "1", True}
        bounded_side = "query" if input_transposed else "target"
        source_side = "target" if input_transposed else "query"
        bounded_id = (
            row["query_occurrence_id"] if bounded_side == "query"
            else row["subject_occurrence_id"]
        )
        source_id = (
            row["query_occurrence_id"] if source_side == "query"
            else row["subject_occurrence_id"]
        )
        bounded = occurrence_by_id.get(bounded_id, {})
        source = occurrence_by_id.get(source_id, {})
        retained = [
            record for record in row.get("_candidate_records", ())
            if record.get("candidate_id") in _value_tokens(row.get("retained_candidate_ids"))
        ]
        flank_pairs = {
            (
                tuple(sorted(_value_tokens(record.get("left_anchor_id")))),
                tuple(sorted(_value_tokens(record.get("right_anchor_id")))),
            )
            for record in retained
        }
        flank_pairs.discard((tuple(), tuple()))
        inferred_pairs = _context_flank_pairs(
            row,
            source=source,
            bounded=bounded,
            candidate_owner=candidate_owner,
            candidate_record=candidate_record,
            occurrence_by_id=occurrence_by_id,
            transcript_paths=transcript_paths,
        )
        if transcript_paths:
            flank_pairs = inferred_pairs
        if len(flank_pairs) != 1:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "unique_same_context_flank_pair_unavailable"
            continue
        left_ids, right_ids = next(iter(flank_pairs))
        if len(left_ids) != 1 or len(right_ids) != 1:
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "competing_same_context_flank_pairs"
            continue

        bounded_copy = occurrence_copy_key(bounded)
        locus = gene_loci.get((bounded.get("species"), bounded.get("gene_copy_id")))
        if (
            locus is None
            or locus.get("contig") != bounded.get("contig")
            or locus.get("strand") != bounded.get("strand")
        ):
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "target_gene_locus_sequence_unavailable"
            continue

        left_id, right_id = left_ids[0], right_ids[0]
        left_record = candidate_record.get(left_id)
        right_record = candidate_record.get(right_id)
        left_owner = candidate_owner.get(left_id)
        right_owner = candidate_owner.get(right_id)
        if not left_record or not right_record or not left_owner or not right_owner:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flank_candidate_coordinates_unavailable"
            continue
        flank_occurrences = (
            _owner_occurrence_for_copy(
                left_owner, occurrence_copy_key(source), occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                right_owner, occurrence_copy_key(source), occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                left_owner, bounded_copy, occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                right_owner, bounded_copy, occurrence_by_id,
            ),
        )
        expected_geometry = (
            (source.get("contig"), source.get("strand")),
            (source.get("contig"), source.get("strand")),
            (bounded.get("contig"), bounded.get("strand")),
            (bounded.get("contig"), bounded.get("strand")),
        )
        if any(
            occurrence is None
            or (occurrence.get("contig"), occurrence.get("strand")) != expected
            for occurrence, expected in zip(flank_occurrences, expected_geometry)
        ):
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flanks_are_not_on_matching_contigs_and_gene_strands"
            continue
        left_blocks = _candidate_blocks_for_copy(
            left_record, left_owner, bounded_copy, occurrence_by_id,
        )
        right_blocks = _candidate_blocks_for_copy(
            right_record, right_owner, bounded_copy, occurrence_by_id,
        )
        search_interval = _interval_between_transcript_flanks(
            left_blocks, right_blocks, bounded.get("strand"),
        )
        if (
            search_interval is None
            or search_interval.length == 0
            or search_interval.start0 < locus["interval"].start0
            or search_interval.end0 > locus["interval"].end0
        ):
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flanks_do_not_define_one_contained_genome_interval"
            continue
        local_search = genome_interval_to_local(
            search_interval, locus["interval"], locus["strand"],
        )
        target_sequence = locus["sequence"][local_search.start0:local_search.end0]
        query_sequence = seqs.get(source_id, "")
        search_metadata = {
            "coordinate_system": "0-based-half-open",
            "contig": locus["contig"],
            "start0": search_interval.start0,
            "end0": search_interval.end0,
            "strand": locus["strand"],
        }
        try:
            candidate_set = anchored_short_alignment(
                query_sequence,
                target_sequence,
                mode="local",
                query_occurrence_id=source_id,
                target_occurrence_id=bounded_id,
                query_transcript_id=source.get("transcript_id"),
                target_transcript_id=bounded.get("transcript_id"),
                left_anchor_id=left_id,
                right_anchor_id=right_id,
                search_interval=search_metadata,
            )
        except AlignmentBackendError as error:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["short_context_route"] = "anchor_bounded_unavailable"
            row["incomplete_reason"] = str(error)
            continue

        records = []
        try:
            bounded_locus = ClosedInterval1(
                int(bounded["start"]), int(bounded["end"]),
            ).to_interval0()
        except (KeyError, TypeError, ValueError):
            bounded_locus = None
        for rank, candidate in enumerate(candidate_set.candidates, start=1):
            record = _candidate_record(
                candidate, rank, candidate.backend, candidate.score_scheme,
            )
            source_genomic = _genomic_blocks0(
                source, record["aligned_blocks"], "query",
            )
            bounded_genomic = tuple(
                local_interval_to_genome(
                    block.target, search_interval, locus["strand"],
                )
                for block in record["aligned_blocks"]
            )
            try:
                if bounded_locus is None:
                    raise ValueError("bounded parent feature is unavailable")
                bounded_local = tuple(
                    genome_interval_to_local(block, bounded_locus, bounded.get("strand"))
                    for block in bounded_genomic
                )
                parent_gap_blocks = _anchor_bounded_gap_blocks(
                    record.get("gap_blocks", ()),
                    search_interval,
                    bounded_locus,
                    bounded.get("strand"),
                )
            except (KeyError, TypeError, ValueError):
                bounded_local = tuple()
                parent_gap_blocks = []
            if (
                len(source_genomic) == len(record["aligned_blocks"])
                and len(bounded_local) == len(record["aligned_blocks"])
            ):
                row_blocks = tuple(
                    CoordinateBlock(source_block.query, target_local)
                    for source_block, target_local in zip(
                        record["aligned_blocks"], bounded_local,
                    )
                )
                if input_transposed:
                    row_blocks = tuple(
                        CoordinateBlock(block.target, block.query) for block in row_blocks
                    )
                    parent_gap_blocks = _transpose_gap_blocks(parent_gap_blocks)
            else:
                row_blocks = tuple()
                parent_gap_blocks = []
            if input_transposed:
                record = _transpose_candidate_record(record)
                record["query_genomic_blocks0"] = bounded_genomic
                record["target_genomic_blocks0"] = source_genomic
            else:
                record["query_genomic_blocks0"] = source_genomic
                record["target_genomic_blocks0"] = bounded_genomic
            record["aligned_blocks"] = row_blocks
            record["gap_blocks"] = parent_gap_blocks
            if row_blocks:
                record["query_start0"] = min(block.query.start0 for block in row_blocks)
                record["query_end0"] = max(block.query.end0 for block in row_blocks)
                record["target_start0"] = min(block.target.start0 for block in row_blocks)
                record["target_end0"] = max(block.target.end0 for block in row_blocks)
            record["candidate_id"] = (
                f"{row['match_id']}.anchor_bounded_candidate_{rank:03d}"
            )
            record["left_anchor_id"] = left_id
            record["right_anchor_id"] = right_id
            record["search_interval"] = search_metadata
            record["search_interval_side"] = bounded_side
            record["source"] = "nucleotide_alignment"
            record["short_sequence_coverage"] = candidate.query_coverage
            record["accepted"] = int(bool(
                row_blocks
                and _candidate_sequence_accepted(
                    record, to_float(row.get("threshold"), 0.0), short_context=True,
                )
            ))
            record["acceptance_threshold"] = row.get("threshold", "NA")
            records.append(record)

        final_candidate_ids = [record["candidate_id"] for record in records]
        for record in records:
            record["hit_count"] = len(records)
            record["alternative_candidate_ids"] = tuple(
                candidate_id for candidate_id in final_candidate_ids
                if candidate_id != record["candidate_id"]
            )

        row["_candidate_records"] = records
        row["short_context_route"] = "anchor_bounded_local"
        row["flanking_anchor_status"] = "ordered_double_sided_homologous_flanks_same_path"
        row["left_anchor_id"] = left_id
        row["right_anchor_id"] = right_id
        row["search_interval"] = _format_contract_value(_public_interval(search_metadata))
        row["search_interval_side"] = bounded_side
        row["local_boundary_range"] = row["search_interval"]
        row["enumeration_complete"] = int(candidate_set.enumeration_complete)
        row["candidate_enumeration_status"] = (
            "complete" if candidate_set.enumeration_complete else "incomplete"
        )
        row["incomplete_reason"] = candidate_set.incomplete_reason or "NA"
        row["hit_count"] = len(records)
        row["ambiguous_hit_count"] = max(0, len(records) - 1)
        accepted = [record for record in records if record.get("accepted") == 1]
        row["match_status"] = "mapped" if accepted else "candidate_low_similarity"
        row["candidate_resolution"] = "unassessed"
        primary = records[0] if records else None
        if primary is not None:
            row["alignment_score"] = f"{to_float(primary.get('identity'), 0.0):.6g}"
            row["coverage_score"] = f"{to_float(primary.get('coverage'), 0.0):.6g}"
            row["sequence_score"] = f"{(0.70 * to_float(primary.get('identity'), 0.0) + 0.30 * to_float(primary.get('coverage'), 0.0)):.6g}"
            row["correspondence_score"] = row["sequence_score"]
            row["candidate_id"] = primary["candidate_id"]
            row["candidate_ids"] = ";".join(record["candidate_id"] for record in records)
            row["alternative_candidate_ids"] = ";".join(
                record["candidate_id"] for record in records[1:]
            ) or "NA"
            row["matched_blocks"] = _format_alignment_blocks(primary["aligned_blocks"])
            row["aligned_blocks"] = row["matched_blocks"]
            row["projected_reference_blocks"] = row["matched_blocks"]
            query_blocks = tuple(
                _coordinate_block0(block).query for block in primary["aligned_blocks"]
            )
            target_blocks = tuple(
                _coordinate_block0(block).target for block in primary["aligned_blocks"]
            )
            if query_blocks:
                public_query = ClosedInterval1.from_interval0(Interval0(
                    min(block.start0 for block in query_blocks),
                    max(block.end0 for block in query_blocks),
                ))
                public_target = ClosedInterval1.from_interval0(Interval0(
                    min(block.start0 for block in target_blocks),
                    max(block.end0 for block in target_blocks),
                ))
                row["query_alignment_start"] = public_query.start
                row["query_alignment_end"] = public_query.end
                row["target_alignment_start"] = public_target.start
                row["target_alignment_end"] = public_target.end
            row["query_genomic_matched_blocks"] = _format_genomic_blocks(
                occurrence_by_id[row["query_occurrence_id"]].get("contig", "NA"),
                occurrence_by_id[row["query_occurrence_id"]].get("strand", "NA"),
                primary.get("query_genomic_blocks0", ()),
            )
            row["subject_genomic_matched_blocks"] = _format_genomic_blocks(
                occurrence_by_id[row["subject_occurrence_id"]].get("contig", "NA"),
                occurrence_by_id[row["subject_occurrence_id"]].get("strand", "NA"),
                primary.get("target_genomic_blocks0", ()),
            )
            row["gap_blocks"] = json.dumps(
                _public_gap_blocks(primary.get("gap_blocks", ())),
                sort_keys=True, separators=(",", ":"),
            )
            row["alignment_cigar"] = primary.get("cigar", "NA")
            row["raw_alignment_score"] = _format_optional_number(primary.get("score"))
            row["raw_score"] = row["raw_alignment_score"]
            row["short_sequence_coverage"] = primary.get("short_sequence_coverage", "NA")
        public_records = [_public_candidate_record(record) for record in records]
        row["candidate_assessments"] = json.dumps(
            public_records, sort_keys=True, separators=(",", ":"),
        )
        row["dna_candidate_assessments"] = row["candidate_assessments"]
        row["alternative_hits"] = json.dumps(
            public_records[1:], sort_keys=True, separators=(",", ":"),
        )
        changed = True
    return changed
