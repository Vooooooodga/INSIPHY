"""observations / junction coordinates: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.coordinates import parse_legacy_blocks
from intraphy.observations.support import EXONIC_ROLES
from intraphy.observations.support import MISSING
from intraphy.observations.support import _is_true
from intraphy.observations.support import _json_records
from intraphy.observations.support import _tokens
from intraphy.storage.tabular import read_tsv
from pathlib import Path


def _read_match_rows(input_dir, output_dir):
    for base in (Path(output_dir), Path(input_dir)):
        rows = read_tsv(base / "segment_matches.tsv", optional=True)
        if rows:
            return rows
    return []


def _occurrence_length(occurrence):
    try:
        return int(occurrence["end"]) - int(occurrence["start"]) + 1
    except (KeyError, TypeError, ValueError):
        return 0


def _has_resolved_protein_projection(row):
    """Protein coordinates are an independent evidence channel, not DNA candidates."""
    return (
        "annotated_CDS_protein" in str(row.get("correspondence_basis", ""))
        and _is_true(row.get("protein_hard_observation_eligible"))
        and _is_true(row.get("protein_position_eligible"))
        and row.get("protein_mapping_status") in {
            "resolved_local", "resolved_joint_terminal_partition",
        }
        and row.get("protein_projected_blocks") not in MISSING
    )


def _position_eligible(row):
    if _has_resolved_protein_projection(row):
        return True
    if row.get("match_status") != "mapped":
        return False
    explicit_edge = row.get("position_edge_eligible")
    explicit_resolution = row.get("candidate_resolution")
    if explicit_edge not in MISSING or explicit_resolution not in MISSING:
        return _is_true(explicit_edge) and explicit_resolution == "resolved"
    legacy_blocks = row.get("matched_blocks", row.get("projected_reference_blocks"))
    legacy_strand = row.get(
        "alignment_strand", row.get("projected_reference_strand", "NA")
    )
    return legacy_blocks not in MISSING and legacy_strand in {"+", "-"}


def _legacy_projected_reference_blocks(row):
    if not _position_eligible(row):
        return []
    try:
        block_field = (
            row.get("protein_projected_blocks")
            if "annotated_CDS_protein" in str(row.get("correspondence_basis", ""))
            else row.get("matched_blocks", row.get("projected_reference_blocks"))
        )
        blocks = parse_legacy_blocks(block_field)
    except (TypeError, ValueError):
        return []
    return [{"query_start": block.query.start0 + 1, "query_end": block.query.end0,
             "target_start": block.target.start0 + 1, "target_end": block.target.end0,
             "strand": row.get(
                 "alignment_strand", row.get("projected_reference_strand", "+")
             )} for block in blocks]


def _candidate_record_blocks(record, row):
    blocks = []
    for block in record.get("aligned_blocks", ()):
        try:
            query_start, query_end, target_start, target_end = map(int, block)
        except (TypeError, ValueError):
            continue
        if min(query_start, query_end, target_start, target_end) < 1:
            continue
        blocks.append({
            "query_start": min(query_start, query_end),
            "query_end": max(query_start, query_end),
            "target_start": min(target_start, target_end),
            "target_end": max(target_start, target_end),
            "strand": record.get(
                "relative_strand", row.get("alignment_strand", "+"),
            ),
        })
    return blocks


def _retained_candidate_reference_blocks(row):
    """Return one projected block set for every retained alignment candidate."""
    if not _position_eligible(row):
        return []
    if _has_resolved_protein_projection(row):
        # Do not project a resolved coding correspondence through stale DNA hits.
        blocks = _legacy_projected_reference_blocks(row)
        for block in blocks:
            block["strand"] = "+"  # coding positions are in transcript direction
        return [blocks]
    encoded_candidates = row.get("candidate_assessments")
    if encoded_candidates is None or (
        isinstance(encoded_candidates, str)
        and encoded_candidates in {"", "NA", ".", "unknown", "unavailable"}
    ):
        blocks = _legacy_projected_reference_blocks(row)
        return [blocks]
    candidates = _json_records(encoded_candidates)
    retained_ids = _tokens(row.get("retained_candidate_ids"))
    records_by_id = {
        str(record.get("candidate_id")): record
        for record in candidates
        if record.get("candidate_id") not in MISSING
    }
    if not retained_ids:
        return [[]]
    block_sets = []
    for candidate_id in sorted(retained_ids):
        record = records_by_id.get(candidate_id, {})
        block_sets.append(_candidate_record_blocks(record, row))
    return block_sets


def _parse_projected_reference_blocks(row):
    """Return a unique retained projection for position-level inference."""
    block_sets = _retained_candidate_reference_blocks(row)
    if not block_sets or not block_sets[0]:
        return []

    def signature(blocks):
        return tuple(sorted(
            (
                block["query_start"], block["query_end"],
                block["target_start"], block["target_end"], block["strand"],
            )
            for block in blocks
        ))

    first_signature = signature(block_sets[0])
    if any(signature(blocks) != first_signature for blocks in block_sets[1:]):
        return []
    return block_sets[0]


def _project_query_base_to_reference(blocks, position):
    coordinates = {
        (
            block["target_end"] - (position - block["query_start"])
            if block["strand"] == "-"
            else block["target_start"] + position - block["query_start"]
        )
        for block in blocks
        if block["query_start"] <= position <= block["query_end"]
    }
    return next(iter(coordinates)) if len(coordinates) == 1 else None


def _project_target_base_to_reference(blocks, position):
    coordinates = {
        (
            block["query_end"] - (position - block["target_start"])
            if block["strand"] == "-"
            else block["query_start"] + position - block["target_start"]
        )
        for block in blocks
        if block["target_start"] <= position <= block["target_end"]
    }
    return next(iter(coordinates)) if len(coordinates) == 1 else None


def _mapped_reference_coordinate(match_rows, query_occurrence, reference_occurrence, position):
    if query_occurrence == reference_occurrence:
        return position
    coordinates = set()
    for row in match_rows:
        blocks = _parse_projected_reference_blocks(row)
        if row.get("query_occurrence_id") == query_occurrence and row.get("subject_occurrence_id") == reference_occurrence:
            coordinate = _project_query_base_to_reference(blocks, position)
        elif row.get("query_occurrence_id") == reference_occurrence and row.get("subject_occurrence_id") == query_occurrence:
            coordinate = _project_target_base_to_reference(blocks, position)
        else:
            continue
        if coordinate is not None:
            coordinates.add(coordinate)
    return next(iter(coordinates)) if len(coordinates) == 1 else None


def _within_exon_boundary(family, element, left_occurrence, right_occurrence,
                          occ_by_id, reference_by_element, match_rows):
    """Compatibility entry for one exact, position-qualified within-element cutpoint."""
    reference = reference_by_element.get(element)
    left_length = _occurrence_length(occ_by_id.get(left_occurrence, {}))
    if not reference or left_length <= 0:
        return None
    donor = _mapped_reference_coordinate(
        match_rows, left_occurrence, reference, left_length)
    acceptor = _mapped_reference_coordinate(
        match_rows, right_occurrence, reference, 1)
    if donor is None or acceptor is None or acceptor <= donor:
        return None
    return {
        "family_id": family,
        "site_id": f"JG_{element}_REF_{reference}_D{donor}_A{acceptor}",
        "element_id": element,
        "reference_occurrence_id": reference,
        "acceptor_reference_occurrence_id": reference,
        "donor_projection": donor,
        "acceptor_projection": acceptor,
        "cutpoint0": donor,
        "left_occurrence_id": left_occurrence,
        "right_occurrence_id": right_occurrence,
        "projection_method": "unique_position_eligible_actual_blocks",
        "site_kind": "within_element_junction",
        "linked_group_id": "NA",
        "position_edge_eligible": "1",
        "unavailable_reason": "NA",
    }


def _continuous_reference_block_covers(match_rows, query_occurrence, reference_occurrence,
                                       left_coordinate, right_coordinate):
    low, high = sorted((left_coordinate, right_coordinate))
    if query_occurrence == reference_occurrence:
        return True
    conclusions = []
    for row in match_rows:
        if row.get("query_occurrence_id") == query_occurrence and row.get("subject_occurrence_id") == reference_occurrence:
            conclusions.extend(
                any(block["target_start"] <= low and high <= block["target_end"]
                    for block in blocks)
                for blocks in _retained_candidate_reference_blocks(row)
            )
        elif row.get("query_occurrence_id") == reference_occurrence and row.get("subject_occurrence_id") == query_occurrence:
            conclusions.extend(
                any(block["query_start"] <= low and high <= block["query_end"]
                    for block in blocks)
                for blocks in _retained_candidate_reference_blocks(row)
            )
    return bool(conclusions) and all(conclusions)


def _genomically_contiguous(left, right):
    if not left or not right or left.get("contig") in MISSING or left.get("contig") != right.get("contig"):
        return False
    if left.get("strand") not in {"+", "-"} or left.get("strand") != right.get("strand"):
        return False
    try:
        left_start, left_end = int(left["start"]), int(left["end"])
        right_start, right_end = int(right["start"]), int(right["end"])
    except (KeyError, TypeError, ValueError):
        return False
    return right_end + 1 == left_start if left["strand"] == "-" else left_end + 1 == right_start


def _phase_status(left_row, right_row, left, right):
    roles = {left.get("role", left_row.get("role", "unknown")),
             right.get("role", right_row.get("role", "unknown"))}
    if "CDS" not in roles:
        return "not_applicable"
    phases = (str(left_row.get("phase", left.get("phase", "."))),
              str(right_row.get("phase", right.get("phase", "."))))
    return "known" if all(phase in {"0", "1", "2"} for phase in phases) else "unknown"


def _explicit_intron_between(path_rows, left_index, right_index, left, right,
                             intron_rows, transcript_id):
    intervening = path_rows[left_index + 1:right_index]
    if any(row.get("role") != "intron" for row in intervening):
        return False  # retained-marker adjacency is not an RNA splice connection
    if intervening:
        return True
    left_features, right_features = _tokens(left.get("source_feature_id")), _tokens(right.get("source_feature_id"))
    return any(row.get("transcript_id") == transcript_id
               and row.get("left_feature_id") in left_features
               and row.get("right_feature_id") in right_features for row in intron_rows)


def _reference_by_element(element_rows, occ_by_id):
    by_element = defaultdict(list)
    for row in element_rows:
        if row.get("membership_call", "core_member") == "core_member":
            by_element[row.get("element_id")].append(row)
    references = {}
    for element, rows in by_element.items():
        explicit = {row.get("reference_occurrence_id") for row in rows
                    if row.get("reference_occurrence_id") not in MISSING}
        if len(explicit) == 1:
            references[element] = next(iter(explicit))
            continue
        exonic = [row.get("occurrence_id") for row in rows
                  if occ_by_id.get(row.get("occurrence_id"), {}).get("role") in EXONIC_ROLES]
        if exonic:
            references[element] = max(sorted(exonic), key=lambda item: _occurrence_length(occ_by_id.get(item, {})))
    return references


def _member_position_eligible(membership, occurrence_id, reference_id):
    if occurrence_id == reference_id:
        return True
    explicit_edge = membership.get("position_edge_eligible")
    explicit_status = membership.get("correspondence_status")
    if explicit_edge in MISSING and explicit_status in MISSING:
        return True
    return (
        _is_true(explicit_edge)
        and explicit_status == "resolved"
        and membership.get("reference_coverage_relation") != "repeated_overlap"
    )


def _boundary_from_pair(family, left_element, right_element, left_id, right_id,
                        occ_by_id, membership_by_occ, references, match_rows, linked_group):
    left_reference, right_reference = references.get(left_element), references.get(right_element)
    if not left_reference or not right_reference:
        return None, "reference_axis_unavailable"
    if not _member_position_eligible(membership_by_occ.get((left_element, left_id), {}), left_id, left_reference):
        return None, "left_position_mapping_unavailable"
    if not _member_position_eligible(membership_by_occ.get((right_element, right_id), {}), right_id, right_reference):
        return None, "right_position_mapping_unavailable"
    donor = _mapped_reference_coordinate(match_rows, left_id, left_reference, _occurrence_length(occ_by_id.get(left_id, {})))
    acceptor = _mapped_reference_coordinate(match_rows, right_id, right_reference, 1)
    if donor is None or acceptor is None:
        return None, "unique_exact_boundary_projection_unavailable"
    if left_element == right_element and (left_reference != right_reference or acceptor <= donor):
        return None, "boundary_not_an_exact_common_axis_cutpoint"
    if left_element == right_element:
        site_id, site_kind, cutpoint = f"JG_{left_element}_REF_{left_reference}_D{donor}_A{acceptor}", "within_element_junction", donor
    else:
        site_id = f"JG_{left_element}_REF_{left_reference}_D{donor}__{right_element}_REF_{right_reference}_A{acceptor}"
        site_kind, cutpoint = "between_element_junction", "NA"
    return {
        "family_id": family, "site_id": site_id,
        "element_id": left_element if left_element == right_element else f"{left_element};{right_element}",
        "reference_occurrence_id": left_reference,
        "acceptor_reference_occurrence_id": right_reference,
        "donor_projection": donor, "acceptor_projection": acceptor, "cutpoint0": cutpoint,
        "left_occurrence_id": left_id, "right_occurrence_id": right_id,
        "projection_method": "unique_position_eligible_actual_blocks",
        "site_kind": site_kind, "linked_group_id": linked_group,
        "position_edge_eligible": "1", "unavailable_reason": "NA",
    }, None
