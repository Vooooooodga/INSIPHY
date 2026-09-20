"""Exact interval membership refinement and its coordinate evidence."""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from insiphy.coordinates import ClosedInterval1, Interval0, parse_legacy_blocks
from insiphy.elements import collect_element_profiles, element_class_for_occurrence, element_role_from_occurrences, membership_call
from insiphy.alignment import global_alignment_stats
from insiphy.io import norm_state, parse_fasta, read_tsv, to_float, uniq, write_tsv
from insiphy.topology import SpeciesTree



CORRESPONDENCE_EVIDENCE_FIELDS = (
    "matched_blocks",
    "query_genomic_matched_blocks",
    "subject_genomic_matched_blocks",
    "query_parent_feature_ids",
    "subject_parent_feature_ids",
    "query_transcript_ids",
    "subject_transcript_ids",
    "candidate_id",
    "alternative_candidate_ids",
    "aligned_blocks",
    "gap_blocks",
    "sequence_kind",
    "backend",
    "backend_version",
    "raw_score",
    "nt_identity",
    "aa_identity",
    "known_aligned_pairs",
    "unknown_aligned_pairs",
    "query_covered_bases",
    "target_covered_bases",
    "query_length",
    "target_length",
    "relative_strand",
    "is_secondary",
    "left_anchor_id",
    "right_anchor_id",
    "search_interval",
    "search_interval_side",
    "short_sequence_coverage",
    "alignment_input_transposed",
    "mapq",
    "mapping_quality",
    "hit_count",
    "ambiguous_hit_count",
    "alternative_hits",
    "candidate_assessments",
    "dna_candidate_assessments",
    "candidate_ids",
    "score_scheme",
    "raw_alignment_score",
    "enumeration_complete",
    "candidate_enumeration_status",
    "incomplete_reason",
    "short_context_route",
    "local_boundary_range",
    "candidate_resolution",
    "retained_candidate_ids",
    "best_path_candidate_ids",
    "chain_best_score",
    "chain_score_delta",
    "chain_configuration",
    "chain_delta_rule",
    "chain_local_mode",
    "chain_status",
    "chain_ambiguity",
    "chain_start_anchor_ids",
    "chain_end_anchor_ids",
    "chain_retained_edges",
    "chain_best_path_count_capped",
    "chain_near_optimal_path_count_capped",
    "flanking_anchor_status",
    "membership_edge_eligible",
    "membership_edge_reason",
    "position_edge_eligible",
    "position_edge_reason",
    "true_absence_eligible",
    "true_absence_evidence_status",
    "true_absence_reason",
    "protein_mapping_status",
    "protein_candidate_evidence_available",
    "protein_membership_eligible",
    "protein_position_eligible",
    "protein_hard_observation_eligible",
    "protein_known_aa_pairs",
    "protein_blosum62_score",
    "protein_gap_fraction",
    "protein_left_anchor_pairs",
    "protein_right_anchor_pairs",
    "protein_left_anchor_score",
    "protein_right_anchor_score",
    "protein_left_anchor_supported",
    "protein_right_anchor_supported",
    "protein_terminal_side",
    "protein_msa_column_start0",
    "protein_msa_column_end0",
    "protein_msa_column_interval",
    "protein_msa_mode",
    "protein_candidate_mapping_count",
    "protein_candidate_coordinate_consensus",
    "protein_candidate_details",
    "protein_competing_occurrences",
    "protein_query_source_features",
    "protein_target_source_features",
)


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


def _membership_match_details(element_rows, occurrences, scored):
    """Replace whole-feature membership rows with resolved local memberships."""

    occurrence_by_id = {row["occurrence_id"]: row for row in occurrences}
    templates_by_occurrence = defaultdict(dict)
    for row in element_rows:
        templates_by_occurrence[row["occurrence_id"]].setdefault(
            row["element_id"], dict(row),
        )
    base_elements_by_occurrence = {
        occurrence_id: set(templates)
        for occurrence_id, templates in templates_by_occurrence.items()
    }
    parent = {}

    def find(item):
        parent.setdefault(item, item)
        if parent[item] != item:
            parent[item] = find(parent[item])
        return parent[item]

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    node_info = {}
    nodes_by_occurrence = defaultdict(list)
    matches_by_occurrence = defaultdict(list)
    for match in scored:
        query_id = match.get("query_occurrence_id")
        subject_id = match.get("subject_occurrence_id")
        matches_by_occurrence[query_id].append((match, "query"))
        matches_by_occurrence[subject_id].append((match, "subject"))
        if not _hard_position_eligible(match):
            continue
        shared_base_elements = (
            base_elements_by_occurrence.get(query_id, set())
            & base_elements_by_occurrence.get(subject_id, set())
        )
        if not shared_base_elements:
            continue
        query_blocks = _parse_genomic_blocks(match.get("query_genomic_matched_blocks"))
        subject_blocks = _parse_genomic_blocks(match.get("subject_genomic_matched_blocks"))
        if not query_blocks or not subject_blocks:
            continue
        for base_element in sorted(shared_base_elements):
            query_node = (match["match_id"], "query", base_element)
            subject_node = (match["match_id"], "subject", base_element)
            node_info[query_node] = {
                "occurrence_id": query_id,
                "base_element": base_element,
                "blocks": query_blocks,
                "match": match,
                "side": "query",
            }
            node_info[subject_node] = {
                "occurrence_id": subject_id,
                "base_element": base_element,
                "blocks": subject_blocks,
                "match": match,
                "side": "subject",
            }
            nodes_by_occurrence[(query_id, base_element)].append(query_node)
            nodes_by_occurrence[(subject_id, base_element)].append(subject_node)
            union(query_node, subject_node)

    for nodes in nodes_by_occurrence.values():
        for index, left in enumerate(nodes):
            for right in nodes[index + 1:]:
                if _blocks_overlap(node_info[left]["blocks"], node_info[right]["blocks"]):
                    union(left, right)

    components = defaultdict(list)
    for node in node_info:
        components[find(node)].append(node)
    components_by_base = defaultdict(list)
    for nodes in components.values():
        bases = {
            node_info[node]["base_element"]
            for node in nodes
        } - {None}
        if len(bases) == 1:
            components_by_base[next(iter(bases))].append(nodes)

    def component_key(nodes):
        coordinates = [
            (contig, interval.start0, interval.end0, node_info[node]["occurrence_id"])
            for node in nodes
            for contig, _strand, interval in node_info[node]["blocks"]
        ]
        return min(coordinates) if coordinates else ("", 0, 0, "")

    localized_rows = []
    rows_by_occurrence = defaultdict(list)
    for base_element, base_components in sorted(components_by_base.items()):
        ordered_components = sorted(base_components, key=component_key)
        for component_index, nodes in enumerate(ordered_components, start=1):
            element_id = (
                base_element if len(ordered_components) == 1
                else f"{base_element}.local_{component_index:03d}"
            )
            by_occurrence = defaultdict(list)
            for node in nodes:
                by_occurrence[node_info[node]["occurrence_id"]].append(node)
            for occurrence_index, (occurrence_id, occurrence_nodes) in enumerate(
                sorted(by_occurrence.items()), start=1,
            ):
                template = dict(
                    templates_by_occurrence[occurrence_id][base_element]
                )
                occurrence = occurrence_by_id.get(occurrence_id, {})
                blocks = sorted({
                    (contig, strand, interval)
                    for node in occurrence_nodes
                    for contig, strand, interval in node_info[node]["blocks"]
                }, key=lambda item: (item[0], item[2].start0, item[2].end0, item[1]))
                block_tokens = []
                actual_blocks = []
                for block_index, (contig, strand, interval) in enumerate(blocks, start=1):
                    public = ClosedInterval1.from_interval0(interval)
                    block_tokens.append(
                        f"{contig}:{public.start}-{public.end}:{strand}"
                    )
                    actual_blocks.append({
                        "block_id": f"{element_id}.{occurrence_id}.block_{block_index:03d}",
                        "target_contig": contig,
                        "target_start": public.start,
                        "target_end": public.end,
                        "target_strand": strand,
                    })
                matches = [node_info[node]["match"] for node in occurrence_nodes]
                candidate_ids = {
                    candidate_id
                    for match in matches
                    for candidate_id in (
                        _tokens(match.get("candidate_ids"))
                        | _tokens(match.get("retained_candidate_ids"))
                    )
                }
                template.update({
                    "element_id": element_id,
                    "parent_element_id": base_element,
                    "member_interval_id": f"{element_id}.{occurrence_id}.member_{occurrence_index:03d}",
                    "matched_blocks": _tagged_membership_blocks(
                        occurrence_nodes, node_info,
                    ),
                    "genomic_matched_blocks": ";".join(block_tokens) or "NA",
                    "actual_matched_blocks": json.dumps(
                        actual_blocks, sort_keys=True, separators=(",", ":"),
                    ),
                    "parent_feature_ids": occurrence.get("source_feature_id", "NA"),
                    "transcript_ids": occurrence.get("transcript_id", "NA"),
                    "path_roles": occurrence.get("path_roles", "unknown"),
                    "coding_roles": occurrence.get("coding_roles", "unknown"),
                    "position_roles": occurrence.get("position_roles", "unknown"),
                    "annotation_source": occurrence.get("annotation_source", "NA"),
                    "original_attributes": occurrence.get("original_attributes", "NA"),
                    "partial_start": occurrence.get("partial_start", "NA"),
                    "partial_end": occurrence.get("partial_end", "NA"),
                    "path_role_records": occurrence.get("path_role_records", "[]"),
                    "repeat_instance_id": "NA",
                    "correspondence_status": "resolved",
                    "candidate_ids": ";".join(sorted(candidate_ids)) or "NA",
                    "membership_edge_eligible": 1,
                    "membership_edge_reason": "resolved_local_matched_subinterval",
                    "position_edge_eligible": 1,
                    "position_edge_reason": "unique_resolved_actual_coordinates",
                    "reference_occurrence_id": "NA",
                    "reference_length": "NA",
                    "reference_length_source": "localized_matched_block_union",
                    "reference_coverage_relation": "localized_actual_blocks",
                    "reference_covered_bases": _block_union_length(blocks),
                    "reference_overlap_bases": 0,
                    "reference_uncovered_bases": "NA",
                    "_interval_records": tuple(blocks),
                })
                localized_rows.append(template)
                rows_by_occurrence[occurrence_id].append(template)

    localized_keys = {
        (row["occurrence_id"], row["parent_element_id"])
        for row in localized_rows
    }
    for occurrence_id, templates in sorted(templates_by_occurrence.items()):
        for base_element, template_source in sorted(templates.items()):
            if (occurrence_id, base_element) in localized_keys:
                continue
            occurrence = occurrence_by_id.get(occurrence_id, {})
            relevant = matches_by_occurrence.get(occurrence_id, ())
            candidate_ids = {
                candidate_id
                for match, _side in relevant
                for candidate_id in (
                    _tokens(match.get("candidate_ids"))
                    | _tokens(match.get("retained_candidate_ids"))
                )
            }
            ambiguous = any(
                match.get("candidate_resolution") == "ambiguous"
                or match.get("match_status") == "candidate_ambiguous"
                or match.get("protein_mapping_status") == "ambiguous_mapping"
                for match, _side in relevant
            )
            template = dict(template_source)
            membership_supported = any(
                _hard_membership_eligible(match) for match, _side in relevant
            )
            template.update({
                "parent_element_id": base_element,
                "member_interval_id": "NA",
                "matched_blocks": "NA",
                "genomic_matched_blocks": "NA",
                "actual_matched_blocks": "[]",
                "parent_feature_ids": occurrence.get("source_feature_id", "NA"),
                "transcript_ids": occurrence.get("transcript_id", "NA"),
                "path_roles": occurrence.get("path_roles", "unknown"),
                "coding_roles": occurrence.get("coding_roles", "unknown"),
                "position_roles": occurrence.get("position_roles", "unknown"),
                "annotation_source": occurrence.get("annotation_source", "NA"),
                "original_attributes": occurrence.get("original_attributes", "NA"),
                "partial_start": occurrence.get("partial_start", "NA"),
                "partial_end": occurrence.get("partial_end", "NA"),
                "path_role_records": occurrence.get("path_role_records", "[]"),
                "repeat_instance_id": "NA",
                "correspondence_status": "ambiguous" if ambiguous else "candidate",
                "candidate_ids": ";".join(sorted(candidate_ids)) or "NA",
                "membership_edge_eligible": int(membership_supported),
                "membership_edge_reason": (
                    "retained_membership_with_unresolved_local_position"
                    if membership_supported else "local_matched_subinterval_unresolved"
                ),
                "position_edge_eligible": 0,
                "position_edge_reason": "actual_local_coordinates_unavailable",
                "reference_occurrence_id": "NA",
                "reference_length": "NA",
                "reference_length_source": "localized_matched_block_union",
                "reference_coverage_relation": "uncovered",
                "reference_covered_bases": 0,
                "reference_overlap_bases": 0,
                "reference_uncovered_bases": "NA",
                "_interval_records": tuple(),
            })
            localized_rows.append(template)
            rows_by_occurrence[occurrence_id].append(template)

    for match in scored:
        if _hard_position_eligible(match):
            continue
        if match.get("match_status") in {
            "candidate_chain_excluded", "candidate_low_similarity", "low_similarity",
        }:
            continue
        if not (
            str(match.get("match_status", "")).startswith("candidate_")
            or match.get("candidate_resolution") in {"candidate", "ambiguous"}
            or _edge_eligible(match, "protein_candidate_evidence_available", False)
        ):
            continue
        for occurrence_id, side in (
            (match.get("query_occurrence_id"), "query"),
            (match.get("subject_occurrence_id"), "subject"),
        ):
            for candidate_blocks, local_candidate_ids in _candidate_block_evidence(
                match, side,
            ):
                for row in rows_by_occurrence.get(occurrence_id, ()):
                    resolved_blocks = row.get("_interval_records", ())
                    if resolved_blocks and not _blocks_overlap(
                        resolved_blocks, candidate_blocks,
                    ):
                        continue
                    row["correspondence_status"] = "ambiguous"
                    if not _hard_membership_eligible(match):
                        row["membership_edge_eligible"] = 0
                    row["position_edge_eligible"] = 0
                    if not _hard_membership_eligible(match):
                        row["membership_edge_reason"] = "overlapping_competing_candidate"
                    row["position_edge_reason"] = "overlapping_competing_candidate"
                    row["candidate_ids"] = ";".join(sorted(
                        _tokens(row.get("candidate_ids"))
                        | local_candidate_ids
                        | _tokens(match.get("retained_candidate_ids"))
                    )) or "NA"

    by_element_copy = defaultdict(list)
    by_element = defaultdict(list)
    for row in localized_rows:
        by_element[row["element_id"]].append(row)
        by_element_copy[
            (row["element_id"], row.get("species", "NA"), row.get("gene_copy_id", "NA"))
        ].append(row)
    for element_id, rows in by_element.items():
        occurrence_ids = {row["occurrence_id"] for row in rows}
        reference_id = max(
            sorted(occurrence_ids),
            key=lambda item: int(occurrence_by_id.get(item, {}).get("end", 0))
            - int(occurrence_by_id.get(item, {}).get("start", 1)) + 1,
        )
        reference_length = (
            int(occurrence_by_id.get(reference_id, {}).get("end", 0))
            - int(occurrence_by_id.get(reference_id, {}).get("start", 1)) + 1
        )
        intervals_by_occurrence = defaultdict(list)
        for match in scored:
            query_id = match.get("query_occurrence_id")
            subject_id = match.get("subject_occurrence_id")
            if reference_id not in {query_id, subject_id}:
                continue
            try:
                blocks = parse_legacy_blocks(match.get("matched_blocks"))
            except (TypeError, ValueError):
                continue
            member_id = subject_id if query_id == reference_id else query_id
            if member_id not in occurrence_ids:
                continue
            intervals_by_occurrence[member_id].extend(
                block.query if query_id == reference_id else block.target
                for block in blocks
            )
        copy_groups = defaultdict(list)
        for row in rows:
            copy_groups[(row.get("species", "NA"), row.get("gene_copy_id", "NA"))].append(row)
        for copy_rows in copy_groups.values():
            intervals = sorted(
                interval
                for row in copy_rows
                for interval in intervals_by_occurrence.get(row["occurrence_id"], ())
            )
            if not intervals:
                continue
            covered = _block_union_length([
                ("reference", "+", interval) for interval in intervals
            ])
            overlap = sum(
                max(0, left.end0 - right.start0)
                for left, right in zip(intervals, intervals[1:])
            )
            relation = (
                "repeated_overlap" if overlap > 0
                else "complementary_complete" if covered >= reference_length and len(copy_rows) > 1
                else "single_partial"
            )
            for row in copy_rows:
                row["reference_occurrence_id"] = reference_id
                row["reference_length"] = reference_length
                row["reference_length_source"] = "alignment_sequence_length"
                row["reference_coverage_relation"] = relation
                row["reference_covered_bases"] = covered
                row["reference_overlap_bases"] = overlap
                row["reference_uncovered_bases"] = max(0, reference_length - covered)
            if relation == "repeated_overlap":
                for index, row in enumerate(
                    sorted(copy_rows, key=lambda item: item["occurrence_id"]), start=1,
                ):
                    row["repeat_instance_id"] = f"{element_id}.repeat_{index:03d}"
    for (element_id, _species, _copy), rows in by_element_copy.items():
        resolved = [row for row in rows if row.get("_interval_records")]
        if len({row["occurrence_id"] for row in resolved}) < 2:
            continue
        locus_groups = []
        for row in sorted(resolved, key=lambda item: item["occurrence_id"]):
            overlapping_groups = [
                group
                for group in locus_groups
                if any(
                    _blocks_overlap(
                        row["_interval_records"], member["_interval_records"],
                    )
                    for member in group
                )
            ]
            if not overlapping_groups:
                locus_groups.append([row])
                continue
            primary = overlapping_groups[0]
            primary.append(row)
            for group in overlapping_groups[1:]:
                primary.extend(group)
                locus_groups.remove(group)
        if len(locus_groups) < 2:
            continue
        locus_groups.sort(key=lambda group: min(
            (
                contig, interval.start0, interval.end0, row["occurrence_id"],
            )
            for row in group
            for contig, _strand, interval in row["_interval_records"]
        ))
        for index, group in enumerate(locus_groups, start=1):
            repeat_id = f"{element_id}.repeat_{index:03d}"
            for row in group:
                row["repeat_instance_id"] = repeat_id
                row["reference_coverage_relation"] = "repeated_overlap"
    for element_id, rows in by_element.items():
        if any(row.get("reference_length_source") == "alignment_sequence_length" for row in rows):
            reference_ids = {
                row.get("reference_occurrence_id") for row in rows
                if row.get("reference_occurrence_id") not in {None, "", "NA"}
            }
            if len(reference_ids) == 1:
                reference_id = next(iter(reference_ids))
                for row in rows:
                    if (
                        row["occurrence_id"] == reference_id
                        and row.get("reference_coverage_relation") != "repeated_overlap"
                    ):
                        row["reference_coverage_relation"] = "reference_axis"
                continue
        reference = max(
            rows,
            key=lambda row: (int(row.get("reference_covered_bases", 0)), row["occurrence_id"]),
        )
        reference_id = reference["occurrence_id"]
        reference_length = reference.get("reference_covered_bases", 0) or "NA"
        for row in rows:
            row["reference_occurrence_id"] = reference_id
            row["reference_length"] = reference_length
            if (
                row["occurrence_id"] == reference_id
                and row.get("reference_coverage_relation") != "repeated_overlap"
            ):
                row["reference_coverage_relation"] = "reference_axis"

    for row in localized_rows:
        row.pop("_interval_records", None)
    element_rows[:] = sorted(
        localized_rows,
        key=lambda row: (
            row.get("family_id", ""), row.get("element_id", ""),
            row.get("species", ""), row.get("gene_copy_id", ""),
            row.get("occurrence_id", ""), row.get("member_interval_id", ""),
        ),
    )
