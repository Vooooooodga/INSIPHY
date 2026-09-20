"""mapping / member refinement: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.coordinates import ClosedInterval1
from intraphy.mapping.member_intervals import _block_union_length
from intraphy.mapping.member_intervals import _blocks_overlap
from intraphy.mapping.member_intervals import _candidate_block_evidence
from intraphy.mapping.member_intervals import _edge_eligible
from intraphy.mapping.member_intervals import _hard_membership_eligible
from intraphy.mapping.member_intervals import _hard_position_eligible
from intraphy.mapping.member_intervals import _parse_genomic_blocks
from intraphy.mapping.member_intervals import _tagged_membership_blocks
from intraphy.mapping.member_intervals import _tokens
from intraphy.mapping.reference_coverage import _assign_reference_coverage
import json


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
    for rows in by_element.values():
        _assign_reference_coverage(rows, scored, occurrence_by_id)

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
