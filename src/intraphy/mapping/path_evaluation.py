"""mapping / path evaluation: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.mapping.candidate_coordinates import _block_signature
from intraphy.mapping.chain_types import DEFAULT_CHAIN_CONFIGURATION
from intraphy.mapping.ordered_paths import genomic_candidate_chain
from intraphy.mapping.ordered_paths import ordered_candidate_chain
from intraphy.mapping.path_membership import _chain_precedes


def _evaluate_candidate_paths(original_intervals, record_by_candidate, owner_by_candidate, groups, transcript_paths, candidate_by_id):
    retained_ids = set()
    best_ids = set()
    chain_metadata = {}
    anchor_status_by_candidate = {}
    left_flank_ids_by_candidate = defaultdict(set)
    right_flank_ids_by_candidate = defaultdict(set)

    def owner_has_resolved_anchor_position(owner):
        if owner.get("enumeration_complete") not in {1, "1", True}:
            return False
        if owner.get("short_context_route") in {
            "feature_bounded_candidate", "anchor_bounded_unavailable",
        }:
            return False
        if (
            "annotated_CDS_protein" in str(owner.get("correspondence_basis", ""))
            and owner.get("protein_hard_observation_eligible") not in {1, "1", True}
        ):
            return False
        signatures = {
            (
                original_intervals[candidate_id][0],
                original_intervals[candidate_id][1],
                tuple(
                    _block_signature(block)
                    for block in record_by_candidate[candidate_id].get("aligned_blocks", ())
                ),
            )
            for record in owner.get("_candidate_records", ())
            for candidate_id in (record.get("candidate_id"),)
            if record.get("accepted") in {1, "1", True}
            and candidate_id in original_intervals
        }
        return len(signatures) == 1

    def fixed_anchor_ids(group):
        anchors = {}
        contexts = {
            membership.context
            for candidate in group
            for membership in candidate.path_memberships
        }
        for context in contexts:
            start_ids, end_ids = set(), set()
            path_candidates = []
            for candidate in group:
                membership = next(
                    (item for item in candidate.path_memberships if item.context == context),
                    None,
                )
                if membership is None or not owner_has_resolved_anchor_position(
                    owner_by_candidate[candidate.candidate_id]
                ):
                    continue
                path_candidates.append((candidate, membership))
            if len({
                (
                    owner_by_candidate[candidate.candidate_id]["query_occurrence_id"],
                    owner_by_candidate[candidate.candidate_id]["subject_occurrence_id"],
                )
                for candidate, _membership in path_candidates
            }) < 2:
                anchors[context] = (start_ids, end_ids)
                continue
            ordered = sorted(
                path_candidates,
                key=lambda item: (
                    item[1].query_order,
                    item[1].target_order,
                    item[0].query.start0,
                    item[0].target.start0,
                    item[0].candidate_id,
                ),
            )
            first_order = (ordered[0][1].query_order, ordered[0][1].target_order)
            last_order = (ordered[-1][1].query_order, ordered[-1][1].target_order)
            if first_order == last_order:
                anchors[context] = (start_ids, end_ids)
                continue
            start_ids.update(
                candidate.candidate_id
                for candidate, membership in ordered
                if (membership.query_order, membership.target_order) == first_order
            )
            end_ids.update(
                candidate.candidate_id
                for candidate, membership in ordered
                if (membership.query_order, membership.target_order) == last_order
            )
            anchors[context] = (start_ids, end_ids)
        return anchors

    physical_retained = set()
    physical_best = set()
    physical_metadata = {}
    for group in groups.values():
        dna_candidates = [c for c in group if c.relative_strand == "+"]
        if dna_candidates:
            exact_dna = genomic_candidate_chain(dna_candidates, 0.0)
            dna_delta = DEFAULT_CHAIN_CONFIGURATION.score_delta(
                exact_dna.best_score, dna_candidates[0].score_scheme)
            dna = genomic_candidate_chain(dna_candidates, dna_delta)
            physical_retained.update(dna.retained_ids)
            physical_best.update(dna.best_path_member_ids)
            for c in dna_candidates:
                physical_metadata[c.candidate_id] = (dna, dna_delta)
        collinear = [
            candidate
            for candidate in group
            if candidate.relative_strand == "+"
            and (not transcript_paths or candidate.path_memberships)
        ]
        for candidate in group:
            if candidate.relative_strand == "-":
                chain_metadata[candidate.candidate_id] = {
                    "status": "noncollinear_candidate",
                    "best_score": "NA",
                    "score_delta": "NA",
                    "local_mode": "NA",
                }
            elif transcript_paths and not candidate.path_memberships:
                chain_metadata[candidate.candidate_id] = {
                    "status": "incompatible_transcript_paths",
                    "best_score": "NA",
                    "score_delta": "NA",
                    "local_mode": "NA",
                }
        if not collinear:
            continue
        context_anchor_ids = fixed_anchor_ids(collinear)
        exact = ordered_candidate_chain(
            collinear,
            0.0,
            context_anchor_ids=context_anchor_ids,
            configuration_name=DEFAULT_CHAIN_CONFIGURATION.name,
        )
        scheme = collinear[0].score_scheme
        delta = DEFAULT_CHAIN_CONFIGURATION.score_delta(
            exact.best_score, scheme,
        )
        result = ordered_candidate_chain(
            collinear,
            delta,
            context_anchor_ids=context_anchor_ids,
            configuration_name=DEFAULT_CHAIN_CONFIGURATION.name,
        )
        retained_ids.update(result.retained_ids)
        best_ids.update(result.best_path_member_ids)
        for candidate in collinear:
            status = "outside_near_optimal_chain"
            if candidate.candidate_id in result.retained_ids:
                status = "retained_near_optimal"
            if candidate.candidate_id in result.best_path_member_ids:
                status = "best_path_member"
            chain_metadata[candidate.candidate_id] = {
                "status": status,
                "best_score": f"{result.best_score:.6g}",
                "score_delta": f"{result.score_delta:.6g}",
                "local_mode": int(result.local_mode),
                "configuration": result.configuration_name,
                "ambiguity": result.ambiguity_status,
                "candidate_ambiguous": candidate.candidate_id in result.ambiguous_ids,
                "start_anchor_ids": ";".join(sorted(result.start_anchor_ids)) or "NA",
                "end_anchor_ids": ";".join(sorted(result.end_anchor_ids)) or "NA",
                "retained_edges": ";".join(
                    f"{left}>{right}" for left, right in sorted(result.retained_edges)
                ) or "NA",
                "best_path_count_capped": result.best_path_count_capped,
                "near_optimal_path_count_capped": result.near_optimal_path_count_capped,
            }
            focus_owner = owner_by_candidate[candidate.candidate_id]

            def independent_structure_unit(peer):
                peer_owner = owner_by_candidate[peer.candidate_id]
                return (
                    peer_owner["query_occurrence_id"] != focus_owner["query_occurrence_id"]
                    and peer_owner["subject_occurrence_id"] != focus_owner["subject_occurrence_id"]
                    and owner_has_resolved_anchor_position(peer_owner)
                    and peer.candidate_id not in result.ambiguous_ids
                )

            saw_one_sided = False
            saw_double_sided = False
            for summary in result.context_summaries:
                if candidate.candidate_id not in summary.retained_ids:
                    continue
                context_candidates = [
                    candidate_by_id[candidate_id]
                    for candidate_id in summary.retained_ids
                ]

                def path_reachable(source_id, target_id):
                    pending = [source_id]
                    seen = set()
                    while pending:
                        current = pending.pop()
                        if current == target_id:
                            return True
                        if current in seen:
                            continue
                        seen.add(current)
                        pending.extend(
                            right for left, right in summary.retained_edges
                            if left == current
                        )
                    return False

                left_anchors = [
                    peer for peer in context_candidates
                    if peer.candidate_id != candidate.candidate_id
                    and independent_structure_unit(peer)
                    and (
                        path_reachable(peer.candidate_id, candidate.candidate_id)
                        if summary.retained_edges else _chain_precedes(peer, candidate)
                    )
                ]
                right_anchors = [
                    peer for peer in context_candidates
                    if peer.candidate_id != candidate.candidate_id
                    and independent_structure_unit(peer)
                    and (
                        path_reachable(candidate.candidate_id, peer.candidate_id)
                        if summary.retained_edges else _chain_precedes(candidate, peer)
                    )
                ]
                if left_anchors:
                    nearest_left = max(
                        left_anchors,
                        key=lambda item: (
                            item.query.end0, item.target.end0, item.candidate_id,
                        ),
                    )
                    left_flank_ids_by_candidate[candidate.candidate_id].add(
                        nearest_left.candidate_id
                    )
                if right_anchors:
                    nearest_right = min(
                        right_anchors,
                        key=lambda item: (
                            item.query.start0, item.target.start0, item.candidate_id,
                        ),
                    )
                    right_flank_ids_by_candidate[candidate.candidate_id].add(
                        nearest_right.candidate_id
                    )
                saw_one_sided |= bool(left_anchors or right_anchors)
                saw_double_sided |= bool(left_anchors and right_anchors)
            if saw_double_sided:
                anchor_status_by_candidate[candidate.candidate_id] = (
                    "ordered_double_sided_homologous_flanks_same_path"
                )
            elif saw_one_sided:
                anchor_status_by_candidate[candidate.candidate_id] = (
                    "ordered_one_sided_independent_homologous_flank"
                )
            else:
                anchor_status_by_candidate[candidate.candidate_id] = (
                    "no_independent_homologous_flanks_on_same_path"
                )
    return physical_retained, physical_metadata, retained_ids, physical_best, best_ids, chain_metadata, left_flank_ids_by_candidate, right_flank_ids_by_candidate, anchor_status_by_candidate
