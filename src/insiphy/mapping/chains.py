"""mapping / chains: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.candidate_chain import ChainCandidate
from insiphy.candidate_chain import ChainPathMembership
from insiphy.candidate_chain import DEFAULT_CHAIN_CONFIGURATION
from insiphy.candidate_chain import ordered_candidate_chain
from insiphy.candidate_chain import genomic_candidate_chain
from insiphy.coordinates import ClosedInterval1
from insiphy.coordinates import Interval0
from insiphy.coordinates import local_interval_to_genome
from insiphy.mapping.candidate_codec import _block_signature
from insiphy.mapping.candidate_codec import _public_candidate_record
from insiphy.mapping.policies import occurrence_copy_key
from insiphy.storage.values import to_float
import json


def _copy_transcription_bounds(occurrences):
    bounds = {}
    grouped = defaultdict(list)
    for occurrence in occurrences:
        grouped[occurrence_copy_key(occurrence)].append(occurrence)
    for key, rows in grouped.items():
        bounds[key] = (
            min(int(row["start"]) for row in rows),
            max(int(row["end"]) for row in rows),
        )
    return bounds


def _candidate_copy_interval(occurrence, record, side, copy_bounds):
    genomic_blocks = tuple(record.get(f"{side}_genomic_blocks0", ()) or ())
    if genomic_blocks:
        genomic = Interval0(
            min(interval.start0 for interval in genomic_blocks),
            max(interval.end0 for interval in genomic_blocks),
        )
    else:
        start0 = record.get(f"{side}_start0")
        end0 = record.get(f"{side}_end0")
        if start0 in {None, "", "NA"} or end0 in {None, "", "NA"}:
            return None
        try:
            locus = ClosedInterval1(
                int(occurrence["start"]), int(occurrence["end"]),
            ).to_interval0()
            genomic = local_interval_to_genome(
                Interval0(int(start0), int(end0)), locus, occurrence.get("strand"),
            )
        except (KeyError, TypeError, ValueError):
            return None
    if genomic.length == 0:
        return None
    copy_start, copy_end = copy_bounds[occurrence_copy_key(occurrence)]
    if occurrence.get("strand") == "-":
        return Interval0(copy_end - genomic.end0, copy_end - genomic.start0)
    return Interval0(genomic.start0 - (copy_start - 1), genomic.end0 - (copy_start - 1))


def _chain_precedes(left, right):
    return left.query.end0 <= right.query.start0 and left.target.end0 <= right.target.start0


def _candidate_path_memberships(
    query,
    subject,
    transcript_paths,
    reverse=False,
    query_transcript_ids=None,
    subject_transcript_ids=None,
):
    def transcript_tokens(value):
        return {
            token
            for token in str(value or "").replace(",", ";").split(";")
            if token and token != "NA"
        }

    allowed_query_transcripts = transcript_tokens(query_transcript_ids)
    allowed_subject_transcripts = transcript_tokens(subject_transcript_ids)
    by_occurrence = defaultdict(list)
    for path in transcript_paths or ():
        by_occurrence[path.get("occurrence_id")].append(path)
    memberships = []
    for query_path in by_occurrence.get(query.get("occurrence_id"), ()):
        if (
            allowed_query_transcripts
            and query_path.get("transcript_id") not in allowed_query_transcripts
        ):
            continue
        for subject_path in by_occurrence.get(subject.get("occurrence_id"), ()):
            if (
                allowed_subject_transcripts
                and subject_path.get("transcript_id") not in allowed_subject_transcripts
            ):
                continue
            if query_path.get("transcript_id") in {None, "", "NA"}:
                continue
            if subject_path.get("transcript_id") in {None, "", "NA"}:
                continue
            try:
                query_order = int(query_path.get("path_rank", query_path.get("transcript_order")))
                subject_order = int(subject_path.get("path_rank", subject_path.get("transcript_order")))
            except (TypeError, ValueError):
                continue
            query_membership = {
                "path_id": "|".join(
                    str(query_path.get(field, query.get(field, "NA")))
                    for field in ("family_id", "species", "gene_copy_id", "transcript_id")
                ),
                "order": query_order,
                "parent_id": query_path.get("occurrence_id", query.get("occurrence_id", "")),
                "contig": query_path.get("contig", query.get("contig", "NA")),
                "strand": query_path.get("strand", query.get("strand", "NA")),
            }
            subject_membership = {
                "path_id": "|".join(
                    str(subject_path.get(field, subject.get(field, "NA")))
                    for field in ("family_id", "species", "gene_copy_id", "transcript_id")
                ),
                "order": subject_order,
                "parent_id": subject_path.get("occurrence_id", subject.get("occurrence_id", "")),
                "contig": subject_path.get("contig", subject.get("contig", "NA")),
                "strand": subject_path.get("strand", subject.get("strand", "NA")),
            }
            if reverse:
                query_membership, subject_membership = subject_membership, query_membership
            if (
                query_membership["strand"] not in {"+", "-"}
                or subject_membership["strand"] not in {"+", "-"}
                or query_membership["contig"] in {None, "", "NA"}
                or subject_membership["contig"] in {None, "", "NA"}
            ):
                continue
            memberships.append(
                ChainPathMembership(
                    query_path_id=query_membership["path_id"],
                    target_path_id=subject_membership["path_id"],
                    query_order=query_membership["order"],
                    target_order=subject_membership["order"],
                    query_contig=query_membership["contig"],
                    target_contig=subject_membership["contig"],
                    query_strand=query_membership["strand"],
                    target_strand=subject_membership["strand"],
                    query_parent_id=query_membership["parent_id"],
                    target_parent_id=subject_membership["parent_id"],
                )
            )
    return tuple(sorted(set(memberships), key=lambda item: (item.context, item.query_order, item.target_order)))


def _apply_ordered_candidate_chains(rows, occurrence_by_id, occurrences, transcript_paths=None):
    for row in rows:
        row["membership_edge_eligible"] = 0
        row["position_edge_eligible"] = 0
        row["_membership_edge_eligible"] = False
        row["_position_edge_eligible"] = False
        row["retained_candidate_ids"] = "NA"
        row["best_path_candidate_ids"] = "NA"
        row["chain_best_path_count_capped"] = 0
        row["chain_near_optimal_path_count_capped"] = 0
        row["chain_status"] = "unassessed"
        row["chain_ambiguity"] = "unassessed"
        row["chain_start_anchor_ids"] = "NA"
        row["chain_end_anchor_ids"] = "NA"
        row["chain_retained_edges"] = "NA"
        row["left_anchor_id"] = "NA"
        row["right_anchor_id"] = "NA"
    copy_bounds = _copy_transcription_bounds(occurrences)
    groups = defaultdict(list)
    owner_by_candidate = {}
    candidate_by_id = {}
    record_by_candidate = {}
    original_intervals = {}
    for row in rows:
        if row.get("match_status") != "mapped":
            continue
        query = occurrence_by_id.get(row.get("query_occurrence_id"))
        subject = occurrence_by_id.get(row.get("subject_occurrence_id"))
        if not query or not subject:
            continue
        query_key = occurrence_copy_key(query)
        subject_key = occurrence_copy_key(subject)
        reverse = subject_key < query_key
        group_copies = tuple(sorted((query_key, subject_key)))
        for record in row.get("_candidate_records", ()):
            if record.get("accepted") not in {1, "1", True}:
                continue
            query_interval = _candidate_copy_interval(
                query, record, "query", copy_bounds,
            )
            target_interval = _candidate_copy_interval(
                subject, record, "target", copy_bounds,
            )
            if query_interval is None or target_interval is None:
                continue
            original_query_interval = query_interval
            original_target_interval = target_interval
            if reverse:
                query_interval, target_interval = target_interval, query_interval
            try:
                raw_score = float(record.get("score"))
            except (TypeError, ValueError):
                continue
            score_scheme = str(record.get("score_scheme") or row.get("score_scheme") or "unspecified")
            alignment_strand = str(record.get("strand") or row.get("alignment_strand") or "+")
            candidate = ChainCandidate(
                candidate_id=record["candidate_id"],
                query=query_interval,
                target=target_interval,
                score=raw_score,
                score_scheme=score_scheme,
                relative_strand=alignment_strand if alignment_strand in {"+", "-"} else "+",
                path_memberships=_candidate_path_memberships(
                    query,
                    subject,
                    transcript_paths,
                    reverse=reverse,
                    query_transcript_ids=record.get(
                        "query_transcript_id", row.get("query_transcript_ids"),
                    ),
                    subject_transcript_ids=record.get(
                        "target_transcript_id", row.get("subject_transcript_ids"),
                    ),
                ),
            )
            owner_by_candidate[candidate.candidate_id] = row
            candidate_by_id[candidate.candidate_id] = candidate
            record_by_candidate[candidate.candidate_id] = record
            original_intervals[candidate.candidate_id] = (
                original_query_interval, original_target_interval,
            )
            groups[(*group_copies, score_scheme)].append(candidate)

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

    for candidate_id in physical_retained:
        candidate = candidate_by_id[candidate_id]
        if candidate.path_memberships:
            continue
        dna, delta = physical_metadata[candidate_id]
        retained_ids.add(candidate_id)
        if candidate_id in physical_best:
            best_ids.add(candidate_id)
        chain_metadata[candidate_id] = {
            "status": "genomic_DNA_only", "best_score": f"{dna.best_score:.6g}",
            "score_delta": f"{delta:.6g}", "local_mode": 1,
            "configuration": dna.configuration_name, "ambiguity": dna.ambiguity_status,
            "candidate_ambiguous": candidate_id in dna.ambiguous_ids,
            "start_anchor_ids": "NA", "end_anchor_ids": "NA",
            "retained_edges": ";".join(f"{a}>{b}" for a,b in sorted(dna.retained_edges)) or "NA",
            "best_path_count_capped": dna.best_path_count_capped,
            "near_optimal_path_count_capped": dna.near_optimal_path_count_capped,
        }

    partners_by_query = defaultdict(set)
    partners_by_subject = defaultdict(set)
    for row in rows:
        own_ids = {record["candidate_id"] for record in row.get("_candidate_records", ())}
        retained = own_ids & retained_ids
        best = own_ids & best_ids
        row["genomic_retained_candidate_ids"] = ";".join(sorted(own_ids & physical_retained)) or "NA"
        row["transcript_retained_candidate_ids"] = ";".join(sorted(
            cid for cid in retained if candidate_by_id[cid].path_memberships)) or "NA"
        row["correspondence_channels"] = (
            "DNA_and_annotated_paths" if row["transcript_retained_candidate_ids"] != "NA"
            else "DNA_only_no_observed_transcript_path" if row["genomic_retained_candidate_ids"] != "NA"
            else "unresolved_correspondence")
        row["retained_candidate_ids"] = ";".join(sorted(retained)) or "NA"
        row["best_path_candidate_ids"] = ";".join(sorted(best)) or "NA"
        metadata = [chain_metadata[candidate_id] for candidate_id in retained if candidate_id in chain_metadata]
        if metadata:
            row["chain_best_score"] = max(
                metadata, key=lambda item: to_float(item["best_score"], float("-inf"))
            )["best_score"]
            row["chain_score_delta"] = max(
                metadata, key=lambda item: to_float(item["score_delta"], float("-inf"))
            )["score_delta"]
            row["chain_configuration"] = metadata[0]["configuration"]
            row["chain_delta_rule"] = DEFAULT_CHAIN_CONFIGURATION.delta_rule
            row["chain_local_mode"] = int(any(item["local_mode"] == 1 for item in metadata))
            row["chain_ambiguity"] = (
                "multiple_near_optimal_chains"
                if any(item["candidate_ambiguous"] for item in metadata)
                else "unique_within_reported_candidates"
            )
            row["chain_best_path_count_capped"] = max(
                item["best_path_count_capped"] for item in metadata
            )
            row["chain_near_optimal_path_count_capped"] = max(
                item["near_optimal_path_count_capped"] for item in metadata
            )
            row["chain_start_anchor_ids"] = ";".join(sorted({
                token
                for item in metadata
                for token in item["start_anchor_ids"].split(";")
                if token != "NA"
            })) or "NA"
            row["chain_end_anchor_ids"] = ";".join(sorted({
                token
                for item in metadata
                for token in item["end_anchor_ids"].split(";")
                if token != "NA"
            })) or "NA"
            row["left_anchor_id"] = ";".join(sorted({
                anchor_id
                for candidate_id in retained
                for anchor_id in left_flank_ids_by_candidate[candidate_id]
            })) or "NA"
            row["right_anchor_id"] = ";".join(sorted({
                anchor_id
                for candidate_id in retained
                for anchor_id in right_flank_ids_by_candidate[candidate_id]
            })) or "NA"
            row["chain_retained_edges"] = ";".join(sorted({
                token
                for item in metadata
                for token in item["retained_edges"].split(";")
                if token != "NA"
            })) or "NA"
            row["chain_status"] = (
                "best_path_member" if best else "retained_near_optimal"
            )
        elif own_ids & set(chain_metadata):
            own_statuses = {
                chain_metadata[candidate_id]["status"]
                for candidate_id in own_ids
                if candidate_id in chain_metadata
            }
            row["chain_status"] = (
                "incompatible_transcript_paths"
                if "incompatible_transcript_paths" in own_statuses
                else "noncollinear_candidate"
                if "noncollinear_candidate" in own_statuses
                else "outside_near_optimal_chain"
            )
        anchor_states = {
            anchor_status_by_candidate[candidate_id]
            for candidate_id in retained
            if candidate_id in anchor_status_by_candidate
        }
        if "ordered_double_sided_homologous_flanks_same_path" in anchor_states:
            row["flanking_anchor_status"] = "ordered_double_sided_homologous_flanks_same_path"
        elif "ordered_one_sided_independent_homologous_flank" in anchor_states:
            row["flanking_anchor_status"] = "ordered_one_sided_independent_homologous_flank"
        else:
            row["flanking_anchor_status"] = "no_independent_homologous_flanks_on_same_path"
        for record in row.get("_candidate_records", ()):
            if record.get("candidate_id") not in retained:
                continue
            candidate_id = record["candidate_id"]
            record["left_anchor_id"] = ";".join(sorted(
                left_flank_ids_by_candidate[candidate_id]
            )) or "NA"
            record["right_anchor_id"] = ";".join(sorted(
                right_flank_ids_by_candidate[candidate_id]
            )) or "NA"
            record["chain_configuration"] = row.get("chain_configuration", DEFAULT_CHAIN_CONFIGURATION.name)
            record["chain_score_delta"] = row.get("chain_score_delta", "NA")
        public_candidates = [
            _public_candidate_record(record)
            for record in row.get("_candidate_records", ())
        ]
        row["candidate_assessments"] = json.dumps(
            public_candidates, sort_keys=True, separators=(",", ":"),
        )
        row["alternative_hits"] = json.dumps(
            public_candidates[1:], sort_keys=True, separators=(",", ":"),
        )
        if retained and row.get("match_status") == "mapped":
            query_id = row["query_occurrence_id"]
            subject_id = row["subject_occurrence_id"]
            query_copy = occurrence_copy_key(occurrence_by_id[query_id])
            subject_copy = occurrence_copy_key(occurrence_by_id[subject_id])
            partners_by_query[(query_id, subject_copy)].add(subject_id)
            partners_by_subject[(subject_id, query_copy)].add(query_id)

    def compatible_partner_projections(shared_id, other_copy, shared_side):
        intervals = []
        for candidate_id in retained_ids:
            owner = owner_by_candidate[candidate_id]
            if candidate_id not in original_intervals:
                continue
            query_id = owner["query_occurrence_id"]
            subject_id = owner["subject_occurrence_id"]
            if shared_side == "query":
                if query_id != shared_id or occurrence_copy_key(occurrence_by_id[subject_id]) != other_copy:
                    continue
                intervals.append(original_intervals[candidate_id][0])
            else:
                if subject_id != shared_id or occurrence_copy_key(occurrence_by_id[query_id]) != other_copy:
                    continue
                intervals.append(original_intervals[candidate_id][1])
        intervals.sort()
        return all(left.end0 <= right.start0 for left, right in zip(intervals, intervals[1:]))

    for row in rows:
        retained = set(str(row.get("retained_candidate_ids", "NA")).split(";")) - {"NA", ""}
        # Direct callers and legacy rows may provide accepted candidate
        # records without the chain annotation pass.  Their membership still
        # carries evidence; coordinate eligibility is decided below from the
        # number and identity of placements.
        if not retained:
            retained = {
                record.get("candidate_id")
                for record in row.get("_candidate_records", ())
                if record.get("candidate_id") and record.get("accepted") in {1, "1", True}
            }
        query_id = row["query_occurrence_id"]
        subject_id = row["subject_occurrence_id"]
        query_copy = occurrence_copy_key(occurrence_by_id[query_id])
        subject_copy = occurrence_copy_key(occurrence_by_id[subject_id])
        query_partners = partners_by_query[(query_id, subject_copy)]
        subject_partners = partners_by_subject[(subject_id, query_copy)]
        query_partner_compatible = (
            len(query_partners) <= 1
            or compatible_partner_projections(query_id, subject_copy, "query")
        )
        subject_partner_compatible = (
            len(subject_partners) <= 1
            or compatible_partner_projections(subject_id, query_copy, "subject")
        )
        placement_signatures = {
            (
                original_intervals[candidate_id][0],
                original_intervals[candidate_id][1],
                tuple(
                    _block_signature(block)
                    for block in record_by_candidate[candidate_id].get("aligned_blocks", ())
                ),
            )
            for candidate_id in retained
            if candidate_id in record_by_candidate and candidate_id in original_intervals
        }
        unique_position = len(placement_signatures) == 1
        enumeration_complete = row.get("enumeration_complete") in {1, "1", True}
        chain_unambiguous = (
            row.get("chain_ambiguity", "unique_within_reported_candidates")
            == "unique_within_reported_candidates"
        )
        double_flanks = (
            row.get("flanking_anchor_status")
            == "ordered_double_sided_homologous_flanks_same_path"
        )
        protein_hard = row.get("protein_hard_observation_eligible") in {1, "1", True}
        protein_position = row.get("protein_position_eligible") in {1, "1", True}
        short_route = row.get("short_context_route")
        short_context_supported = short_route not in {
            "feature_bounded_candidate", "bounded_local", "anchor_bounded_local",
            "anchor_bounded_unavailable",
        } or (
            short_route == "anchor_bounded_local" and double_flanks
        )
        chain_membership = bool(
            retained
            and query_partner_compatible
            and subject_partner_compatible
            and short_context_supported
        )
        protein_basis = (
            "annotated_CDS_protein"
            in str(row.get("correspondence_basis", ""))
        )
        membership_eligible = bool(
            row.get("match_status") == "mapped"
            and (chain_membership or protein_hard)
            and (not protein_basis or protein_hard)
        )
        if row.get("match_status") != "mapped":
            membership_reason = "sequence_correspondence_not_accepted"
        elif protein_basis and not protein_hard:
            membership_reason = "protein_candidate_without_hard_coordinates"
        elif protein_hard and not chain_membership:
            membership_reason = "resolved_annotated_CDS_protein_membership"
        elif not retained:
            membership_reason = "no_retained_accepted_candidate"
        elif not unique_position:
            membership_reason = "multiple_accepted_retained_coordinate_placements"
        elif not chain_unambiguous:
            membership_reason = "multiple_near_optimal_candidate_chains"
        elif not query_partner_compatible or not subject_partner_compatible:
            membership_reason = "overlapping_partner_projections"
        elif not short_context_supported:
            membership_reason = "short_context_without_same_path_independent_double_flanks"
        else:
            membership_reason = "retained_sequence_membership"

        position_eligible = bool(
            membership_eligible
            and enumeration_complete
            and unique_position
            and (len(retained) >= 1 or (protein_basis and protein_position))
            and (
                not protein_basis
                or protein_position
            )
        )
        if not membership_eligible:
            position_reason = membership_reason
        elif not enumeration_complete:
            position_reason = "candidate_enumeration_unassessed_or_incomplete"
        elif not unique_position:
            position_reason = "multiple_accepted_retained_coordinate_placements"
        elif (
            protein_basis
            and not protein_position
        ):
            position_reason = "protein_membership_without_resolved_coordinates"
        else:
            position_reason = "unique_resolved_actual_coordinates"

        row["membership_edge_eligible"] = int(membership_eligible)
        row["membership_edge_reason"] = membership_reason
        row["position_edge_eligible"] = int(position_eligible)
        row["position_edge_reason"] = position_reason
        if row.get("match_status") != "mapped":
            row["candidate_resolution"] = "candidate"
        elif protein_hard and not retained:
            row["candidate_resolution"] = "ambiguous"
        elif not retained:
            row["candidate_resolution"] = "excluded"
        elif not position_eligible:
            row["candidate_resolution"] = "ambiguous"
        else:
            row["candidate_resolution"] = "resolved"
        row["_membership_edge_eligible"] = membership_eligible
        row["_position_edge_eligible"] = position_eligible
        if row.get("match_status") == "mapped" and not membership_eligible:
            if not retained:
                row["match_status"] = "candidate_chain_excluded"
            elif not short_context_supported:
                row["match_status"] = "candidate_unanchored"
            else:
                row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
        elif row.get("match_status") == "mapped" and not enumeration_complete:
            row["match_status"] = "candidate_search_incomplete"
            row["candidate_resolution"] = "ambiguous"

        row["true_absence_eligible"] = 0
        if row.get("alignment_backend") == "genomic_overlap":
            row["true_absence_evidence_status"] = "not_applicable"
            row["true_absence_reason"] = "same_locus_annotation_overlap_is_not_deletion_evidence"
            continue
        absence_reasons = []
        if not double_flanks:
            absence_reasons.append("same_path_independent_double_flanks_not_established")
        if not enumeration_complete:
            absence_reasons.append("acceptable_alternative_alignment_set_unassessed_or_incomplete")
        if not unique_position:
            absence_reasons.append("acceptable_alternatives_do_not_define_one_position")
        absence_reasons.extend(
            [
                "anchor_interval_sequence_not_extracted",
                "assembly_continuity_unassessed",
                "ambiguous_base_status_unassessed",
                "query_only_deletion_gap_unassessed",
                "alternative_alignment_concordance_unassessed",
            ]
        )
        row["true_absence_evidence_status"] = (
            "evidence_candidate" if double_flanks else "insufficient_evidence"
        )
        row["true_absence_reason"] = ";".join(absence_reasons)
