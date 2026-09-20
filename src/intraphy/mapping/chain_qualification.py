"""mapping / chain qualification: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.mapping.candidate_coordinates import _block_signature
from intraphy.mapping.candidate_serialization import _public_candidate_record
from intraphy.mapping.chain_types import DEFAULT_CHAIN_CONFIGURATION
from intraphy.mapping.policies import occurrence_copy_key
from intraphy.storage.values import to_float
import json


def _qualify_chain_edges(rows, retained_ids, best_ids, physical_retained, candidate_by_id, chain_metadata, left_flank_ids_by_candidate, right_flank_ids_by_candidate, anchor_status_by_candidate, occurrence_by_id, original_intervals, owner_by_candidate, record_by_candidate):
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
