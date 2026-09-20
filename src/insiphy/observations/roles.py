"""observations / roles: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.coordinates import ClosedInterval1
from insiphy.observations.support import EXONIC_ROLES
from insiphy.observations.support import KNOWN_NONEXONIC_ROLES
from insiphy.observations.support import MISSING
from insiphy.observations.support import _blocks_overlap_status
from insiphy.observations.support import _completion_blocks
from insiphy.observations.support import _completion_presence
from insiphy.observations.support import _completion_role_prediction_supported
from insiphy.observations.support import _finalize_site_metadata
from insiphy.observations.support import _is_true
from insiphy.observations.support import _json_records
from insiphy.observations.support import _membership_blocks
from insiphy.observations.support import _observation_row
from insiphy.observations.support import _presence_state
from insiphy.observations.support import _role_from_transcript_paths
from insiphy.observations.support import _tokens


def _element_site_rows(occurrences, element_rows, completion_rows, valid_families,
                       species_by_family, transcript_paths=None, annotation_view="repertoire"):
    if annotation_view not in {"repertoire", "canonical"}:
        raise ValueError("annotation_view must be 'repertoire' or 'canonical'")
    transcript_paths = transcript_paths or []
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    mapped = defaultdict(lambda: defaultdict(list))
    family_by_element = {}
    elements_by_homology = defaultdict(set)
    elements_by_source_occurrence = defaultdict(set)
    for membership in element_rows:
        if membership.get("element_class") not in {"exon_like", "candidate_source", "absent"}:
            continue
        occurrence = occ_by_id.get(membership.get("occurrence_id", ""))
        family = membership.get("family_id") or (occurrence or {}).get("family_id")
        species = membership.get("species") or (occurrence or {}).get("species")
        if family not in valid_families or not species:
            continue
        element = membership["element_id"]
        family_by_element[element] = family
        if membership.get("homology_id") not in MISSING:
            homology_id = membership["homology_id"]
            elements_by_homology[homology_id].add(element)
            elements_by_source_occurrence[
                (homology_id, membership.get("occurrence_id", "NA"))
            ].add(element)
        mapped[element][species].append((occurrence, membership))
    completion_by_element = defaultdict(lambda: defaultdict(list))
    for completion in completion_rows:
        homology_id = completion.get("homology_id")
        source_occurrence = completion.get("source_parent_occurrence_id", "NA")
        elements = elements_by_source_occurrence.get((homology_id, source_occurrence))
        if not elements:
            elements = elements_by_homology.get(homology_id, set())
        elements = {
            element for element in elements
            if family_by_element.get(element) == completion.get("family_id")
        }
        if not elements or completion.get("family_id") not in valid_families or not completion.get("species"):
            continue
        assignment = "unique" if len(elements) == 1 else "ambiguous"
        for element in sorted(elements):
            assigned = dict(completion)
            assigned["completion_element_assignment"] = assignment
            if assignment == "ambiguous":
                assigned["candidate_resolution_status"] = "ambiguous"
                assigned["correspondence_status"] = "ambiguous"
            completion_by_element[element][completion["species"]].append(assigned)

    output = []
    transcript_scope = ("explicit_canonical_transcript" if annotation_view == "canonical"
                        else "annotated_transcript_repertoire")
    for element, family in sorted(family_by_element.items()):
        member_count = sum(len(values) for values in mapped[element].values())
        for species in sorted(species_by_family[family]):
            presence_calls, role_calls = set(), []
            evidence, presence_evidence = set(), set()
            parent_ids, interval_ids, evidence_ids = set(), set(), set()
            local_blocks, completeness = [], set()
            legacy_role_occurrences = []
            unknown_role = False
            established_exonic_path = False
            for occurrence, membership in mapped[element].get(species, []):
                evidence_ids.update(_tokens(membership.get("candidate_ids")))
                parent_ids.update(_tokens(membership.get("parent_feature_ids")))
                blocks, block_status = _membership_blocks(membership)
                local_blocks.extend(blocks)
                interval_ids.update(block["block_id"] for block in blocks)
                if not blocks:
                    completeness.add(block_status)
                call = membership.get("membership_call", "core_member")
                eligible = _is_true(membership.get("membership_edge_eligible"))
                legacy = membership.get("membership_edge_eligible") in MISSING
                singleton = call != "core_member" and member_count == 1
                if not legacy and not eligible:
                    evidence.add("ambiguous_correspondence")
                    evidence.add("own_annotation_not_used_as_crossspecies_EG_state")
                    unknown_role = True
                    presence_evidence.add("crossspecies_membership_uncertain")
                    continue
                if legacy and call != "core_member":
                    evidence.add("ambiguous_correspondence")
                    if not singleton:
                        unknown_role = True
                        presence_evidence.add("crossspecies_membership_uncertain")
                        evidence.add("own_annotation_not_used_as_crossspecies_EG_state")
                        continue
                position_status = membership.get("correspondence_status")
                position_unresolved = position_status not in MISSING and position_status != "resolved"
                if position_unresolved:
                    unknown_role = True
                    evidence.add("local_position_correspondence_unresolved")
                if not occurrence:
                    unknown_role = True
                    continue
                parent_ids.update(_tokens(occurrence.get("source_feature_id")))
                presence = _presence_state(occurrence.get("presence_status"))
                if presence in {"present", "absent"}:
                    presence_calls.add(presence)
                    presence_evidence.add("homologous_dna_observed" if presence == "present"
                                          else "explicit_sequence_absence")
                if presence != "present":
                    continue
                if position_unresolved:
                    continue
                role = occurrence.get("role", "unknown")
                legacy_role_occurrences.append((role, occurrence))
                if transcript_paths:
                    path_role, path_reason, path_ids, path_calls = _role_from_transcript_paths(
                        family, species, occurrence.get("gene_copy_id", "NA"), blocks,
                        transcript_paths, occ_by_id, annotation_view)
                    evidence.add(path_reason)
                    evidence.update(
                        f"path_role:{transcript_id}:{call}"
                        for transcript_id, call in sorted(path_calls.items())
                    )
                    if {"exonic", "not_exonic"} <= set(path_calls.values()):
                        evidence.add("alternative_usage")
                    evidence_ids.update(path_ids)
                    if path_role != "unknown":
                        role_calls.append(path_role)
                        established_exonic_path |= path_role == "exonic"
                    elif annotation_view == "canonical":
                        unknown_role = True
                        evidence.add("canonical_transcript_role_unresolved")
                    elif role in KNOWN_NONEXONIC_ROLES and not blocks:
                        role_calls.append("not_exonic")
                        evidence.add("explicit_nonexonic_source_role")
                        evidence.add("homologous_non_exonic_sequence")
                    else:
                        unknown_role = True
                elif annotation_view == "canonical":
                    unknown_role = True
                    evidence.add("canonical_transcript_path_unavailable")
                elif role in EXONIC_ROLES:
                    role_calls.append("exonic")
                    evidence.add("annotated_exon")
                elif role in KNOWN_NONEXONIC_ROLES:
                    role_calls.append("not_exonic")
                    evidence.add("homologous_non_exonic_sequence")
                else:
                    unknown_role = True
                    evidence.add("homologous_sequence_role_unknown")

            role_conflict = False
            prediction_without_blocks = False
            for completion in completion_by_element[element].get(species, []):
                mapped_copies = {
                    occurrence.get("gene_copy_id")
                    for occurrence, _membership in mapped[element].get(species, [])
                    if occurrence
                }
                if (
                    mapped_copies
                    and completion.get("gene_copy_id") not in MISSING
                    and completion.get("gene_copy_id") not in mapped_copies
                ):
                    continue
                if completion.get("completion_element_assignment") == "ambiguous":
                    evidence.add("completion_element_assignment_ambiguous")
                    presence_evidence.add("completion_element_assignment_ambiguous")
                presence, reason = _completion_presence(completion)
                presence_evidence.add(reason)
                evidence.add(reason)
                evidence_ids.update(_tokens(completion.get("evidence_id")))
                parent_ids.update(_tokens(completion.get("target_parent_occurrence_ids")))
                if presence:
                    presence_calls.add(presence)
                completion_local = _completion_blocks(completion, "dna_aligned_blocks")
                if completion_local:
                    local_blocks.extend(completion_local)
                    interval_ids.update(block["block_id"] for block in completion_local)
                supplied = [record for record in _json_records(completion.get("supplied_annotation_overlaps"))
                            if record.get("strand_relation", "sense") == "sense"]
                supplied_roles = {record.get("role") for record in supplied}
                if completion_local and transcript_paths:
                    path_role, path_reason, path_ids, path_calls = _role_from_transcript_paths(
                        family, species, completion.get("gene_copy_id", "NA"), completion_local,
                        transcript_paths, occ_by_id, annotation_view)
                    evidence.add(path_reason)
                    evidence.update(
                        f"path_role:{transcript_id}:{call}"
                        for transcript_id, call in sorted(path_calls.items())
                    )
                    if {"exonic", "not_exonic"} <= set(path_calls.values()):
                        evidence.add("alternative_usage")
                    evidence_ids.update(path_ids)
                    if path_role != "unknown":
                        role_calls.append(path_role)
                        established_exonic_path |= path_role == "exonic"
                    else:
                        unknown_role = True
                elif annotation_view == "canonical":
                    unknown_role = True
                    evidence.add("completion_not_traceable_to_canonical_transcript")
                elif any(role in EXONIC_ROLES for role in supplied_roles):
                    role_calls.append("exonic")
                    evidence.add("supplied_repertoire_exonic_overlap")
                elif (supplied and supplied_roles <= KNOWN_NONEXONIC_ROLES
                      and completion_local
                      and all(_is_true(record.get("contains_block")) for record in supplied)):
                    role_calls.append("not_exonic")
                    evidence.add("supplied_annotation_conditional_nonexonic")
                elif supplied:
                    unknown_role = True
                    evidence.add("supplied_annotation_overlap_unresolved")
                prediction_supported = _completion_role_prediction_supported(completion)
                predicted = _completion_blocks(completion, "predicted_role_blocks") if prediction_supported else []
                conflicts = _completion_blocks(completion, "annotation_conflict_blocks") if prediction_supported else []
                if predicted or conflicts:
                    comparison_blocks = local_blocks
                    predicted_overlap = _blocks_overlap_status(comparison_blocks, predicted)
                    conflict_overlap = _blocks_overlap_status(comparison_blocks, conflicts)
                    if comparison_blocks and (conflict_overlap is True
                                         or (predicted_overlap is True and "not_exonic" in role_calls)):
                        role_conflict = True
                        evidence.add("local_annotation_prediction_conflict")
                        if "not_exonic" in role_calls:
                            evidence.add("nonexonic_annotation_conflicts_with_exon_prediction")
                    elif ((predicted_overlap is None or conflict_overlap is None)
                          and "not_exonic" in role_calls):
                        unknown_role = True
                        evidence.add("local_overlap_strand_unavailable")
                    elif not comparison_blocks:
                        prediction_without_blocks = True
                        unknown_role = True
                        evidence.add("prediction_overlap_unavailable_without_actual_membership_blocks")
                if prediction_supported:
                    evidence.add("predicted_exonic_role_candidate")
                    evidence.add("sequence_supported_exon_completion_role_unresolved")
                    if not role_calls:
                        unknown_role = True

            if not transcript_paths:
                exonic_occurrences = [
                    occurrence for role, occurrence in legacy_role_occurrences
                    if role in EXONIC_ROLES and _presence_state(occurrence.get("presence_status")) == "present"
                ]
                nonexonic_occurrences = [
                    occurrence for role, occurrence in legacy_role_occurrences
                    if role in KNOWN_NONEXONIC_ROLES and _presence_state(occurrence.get("presence_status")) == "present"
                ]
                if exonic_occurrences and nonexonic_occurrences:
                    def same_local_locus(left, right):
                        if (
                            left.get("gene_copy_id") != right.get("gene_copy_id")
                            or left.get("contig") != right.get("contig")
                            or left.get("strand") != right.get("strand")
                        ):
                            return False
                        try:
                            left_interval = ClosedInterval1(
                                int(left["start"]), int(left["end"]),
                            ).to_interval0()
                            right_interval = ClosedInterval1(
                                int(right["start"]), int(right["end"]),
                            ).to_interval0()
                        except (KeyError, TypeError, ValueError):
                            return False
                        return left_interval.overlaps(right_interval)

                    alternative_usage = any(
                        same_local_locus(exonic, nonexonic)
                        for exonic in exonic_occurrences
                        for nonexonic in nonexonic_occurrences
                    )
                    if alternative_usage:
                        evidence.add("alternative_usage")
                        established_exonic_path = True
                    else:
                        unknown_role = True
                        role_conflict = True
                        evidence.add("local_role_assignment_incompatible")

            presence_state = next(iter(presence_calls)) if len(presence_calls) == 1 else "unknown"
            if len(presence_calls) > 1:
                presence_evidence.add("conflicting_presence_states")
            if presence_state != "present":
                role_conflict = False
            output.append(_observation_row(
                family=family, layer="exon_presence", site_id=element, species=species,
                state=presence_state, state_0="absent", state_1="present",
                evidence=presence_evidence,
                applicability="applicable" if presence_state != "unknown" else "undetermined",
                reason=("homologous_dna_presence_resolved" if presence_state != "unknown"
                        else "homologous_dna_presence_unresolved"), transcript_scope="gene_locus",
                parent_ids=parent_ids, interval_ids=interval_ids, evidence_ids=evidence_ids,
                discovery_rule="homologous_dna_element_correspondence",
                site_kind="homologous_dna_element_presence",
                confidence="medium" if presence_state != "unknown" else "low",
                conclusion="observed" if presence_state != "unknown" else "unknown",
                annotation_completeness=";".join(sorted(completeness)) or "reported"))

            if presence_state == "absent":
                role_state, applicability, reason = "unknown", "inapplicable", "homologous_dna_absent"
                evidence.add("sequence_absent_role_not_applicable")
            elif presence_state != "present":
                role_state, applicability, reason = "unknown", "undetermined", "dna_presence_unresolved"
            elif established_exonic_path:
                role_state, applicability = "exonic", "applicable"
                reason = (
                    "canonical_path_exonic"
                    if annotation_view == "canonical"
                    else "repertoire_any_exonic_path"
                )
                if "not_exonic" in role_calls:
                    evidence.add("alternative_usage")
            elif role_conflict:
                role_state, applicability, reason = "unknown", "undetermined", "local_annotation_conflict"
            elif "exonic" in role_calls:
                role_state, applicability = "exonic", "applicable"
                reason = (
                    "canonical_annotation_exonic"
                    if annotation_view == "canonical"
                    else "repertoire_annotation_exonic"
                )
            elif role_calls and set(role_calls) == {"not_exonic"} and not unknown_role:
                role_state, applicability, reason = "not_exonic", "applicable", "annotation_conditional_nonexonic"
            else:
                role_state, applicability = "unknown", "undetermined"
                reason = ("actual_local_blocks_unavailable" if prediction_without_blocks
                          else "transcript_role_unresolved")
            output.append(_observation_row(
                family=family, layer="exon_role", site_id=element, species=species,
                state=role_state, state_0="not_exonic", state_1="exonic", evidence=evidence,
                applicability=applicability, reason=reason, transcript_scope=transcript_scope,
                parent_ids=parent_ids, interval_ids=interval_ids, evidence_ids=evidence_ids,
                discovery_rule="local_role_on_homologous_dna",
                site_kind="homologous_dna_element_role",
                confidence="medium" if role_state != "unknown" else "low",
                conclusion=("annotation_conditional" if role_state == "not_exonic"
                            else "observed" if role_state == "exonic"
                            else "role_conflict_unresolved" if role_conflict else "unknown"),
                annotation_completeness=";".join(sorted(completeness)) or "reported"))
    return _finalize_site_metadata(output)
