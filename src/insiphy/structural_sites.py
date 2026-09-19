"""Construct single-copy structural observations from local correspondence evidence."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from .coordinates import ClosedInterval1, Interval0, parse_legacy_blocks
from .elements import EXON_LIKE_ROLES, NONCODING_ROLES
from .io import (
    STRUCTURAL_SITE_ANNOTATION_VIEWS,
    STRUCTURAL_SITE_SCHEMA_VERSION,
    STRUCTURAL_SITE_SCHEMA_VERSIONS,
    normalize_structural_site_row,
    norm_state,
    read_tsv,
    structural_site_annotation_view,
    write_tsv,
)
from .tree import SpeciesTree


EXONIC_ROLES = set(EXON_LIKE_ROLES)
KNOWN_NONEXONIC_ROLES = set(NONCODING_ROLES)
MISSING = {None, "", "NA", ".", "unknown", "unavailable"}


def _single_copy_families(occurrences):
    copies = defaultdict(set)
    species_by_family = defaultdict(set)
    for row in occurrences:
        key = (row.get("family_id", "NA"), row.get("species", "NA"))
        copies[key].add(row.get("gene_copy_id", "NA"))
        species_by_family[key[0]].add(key[1])
    valid = set(species_by_family)
    excluded = []
    for (family, species), gene_copies in sorted(copies.items()):
        if len(gene_copies) > 1:
            valid.discard(family)
            excluded.append({
                "family_id": family, "species": species,
                "copy_count": len(gene_copies),
                "gene_copy_ids": ";".join(sorted(gene_copies)),
                "reason": "multiple_gene_copies_in_single_copy_mode",
            })
    return valid, species_by_family, excluded


def _tokens(value):
    return {token for token in str(value or "").replace(",", ";").split(";")
            if token and token not in {"NA", "."}}


def _is_true(value):
    return value in {True, 1, "1", "true", "True", "yes"}


def _presence_state(value):
    if value in {0, "0"}:
        return "absent"
    if value in {1, "1"}:
        return "present"
    state = norm_state(value)
    return state


def _json_records(value):
    if value is None or value == "" or value == "NA" or value == ".":
        return []
    try:
        records = json.loads(value) if isinstance(value, str) else value
    except (TypeError, ValueError, json.JSONDecodeError):
        return []
    return records if isinstance(records, list) else []


def _record_interval(record, prefix="target"):
    try:
        start, end = int(record[f"{prefix}_start"]), int(record[f"{prefix}_end"])
    except (KeyError, TypeError, ValueError):
        return None
    start, end = min(start, end), max(start, end)
    if start < 1:
        return None
    return {
        "block_id": str(record.get("block_id", "NA")),
        "contig": record.get(f"{prefix}_contig", record.get("contig", "NA")),
        "strand": record.get(f"{prefix}_strand", record.get("strand", "NA")),
        "interval": ClosedInterval1(start, end).to_interval0(),
    }


def _parse_genomic_blocks(value):
    blocks = []
    for token in str(value or "").split(";"):
        if not token or token in {"NA", "."}:
            continue
        try:
            contig, bounds, strand = token.rsplit(":", 2)
            start, end = (int(part) for part in bounds.split("-", 1))
            interval = ClosedInterval1(min(start, end), max(start, end)).to_interval0()
        except (TypeError, ValueError):
            continue
        blocks.append({"block_id": token, "contig": contig, "strand": strand, "interval": interval})
    return blocks


def _membership_blocks(membership):
    blocks = _parse_genomic_blocks(membership.get("genomic_matched_blocks"))
    if blocks:
        return blocks, "actual_genomic_matched_blocks"
    blocks = [block for record in _json_records(membership.get("actual_matched_blocks"))
              if (block := _record_interval(record))]
    return (blocks, "actual_matched_blocks") if blocks else ([], "actual_local_blocks_unavailable")


def _legacy_occurrence_blocks(occurrence):
    try:
        interval = ClosedInterval1(
            int(occurrence["start"]), int(occurrence["end"]),
        ).to_interval0()
    except (KeyError, TypeError, ValueError):
        return []
    return [{
        "block_id": f"legacy_parent:{occurrence.get('occurrence_id', 'NA')}",
        "contig": occurrence.get("contig", "NA"),
        "strand": occurrence.get("strand", "NA"),
        "interval": interval,
    }]


def _completion_blocks(completion, field):
    blocks = []
    for record in _json_records(completion.get(field)):
        block = _record_interval(record)
        if block:
            blocks.append(block)
            continue
        try:
            start, end = int(record["overlap_start"]), int(record["overlap_end"])
        except (KeyError, TypeError, ValueError):
            continue
        blocks.append({
            "block_id": str(record.get("block_id", "NA")),
            "contig": record.get("target_contig", completion.get("contig", "NA")),
            "strand": record.get("target_strand", completion.get("strand", "NA")),
            "interval": ClosedInterval1(min(start, end), max(start, end)).to_interval0(),
        })
    if blocks or field != "predicted_role_blocks":
        return blocks
    interval = str(completion.get("interval", ""))
    try:
        contig, bounds, strand = interval.rsplit(":", 2)
        start, end = (int(value) for value in bounds.split("-", 1))
        public = ClosedInterval1(min(start, end), max(start, end))
    except (TypeError, ValueError):
        return []
    return [{
        "block_id": f"legacy_completion:{interval}",
        "contig": contig,
        "strand": strand,
        "interval": public.to_interval0(),
    }]


def _blocks_overlap_status(left, right):
    strand_unavailable = False
    for first in left:
        for second in right:
            if (first["contig"] in {"NA", ".", ""}
                    or first["contig"] != second["contig"]
                    or not first["interval"].overlaps(second["interval"])):
                continue
            first_strand, second_strand = first.get("strand"), second.get("strand")
            if first_strand not in {"+", "-"} or second_strand not in {"+", "-"}:
                strand_unavailable = True
            elif first_strand == second_strand:
                return True
    return None if strand_unavailable else False


def _blocks_overlap(left, right):
    return _blocks_overlap_status(left, right) is True


def _completion_presence(row):
    """Resolve homologous DNA presence without promoting a predicted exon role."""
    statuses = []
    for field in ("candidate_resolution_status", "correspondence_status"):
        value = row.get(field)
        status = str(value or "").strip()
        if status not in {"", "NA", "."}:
            statuses.append(status)
    resolved = bool(statuses) and all(status == "resolved" for status in statuses)
    if not resolved:
        return None, "completion_correspondence_unresolved"

    explicit = _presence_state(row.get("homologous_dna_presence"))
    evidence = row.get("homologous_dna_evidence", "NA")
    if explicit == "absent":
        if evidence == "ordered_flank_deletion":
            return "absent", "ordered_flank_deletion"
        return None, "sequence_absence_not_established"
    if explicit == "present":
        if evidence not in MISSING:
            return "present", evidence
        return None, "sequence_presence_evidence_unavailable"
    call = row.get("completion_call")
    if call == "supports_true_absence" and row.get("absence_evidence") == "ordered_flank_deletion":
        return "absent", "ordered_flank_deletion"
    return None, "no_resolved_sequence_evidence"


def _completion_role_prediction_supported(row):
    """Require resolved exon-prediction evidence before testing a local role conflict."""
    call = str(row.get("completion_call", ""))
    has_predicted_blocks = bool(_json_records(row.get("predicted_role_blocks")))
    has_conflict_blocks = bool(_json_records(row.get("annotation_conflict_blocks")))
    if call not in {"predicted_exon_candidate", "boundary_conflict_candidate"} and not (has_predicted_blocks or has_conflict_blocks):
        return False
    if (call in {"predicted_exon_candidate", "boundary_conflict_candidate"} or has_predicted_blocks) and row.get("predicted_role") not in EXONIC_ROLES:
        return False
    if row.get("candidate_resolution_status") not in {None, "", "NA", ".", "resolved"}:
        return False
    if row.get("correspondence_status") not in {None, "", "NA", ".", "resolved"}:
        return False
    if str(row.get("confidence_flag", "")).lower() in {"low", "unknown", "unresolved"}:
        return False
    return True


def _exon_prediction_overlaps_observation(completion, occurrence_or_blocks):
    predicted = (_completion_blocks(completion, "predicted_role_blocks")
                 or _completion_blocks(completion, "annotation_conflict_blocks"))
    observed = occurrence_or_blocks if isinstance(occurrence_or_blocks, list) else []
    return bool(predicted and observed and _blocks_overlap(predicted, observed))


def _path_complete(path_rows, occ_by_id):
    for row in path_rows:
        occurrence = occ_by_id.get(row.get("occurrence_id", ""), {})
        if row.get("partial_start") in {1, "1", True, "true", "True", "yes"}:
            return False
        if row.get("partial_end") in {1, "1", True, "true", "True", "yes"}:
            return False
        text = ";".join(str(value).lower() for value in (
            row.get("path_status", ""), row.get("boundary_class", ""),
            row.get("annotation_completeness", ""), occurrence.get("boundary_class", "")))
        if any(term in text for term in ("partial", "truncated", "incomplete")):
            return False
    return True


def _classify_path_blocks(path_rows, blocks, occ_by_id):
    if not blocks:
        return "uncovered"
    exonic = []
    for row in path_rows:
        occurrence = occ_by_id.get(row.get("occurrence_id", ""), row)
        path_role = row.get("path_role")
        row_role = row.get("role", occurrence.get("role", "unknown"))
        is_exonic = path_role == "exonic" or (
            path_role in MISSING and row_role in EXONIC_ROLES
        )
        if not is_exonic:
            continue
        try:
            interval = ClosedInterval1(int(occurrence["start"]), int(occurrence["end"])).to_interval0()
        except (KeyError, TypeError, ValueError):
            continue
        exonic.append((occurrence.get("contig", row.get("contig", "NA")),
                       occurrence.get("strand", row.get("strand", "NA")), interval))
    if not exonic or len({(contig, strand) for contig, strand, _interval in exonic}) != 1:
        return "uncovered"
    contig, strand = exonic[0][:2]
    span = Interval0(min(interval.start0 for _contig, _strand, interval in exonic),
                     max(interval.end0 for _contig, _strand, interval in exonic))
    calls = []
    for block in blocks:
        if block["contig"] != contig or block.get("strand") not in {"+", "-"}:
            calls.append("uncovered")
            continue
        if block["strand"] != strand:
            calls.append("conflicting")
            continue
        containing = any(interval.start0 <= block["interval"].start0
                         and block["interval"].end0 <= interval.end0
                         for _contig, _strand, interval in exonic)
        overlapping = any(interval.overlaps(block["interval"])
                          for _contig, _strand, interval in exonic)
        if containing:
            calls.append("exonic")
        elif overlapping:
            calls.append("conflicting")
        elif span.start0 <= block["interval"].start0 and block["interval"].end0 <= span.end0:
            calls.append("not_exonic")
        else:
            calls.append("uncovered")
    return calls[0] if calls and len(set(calls)) == 1 else "conflicting"


def _role_from_transcript_paths(family, species, gene_copy, blocks, transcript_paths,
                                occ_by_id, annotation_view):
    relevant = [row for row in transcript_paths if row.get("family_id") == family
                and row.get("species") == species and row.get("gene_copy_id") == gene_copy]
    if annotation_view == "canonical":
        relevant = [row for row in relevant if row.get("path_status") == "canonical_transcript_path"]
        if not relevant:
            return "unknown", "canonical_transcript_not_explicitly_recorded", set(), {}
    by_transcript = defaultdict(list)
    for row in relevant:
        by_transcript[row.get("transcript_id", "NA")].append(row)
    if not by_transcript:
        return "unknown", "transcript_paths_unavailable", set(), {}
    calls = {}
    for transcript_id, path_rows in by_transcript.items():
        calls[transcript_id] = (_classify_path_blocks(path_rows, blocks, occ_by_id)
                                if _path_complete(path_rows, occ_by_id) else "partial")
    if annotation_view == "repertoire" and any(call == "exonic" for call in calls.values()):
        return "exonic", "repertoire_any_exonic_path", set(calls), calls
    if annotation_view == "canonical" and set(calls.values()) == {"exonic"}:
        return "exonic", "canonical_path_exonic", set(calls), calls
    if any(call == "conflicting" for call in calls.values()):
        return "unknown", "local_path_role_conflict", set(calls), calls
    if any(call in {"partial", "uncovered", "conflicting"} for call in calls.values()):
        return "unknown", "path_partial_uncovered_or_conflicting", set(calls), calls
    if calls and set(calls.values()) == {"not_exonic"}:
        return "not_exonic", "annotation_conditional_all_complete_paths_nonexonic", set(calls), calls
    return "unknown", "transcript_role_unresolved", set(calls), calls


def _observation_row(*, family, layer, site_id, species, state, state_0, state_1,
                     evidence, applicability, reason, transcript_scope,
                     parent_ids=(), interval_ids=(), evidence_ids=(), discovery_rule,
                     linked_group="NA", site_kind, confidence="unassessed",
                     conclusion="unassessed", annotation_completeness="unassessed"):
    annotation_view = (
        "view_independent"
        if layer == "exon_presence"
        else "canonical"
        if transcript_scope == "explicit_canonical_transcript"
        else "repertoire"
    )
    return {
        "family_id": family, "layer": layer, "site_id": site_id, "species": species,
        "state": state, "state_0": state_0, "state_1": state_1,
        "evidence": ";".join(sorted(set(evidence))) or "no_observation",
        "applicability": applicability, "observation_reason": reason,
        "transcript_scope": transcript_scope,
        "parent_feature_ids": ";".join(sorted(set(parent_ids))) or "NA",
        "member_interval_ids": ";".join(sorted(set(interval_ids))) or "NA",
        "evidence_ids": ";".join(sorted(set(evidence_ids))) or "NA",
        "discovery_rule": discovery_rule, "discovery_species": "NA",
        "observation_mask": "missing" if state == "unknown" else "observed",
        "linked_group_id": linked_group, "site_kind": site_kind,
        "observation_source": "genome_and_supplied_annotation",
        "annotation_completeness": annotation_completeness,
        "confidence_flag": confidence, "conclusion_flag": conclusion,
        "annotation_view": annotation_view,
    }


def _finalize_site_metadata(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["family_id"], row["layer"], row["site_id"])].append(row)
    for site_rows in grouped.values():
        discovered = sorted(row["species"] for row in site_rows if row["state"] == row["state_1"])
        for row in site_rows:
            if row.get("discovery_species") in MISSING:
                row["discovery_species"] = ";".join(discovered) or "NA"
    return rows


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


def _position_eligible(row):
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
    if any(row.get("role") == "intron" for row in path_rows[left_index + 1:right_index]):
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


def _junction_site_rows(input_dir, output_dir, occurrences, element_rows, valid_families,
                        species_by_family, element_site_rows=None,
                        annotation_view="repertoire"):
    if annotation_view not in {"repertoire", "canonical"}:
        raise ValueError("annotation_view must be 'repertoire' or 'canonical'")
    paths = read_tsv(Path(input_dir) / "transcript_paths.tsv", optional=True)
    if not paths:
        return []
    if annotation_view == "canonical":
        paths = [
            row for row in paths
            if row.get("path_status") == "canonical_transcript_path"
        ]
        if not paths:
            return []
    intron_rows = read_tsv(Path(input_dir) / "intron_sites.tsv", optional=True)
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    match_rows = _read_match_rows(input_dir, output_dir)
    membership_by_occ, elements_by_occ, coverage_by_occ = {}, defaultdict(list), {}
    for row in element_rows:
        if row.get("membership_call", "core_member") != "core_member":
            continue
        element, occurrence = row.get("element_id"), row.get("occurrence_id")
        membership_by_occ[(element, occurrence)] = row
        elements_by_occ[occurrence].append(element)
        coverage_by_occ[occurrence] = row.get("reference_coverage_relation", "uncovered")
    references = _reference_by_element(element_rows, occ_by_id)
    path_groups = defaultdict(list)
    for row in paths:
        occurrence = occ_by_id.get(row.get("occurrence_id", ""), row)
        family = occurrence.get("family_id", row.get("family_id"))
        if family in valid_families:
            key = (family, occurrence.get("species", row.get("species")),
                   occurrence.get("gene_copy_id", row.get("gene_copy_id")), row.get("transcript_id", "NA"))
            path_groups[key].append(row)

    observations = defaultdict(lambda: defaultdict(set))
    evidence = defaultdict(lambda: defaultdict(set))
    transcript_evidence = defaultdict(lambda: defaultdict(set))
    parent_evidence = defaultdict(lambda: defaultdict(set))
    site_info, boundary_rows = {}, []
    linked_groups_by_site = defaultdict(set)
    paths_by_element_species = defaultdict(lambda: defaultdict(list))
    for (family, species, gene_copy, transcript_id), path_rows in sorted(path_groups.items()):
        path_rows.sort(key=lambda row: int(row.get("path_rank", "0") or 0))
        entries = []
        for index, row in enumerate(path_rows):
            occurrence = occ_by_id.get(row.get("occurrence_id", ""), row)
            if occurrence.get("role", row.get("role", "unknown")) not in EXONIC_ROLES:
                continue
            elements = sorted(set(elements_by_occ.get(row.get("occurrence_id", ""), [])))
            if len(elements) == 1:
                entries.append((index, row, occurrence, elements[0]))
                paths_by_element_species[(family, elements[0])][species].append((transcript_id, row["occurrence_id"]))
        counts = defaultdict(int)
        for _index, _row, _occurrence, element in entries:
            counts[element] += 1
        linked = {
            element: (
                f"LG_{family}_{element}_{references.get(element, 'NA')}"
                f"_SP_{species}_GC_{gene_copy}_TX_{transcript_id}"
            )
            for element, count in counts.items()
            if count > 2
        }
        for left_entry, right_entry in zip(entries, entries[1:]):
            left_index, left_row, left, left_element = left_entry
            right_index, right_row, right, right_element = right_entry
            linked_group = linked.get(left_element, "NA") if left_element == right_element else "NA"
            boundary, unavailable = _boundary_from_pair(
                family, left_element, right_element, left.get("occurrence_id"), right.get("occurrence_id"),
                occ_by_id, membership_by_occ, references, match_rows, linked_group)
            phase_status = _phase_status(left_row, right_row, left, right)
            boundary_rows.append({
                "family_id": family, "species": species, "gene_copy_id": gene_copy,
                "transcript_id": transcript_id, "site_id": boundary["site_id"] if boundary else "NA",
                "element_id": left_element if left_element == right_element else f"{left_element};{right_element}",
                "reference_occurrence_id": references.get(left_element, "NA"),
                "acceptor_reference_occurrence_id": references.get(right_element, "NA"),
                "donor_projection": boundary["donor_projection"] if boundary else "NA",
                "acceptor_projection": boundary["acceptor_projection"] if boundary else "NA",
                "cutpoint0": boundary["cutpoint0"] if boundary else "NA",
                "projection_method": boundary["projection_method"] if boundary else "unavailable",
                "left_occurrence_id": left.get("occurrence_id", "NA"),
                "right_occurrence_id": right.get("occurrence_id", "NA"),
                "strand": left.get("strand", "NA"),
                "site_kind": boundary["site_kind"] if boundary else "unresolved_junction",
                "position_edge_eligible": "1" if boundary else "0", "phase_status": phase_status,
                "linked_group_id": linked_group,
                "unavailable_reason": unavailable or ("phase_unknown" if phase_status == "unknown" else "NA"),
            })
            if boundary is None:
                continue
            site_id, site_info[boundary["site_id"]] = boundary["site_id"], boundary
            if linked_group != "NA":
                linked_groups_by_site[site_id].add(linked_group)
            parent_evidence[site_id][species].update({
                left.get("occurrence_id", "NA"), right.get("occurrence_id", "NA")})
            has_intron = _explicit_intron_between(path_rows, left_index, right_index, left, right,
                                                  intron_rows, transcript_id)
            transcript_evidence[site_id][species].add(transcript_id)
            if phase_status == "unknown":
                evidence[site_id][species].add("junction_phase_unknown")
            if has_intron:
                observations[site_id][species].add("present")
                evidence[site_id][species].add("explicit_annotated_intron_at_exact_cutpoint")
            elif _genomically_contiguous(left, right):
                observations[site_id][species].add("absent")
                evidence[site_id][species].add("single_transcript_continuously_spans_exact_cutpoint")
            else:
                evidence[site_id][species].add("no_explicit_intron_and_no_continuous_cutpoint_span")
                evidence[site_id][species].add("missing_intron_record_with_genomic_gap")

    for site_id, info in site_info.items():
        if info["site_kind"] != "within_element_junction":
            continue
        family, element, reference = info["family_id"], info["element_id"], info["reference_occurrence_id"]
        donor, acceptor = int(info["donor_projection"]), int(info["acceptor_projection"])
        for species, pairs in paths_by_element_species[(family, element)].items():
            by_transcript = defaultdict(set)
            for transcript_id, occurrence_id in pairs:
                by_transcript[transcript_id].add(occurrence_id)
            for transcript_id, occurrence_ids in by_transcript.items():
                if len(occurrence_ids) != 1:
                    continue
                occurrence_id = next(iter(occurrence_ids))
                if coverage_by_occ.get(occurrence_id) == "repeated_overlap":
                    evidence[site_id][species].add("repeated_overlap_position_ambiguous")
                    continue
                membership = membership_by_occ.get((element, occurrence_id), {})
                if not _member_position_eligible(membership, occurrence_id, reference):
                    evidence[site_id][species].add("position_mapping_unavailable")
                elif _continuous_reference_block_covers(match_rows, occurrence_id, reference, donor, acceptor):
                    observations[site_id][species].add("absent")
                    evidence[site_id][species].add("single_transcript_continuous_actual_block_spans_cutpoint")
                    evidence[site_id][species].add(
                        "single_transcript_continuous_alignment_block_spans_both_boundary_anchors"
                    )
                    transcript_evidence[site_id][species].add(transcript_id)
                    parent_evidence[site_id][species].add(occurrence_id)
                else:
                    evidence[site_id][species].add("actual_block_does_not_span_both_cutpoint_sides")

    presence_lookup = {(row.get("site_id"), row.get("species")): row.get("state")
                       for row in element_site_rows or [] if row.get("layer") == "exon_presence"}
    rows = []
    transcript_scope = (
        "explicit_canonical_transcript"
        if annotation_view == "canonical"
        else "annotated_transcript_repertoire"
    )
    for site_id, info in sorted(site_info.items()):
        family, elements = info["family_id"], info["element_id"].split(";")
        for species in sorted(species_by_family[family]):
            values, ev = observations[site_id].get(species, set()), set(evidence[site_id].get(species, set()))
            if any(presence_lookup.get((element, species)) == "absent" for element in elements):
                state, applicability, reason = "unknown", "inapplicable", "homologous_dna_required_for_junction_absent"
                ev.add("sequence_absent_junction_not_applicable")
            elif "present" in values:
                reason = (
                    "canonical_transcript_contains_annotated_junction"
                    if annotation_view == "canonical"
                    else "repertoire_contains_annotated_junction"
                )
                state, applicability = "present", "applicable"
                if "absent" in values:
                    ev.add("alternative_transcript_without_junction")
                    ev.add("alternative_transcript_without_this_junction")
                    ev.add("reference_unsplit_exon_spans_boundary")
            elif values == {"absent"}:
                state, applicability, reason = "absent", "applicable", "continuous_path_spans_exact_cutpoint"
            else:
                state, applicability, reason = "unknown", "undetermined", "junction_position_or_path_coverage_unresolved"
            rows.append(_observation_row(
                family=family, layer="splice_junction", site_id=site_id, species=species,
                state=state, state_0="absent", state_1="present", evidence=ev,
                applicability=applicability, reason=reason,
                transcript_scope=transcript_scope,
                parent_ids=(parent_evidence[site_id].get(species, set())
                            or {info["left_occurrence_id"], info["right_occurrence_id"]}),
                interval_ids={f"cut0:{info['cutpoint0']}"} if info["cutpoint0"] != "NA" else set(),
                evidence_ids=transcript_evidence[site_id].get(species, set()),
                discovery_rule="unique_exact_position_correspondence",
                linked_group=(";".join(sorted(linked_groups_by_site.get(site_id, set()))) or "NA"),
                site_kind=info["site_kind"],
                confidence="medium" if state != "unknown" else "low",
                conclusion=(
                    "repertoire_present"
                    if state == "present" and annotation_view == "repertoire"
                    else "observed" if state != "unknown" else "unknown"
                )))
    write_tsv(Path(output_dir) / "splice_boundary_correspondence.tsv", boundary_rows, [
        "family_id", "species", "gene_copy_id", "transcript_id", "site_id", "element_id",
        "reference_occurrence_id", "acceptor_reference_occurrence_id", "donor_projection",
        "acceptor_projection", "cutpoint0", "projection_method", "left_occurrence_id",
        "right_occurrence_id", "strand", "site_kind", "position_edge_eligible", "phase_status",
        "linked_group_id", "unavailable_reason"])
    return _finalize_site_metadata(rows)


def _structural_site_universe(input_dir, output_dir):
    for base in (Path(input_dir), Path(output_dir)):
        rows = read_tsv(base / "structural_site_universe.tsv", optional=True)
        if rows:
            return rows
    return []


def _state_labels_for_layer(layer):
    return ("not_exonic", "exonic") if layer == "exon_role" else ("absent", "present")


def _merge_catalogue_observations(
    rows, catalogue_rows, valid_families, species_by_family, annotation_view
):
    if annotation_view not in {"canonical", "repertoire"}:
        raise ValueError("annotation_view must be 'repertoire' or 'canonical'")
    merged = {(row["family_id"], row["layer"], row["site_id"], row["species"]): dict(row)
              for row in rows}
    for catalogue in catalogue_rows:
        family, layer = catalogue.get("family_id"), catalogue.get("layer")
        site_id, species = catalogue.get("site_id"), catalogue.get("species")
        state = norm_state(catalogue.get("state"))
        if not all(value and value != "NA" for value in (family, layer, site_id)):
            continue
        catalogue_view = structural_site_annotation_view(catalogue)
        if layer == "exon_presence":
            if catalogue_view != "view_independent":
                raise ValueError("structural-site presence catalogue entries are view-independent")
        elif catalogue_view not in {"canonical", "repertoire"}:
            raise ValueError(
                "structural-site role and junction catalogue entries require an explicit "
                "or inferable canonical/repertoire view"
            )
        elif catalogue_view != annotation_view:
            raise ValueError(
                "structural-site catalogue annotation view conflicts with this analysis: "
                f"catalogue={catalogue_view}, analysis={annotation_view}"
            )
        if not species or species == "NA":
            continue
        if state == "unknown" and catalogue.get("state") in MISSING:
            continue
        if family not in valid_families or species not in species_by_family.get(family, set()):
            continue
        state_0, state_1 = catalogue.get("state_0"), catalogue.get("state_1")
        if state_0 in MISSING or state_1 in MISSING:
            state_0, state_1 = _state_labels_for_layer(layer)
        discovery_rule = catalogue.get("discovery_rule")
        if discovery_rule in MISSING:
            discovery_rule = "independent_catalogue"
        row = _observation_row(
            family=family, layer=layer, site_id=site_id, species=species, state=state,
            state_0=state_0, state_1=state_1,
            evidence={catalogue.get("evidence", "structural_site_universe_explicit_observation")},
            applicability=catalogue.get("applicability", "applicable" if state != "unknown" else "undetermined"),
            reason=catalogue.get("observation_reason", "explicit_catalogue_observation"),
            transcript_scope=catalogue.get("transcript_scope", "catalogue_defined"),
            parent_ids=_tokens(catalogue.get("parent_feature_ids")),
            interval_ids=_tokens(catalogue.get("member_interval_ids")),
            evidence_ids=_tokens(catalogue.get("evidence_ids")),
            discovery_rule=discovery_rule,
            linked_group=catalogue.get("linked_group_id", "NA"),
            site_kind=catalogue.get("site_kind", "catalogue_explicit_observation"),
            confidence=catalogue.get("confidence_flag", "curated"),
            conclusion=catalogue.get("conclusion_flag", "explicit_catalogue_observation"),
            annotation_completeness=catalogue.get("annotation_completeness", "catalogue_reported"))
        row["observation_source"] = catalogue.get("observation_source", "independent_structural_site_catalogue")
        row["discovery_species"] = catalogue.get("discovery_species", "NA")
        row["annotation_view"] = catalogue_view
        catalogue_mask = catalogue.get("observation_mask")
        if catalogue_mask not in {None, ""}:
            row["observation_mask"] = catalogue_mask
        merged[(family, layer, site_id, species)] = row
    return _finalize_site_metadata([merged[key] for key in sorted(merged)])


def _complete_tree_tip_observations(rows, tip_labels):
    """Materialize unknown observations instead of letting inference invent missing tips."""

    by_site = defaultdict(list)
    for row in rows:
        by_site[(row["family_id"], row["layer"], row["site_id"])].append(row)

    completed = list(rows)
    for site_key, site_rows in sorted(by_site.items()):
        template = site_rows[0]
        observed_species = {row["species"] for row in site_rows}
        for species in sorted(set(tip_labels) - observed_species):
            missing = dict(template)
            missing.update({
                "species": species,
                "state": "unknown",
                "evidence": "no_structural_observation_for_tree_tip",
                "applicability": "undetermined",
                "observation_reason": "tree_tip_missing_from_structural_evidence",
                "parent_feature_ids": "NA",
                "member_interval_ids": "NA",
                "evidence_ids": "NA",
                "observation_mask": "missing",
                "observation_source": "no_observation_for_tree_tip",
                "annotation_completeness": "tree_tip_not_observed",
                "confidence_flag": "low",
                "conclusion_flag": "unknown",
            })
            completed.append(missing)
    return completed


def build_structural_site_matrix(input_dir, output_dir, annotation_view="repertoire"):
    occurrences = read_tsv(Path(input_dir) / "segment_occurrences.tsv",
                           ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"])
    valid, species_by_family, excluded = _single_copy_families(occurrences)
    element_rows = read_tsv(Path(output_dir) / "element_correspondence.tsv", optional=True)
    if not element_rows:
        raise SystemExit("element_correspondence.tsv is required before single-copy phylogenetic analysis")
    completion_rows = read_tsv(Path(output_dir) / "annotation_completion_candidates.tsv", optional=True)
    transcript_paths = read_tsv(Path(input_dir) / "transcript_paths.tsv", optional=True)
    rows = _element_site_rows(occurrences, element_rows, completion_rows, valid, species_by_family,
                              transcript_paths=transcript_paths, annotation_view=annotation_view)
    rows.extend(_junction_site_rows(
        input_dir,
        output_dir,
        occurrences,
        element_rows,
        valid,
        species_by_family,
        element_site_rows=rows,
        annotation_view=annotation_view,
    ))
    universe_rows = _structural_site_universe(input_dir, output_dir)
    if universe_rows:
        rows = _merge_catalogue_observations(
            rows, universe_rows, valid, species_by_family, annotation_view
        )
    tree_path = Path(input_dir) / "species_tree.tsv"
    if tree_path.exists():
        tree = SpeciesTree(read_tsv(tree_path, ["node_id", "parent_id", "label"]))
        rows = _complete_tree_tip_observations(rows, tree.leaf_by_label)
    for row in rows:
        expected_view = "view_independent" if row.get("layer") == "exon_presence" else annotation_view
        actual_view = structural_site_annotation_view(row)
        if actual_view != expected_view:
            raise ValueError(
                f"structural observation view conflicts with analysis: "
                f"site={row.get('site_id')}, row={actual_view}, analysis={expected_view}"
            )
    return [normalize_structural_site_row(row) for row in rows], excluded


def structural_matrix_annotation_view(rows):
    views = set()
    for row in rows:
        schema_version = str(row.get("schema_version", "")).strip()
        if schema_version and schema_version not in STRUCTURAL_SITE_SCHEMA_VERSIONS:
            raise ValueError(f"unsupported structural-site schema_version: {schema_version!r}")
        if schema_version == STRUCTURAL_SITE_SCHEMA_VERSION and not str(
            row.get("annotation_view", "")
        ).strip():
            raise ValueError("schema-v3 structural-site rows require annotation_view")
        try:
            view = structural_site_annotation_view(row)
        except ValueError as exc:
            raise ValueError(f"invalid structural-site annotation view: {exc}") from exc
        if view not in STRUCTURAL_SITE_ANNOTATION_VIEWS:
            raise ValueError(
                "Frozen role/junction observations require an explicit annotation_view; "
                "legacy rows must have an inferable transcript_scope"
            )
        views.add(view)

    model_views = views - {"view_independent"}
    if len(model_views) > 1:
        raise ValueError("Structural-site matrix mixes canonical and repertoire observations")
    if model_views:
        return next(iter(model_views))
    return "view_independent"
