"""observations / support: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0
from intraphy.elements import EXON_LIKE_ROLES
from intraphy.elements import NONCODING_ROLES
from intraphy.storage.values import norm_state
import json


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
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(f"malformed structural coordinate records: {value!r}") from exc
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
