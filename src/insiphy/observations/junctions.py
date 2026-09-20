"""observations / junctions: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.coordinates import parse_legacy_blocks
from insiphy.observations.support import EXONIC_ROLES
from insiphy.observations.support import MISSING
from insiphy.observations.support import _finalize_site_metadata
from insiphy.observations.support import _is_true
from insiphy.observations.support import _json_records
from insiphy.observations.support import _observation_row
from insiphy.observations.support import _tokens
from insiphy.storage.tabular import read_tsv
from insiphy.storage.tabular import write_tsv
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
            intervening = path_rows[left_index + 1:right_index]
            if any(occ_by_id.get(row.get("occurrence_id"), row).get("role", row.get("role")) != "intron"
                   for row in intervening):
                boundary, unavailable = None, "unresolved_intervening_annotated_feature"
            else:
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
