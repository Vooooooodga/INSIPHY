"""Structural site matrix construction for single-copy INSIPHY analyses."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from .elements import EXON_LIKE_ROLES
from .io import norm_state, read_tsv, to_float, write_tsv


EXONIC_ROLES = set(EXON_LIKE_ROLES)
EXON_COMPLETION_CALLS = {
    "hidden_segment_candidate",
    "shifted_splice_site_candidate",
    "joined_exon_candidate",
    "hidden_segment_with_frame_disruption",
}
HOMOLOGOUS_SEQUENCE_STATUSES = {"homologous_sequence_candidate"}


def _single_copy_families(occurrences):
    copies = defaultdict(set)
    species_by_family = defaultdict(set)
    for row in occurrences:
        key = (row.get("family_id", "NA"), row.get("species", "NA"))
        copies[key].add(row.get("gene_copy_id", "NA"))
        species_by_family[key[0]].add(key[1])
    excluded = []
    valid = set(species_by_family)
    for (family, species), gene_copies in sorted(copies.items()):
        if len(gene_copies) > 1:
            valid.discard(family)
            excluded.append(
                {
                    "family_id": family,
                    "species": species,
                    "copy_count": len(gene_copies),
                    "gene_copy_ids": ";".join(sorted(gene_copies)),
                    "reason": "multiple_gene_copies_in_single_copy_mode",
                }
            )
    return valid, species_by_family, excluded


def _completion_presence(row):
    call = row.get("completion_call")
    status = row.get("evidence_status")
    inferred_role = row.get("inferred_role", "unknown")
    if call == "supports_true_absence":
        return "absent", "sequence_supported_true_absence", None, None
    if (call == "supports_annotation" or status == "supports_annotation") and inferred_role in EXONIC_ROLES:
        return "present", "annotated_exon_sequence_present", "exonic", "annotated_exon"
    if call in EXON_COMPLETION_CALLS and inferred_role in EXONIC_ROLES:
        return "present", "sequence_supported_hidden_exon", "exonic", "sequence_supported_exon_completion"
    if call == "boundary_conflict_candidate":
        return "present", "protein_projected_sequence_boundary_conflict", None, "exon_boundary_conflict"
    if status in HOMOLOGOUS_SEQUENCE_STATUSES or call == "homologous_sequence_candidate":
        return "present", "homologous_sequence_candidate", None, "homologous_sequence_role_unknown"
    return None, None, None, None


def _element_site_rows(
    occurrences,
    element_rows,
    completion_rows,
    valid_families,
    species_by_family,
):
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    mapped = defaultdict(lambda: defaultdict(list))
    family_by_element = {}
    element_by_homology = {}
    for row in element_rows:
        if row.get("element_class") != "exon_like":
            continue
        occ = occ_by_id.get(row.get("occurrence_id", ""))
        if not occ or occ.get("family_id") not in valid_families:
            continue
        family = occ["family_id"]
        element = row["element_id"]
        family_by_element[element] = family
        element_by_homology[row.get("homology_id", "")] = element
        mapped[element][occ["species"]].append((occ, row))

    completion_by_element = defaultdict(lambda: defaultdict(list))
    for row in completion_rows:
        element = element_by_homology.get(row.get("homology_id", ""))
        family = row.get("family_id")
        species = row.get("species")
        if element and family in valid_families and species:
            completion_by_element[element][species].append(row)

    rows = []
    for element, family in sorted(family_by_element.items()):
        for species in sorted(species_by_family[family]):
            observations = mapped[element].get(species, [])
            presence_states = set()
            presence_evidence = set()
            states = set()
            evidence = set()
            confidence_flags = set()
            conclusion_flags = set()
            for occ, membership in observations:
                if membership.get("membership_call", "core_member") != "core_member":
                    presence_evidence.add("ambiguous_correspondence")
                    evidence.add("ambiguous_correspondence")
                    confidence_flags.add("low")
                    continue
                presence = norm_state(occ.get("presence_status"))
                if presence == "unknown":
                    presence_evidence.add("uncertain_sequence_or_annotation")
                    evidence.add("uncertain_sequence_or_annotation")
                    confidence_flags.add("low")
                    continue
                presence_states.add(presence)
                presence_evidence.add(
                    "homologous_sequence_observed" if presence == "present" else "explicit_sequence_absence"
                )
                confidence_flags.add("high" if presence == "present" else "medium")
                if presence == "absent":
                    conclusion_flags.add("explicit_absence")
                    continue
                role = occ.get("role", "unknown")
                completion = [
                    row
                    for row in completion_by_element[element].get(species, [])
                    if row.get("gene_copy_id") == occ.get("gene_copy_id")
                ]
                completed_exon = any(
                    row.get("completion_call") in EXON_COMPLETION_CALLS
                    and row.get("inferred_role") in EXONIC_ROLES
                    for row in completion
                )
                if completed_exon:
                    states.add("exonic")
                    evidence.add("sequence_supported_exon_completion")
                    conclusion_flags.add("predicted_exon")
                else:
                    states.add("exonic" if role in EXONIC_ROLES else "not_exonic")
                    evidence.add("annotated_exon" if role in EXONIC_ROLES else "homologous_non_exonic_sequence")
                    conclusion_flags.add("annotated_or_mapped_role")
            for completion in completion_by_element[element].get(species, []):
                presence, presence_ev, role_state, role_ev = _completion_presence(completion)
                if presence:
                    presence_states.add(presence)
                    presence_evidence.add(presence_ev)
                    confidence_flags.add(completion.get("confidence_flag", "medium"))
                    conclusion_flags.add(completion.get("evidence_conclusion", "sequence_evidence"))
                if role_state:
                    states.add(role_state)
                if role_ev:
                    evidence.add(role_ev)
            presence_state = next(iter(presence_states)) if len(presence_states) == 1 else "unknown"
            if len(presence_states) > 1:
                presence_evidence.add("conflicting_presence_states")
                conclusion_flags.add("conflict_unknown")
            rows.append(
                {
                    "family_id": family,
                    "layer": "exon_presence",
                    "site_id": element,
                    "species": species,
                    "state": presence_state,
                    "state_0": "absent",
                    "state_1": "present",
                    "evidence": ";".join(sorted(presence_evidence)) or "no_observation",
                    "site_kind": "element_presence",
                    "confidence_flag": ";".join(sorted(confidence_flags)) or "unknown",
                    "conclusion_flag": ";".join(sorted(conclusion_flags)) or "unknown",
                }
            )
            state = next(iter(states)) if len(states) == 1 else "unknown"
            if len(states) > 1:
                evidence.add("conflicting_transcript_states")
                conclusion_flags.add("conflict_unknown")
            rows.append(
                {
                    "family_id": family,
                    "layer": "exon_role",
                    "site_id": element,
                    "species": species,
                    "state": state,
                    "state_0": "not_exonic",
                    "state_1": "exonic",
                    "evidence": ";".join(sorted(evidence)) or "no_observation",
                    "site_kind": "element_role",
                    "confidence_flag": ";".join(sorted(confidence_flags)) or "unknown",
                    "conclusion_flag": ";".join(sorted(conclusion_flags)) or "unknown",
                }
            )
    return rows


def _read_match_rows(input_dir, output_dir):
    for base in (Path(output_dir), Path(input_dir)):
        path = base / "segment_matches.tsv"
        rows = read_tsv(path, optional=True)
        if rows:
            return rows
    return []


def _occurrence_length(occurrence):
    try:
        return int(occurrence["end"]) - int(occurrence["start"]) + 1
    except (KeyError, TypeError, ValueError):
        return 0


def _parse_projected_reference_blocks(row):
    text = row.get("projected_reference_blocks", "")
    if not text or text == "NA":
        return []
    blocks = []
    strand = row.get("projected_reference_strand") or row.get("strand") or row.get("alignment_strand") or "+"
    for token in text.split(";"):
        if not token or ":" not in token:
            continue
        query_part, target_part = token.split(":", 1)
        if "-" not in query_part or "-" not in target_part:
            continue
        try:
            query_start, query_end = (int(value) for value in query_part.split("-", 1))
            target_start, target_end = (int(value) for value in target_part.split("-", 1))
        except ValueError:
            continue
        blocks.append(
            {
                "query_start": min(query_start, query_end),
                "query_end": max(query_start, query_end),
                "target_start": target_start,
                "target_end": target_end,
                "target_min": min(target_start, target_end),
                "target_max": max(target_start, target_end),
                "strand": strand,
            }
        )
    return blocks


def _project_query_base_to_reference(blocks, query_position):
    for block in blocks:
        if block["strand"] == "-" or block["target_start"] > block["target_end"]:
            continue
        if block["query_start"] <= query_position <= block["query_end"]:
            offset = query_position - block["query_start"]
            return block["target_min"] + offset
    return None


def _project_target_base_to_reference(blocks, target_position):
    for block in blocks:
        if block["strand"] == "-" or block["target_start"] > block["target_end"]:
            continue
        if block["target_min"] <= target_position <= block["target_max"]:
            offset = target_position - block["target_min"]
            return block["query_start"] + offset
    return None


def _mapped_reference_coordinate(match_rows, query_occurrence, reference_occurrence, query_position):
    if query_occurrence == reference_occurrence:
        return query_position
    coordinates = set()
    for row in match_rows:
        if row.get("match_status") != "mapped":
            continue
        query = row.get("query_occurrence_id")
        subject = row.get("subject_occurrence_id")
        blocks = _parse_projected_reference_blocks(row)
        if not blocks:
            continue
        if query == query_occurrence and subject == reference_occurrence:
            coordinate = _project_query_base_to_reference(blocks, query_position)
        elif query == reference_occurrence and subject == query_occurrence:
            coordinate = _project_target_base_to_reference(blocks, query_position)
        else:
            continue
        if coordinate is not None:
            coordinates.add(coordinate)
    if len(coordinates) == 1:
        return next(iter(coordinates))
    return None


def _continuous_reference_block_covers(match_rows, query_occurrence, reference_occurrence, left_coordinate, right_coordinate):
    low = min(left_coordinate, right_coordinate)
    high = max(left_coordinate, right_coordinate)
    for row in match_rows:
        if row.get("match_status") != "mapped":
            continue
        query = row.get("query_occurrence_id")
        subject = row.get("subject_occurrence_id")
        blocks = _parse_projected_reference_blocks(row)
        if not blocks:
            continue
        if query == query_occurrence and subject == reference_occurrence:
            for block in blocks:
                if block["strand"] == "-" or block["target_start"] > block["target_end"]:
                    continue
                if block["target_min"] <= low and block["target_max"] >= high:
                    return True
        elif query == reference_occurrence and subject == query_occurrence:
            for block in blocks:
                if block["strand"] == "-" or block["target_start"] > block["target_end"]:
                    continue
                if block["query_start"] <= low and block["query_end"] >= high:
                    return True
    return False


def _within_exon_boundary(
    family,
    element,
    left_occurrence,
    right_occurrence,
    occ_by_id,
    reference_by_element,
    match_rows,
):
    left_length = _occurrence_length(occ_by_id.get(left_occurrence, {}))
    if left_length <= 0:
        return None
    reference = reference_by_element.get(element)
    if not reference or reference in {left_occurrence, right_occurrence}:
        return None
    donor = _mapped_reference_coordinate(match_rows, left_occurrence, reference, left_length)
    acceptor = _mapped_reference_coordinate(match_rows, right_occurrence, reference, 1)
    if donor is None or acceptor is None or donor >= acceptor:
        return None
    site_id = f"JG_{element}_REF_{reference}_D{donor}_A{acceptor}"
    return {
        "family_id": family,
        "site_id": site_id,
        "element_id": element,
        "reference_occurrence_id": reference,
        "acceptor_reference_occurrence_id": reference,
        "donor_projection": donor,
        "acceptor_projection": acceptor,
        "left_occurrence_id": left_occurrence,
        "right_occurrence_id": right_occurrence,
        "projection_method": "projected_reference_blocks_exact_boundary_bases",
        "site_kind": "within_exon_boundary",
    }


def _between_exon_boundary(
    family,
    left_element,
    right_element,
    left_occurrence,
    right_occurrence,
    occ_by_id,
    reference_by_element,
    match_rows,
):
    left_length = _occurrence_length(occ_by_id.get(left_occurrence, {}))
    if left_length <= 0:
        return None
    left_reference = reference_by_element.get(left_element)
    right_reference = reference_by_element.get(right_element)
    if not left_reference or not right_reference:
        return None
    donor = _mapped_reference_coordinate(match_rows, left_occurrence, left_reference, left_length)
    acceptor = _mapped_reference_coordinate(match_rows, right_occurrence, right_reference, 1)
    if donor is None or acceptor is None:
        return None
    site_id = (
        f"JG_{left_element}_REF_{left_reference}_D{donor}"
        f"__{right_element}_REF_{right_reference}_A{acceptor}"
    )
    return {
        "family_id": family,
        "site_id": site_id,
        "element_id": f"{left_element};{right_element}",
        "reference_occurrence_id": left_reference,
        "acceptor_reference_occurrence_id": right_reference,
        "donor_projection": donor,
        "acceptor_projection": acceptor,
        "left_occurrence_id": left_occurrence,
        "right_occurrence_id": right_occurrence,
        "projection_method": "projected_reference_blocks_exact_boundary_bases",
        "site_kind": "between_exon_boundary",
    }


def _junction_site_rows(input_dir, output_dir, occurrences, element_rows, valid_families, species_by_family):
    paths = read_tsv(Path(input_dir) / "transcript_paths.tsv", optional=True)
    if not paths:
        return []
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    match_rows = _read_match_rows(input_dir, output_dir)
    elements_by_occ = defaultdict(list)
    element_occurrences = defaultdict(set)
    for row in element_rows:
        if row.get("element_class") == "exon_like" and row.get("membership_call", "core_member") == "core_member":
            elements_by_occ[row["occurrence_id"]].append(row["element_id"])
            element_occurrences[row["element_id"]].add(row["occurrence_id"])
    reference_by_element = {}
    for element, occurrences_for_element in element_occurrences.items():
        candidates = [
            occurrence_id
            for occurrence_id in occurrences_for_element
            if occ_by_id.get(occurrence_id, {}).get("role") in EXONIC_ROLES
        ]
        if candidates:
            reference_by_element[element] = max(
                sorted(candidates),
                key=lambda occurrence_id: _occurrence_length(occ_by_id.get(occurrence_id, {})),
            )

    path_groups = defaultdict(list)
    for row in paths:
        occ = occ_by_id.get(row.get("occurrence_id", ""))
        if not occ or occ.get("family_id") not in valid_families:
            continue
        key = (occ["family_id"], occ["species"], occ["gene_copy_id"], row.get("transcript_id", "NA"))
        path_groups[key].append(row)

    observations = defaultdict(lambda: defaultdict(set))
    evidence = defaultdict(lambda: defaultdict(set))
    site_kind = {}
    family_by_site = {}
    element_counts = defaultdict(lambda: defaultdict(list))
    element_occ_by_species = defaultdict(lambda: defaultdict(list))
    boundary_rows = []
    within_boundaries = {}

    for (family, species, _copy, _transcript), path_rows in sorted(path_groups.items()):
        path_rows.sort(key=lambda row: int(row.get("path_rank", "0")))
        previous = None
        previous_occurrence = None
        previous_path = None
        intron_between = False
        path_elements = []
        for row in path_rows:
            occ_id = row.get("occurrence_id", "")
            role = occ_by_id.get(occ_id, {}).get("role", row.get("role", "unknown"))
            if role == "intron":
                if previous is not None:
                    intron_between = True
                continue
            current_elements = sorted(set(elements_by_occ.get(occ_id, [])))
            if not current_elements:
                previous = None
                previous_occurrence = None
                previous_path = None
                intron_between = False
                continue
            current = current_elements[0]
            path_elements.append(current)
            element_occ_by_species[(family, current)][species].append(occ_id)
            if previous is not None:
                if previous == current:
                    boundary = _within_exon_boundary(
                        family, current, previous_occurrence, occ_id, occ_by_id,
                        reference_by_element, match_rows,
                    )
                    if boundary is None:
                        previous = current
                        previous_occurrence = occ_id
                        previous_path = row
                        intron_between = False
                        continue
                    site_id = boundary["site_id"]
                    boundary_rows.append(
                        {
                            **boundary,
                            "species": species,
                            "gene_copy_id": _copy,
                            "transcript_id": _transcript,
                        }
                    )
                    within_boundaries[site_id] = boundary
                    site_kind[site_id] = "within_exon_boundary"
                    family_by_site[site_id] = family
                    if intron_between:
                        observations[site_id][species].add("present")
                        evidence[site_id][species].add("annotated_intron_splits_exon_fragments_projecting_to_ordered_reference_parts")
                    else:
                        evidence[site_id][species].add("split_fragments_project_to_reference_without_annotated_intron")
                else:
                    boundary = _between_exon_boundary(
                        family, previous, current, previous_occurrence, occ_id,
                        occ_by_id, reference_by_element, match_rows,
                    )
                    if boundary is None:
                        previous = current
                        previous_occurrence = occ_id
                        previous_path = row
                        intron_between = False
                        continue
                    site_id = boundary["site_id"]
                    site_kind[site_id] = "between_exon_boundary"
                    family_by_site[site_id] = family
                    observations[site_id][species].add("present" if intron_between else "absent")
                    evidence[site_id][species].add(
                        "annotated_intron_between_homologous_exons" if intron_between else "direct_exonic_adjacency"
                    )
                    boundary_rows.append(
                        {
                            "family_id": family,
                            "species": species,
                            "gene_copy_id": _copy,
                            "transcript_id": _transcript,
                            **boundary,
                        }
                    )
            previous = current
            previous_occurrence = occ_id
            previous_path = row
            intron_between = False
        for element in set(path_elements):
            element_counts[(family, element)][species].append(path_elements.count(element))

    for site_id, boundary in within_boundaries.items():
        family = boundary["family_id"]
        element = boundary["element_id"]
        reference = boundary["reference_occurrence_id"]
        donor = int(boundary["donor_projection"])
        acceptor = int(boundary["acceptor_projection"])
        for species, counts in element_counts[(family, element)].items():
            if not counts or set(counts) != {1}:
                continue
            occurrences_in_species = element_occ_by_species[(family, element)].get(species, [])
            if len(set(occurrences_in_species)) != 1:
                evidence[site_id][species].add("single_element_count_without_unique_occurrence")
                continue
            occurrence = occurrences_in_species[0]
            if occurrence == reference:
                observations[site_id][species].add("absent")
                evidence[site_id][species].add("reference_unsplit_exon_spans_boundary")
                continue
            if _continuous_reference_block_covers(match_rows, occurrence, reference, donor, acceptor):
                observations[site_id][species].add("absent")
                evidence[site_id][species].add("single_continuous_alignment_block_spans_both_boundary_anchors")
            else:
                evidence[site_id][species].add("boundary_not_covered_by_continuous_mapped_exon")

    rows = []
    for site_id, family in sorted(family_by_site.items()):
        for species in sorted(species_by_family[family]):
            values = observations[site_id].get(species, set())
            state = next(iter(values)) if len(values) == 1 else "unknown"
            ev = set(evidence[site_id].get(species, set()))
            conclusion = "observed" if state in {"present", "absent"} else "unknown"
            confidence = "medium" if state in {"present", "absent"} else "low"
            if len(values) > 1:
                ev.add("explicit_transcript_conflict")
                conclusion = "conflict_unknown"
                confidence = "low"
            rows.append(
                {
                    "family_id": family,
                    "layer": "splice_junction",
                    "site_id": site_id,
                    "species": species,
                    "state": state,
                    "state_0": "absent",
                    "state_1": "present",
                    "evidence": ";".join(sorted(ev)) or "not_covered",
                    "site_kind": site_kind.get(site_id, "unknown_boundary"),
                    "confidence_flag": confidence,
                    "conclusion_flag": conclusion,
                }
            )
    write_tsv(
        Path(output_dir) / "splice_boundary_correspondence.tsv",
        boundary_rows,
        [
            "family_id", "species", "gene_copy_id", "transcript_id", "site_id", "element_id",
            "reference_occurrence_id", "acceptor_reference_occurrence_id", "donor_projection", "acceptor_projection",
            "projection_method", "left_occurrence_id", "right_occurrence_id", "site_kind", "strand",
        ],
    )
    return rows


def _structural_site_universe(input_dir, output_dir):
    for base in (Path(input_dir), Path(output_dir)):
        rows = read_tsv(base / "structural_site_universe.tsv", optional=True)
        if rows:
            return rows
    return []


def _state_labels_for_layer(layer):
    if layer == "exon_role":
        return "not_exonic", "exonic"
    return "absent", "present"


def _merge_catalogue_observations(rows, catalogue_rows, valid_families, species_by_family):
    """Merge curated states; site-only declarations carry no observations."""
    merged = {
        (row["family_id"], row["layer"], row["site_id"], row["species"]): dict(row)
        for row in rows
    }
    for row in catalogue_rows:
        family = row.get("family_id")
        layer = row.get("layer")
        site_id = row.get("site_id")
        species = row.get("species")
        state = row.get("state")
        if not all(value and value != "NA" for value in (family, layer, site_id, species, state)):
            continue
        if family not in valid_families or species not in species_by_family.get(family, set()):
            continue
        state_0 = row.get("state_0") or "NA"
        state_1 = row.get("state_1") or "NA"
        if state_0 == "NA" or state_1 == "NA":
            state_0, state_1 = _state_labels_for_layer(layer)
        merged[(family, layer, site_id, species)] = {
            "family_id": family,
            "layer": layer,
            "site_id": site_id,
            "species": species,
            "state": state,
            "state_0": state_0,
            "state_1": state_1,
            "evidence": row.get("evidence", "structural_site_universe_explicit_observation"),
            "site_kind": row.get("site_kind", "catalogue_explicit_observation"),
            "confidence_flag": row.get("confidence_flag", "curated"),
            "conclusion_flag": row.get("conclusion_flag", "explicit_catalogue_observation"),
        }
    return [merged[key] for key in sorted(merged)]


def build_structural_site_matrix(input_dir, output_dir):
    occurrences = read_tsv(
        Path(input_dir) / "segment_occurrences.tsv",
        ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"],
    )
    valid, species_by_family, excluded = _single_copy_families(occurrences)
    universe_rows = _structural_site_universe(input_dir, output_dir)
    element_rows = read_tsv(Path(output_dir) / "element_correspondence.tsv", optional=True)
    if not element_rows:
        raise SystemExit("element_correspondence.tsv is required before single-copy phylogenetic analysis")
    completion_rows = read_tsv(Path(output_dir) / "annotation_completion_candidates.tsv", optional=True)
    rows = _element_site_rows(
        occurrences,
        element_rows,
        completion_rows,
        valid,
        species_by_family,
    )
    rows.extend(_junction_site_rows(input_dir, output_dir, occurrences, element_rows, valid, species_by_family))
    if universe_rows:
        rows = _merge_catalogue_observations(rows, universe_rows, valid, species_by_family)
    return rows, excluded
