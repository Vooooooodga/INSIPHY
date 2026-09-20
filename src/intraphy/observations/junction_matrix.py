"""observations / junction matrix: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.observations.junction_coordinates import _boundary_from_pair
from intraphy.observations.junction_coordinates import _continuous_reference_block_covers
from intraphy.observations.junction_coordinates import _explicit_intron_between
from intraphy.observations.junction_coordinates import _genomically_contiguous
from intraphy.observations.junction_coordinates import _member_position_eligible
from intraphy.observations.junction_coordinates import _phase_status
from intraphy.observations.junction_coordinates import _read_match_rows
from intraphy.observations.junction_coordinates import _reference_by_element
from intraphy.observations.support import EXONIC_ROLES
from intraphy.observations.support import _finalize_site_metadata
from intraphy.observations.support import _observation_row
from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from pathlib import Path


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

    presence_lookup = {(row.get("family_id"), row.get("site_id"), row.get("species")): row.get("state")
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
            if any(presence_lookup.get((family, element, species)) == "absent" for element in elements):
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
