"""Structural character identities and coordinate evidence, independent of inference."""
from collections import defaultdict
import json
from pathlib import Path
from urllib.parse import quote

from .schema import structural_site_observed_state
from ..storage.tabular import read_tsv, write_tsv

KEY = ("family_id", "layer", "site_id")
CATALOG_FIELDS = [*KEY, "character_id", "count_unit_id", "biological_character",
                  "state_0", "state_1", "annotation_view", "site_kind",
                  "dependency_groups", "discovery_rule", "coordinates_available",
                  "count_interpretation", "catalogue_version"]
COORDINATE_FIELDS = [*KEY, "species", "gene_copy_id", "occurrence_id", "contig",
                     "start", "end", "strand", "coordinate_type", "element_id"]
DESCRIPTIONS = {
    "exon_presence": "homologous_sequence_presence",
    "exon_role": "annotation_conditional_exon_identity",
    "splice_junction": "splice_junction_at_corresponding_position",
}


def site_key(row):
    return tuple(str(row[field]) for field in KEY)


def character_id(key):
    """Unambiguous identifier; no evidence score or biological direction is encoded."""
    return "/".join(quote(str(value), safe="") for value in key)


def tokens(value):
    return {x for x in str(value or "").split(";") if x not in {"", "NA", ".", "unknown"}}


def _member_coordinates(members, keys):
    for member in members:
        blocks = member.get("actual_matched_blocks", "NA")
        records = json.loads(blocks) if blocks not in {"", "NA"} else []
        for layer in ("exon_presence", "exon_role"):
            key = (member.get("family_id"), layer, member.get("element_id"))
            if key not in keys:
                continue
            for block in records:
                if not all(block.get(f"target_{x}") is not None for x in ("start", "end", "contig", "strand")):
                    continue
                yield dict(zip(KEY, key), species=member["species"],
                           gene_copy_id=member.get("gene_copy_id", "NA"),
                           occurrence_id=member.get("occurrence_id", "NA"),
                           contig=block["target_contig"], start=int(block["target_start"]),
                           end=int(block["target_end"]), strand=block["target_strand"],
                           coordinate_type="genomic_aligned_interval_1_based_closed",
                           element_id=member["element_id"])


def _junction_coordinates(boundaries, keys):
    for boundary in boundaries:
        key = (boundary.get("family_id"), "splice_junction", boundary.get("site_id"))
        if key not in keys or boundary.get("position_edge_eligible") != "1":
            continue
        for side, position, element in (
            ("reference_occurrence_id", "donor_projection", "donor"),
            ("acceptor_reference_occurrence_id", "acceptor_projection", "acceptor"),
        ):
            if boundary.get(position) in {None, "", "NA"}:
                continue
            yield dict(zip(KEY, key), species=boundary["species"],
                       gene_copy_id=boundary.get("gene_copy_id", "NA"),
                       occurrence_id=boundary.get(side, "NA"), contig="NA",
                       start=boundary[position], end=boundary[position],
                       strand=boundary.get("strand", "NA"),
                       coordinate_type=f"{element}_reference_sequence_position_1_based",
                       element_id=boundary.get("element_id", "NA"))


def coordinate_evidence(rows, evidence_dir, frozen=False):
    """Frozen analyses use only their own archived evidence, never current annotations."""
    keys = {site_key(row) for row in rows}
    directory = Path(evidence_dir)
    saved = directory / "character_coordinates.tsv"
    if frozen:
        return [row for row in read_tsv(saved, optional=True) if site_key(row) in keys]
    members = read_tsv(directory / "element_correspondence.tsv", optional=True)
    boundaries = read_tsv(directory / "splice_boundary_correspondence.tsv", optional=True)
    values = [*_member_coordinates(members, keys), *_junction_coordinates(boundaries, keys)]
    unique = {tuple(str(row.get(field, "NA")) for field in COORDINATE_FIELDS): row for row in values}
    return [unique[key] for key in sorted(unique)]


def add_interval_dependencies(rows, coordinates):
    """Mark overlapping sequence characters without merging their states or counts.

    Shared coordinates are a conservative dependence diagnostic, not evidence
    for a shared mutation. Transcript membership alone never creates a group.
    """
    by_context = defaultdict(list)
    for row in coordinates:
        if row["coordinate_type"] == "genomic_aligned_interval_1_based_closed":
            context = (row["family_id"], row["layer"], row["species"],
                       row["gene_copy_id"], row["contig"], row["strand"])
            by_context[context].append(row)
    extra = defaultdict(set)
    for context, intervals in sorted(by_context.items()):
        ordered = sorted(intervals, key=lambda row: (int(row["start"]), int(row["end"]), row["site_id"]))
        for index, left in enumerate(ordered):
            for right in ordered[index + 1:]:
                if int(right["start"]) > int(left["end"]):
                    break
                if left["site_id"] == right["site_id"]:
                    continue
                group = "overlap:" + character_id((*context,
                    max(int(left["start"]), int(right["start"])),
                    min(int(left["end"]), int(right["end"]))))
                extra[site_key(left)].add(group)
                extra[site_key(right)].add(group)
    return [dict(row, linked_group_id=";".join(sorted(
        tokens(row.get("linked_group_id")) | extra[site_key(row)])) or "NA") for row in rows]


def validate_sequence_applicability(rows):
    """Reject an explicitly observed exon role when its homologous DNA is absent."""
    presence = {(r["family_id"], r["site_id"], r["species"]): structural_site_observed_state(r)
                for r in rows if r["layer"] == "exon_presence"}
    for row in rows:
        key = (row["family_id"], row["site_id"], row["species"])
        if row["layer"] == "exon_role" and presence.get(key) == "absent":
            if structural_site_observed_state(row) != "unknown":
                raise ValueError("Exon role is inapplicable where homologous DNA is absent: " + "/".join(key))


def write_character_catalogue(output_dir, rows, coordinates):
    grouped = defaultdict(list)
    available = {site_key(row) for row in coordinates}
    for row in rows:
        grouped[site_key(row)].append(row)
    catalogue, dependencies = [], defaultdict(set)
    for key, observations in sorted(grouped.items()):
        first = observations[0]
        groups = set().union(*(tokens(row.get("linked_group_id")) for row in observations))
        for group in groups:
            dependencies[(key[0], group)].add(key)
        catalogue.append(dict(zip(KEY, key), character_id=character_id(key),
            count_unit_id=character_id(key), biological_character=DESCRIPTIONS[key[1]],
            state_0=first["state_0"], state_1=first["state_1"],
            annotation_view=first["annotation_view"], site_kind=first.get("site_kind", "unspecified"),
            dependency_groups=";".join(sorted(groups)) or "NA",
            discovery_rule=first.get("discovery_rule", "NA"),
            coordinates_available=int(key in available),
            count_interpretation="character_state_transitions_not_mutation_events", catalogue_version=1))
    dependency_rows = [{"family_id": family, "dependency_group": group,
        "character_ids": ";".join(character_id(key) for key in sorted(keys)),
        "n_characters": len(keys),
        "basis": "overlapping_aligned_intervals" if group.startswith("overlap:") else "declared_linked_characters",
        "shared_mutation_inferred": "false"}
        for (family, group), keys in sorted(dependencies.items())]
    directory = Path(output_dir)
    write_tsv(directory / "character_catalogue.tsv", catalogue, CATALOG_FIELDS)
    write_tsv(directory / "character_coordinates.tsv", coordinates, COORDINATE_FIELDS)
    write_tsv(directory / "character_dependencies.tsv", dependency_rows,
              ["family_id", "dependency_group", "character_ids", "n_characters", "basis", "shared_mutation_inferred"])
    return catalogue
