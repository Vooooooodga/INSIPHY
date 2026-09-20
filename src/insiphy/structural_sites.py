"""structural_sites: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from insiphy.observations.catalogue import _complete_tree_tip_observations
from insiphy.observations.catalogue import _merge_catalogue_observations
from insiphy.observations.catalogue import _structural_site_universe
from insiphy.observations.junctions import _junction_site_rows
from insiphy.observations.roles import _element_site_rows
from insiphy.observations.schema import STRUCTURAL_SITE_ANNOTATION_VIEWS
from insiphy.observations.schema import STRUCTURAL_SITE_SCHEMA_VERSION
from insiphy.observations.schema import STRUCTURAL_SITE_SCHEMA_VERSIONS
from insiphy.observations.schema import normalize_structural_site_row
from insiphy.observations.schema import structural_site_annotation_view
from insiphy.observations.support import _single_copy_families
from insiphy.storage.tabular import read_tsv
from insiphy.topology import SpeciesTree
from pathlib import Path


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


# Backward-compatible symbol exports; no alternate implementations.
from insiphy.observations.catalogue import (
    _structural_site_universe,
    _state_labels_for_layer,
    _merge_catalogue_observations,
    _complete_tree_tip_observations,
)
from insiphy.observations.junctions import (
    _read_match_rows,
    _occurrence_length,
    _position_eligible,
    _legacy_projected_reference_blocks,
    _candidate_record_blocks,
    _retained_candidate_reference_blocks,
    _parse_projected_reference_blocks,
    _project_query_base_to_reference,
    _project_target_base_to_reference,
    _mapped_reference_coordinate,
    _within_exon_boundary,
    _continuous_reference_block_covers,
    _genomically_contiguous,
    _phase_status,
    _explicit_intron_between,
    _reference_by_element,
    _member_position_eligible,
    _boundary_from_pair,
    _junction_site_rows,
)
from insiphy.observations.roles import (
    _element_site_rows,
)
from insiphy.observations.support import (
    EXONIC_ROLES,
    KNOWN_NONEXONIC_ROLES,
    MISSING,
    _single_copy_families,
    _tokens,
    _is_true,
    _presence_state,
    _json_records,
    _record_interval,
    _parse_genomic_blocks,
    _membership_blocks,
    _legacy_occurrence_blocks,
    _completion_blocks,
    _blocks_overlap_status,
    _blocks_overlap,
    _completion_presence,
    _completion_role_prediction_supported,
    _exon_prediction_overlaps_observation,
    _path_complete,
    _classify_path_blocks,
    _role_from_transcript_paths,
    _observation_row,
    _finalize_site_metadata,
)
import json
