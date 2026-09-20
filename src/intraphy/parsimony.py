"""parsimony: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from intraphy.inference.parsimony_export import _branch_support_kind
from intraphy.inference.parsimony_export import _conditional_scope
from intraphy.inference.parsimony_export import _encode_sites
from intraphy.inference.parsimony_export import _event_type
from intraphy.inference.parsimony_export import _inference_availability
from intraphy.inference.parsimony_export import _load_site_matrix
from intraphy.inference.parsimony_export import _node_status
from intraphy.inference.parsimony_export import _node_support_kind
from intraphy.inference.parsimony_export import _observation_summary
from intraphy.inference.parsimony_export import _pattern_class
from intraphy.inference.parsimony_export import _placement_status
from intraphy.inference.parsimony_export import _site_key
from intraphy.inference.parsimony_export import _site_support_kind
from intraphy.inference.parsimony_export import _state_name
from intraphy.inference.parsimony_export import _structural_relation
from intraphy.inference.parsimony_export import _tree_rows_for_topology
from intraphy.inference.parsimony_export import _within_exon_boundary_sites
from intraphy.inference.sankoff import PARSIMONY_INFERENCE_METHOD
from intraphy.inference.sankoff import _fmt
from intraphy.inference.sankoff import _parsimony_tables
from intraphy.observations.schema import STRUCTURAL_SITE_SCHEMA_VERSION
from intraphy.observations.schema import validate_structural_site_tip_rows
from intraphy.run_result import record_run_result
from intraphy.storage.tabular import write_tsv
from intraphy.structural_sites import structural_matrix_annotation_view
from intraphy.topology import SpeciesTree as RootedTree
from pathlib import Path
import json


def infer_single_copy_parsimony(
    input_dir,
    output_dir,
    threads=1,
    structural_site_matrix_path=None,
    annotation_view="repertoire",
    *, observation_matrix=None, analysis_range="all", min_callable_fraction=0.70,
):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    threads = max(1, int(threads))
    # This obsolete artifact must never be mistaken for a current inference.
    (output_dir / "compound_structural_events.tsv").unlink(missing_ok=True)
    if observation_matrix is None:
        site_rows, excluded_rows, matrix_source, matrix_mode = _load_site_matrix(
            input_dir, output_dir, structural_site_matrix_path, annotation_view,
            analysis_range=analysis_range, min_callable_fraction=min_callable_fraction,
        )
    else:
        site_rows, excluded_rows, matrix_source, matrix_mode = observation_matrix.as_legacy_tuple()
    matrix_annotation_view = structural_matrix_annotation_view(site_rows)
    effective_annotation_view = (
        annotation_view if matrix_annotation_view == "view_independent" else matrix_annotation_view
    )
    if matrix_annotation_view in {"canonical", "repertoire"} and matrix_annotation_view != annotation_view:
        raise SystemExit(
            "requested annotation view does not match the frozen structural-site matrix: "
            f"requested={annotation_view}, matrix={matrix_annotation_view}"
        )
    write_tsv(
        output_dir / "excluded_families.tsv",
        excluded_rows,
        ["family_id", "species", "copy_count", "gene_copy_ids", "reason"],
    )

    tree = RootedTree(_tree_rows_for_topology(input_dir / "species_tree.tsv"))
    validate_structural_site_tip_rows(site_rows, tree.leaf_by_label)

    matrix_versions = sorted({row.get("schema_version", "NA") for row in site_rows})
    matrix_schema_version = ";".join(matrix_versions) if matrix_versions else STRUCTURAL_SITE_SCHEMA_VERSION
    within_exon_boundary = _within_exon_boundary_sites(site_rows, output_dir)
    metadata_by_site = defaultdict(lambda: defaultdict(set))
    for row in site_rows:
        key = _site_key(row)
        for field in ("discovery_rule", "observation_mask", "linked_group_id", "site_kind"):
            metadata_by_site[key][field].add(str(row.get(field, "NA")))
    node_rows = []
    branch_rows = []
    summary_rows = []
    cached = {}
    histories = []
    from intraphy.inference.history import history_rows
    for family_id, layer, site_id, state_0, state_1, observations in _encode_sites(site_rows, tree):
        pattern_key = tuple((label, observations[label]) for label in sorted(tree.leaf_by_label))
        cache_key = (layer, state_0, state_1, pattern_key)
        if cache_key not in cached:
            cached[cache_key] = _parsimony_tables(tree, observations)
        minimum, node_states, branch_pairs = cached[cache_key]
        pattern_class = _pattern_class(observations)
        inference_status, unavailable_reason = _inference_availability(pattern_class)
        count_0, count_1, missing, summary = _observation_summary(observations)
        site_key = (family_id, layer, site_id)
        histories.extend(history_rows(site_key, tree, observations, minimum, node_states, branch_pairs))
        metadata = {
            field: ";".join(sorted(values))
            for field, values in metadata_by_site.get(site_key, {}).items()
        }
        summary_rows.append(
            {
                "family_id": family_id,
                "layer": layer,
                "site_id": site_id,
                "state_0": state_0,
                "state_1": state_1,
                "observed_state_0_count": count_0,
                "observed_state_1_count": count_1,
                "missing_count": missing,
                "observed_tip_count": count_0 + count_1,
                "pattern_class": pattern_class,
                "inference_method": PARSIMONY_INFERENCE_METHOD,
                "support_kind": _site_support_kind(pattern_class),
                "conditional_scope": _conditional_scope(layer, effective_annotation_view),
                "inference_status": inference_status,
                "unavailable_reason": unavailable_reason,
                "min_changes": _fmt(minimum),
                "compressed_pattern_key": ";".join(f"{label}:{value}" for label, value in pattern_key),
                "observation_summary": summary,
                "discovery_rule": metadata.get("discovery_rule", "NA"),
                "observation_mask": metadata.get("observation_mask", "NA"),
                "linked_group_id": metadata.get("linked_group_id", "NA"),
                "site_kind": metadata.get("site_kind", "NA"),
                "matrix_schema_version": matrix_schema_version,
            }
        )
        for node in tree.preorder():
            optimal_states = node_states[node]
            node_rows.append(
                {
                    "family_id": family_id,
                    "layer": layer,
                    "site_id": site_id,
                    "node_id": node,
                    "node_label": tree.label[node],
                    "state": ";".join(_state_name(state, state_0, state_1) for state in sorted(optimal_states)),
                    "state_indices": ";".join(str(state) for state in sorted(optimal_states)),
                    "placement_status": _node_status(pattern_class, optimal_states),
                    "inference_method": PARSIMONY_INFERENCE_METHOD,
                    "support_kind": _node_support_kind(pattern_class, optimal_states),
                    "conditional_scope": _conditional_scope(layer, effective_annotation_view),
                    "inference_status": inference_status,
                    "unavailable_reason": unavailable_reason,
                    "min_changes": _fmt(minimum),
                    "observation_summary": summary,
                    "linked_group_id": metadata.get("linked_group_id", "NA"),
                    "site_kind": metadata.get("site_kind", "NA"),
                }
            )
        for parent, child in tree.edges():
            pairs = sorted(branch_pairs[(parent, child)])
            for parent_state, child_state in pairs:
                placement_status = _placement_status(
                    pattern_class, (parent_state, child_state), pairs
                )
                branch_rows.append(
                    {
                        "family_id": family_id,
                        "layer": layer,
                        "site_id": site_id,
                        "parent_id": parent,
                        "child_id": child,
                        "parent_label": tree.label[parent],
                        "child_label": tree.label[child],
                        "branch_scope": f"{tree.label[parent]}->{tree.label[child]}",
                        "event_type": _event_type(layer, site_key, parent_state, child_state, within_exon_boundary),
                        "structural_relation": _structural_relation(
                            layer,
                            site_key,
                            parent_state,
                            child_state,
                            within_exon_boundary,
                        ),
                        "placement_status": placement_status,
                        "inference_method": PARSIMONY_INFERENCE_METHOD,
                        "support_kind": _branch_support_kind(layer, placement_status),
                        "conditional_scope": _conditional_scope(layer, effective_annotation_view),
                        "parent_state": _state_name(parent_state, state_0, state_1),
                        "child_state": _state_name(child_state, state_0, state_1),
                        "endpoint_pairs": f"{_state_name(parent_state, state_0, state_1)}->{_state_name(child_state, state_0, state_1)}",
                        "parent_state_index": parent_state,
                        "child_state_index": child_state,
                        "optimal_endpoint_pair_count": len(pairs),
                        "inference_status": inference_status,
                        "unavailable_reason": unavailable_reason,
                        "min_changes": _fmt(minimum),
                        "observation_summary": summary,
                        "linked_group_id": metadata.get("linked_group_id", "NA"),
                        "site_kind": metadata.get("site_kind", "NA"),
                    }
                )

    from intraphy.inference.change_summary import annotate_changes, write_change_summary
    annotate_changes(branch_rows)
    write_change_summary(output_dir, summary_rows, branch_rows, histories, node_rows)

    write_tsv(
        output_dir / "node_structural_states.tsv",
        node_rows,
        [
            "family_id", "layer", "site_id", "node_id", "node_label", "state", "state_indices",
            "placement_status", "inference_method", "support_kind", "conditional_scope",
            "inference_status", "unavailable_reason", "min_changes", "observation_summary",
            "linked_group_id", "site_kind",
        ],
    )
    write_tsv(
        output_dir / "branch_structural_events.tsv",
        branch_rows,
        [
            "family_id", "layer", "site_id", "parent_id", "child_id", "parent_label", "child_label",
            "character_id", "count_unit_id", "event_id", "count_interpretation",
            "branch_scope", "event_type", "structural_relation", "placement_status",
            "parent_state", "child_state",
            "inference_method", "support_kind", "conditional_scope", "endpoint_pairs",
            "parent_state_index", "child_state_index", "optimal_endpoint_pair_count",
            "inference_status", "unavailable_reason", "min_changes", "observation_summary",
            "linked_group_id", "site_kind",
        ],
    )
    write_tsv(
        output_dir / "structural_site_summary.tsv",
        summary_rows,
        [
            "family_id", "layer", "site_id", "state_0", "state_1", "observed_state_0_count",
            "observed_state_1_count", "missing_count", "observed_tip_count", "pattern_class",
            "inference_method", "support_kind", "conditional_scope", "inference_status",
            "unavailable_reason", "min_changes", "compressed_pattern_key", "observation_summary",
            "discovery_rule", "observation_mask", "linked_group_id", "site_kind",
            "matrix_schema_version",
        ],
    )
    write_tsv(
        output_dir / "phylogeny_scope.tsv",
        [{
            "scope": "single_copy_structural_sites",
            "tree_file": "species_tree.tsv",
            "tree_scope": "species_tree",
            "layers": "exon_presence;exon_role;splice_junction",
            "matrix_source": str(matrix_source),
            "matrix_mode": matrix_mode,
            "matrix_schema_version": matrix_schema_version,
            "annotation_view": effective_annotation_view,
            "inference_method": PARSIMONY_INFERENCE_METHOD,
            "support_kind": "all_optimal_histories",
            "conditional_scope": "layer_specific;see_output_rows",
            "note": "Qualitative equal-cost parsimony on the supplied rooted topology; missing observations allow both states.",
        }],
        [
            "scope", "tree_file", "tree_scope", "layers", "matrix_source", "matrix_mode",
            "matrix_schema_version", "annotation_view", "inference_method", "support_kind",
            "conditional_scope", "note",
        ],
    )
    parameters = {
        "analysis_scope": "single-copy",
        "scope": "single_copy_structural_sites",
        "model": "parsimony",
        "rootedtree": str(input_dir / "species_tree.tsv"),
        "tree_file": str(input_dir / "species_tree.tsv"),
        "tree_scope": "species_tree",
        "threads": threads,
        "structural_site_matrix": str(matrix_source),
        "structural_site_matrix_mode": matrix_mode,
        "structural_site_matrix_schema_version": matrix_schema_version,
        "annotation_view": effective_annotation_view,
        "cost_assumptions": {
            "no_change_cost": 0.0,
            "gain_cost": 1.0,
            "loss_cost": 1.0,
            "missing_observation": "permits_both_states",
            "branch_lengths": "ignored_for_equal_cost_qualitative_parsimony",
            "topology": "supplied_rooted_tree",
        },
        "outputs": [
            "structural_site_matrix.tsv",
            "node_structural_states.tsv",
            "branch_structural_events.tsv",
            "character_catalogue.tsv",
            "character_coordinates.tsv",
            "character_dependencies.tsv",
            "gene_change_summary.tsv",
            "minimum_change_history.tsv",
            "ancestral_state_consistency.tsv",
            "structural_site_summary.tsv",
            "phylogeny_scope.tsv",
            "excluded_families.tsv",
        ],
        "excluded_family_count": len({row.get("family_id", "NA") for row in excluded_rows}),
        "structural_site_count": len(summary_rows),
        "mutation_event_count": "not_estimated",
    }
    (output_dir / "run_parameters.json").write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")
    record_run_result(output_dir, input_dir, "parsimony", "single-copy")
    return node_rows, branch_rows, summary_rows


# Public entry points; implementations have a single owner.
from intraphy.inference.parsimony_export import (
    _tree_rows_for_topology,
    _site_key,
    _state_index,
    _state_name,
    _pattern_class,
    _observation_summary,
    _within_exon_boundary_sites,
    _event_type,
    _structural_relation,
    _placement_status,
    _branch_support_kind,
    _node_status,
    _node_support_kind,
    _site_support_kind,
    _conditional_scope,
    _inference_availability,
    _masked_state,
    _encode_sites,
    _load_site_matrix,
)
from intraphy.inference.sankoff import (
    STATES,
    PARSIMONY_INFERENCE_METHOD,
    _fmt,
    _cost,
    _minimize,
    _is_global_optimum,
    _parsimony_tables,
)
import json
import math

from intraphy.structural_sites import build_structural_site_matrix
