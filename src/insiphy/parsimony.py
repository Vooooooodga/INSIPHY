"""Qualitative parsimony reconstruction for single-copy structural sites."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

from .io import (
    STRUCTURAL_SITE_SCHEMA_VERSION,
    read_structural_site_matrix,
    read_tsv,
    structural_site_observed_state,
    validate_structural_site_tip_rows,
    write_structural_site_matrix,
    write_tsv,
)
from .structural_sites import build_structural_site_matrix, structural_matrix_annotation_view
from .tree import SpeciesTree as RootedTree


STATES = (0, 1)
PARSIMONY_INFERENCE_METHOD = "equal_cost_maximum_parsimony"


def _fmt(value):
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):.8g}"


def _tree_rows_for_topology(path):
    rows = read_tsv(path, ["node_id", "parent_id", "label"])
    out = []
    for row in rows:
        item = dict(row)
        item["branch_length"] = 0.0 if not row.get("parent_id") else 1.0
        out.append(item)
    return out


def _site_key(row):
    return (row.get("family_id", "NA"), row.get("layer", "NA"), row.get("site_id", "NA"))


def _state_index(state, state_0, state_1):
    if state == state_0:
        return 0
    if state == state_1:
        return 1
    return "unknown"


def _state_name(index, state_0, state_1):
    return state_0 if index == 0 else state_1


def _cost(src, dst):
    return 0.0 if src == dst else 1.0


def _minimize(values):
    best = min(values.values())
    return best, {state for state, score in values.items() if abs(score - best) <= 1e-9}


def _is_global_optimum(score, minimum):
    return math.isfinite(score) and abs(score - minimum) <= 1e-9


def _parsimony_tables(tree, observations):
    allowed = {}
    for label in tree.leaf_by_label:
        value = observations.get(label, "unknown")
        allowed[tree.leaf_by_label[label]] = set(STATES) if value == "unknown" else {int(value)}

    inside = {}
    sibling_costs = {}
    for node in tree.postorder():
        inside[node] = {}
        if not tree.children.get(node):
            node_allowed = allowed.get(node, set(STATES))
            for state in STATES:
                inside[node][state] = 0.0 if state in node_allowed else math.inf
            continue

        children = list(tree.children[node])
        child_messages = []
        for child in children:
            message = {}
            for state in STATES:
                values = {
                    child_state: inside[child][child_state] + _cost(state, child_state)
                    for child_state in STATES
                }
                message[state], _best_states = _minimize(values)
            child_messages.append(message)

        prefix = [{state: 0.0 for state in STATES}]
        for message in child_messages:
            prefix.append({state: prefix[-1][state] + message[state] for state in STATES})
        suffix = [{state: 0.0 for state in STATES} for _index in range(len(children) + 1)]
        for index in range(len(children) - 1, -1, -1):
            suffix[index] = {
                state: suffix[index + 1][state] + child_messages[index][state]
                for state in STATES
            }

        inside[node] = dict(prefix[-1])
        for index, child in enumerate(children):
            sibling_costs[(node, child)] = {
                state: prefix[index][state] + suffix[index + 1][state]
                for state in STATES
            }

    outside = {tree.root: {0: 0.0, 1: 0.0}}
    for parent in tree.preorder():
        for child in tree.children.get(parent, []):
            outside[child] = {}
            for child_state in STATES:
                candidates = []
                for parent_state in STATES:
                    candidates.append(
                        outside[parent][parent_state]
                        + sibling_costs[(parent, child)][parent_state]
                        + _cost(parent_state, child_state)
                    )
                outside[child][child_state] = min(candidates)

    minimum = min(inside[tree.root].values())
    node_states = {}
    for node in tree.preorder():
        values = {state: outside[node][state] + inside[node][state] for state in STATES}
        node_states[node] = {
            state for state, score in values.items() if _is_global_optimum(score, minimum)
        }

    branch_pairs = {}
    for parent, child in tree.edges():
        pairs = set()
        for parent_state in STATES:
            sibling_score = sibling_costs[(parent, child)][parent_state]
            for child_state in STATES:
                score = (
                    outside[parent][parent_state]
                    + sibling_score
                    + _cost(parent_state, child_state)
                    + inside[child][child_state]
                )
                if _is_global_optimum(score, minimum):
                    pairs.add((parent_state, child_state))
        branch_pairs[(parent, child)] = pairs
    return minimum, node_states, branch_pairs


def _pattern_class(observations):
    observed = [value for value in observations.values() if value in STATES]
    if not observed:
        return "all_missing"
    if len(set(observed)) < 2:
        return "no_observed_contrast"
    return "observed_contrast"


def _observation_summary(observations):
    count_0 = sum(1 for value in observations.values() if value == 0)
    count_1 = sum(1 for value in observations.values() if value == 1)
    missing = sum(1 for value in observations.values() if value == "unknown")
    return count_0, count_1, missing, f"state0={count_0};state1={count_1};missing={missing}"


def _within_exon_boundary_sites(site_rows, output_dir):
    return {
        _site_key(row)
        for row in site_rows
        if row.get("site_kind") in {"within_element_junction", "within_exon_boundary"}
    }


def _event_type(layer, site_key, parent_state, child_state, within_exon_boundary):
    if parent_state == child_state:
        return "no_change"
    if layer == "splice_junction":
        return "intron_gain" if (parent_state, child_state) == (0, 1) else "intron_loss"
    if layer == "exon_presence":
        return "sequence_gain" if (parent_state, child_state) == (0, 1) else "sequence_loss"
    if layer == "exon_role":
        return "exon_role_gain" if (parent_state, child_state) == (0, 1) else "exon_role_loss"
    return "state_gain" if (parent_state, child_state) == (0, 1) else "state_loss"


def _structural_relation(layer, site_key, parent_state, child_state, within_exon_boundary):
    if layer != "splice_junction" or site_key not in within_exon_boundary:
        return "NA"
    if (parent_state, child_state) == (0, 1):
        return "exon_split"
    if (parent_state, child_state) == (1, 0):
        return "exon_fusion"
    return "NA"


def _placement_status(pattern_class, pair, endpoint_pairs):
    if pattern_class == "all_missing":
        return "uninformative"
    if pair[0] == pair[1]:
        return "no_change"
    return "required" if all(candidate == pair for candidate in endpoint_pairs) else "possible"


def _branch_support_kind(layer, placement_status):
    if layer == "splice_junction" and placement_status == "possible":
        return "possible_non_joint"
    return placement_status


def _node_status(pattern_class, optimal_states):
    if pattern_class == "all_missing":
        return "uninformative"
    if len(optimal_states) == 1:
        return "in_all_optimal_reconstructions"
    return "in_some_optimal_reconstructions"


def _node_support_kind(pattern_class, optimal_states):
    if pattern_class == "all_missing":
        return "uninformative"
    return "required" if len(optimal_states) == 1 else "possible"


def _site_support_kind(pattern_class):
    if pattern_class == "all_missing":
        return "uninformative"
    if pattern_class == "no_observed_contrast":
        return "no_change"
    return "required"


def _conditional_scope(layer, annotation_view):
    if layer == "exon_presence":
        return "homologous_sequence_unit_presence"
    if layer == "exon_role":
        return f"homologous_sequence_present;supplied_annotation_view={annotation_view}"
    if layer == "splice_junction":
        return f"homologous_sequence_present;supplied_transcript_view={annotation_view}"
    return "layer_defined_binary_state"


def _inference_availability(pattern_class):
    if pattern_class == "all_missing":
        return "not_analyzable", "no_observed_states"
    if pattern_class == "no_observed_contrast":
        return "observed_conserved", "event_direction_not_applicable"
    return "analyzed", "NA"


def _masked_state(row):
    return structural_site_observed_state(row)


def _linked_junction_groups(site_rows):
    groups = defaultdict(set)
    for row in site_rows:
        if row.get("layer") != "splice_junction":
            continue
        linked_group_ids = {
            token
            for token in str(row.get("linked_group_id", "NA")).split(";")
            if token not in {"", "NA"}
            and all(marker in token for marker in ("_SP_", "_GC_", "_TX_"))
        }
        for linked_group_id in linked_group_ids:
            groups[(row.get("family_id", "NA"), linked_group_id)].add(
                row.get("site_id", "NA")
            )
    return {
        key: tuple(sorted(site_ids))
        for key, site_ids in groups.items()
        if len(site_ids) >= 2
    }


def _compound_event_rows(
    site_rows,
    branch_rows,
    tree,
    annotation_view,
    matrix_schema_version,
):
    """Summarize linked junction changes without collapsing their source rows."""

    linked_groups = _linked_junction_groups(site_rows)
    events_by_site_branch_pair = {
        (
            row["family_id"],
            row["site_id"],
            row["parent_id"],
            row["child_id"],
            int(row["parent_state_index"]),
            int(row["child_state_index"]),
        ): row
        for row in branch_rows
        if row["layer"] == "splice_junction"
    }
    compound_rows = []
    for (family_id, linked_group_id), site_ids in sorted(linked_groups.items()):
        junction_count = len(site_ids)
        segment_count = junction_count + 1
        for parent, child in tree.edges():
            for parent_state, child_state, direction in (
                (0, 1, "intron_gain"),
                (1, 0, "intron_loss"),
            ):
                components = []
                for site_id in site_ids:
                    component = events_by_site_branch_pair.get(
                        (family_id, site_id, parent, child, parent_state, child_state)
                    )
                    if component is None or component["placement_status"] not in {
                        "required",
                        "possible",
                    }:
                        components = []
                        break
                    components.append(component)
                if len(components) != junction_count:
                    continue

                support_kinds = [row["placement_status"] for row in components]
                all_required = all(kind == "required" for kind in support_kinds)
                pattern_status = (
                    "supported_compound_pattern"
                    if all_required
                    else "candidate_compound_pattern"
                )
                support_kind = "required" if all_required else "possible_non_joint"
                if direction == "intron_gain":
                    compound_event_type = f"one_to_{segment_count}_split"
                else:
                    compound_event_type = f"{segment_count}_to_one_fusion"
                compound_rows.append(
                    {
                        "family_id": family_id,
                        "linked_group_id": linked_group_id,
                        "layer": "splice_junction",
                        "parent_id": parent,
                        "child_id": child,
                        "parent_label": tree.label[parent],
                        "child_label": tree.label[child],
                        "branch_scope": f"{tree.label[parent]}->{tree.label[child]}",
                        "compound_event_type": compound_event_type,
                        "pattern_status": pattern_status,
                        "support_kind": support_kind,
                        "change_direction": direction,
                        "junction_count": junction_count,
                        "segment_count": segment_count,
                        "component_site_ids": ";".join(site_ids),
                        "component_event_types": ";".join(
                            f"{row['site_id']}:{row['event_type']}" for row in components
                        ),
                        "component_support_kinds": ";".join(
                            f"{row['site_id']}:{row['placement_status']}" for row in components
                        ),
                        "joint_support": (
                            "all_components_required_on_same_branch_and_direction"
                            if all_required
                            else "not_established_from_marginal_site_support"
                        ),
                        "inference_method": PARSIMONY_INFERENCE_METHOD,
                        "conditional_scope": _conditional_scope(
                            "splice_junction", annotation_view
                        ),
                        "unavailable_reason": (
                            "NA"
                            if all_required
                            else "joint_history_not_inferred_from_possible_components"
                        ),
                        "matrix_schema_version": matrix_schema_version,
                        "annotation_view": annotation_view,
                    }
                )
    return compound_rows


def _encode_sites(site_rows, tree):
    grouped = defaultdict(dict)
    labels = {}
    for row in site_rows:
        key = _site_key(row)
        grouped[key][row.get("species", "NA")] = _masked_state(row)
        labels[key] = (row.get("state_0", "state_0"), row.get("state_1", "state_1"))

    encoded = []
    for (family_id, layer, site_id), species_states in sorted(grouped.items()):
        state_0, state_1 = labels[(family_id, layer, site_id)]
        observations = {
            label: _state_index(species_states[label], state_0, state_1)
            for label in tree.leaf_by_label
        }
        encoded.append((family_id, layer, site_id, state_0, state_1, observations))
    return encoded


def _load_site_matrix(
    input_dir,
    output_dir,
    structural_site_matrix_path,
    annotation_view,
):
    if structural_site_matrix_path is not None:
        source = Path(structural_site_matrix_path)
        site_rows = read_structural_site_matrix(source)
        excluded_rows = read_tsv(source.parent / "excluded_families.tsv", optional=True)
        output_matrix = Path(output_dir) / "structural_site_matrix.tsv"
        if source.resolve() != output_matrix.resolve():
            write_structural_site_matrix(output_matrix, site_rows)
        return site_rows, excluded_rows, source, "frozen_matrix"
    site_rows, excluded_rows = build_structural_site_matrix(
        input_dir,
        output_dir,
        annotation_view=annotation_view,
    )
    source = Path(output_dir) / "structural_site_matrix.tsv"
    write_structural_site_matrix(source, site_rows)
    return read_structural_site_matrix(source), excluded_rows, source, "prepared_in_analysis"


def infer_single_copy_parsimony(
    input_dir,
    output_dir,
    threads=1,
    structural_site_matrix_path=None,
    annotation_view="repertoire",
):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    threads = max(1, int(threads))
    site_rows, excluded_rows, matrix_source, matrix_mode = _load_site_matrix(
        input_dir,
        output_dir,
        structural_site_matrix_path,
        annotation_view,
    )
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

    compound_rows = _compound_event_rows(
        site_rows,
        branch_rows,
        tree,
        effective_annotation_view,
        matrix_schema_version,
    )

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
        output_dir / "compound_structural_events.tsv",
        compound_rows,
        [
            "family_id", "linked_group_id", "layer", "parent_id", "child_id",
            "parent_label", "child_label", "branch_scope", "compound_event_type",
            "pattern_status", "support_kind", "change_direction", "junction_count",
            "segment_count", "component_site_ids", "component_event_types",
            "component_support_kinds", "joint_support", "inference_method",
            "conditional_scope", "unavailable_reason", "matrix_schema_version",
            "annotation_view",
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
            "compound_structural_events.tsv",
            "structural_site_summary.tsv",
            "phylogeny_scope.tsv",
            "excluded_families.tsv",
        ],
        "excluded_family_count": len({row.get("family_id", "NA") for row in excluded_rows}),
        "structural_site_count": len(summary_rows),
        "compound_structural_event_count": len(compound_rows),
    }
    (output_dir / "run_parameters.json").write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")
    return node_rows, branch_rows, summary_rows
