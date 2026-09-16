"""Qualitative parsimony reconstruction for single-copy structural sites."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

from .io import read_tsv, write_tsv
from .structural_sites import build_structural_site_matrix
from .tree import SpeciesTree as RootedTree


STATES = (0, 1)


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


def _parsimony_tables(tree, observations):
    allowed = {}
    for label in tree.leaf_by_label:
        value = observations.get(label, "unknown")
        allowed[tree.leaf_by_label[label]] = set(STATES) if value == "unknown" else {int(value)}

    inside = {}
    child_terms = {}
    for node in tree.postorder():
        inside[node] = {}
        child_terms[node] = {}
        if not tree.children.get(node):
            node_allowed = allowed.get(node, set(STATES))
            for state in STATES:
                inside[node][state] = 0.0 if state in node_allowed else math.inf
            continue
        for state in STATES:
            total = 0.0
            child_terms[node][state] = {}
            for child in tree.children[node]:
                values = {
                    child_state: inside[child][child_state] + _cost(state, child_state)
                    for child_state in STATES
                }
                best, best_states = _minimize(values)
                total += best
                child_terms[node][state][child] = (best, best_states)
            inside[node][state] = total

    outside = {tree.root: {0: 0.0, 1: 0.0}}
    for parent in tree.preorder():
        for child in tree.children.get(parent, []):
            outside[child] = {}
            for child_state in STATES:
                candidates = []
                for parent_state in STATES:
                    siblings = inside[parent][parent_state] - child_terms[parent][parent_state][child][0]
                    candidates.append(outside[parent][parent_state] + siblings + _cost(parent_state, child_state))
                outside[child][child_state] = min(candidates)

    minimum = min(inside[tree.root].values())
    node_states = {}
    for node in tree.preorder():
        values = {state: outside[node][state] + inside[node][state] for state in STATES}
        best = min(values.values())
        node_states[node] = {state for state, score in values.items() if abs(score - best) <= 1e-9}

    branch_pairs = {}
    for parent, child in tree.edges():
        pairs = set()
        for parent_state in STATES:
            sibling_score = inside[parent][parent_state] - child_terms[parent][parent_state][child][0]
            for child_state in STATES:
                score = (
                    outside[parent][parent_state]
                    + sibling_score
                    + _cost(parent_state, child_state)
                    + inside[child][child_state]
                )
                if abs(score - minimum) <= 1e-9:
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


def _within_exon_boundary_sites(output_dir):
    rows = read_tsv(Path(output_dir) / "splice_boundary_correspondence.tsv", optional=True)
    return {
        (row.get("family_id", "NA"), "splice_junction", row.get("site_id", "NA"))
        for row in rows
        if row.get("site_kind") == "within_exon_boundary"
    }


def _event_type(layer, site_key, parent_state, child_state, within_exon_boundary):
    if parent_state == child_state:
        return "no_change"
    if layer == "splice_junction":
        if site_key in within_exon_boundary:
            return "split" if (parent_state, child_state) == (0, 1) else "fusion"
        return "intron_gain" if (parent_state, child_state) == (0, 1) else "intron_loss"
    if layer == "exon_presence":
        return "exon_gain" if (parent_state, child_state) == (0, 1) else "exon_loss"
    if layer == "exon_role":
        return "exonic_role_gain" if (parent_state, child_state) == (0, 1) else "exonic_role_loss"
    return "state_gain" if (parent_state, child_state) == (0, 1) else "state_loss"


def _placement_status(pattern_class, pair, pair_count):
    if pattern_class == "all_missing":
        return "uninformative"
    if pair[0] == pair[1]:
        return "no_change"
    return "required" if pair_count == 1 else "possible"


def _node_status(pattern_class, optimal_states):
    if pattern_class == "all_missing":
        return "uninformative"
    if len(optimal_states) == 1:
        return "in_all_optimal_reconstructions"
    return "in_some_optimal_reconstructions"


def _encode_sites(site_rows, tree):
    grouped = defaultdict(lambda: defaultdict(dict))
    labels = {}
    for row in site_rows:
        key = (row.get("family_id", "NA"), row.get("layer", "NA"))
        grouped[key][row.get("site_id", "NA")][row.get("species", "NA")] = row.get("state", "unknown")
        labels[key] = (row.get("state_0", "state_0"), row.get("state_1", "state_1"))

    encoded = []
    for (family_id, layer), sites in sorted(grouped.items()):
        state_0, state_1 = labels[(family_id, layer)]
        for site_id, species_states in sorted(sites.items()):
            observations = {}
            for label in tree.leaf_by_label:
                observations[label] = _state_index(species_states.get(label, "unknown"), state_0, state_1)
            encoded.append((family_id, layer, site_id, state_0, state_1, observations))
    return encoded


def infer_single_copy_parsimony(input_dir, output_dir, threads=1):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    threads = max(1, int(threads))
    site_rows, excluded_rows = build_structural_site_matrix(input_dir, output_dir)
    write_tsv(
        output_dir / "excluded_families.tsv",
        excluded_rows,
        ["family_id", "species", "copy_count", "gene_copy_ids", "reason"],
    )
    write_tsv(
        output_dir / "structural_site_matrix.tsv",
        site_rows,
        ["family_id", "layer", "site_id", "species", "state", "state_0", "state_1", "evidence"],
    )

    tree = RootedTree(_tree_rows_for_topology(input_dir / "species_tree.tsv"))
    matrix_species = {row.get("species", "NA") for row in site_rows}
    missing_tips = sorted(matrix_species - set(tree.leaf_by_label))
    if missing_tips:
        raise SystemExit("species_tree.tsv lacks structural-matrix species: " + ", ".join(missing_tips))

    within_exon_boundary = _within_exon_boundary_sites(output_dir)
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
        count_0, count_1, missing, summary = _observation_summary(observations)
        site_key = (family_id, layer, site_id)
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
                "min_changes": _fmt(minimum),
                "compressed_pattern_key": ";".join(f"{label}:{value}" for label, value in pattern_key),
                "observation_summary": summary,
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
                    "min_changes": _fmt(minimum),
                    "observation_summary": summary,
                }
            )
        for parent, child in tree.edges():
            pairs = sorted(branch_pairs[(parent, child)])
            for parent_state, child_state in pairs:
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
                        "placement_status": _placement_status(pattern_class, (parent_state, child_state), len(pairs)),
                        "parent_state": _state_name(parent_state, state_0, state_1),
                        "child_state": _state_name(child_state, state_0, state_1),
                        "endpoint_pairs": f"{_state_name(parent_state, state_0, state_1)}->{_state_name(child_state, state_0, state_1)}",
                        "parent_state_index": parent_state,
                        "child_state_index": child_state,
                        "min_changes": _fmt(minimum),
                        "observation_summary": summary,
                    }
                )

    write_tsv(
        output_dir / "node_structural_states.tsv",
        node_rows,
        [
            "family_id", "layer", "site_id", "node_id", "node_label", "state", "state_indices",
            "placement_status", "min_changes", "observation_summary",
        ],
    )
    write_tsv(
        output_dir / "branch_structural_events.tsv",
        branch_rows,
        [
            "family_id", "layer", "site_id", "parent_id", "child_id", "parent_label", "child_label",
            "branch_scope", "event_type", "placement_status", "parent_state", "child_state",
            "endpoint_pairs", "parent_state_index", "child_state_index", "min_changes", "observation_summary",
        ],
    )
    write_tsv(
        output_dir / "structural_site_summary.tsv",
        summary_rows,
        [
            "family_id", "layer", "site_id", "state_0", "state_1", "observed_state_0_count",
            "observed_state_1_count", "missing_count", "observed_tip_count", "pattern_class",
            "min_changes", "compressed_pattern_key", "observation_summary",
        ],
    )
    write_tsv(
        output_dir / "phylogeny_scope.tsv",
        [{
            "scope": "single_copy_structural_sites",
            "tree_file": "species_tree.tsv",
            "tree_scope": "species_tree",
            "layers": "exon_presence;exon_role;splice_junction",
            "note": "Qualitative equal-cost parsimony on the supplied rooted topology; missing observations allow both states.",
        }],
        ["scope", "tree_file", "tree_scope", "layers", "note"],
    )
    parameters = {
        "analysis_scope": "single-copy",
        "scope": "single_copy_structural_sites",
        "model": "parsimony",
        "rootedtree": str(input_dir / "species_tree.tsv"),
        "tree_file": str(input_dir / "species_tree.tsv"),
        "tree_scope": "species_tree",
        "threads": threads,
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
            "structural_site_summary.tsv",
            "phylogeny_scope.tsv",
            "excluded_families.tsv",
        ],
        "excluded_family_count": len({row.get("family_id", "NA") for row in excluded_rows}),
        "structural_site_count": len(summary_rows),
    }
    (output_dir / "run_parameters.json").write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")
    return node_rows, branch_rows, summary_rows
