"""inference / parsimony_export: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from intraphy.inference.sankoff import PARSIMONY_INFERENCE_METHOD
from intraphy.inference.sankoff import STATES
from intraphy.observations.matrix import load_site_matrix
from intraphy.observations.schema import structural_site_observed_state
from intraphy.storage.tabular import read_tsv


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


def _load_site_matrix(input_dir, output_dir, structural_site_matrix_path, annotation_view,
                      *, analysis_range="all", min_callable_fraction=0.70):
    return load_site_matrix(input_dir, output_dir, structural_site_matrix_path, annotation_view,
                           analysis_range=analysis_range, min_callable_fraction=min_callable_fraction)
