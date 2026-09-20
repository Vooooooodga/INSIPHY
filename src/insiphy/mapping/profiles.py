"""Membership scoring and observed pair-level correspondence profiles."""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from insiphy.coordinates import ClosedInterval1, Interval0, parse_legacy_blocks
from insiphy.elements import collect_element_profiles, element_class_for_occurrence, element_role_from_occurrences, membership_call
from insiphy.alignment import global_alignment_stats
from insiphy.io import norm_state, parse_fasta, read_tsv, to_float, uniq, write_tsv
from insiphy.topology import SpeciesTree



def match_total(row):
    if row.get("correspondence_score") not in ("", "NA", None):
        return to_float(row["correspondence_score"])
    if row.get("total_score") not in ("", "NA", None):
        return to_float(row.get("total_score"))
    boundary = to_float(row.get("boundary_score"), None)
    if boundary is None:
        boundary = 0.5 * (to_float(row.get("left_boundary_score"), 0.5) + to_float(row.get("right_boundary_score"), 0.5))
    return (
        0.34 * to_float(row.get("alignment_score"))
        + 0.14 * to_float(row.get("coverage_score"), to_float(row.get("size_ratio"), 1.0))
        + 0.10 * to_float(row.get("left_context_score"), 0.5)
        + 0.10 * to_float(row.get("right_context_score"), 0.5)
        + 0.10 * boundary
        + 0.08 * to_float(row.get("phase_score"), 0.5)
        + 0.06 * to_float(row.get("order_score"), 0.5)
        + 0.04 * to_float(row.get("strand_score"), 0.5)
        + 0.04 * to_float(row.get("splice_score"), 0.5)
    )


def normalized_components(row):
    boundary = to_float(row.get("boundary_score"), None)
    if boundary is None:
        boundary = 0.5 * (to_float(row.get("left_boundary_score"), 0.5) + to_float(row.get("right_boundary_score"), 0.5))
    return {
        "alignment_score": to_float(row.get("alignment_score")),
        "coverage_score": to_float(row.get("coverage_score"), to_float(row.get("size_ratio"), 1.0)),
        "sequence_score": to_float(row.get("sequence_score"), None),
        "structural_context_score": to_float(row.get("structural_context_score"), None),
        "left_context_score": to_float(row.get("left_context_score"), 0.5),
        "right_context_score": to_float(row.get("right_context_score"), 0.5),
        "boundary_score": boundary,
        "phase_score": to_float(row.get("phase_score"), 0.5),
        "order_score": to_float(row.get("order_score"), 0.5),
        "strand_score": to_float(row.get("strand_score"), 0.5),
        "splice_score": to_float(row.get("splice_score"), 0.5),
    }


def tree_distances(input_dir):
    """Prioritize correspondence by supplied lengths, or topology if any are missing.

    Unit edges apply to this ordering only; CTMC branch lengths are unchanged.
    """
    rows = read_tsv(f"{input_dir}/species_tree.tsv", optional=True)
    if not rows:
        return {}
    tree = SpeciesTree(rows)
    use_topology = any(tree.length[node] is None for node in tree.parent if node != tree.root)
    depth = {tree.root: 0.0}
    for node in tree.preorder():
        for child in tree.children.get(node, []):
            depth[child] = depth[node] + (1.0 if use_topology else tree.branch_length(child))

    ancestors = {}
    for node in tree.parent:
        vals = []
        cur = node
        while cur:
            vals.append(cur)
            cur = tree.parent.get(cur)
        ancestors[node] = vals

    distances = {}
    labels = sorted(tree.leaf_by_label)
    for left_label in labels:
        for right_label in labels:
            left = tree.leaf_by_label[left_label]
            right = tree.leaf_by_label[right_label]
            right_ancestors = set(ancestors[right])
            lca = next(node for node in ancestors[left] if node in right_ancestors)
            distances[(left_label, right_label)] = depth[left] + depth[right] - 2 * depth[lca]
    return distances


def progressive_tier(distance):
    if distance == "NA":
        return "tree_not_available"
    value = to_float(distance, 0.0)
    if value <= 2.0:
        return "nearest_species_support"
    if value <= 4.0:
        return "within_clade_support"
    return "deep_tree_support"


def reciprocal_status(scored):
    best_subject_by_query = defaultdict(list)
    best_query_by_subject = defaultdict(list)
    for row in scored:
        score = to_float(row["total_score"])
        q = row["query_occurrence_id"]
        s = row["subject_occurrence_id"]
        best_subject_by_query[q].append((score, s, row["match_id"]))
        best_query_by_subject[s].append((score, q, row["match_id"]))
        best_subject_by_query[s].append((score, q, row["match_id"]))
        best_query_by_subject[q].append((score, s, row["match_id"]))

    query_best = {}
    for query, vals in best_subject_by_query.items():
        best = max(score for score, _node, _mid in vals)
        query_best[query] = {_node for score, _node, _mid in vals if abs(score - best) < 1e-12}
    subject_best = {}
    for subject, vals in best_query_by_subject.items():
        best = max(score for score, _node, _mid in vals)
        subject_best[subject] = {_node for score, _node, _mid in vals if abs(score - best) < 1e-12}

    out = {}
    for row in scored:
        q = row["query_occurrence_id"]
        s = row["subject_occurrence_id"]
        if s in query_best.get(q, set()) and q in subject_best.get(s, set()):
            out[row["match_id"]] = "reciprocal_best"
        elif s in query_best.get(q, set()):
            out[row["match_id"]] = "query_best_only"
        else:
            out[row["match_id"]] = "not_reciprocal_best"
    return out
