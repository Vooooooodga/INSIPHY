"""mapping / policies: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.mapping.fields import CorrespondenceCriteria
from intraphy.mapping.fields import DEFAULT_CORRESPONDENCE_CRITERIA
from intraphy.mapping.fields import EXON_LIKE_ROLES
from intraphy.mapping.fields import STRUCTURAL_ROLES
from intraphy.mapping.fields import UNKNOWN_SOURCE_LABELS
from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from intraphy.storage.values import to_float
from pathlib import Path
import json
import math


def occurrence_copy_key(row):
    return (row.get("family_id", ""), row.get("species", ""), row.get("gene_copy_id", ""))


def roles_compatible(left, right):
    left_role = left.get("role", "")
    right_role = right.get("role", "")
    if left_role == right_role:
        return True
    if left_role in EXON_LIKE_ROLES and right_role in EXON_LIKE_ROLES:
        return True
    return left_role in STRUCTURAL_ROLES and right_role in STRUCTURAL_ROLES


def sequence_supported_mapping(
    evidence,
    threshold,
    criteria=DEFAULT_CORRESPONDENCE_CRITERIA,
):
    identity = evidence["alignment_score"]
    coverage = evidence["coverage_score"]
    if identity >= threshold and coverage >= criteria.regular_min_coverage:
        return True
    if (
        identity >= threshold + criteria.high_identity_offset
        and coverage >= criteria.high_identity_min_coverage
    ):
        return True
    return False


def _candidate_sequence_accepted(
    record,
    threshold,
    short_context=False,
    criteria=DEFAULT_CORRESPONDENCE_CRITERIA,
    *, anchored_microexon=False,
):
    identity = to_float(record.get("identity"), 0.0)
    coverage = to_float(record.get("coverage"), 0.0)
    accepted = (
        (identity >= threshold and coverage >= criteria.regular_min_coverage)
        or (
            identity >= threshold + criteria.high_identity_offset
            and coverage >= criteria.high_identity_min_coverage
        )
    )
    if short_context:
        pairs = int(to_float(record.get("known_aligned_pairs", record.get("aligned_pairs", 0)), 0.0))
        short_coverage = to_float(record.get("short_sequence_coverage", record.get("query_coverage")), 0.0)
        if anchored_microexon and 3 <= pairs < criteria.short_min_aligned_pairs:
            # This is a conservative, testable candidate rule, not a homology
            # probability. Reading frame is an output, not an acceptance filter.
            return bool(accepted and identity >= max(0.90, threshold)
                        and short_coverage >= 0.95)
        accepted = bool(
            accepted
            and identity >= max(criteria.short_min_identity, threshold)
            and to_float(
                record.get("short_sequence_coverage", record.get("query_coverage")),
                0.0,
            )
            >= criteria.short_min_query_coverage
            and int(to_float(
                record.get("known_aligned_pairs", record.get("aligned_pairs", 0)),
                0.0,
            ))
            >= criteria.short_min_aligned_pairs
        )
    return accepted


def assess_short_candidate_thresholds(
    segment_matches,
    output_path,
    short_identity_thresholds=(0.60, 0.70, 0.80),
    short_query_coverage_thresholds=(0.60, 0.80),
):
    identities = sorted({float(value) for value in short_identity_thresholds})
    coverages = sorted({float(value) for value in short_query_coverage_thresholds})
    if not identities or not coverages:
        raise ValueError("at least one identity and coverage threshold is required")
    if any(value < 0.0 or value > 1.0 for value in identities + coverages):
        raise ValueError("identity and coverage thresholds must be within [0, 1]")

    rows = []
    for match in read_tsv(segment_matches):
        if match.get("short_context_route") not in {
            "feature_bounded_candidate", "bounded_local", "anchor_bounded_local",
        }:
            continue
        encoded = match.get("dna_candidate_assessments")
        if encoded in {None, "", "NA"}:
            continue
        candidates = json.loads(encoded)
        if not isinstance(candidates, list):
            raise ValueError("dna_candidate_assessments must encode a JSON list")
        for candidate in candidates:
            source = candidate.get("source", "")
            score_scheme = candidate.get("score_scheme", "")
            if source != "nucleotide_alignment" or not str(score_scheme).startswith("nt_"):
                raise ValueError(
                    "dna_candidate_assessments contains a non-nucleotide candidate: "
                    f"source={source!r}, score_scheme={score_scheme!r}"
                )
            saved_pair_threshold = to_float(candidate.get("acceptance_threshold"), None)
            if saved_pair_threshold is None:
                raise ValueError(
                    "DNA candidate is missing its saved pair-specific acceptance_threshold"
                )
            for identity in identities:
                for coverage in coverages:
                    criteria = CorrespondenceCriteria(
                        short_min_identity=identity,
                        short_min_query_coverage=coverage,
                    )
                    effective_identity_cutoff = max(saved_pair_threshold, identity)
                    rows.append({
                        "match_id": match.get("match_id", "NA"),
                        "query_occurrence_id": match.get("query_occurrence_id", "NA"),
                        "subject_occurrence_id": match.get("subject_occurrence_id", "NA"),
                        "candidate_id": candidate.get("candidate_id", "NA"),
                        "saved_pair_threshold": f"{saved_pair_threshold:.6g}",
                        "short_identity_threshold": f"{identity:.6g}",
                        "effective_identity_cutoff": f"{effective_identity_cutoff:.6g}",
                        "query_coverage_threshold": f"{coverage:.6g}",
                        "candidate_identity": candidate.get("identity", "NA"),
                        "candidate_query_coverage": candidate.get(
                            "short_sequence_coverage", candidate.get("query_coverage", "NA")
                        ),
                        "candidate_aligned_pairs": candidate.get("aligned_pairs", "NA"),
                        "acceptance": int(_candidate_sequence_accepted(
                            candidate,
                            saved_pair_threshold,
                            short_context=True,
                            criteria=criteria,
                        )),
                        "interpretation": "short_DNA_acceptance_rule_sensitivity_only",
                    })
    write_tsv(
        output_path,
        rows,
        [
            "match_id", "query_occurrence_id", "subject_occurrence_id",
            "candidate_id", "saved_pair_threshold", "short_identity_threshold",
            "effective_identity_cutoff", "query_coverage_threshold",
            "candidate_identity", "candidate_query_coverage",
            "candidate_aligned_pairs", "acceptance", "interpretation",
        ],
    )
    return rows


def split_source_labels(value):
    labels = []
    for part in str(value or "").replace("|", ";").replace(",", ";").split(";"):
        label = part.strip()
        if label and label not in UNKNOWN_SOURCE_LABELS:
            labels.append(label)
    return sorted(set(labels))


def known_source_labels(row):
    if row.get("copy_role") == "derived":
        return []
    return split_source_labels(row.get("source_label"))


def inferred_source_label(row, support):
    own = split_source_labels(row.get("source_label"))
    if own and row.get("copy_role") != "derived":
        return ";".join(own)
    if not support:
        return "unknown_source"
    best = max(support.values())
    top = sorted(source for source, score in support.items() if score >= best * 0.90)
    return ";".join(top) if top else "unknown_source"


def graph_components(nodes, edges, occurrence_by_id=None, species_distances=None, same_copy_compatible=None, cross_copy_compatible=None):
    if occurrence_by_id is None:
        import networkx as nx
        graph = nx.Graph()
        graph.add_nodes_from(nodes)
        graph.add_weighted_edges_from(edges)
        return [sorted(component) for component in nx.connected_components(graph)]
    parent = {node: node for node in nodes}
    members = {node: {node} for node in nodes}
    direct = {frozenset((left, right)) for left, right, _score in edges}
    same_copy_compatible = same_copy_compatible or set()
    cross_copy_compatible = cross_copy_compatible or set()

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def compatible(left_root, right_root):
        for left_id in members[left_root]:
            for right_id in members[right_root]:
                left = occurrence_by_id.get(left_id, {})
                right = occurrence_by_id.get(right_id, {})
                if occurrence_copy_key(left) == occurrence_copy_key(right):
                    if frozenset((left_id, right_id)) not in same_copy_compatible:
                        return False
                elif (
                    frozenset((left_id, right_id)) not in direct
                    and frozenset((left_id, right_id)) not in cross_copy_compatible
                ):
                    return False
        return True

    species_distances = species_distances or {}
    def progressive_key(edge):
        left, right, score = edge
        left_species = occurrence_by_id.get(left, {}).get("species")
        right_species = occurrence_by_id.get(right, {}).get("species")
        return (species_distances.get((left_species, right_species), math.inf), -score, left, right)

    for left, right, _score in sorted(edges, key=progressive_key):
        left_root, right_root = find(left), find(right)
        if left_root == right_root or not compatible(left_root, right_root):
            continue
        if len(members[left_root]) < len(members[right_root]):
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        members[left_root].update(members.pop(right_root))
    return [sorted(component) for component in members.values()]


def _species_tree_distances(path):
    """Prioritize correspondence by supplied lengths, or topology if any are missing.

    Unit edges apply to this ordering only; CTMC branch lengths are unchanged.
    """
    if not Path(path).exists():
        return {}
    from intraphy.topology import SpeciesTree
    tree = SpeciesTree(read_tsv(path))
    use_topology = any(tree.length[node] is None for node in tree.parent if node != tree.root)
    root_dist = {tree.root: 0.0}
    ancestors = {}
    for node in tree.preorder():
        ancestors[node] = [node] + (ancestors.get(tree.parent.get(node), []))
        for child in tree.children.get(node, []):
            root_dist[child] = root_dist[node] + (1.0 if use_topology else tree.branch_length(child))
    result = {}
    for left_label, left_node in tree.leaf_by_label.items():
        for right_label, right_node in tree.leaf_by_label.items():
            right_ancestors = set(ancestors[right_node])
            lca = next(node for node in ancestors[left_node] if node in right_ancestors)
            result[(left_label, right_label)] = root_dist[left_node] + root_dist[right_node] - 2 * root_dist[lca]
    return result
