"""Canonical branches for configuration inference, with original-edge provenance.

An unobserved degree-two node must not introduce another source-origin opportunity.
Roots remain fixed. Explicit foreground changes remain meaningful branch boundaries.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from ..topology import SpeciesTree


@dataclass(frozen=True)
class TreeContext:
    tree: SpeciesTree
    foreground: frozenset[str]
    collapsed: tuple[str, ...]
    edge_sources: dict[str, tuple[str, ...]]

    def diagnostics(self):
        return {"collapsed_unobserved_unary_nodes": list(self.collapsed),
                "normalized_edge_to_input_children": self.edge_sources,
                "root_preserved": True,
                "foreground_children": sorted(self.foreground),
                "origin_process": "finite_branch_opportunities_not_continuous_immigration"}


def tree_rows(tree):
    return [{"node_id": n, "parent_id": tree.parent[n], "label": tree.label[n],
             "branch_length": tree.length[n] if tree.length[n] is not None else "NA"}
            for n in tree.preorder()]


def canonical_tree(tree, foreground=frozenset()):
    foreground = frozenset(foreground)
    if not foreground <= set(tree.parent)-{tree.root}:
        raise ValueError("Foreground contains unknown/root nodes")
    parent, length = dict(tree.parent), dict(tree.length)
    children = {n: list(tree.children.get(n, ())) for n in parent}
    paths = {n: (n,) for n in parent if n != tree.root}
    flags = {n: n in foreground for n in paths}
    collapsed = []
    # Postorder makes concatenation associative for arbitrarily subdivided edges.
    for node in tree.postorder():
        if node == tree.root or node not in parent or len(children[node]) != 1:
            continue
        child = children[node][0]
        if flags[node] != flags[child]:
            continue  # A specified rate change is not a redundant representation.
        ancestor = parent[node]
        parent[child] = ancestor
        children[ancestor] = [child if c == node else c for c in children[ancestor]]
        path = paths[node]+paths[child]
        paths[child] = path
        values = [tree.length[n] for n in path]
        length[child] = math.fsum(values) if all(x is not None for x in values) else None
        collapsed.append(node)
        del parent[node], children[node], length[node], paths[node], flags[node]
    rows = [{"node_id": n, "parent_id": parent[n], "label": tree.label[n],
             "branch_length": length[n]} for n in tree.preorder() if n in parent]
    result = SpeciesTree(rows)
    return TreeContext(result, frozenset(n for n, flag in flags.items() if flag),
                       tuple(sorted(collapsed)), paths)
