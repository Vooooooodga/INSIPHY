"""Validated rooted topology; no model, filesystem discovery or legacy mathematics."""
import math
from collections import Counter, defaultdict

class SpeciesTree:
    def __init__(self, rows):
        self.parent = {}
        self.children = defaultdict(list)
        self.label = {}
        self.length = {}
        seen_nodes = set()
        for row in rows:
            node = str(row.get("node_id", "")).strip()
            if not node:
                raise SystemExit("species_tree.tsv contains an empty node_id")
            if node in seen_nodes:
                raise SystemExit(f"species_tree.tsv contains duplicate node_id: {node}")
            seen_nodes.add(node)
            parent = row.get("parent_id", "")
            parent = str(parent).strip() if parent is not None else ""
            self.parent[node] = parent
            self.label[node] = row.get("label", node) or node
            self.length[node] = self._parse_branch_length(row, node)
            if parent:
                self.children[parent].append(node)
        if not self.parent:
            raise SystemExit("species_tree.tsv contains no nodes")
        missing_parents = sorted(
            {parent for parent in self.parent.values() if parent and parent not in self.parent}
        )
        if missing_parents:
            raise SystemExit("species_tree.tsv contains missing parent nodes: " + ", ".join(missing_parents))
        roots = [node for node, parent in self.parent.items() if not parent]
        if len(roots) != 1:
            raise SystemExit("species_tree.tsv must contain exactly one root")
        self.root = roots[0]
        self._validate_reachable_acyclic()
        self.leaves = [node for node in self.parent if not self.children.get(node)]
        leaf_labels = [self.label[node] for node in self.leaves]
        duplicated = sorted(label for label, count in Counter(leaf_labels).items() if count > 1)
        if duplicated:
            raise SystemExit("species_tree.tsv contains duplicate leaf labels: " + ", ".join(duplicated))
        self.leaf_by_label = {self.label[node]: node for node in self.leaves}

    @staticmethod
    def _raw_branch_length(row):
        for key in ("branch_length", "length", "distance"):
            if key in row:
                value = row.get(key)
                if value is None:
                    return None
                text = str(value).strip()
                if text:
                    return text
        return None

    @classmethod
    def _parse_branch_length(cls, row, node):
        raw = cls._raw_branch_length(row)
        if raw is None or str(raw).strip().upper() == "NA":
            return None
        try:
            length = float(raw)
        except ValueError:
            raise SystemExit(f"species_tree.tsv has a non-numeric branch length for node {node}: {raw}")
        if not math.isfinite(length):
            raise SystemExit(f"species_tree.tsv has a non-finite branch length for node {node}")
        if length < 0:
            raise SystemExit(f"species_tree.tsv has a negative branch length for node {node}")
        return length

    def _validate_reachable_acyclic(self):
        visiting = set()
        visited = set()

        def visit(node):
            if node in visiting:
                raise SystemExit("species_tree.tsv contains a cycle")
            if node in visited:
                return
            visiting.add(node)
            for child in self.children.get(node, []):
                visit(child)
            visiting.remove(node)
            visited.add(node)

        visit(self.root)
        unreachable = sorted(set(self.parent) - visited)
        if unreachable:
            raise SystemExit("species_tree.tsv contains nodes unreachable from the root: " + ", ".join(unreachable))

    def postorder(self):
        order = []

        def visit(node):
            for child in self.children.get(node, []):
                visit(child)
            order.append(node)

        visit(self.root)
        return order

    def preorder(self):
        order = []

        def visit(node):
            order.append(node)
            for child in self.children.get(node, []):
                visit(child)

        visit(self.root)
        return order

    def edges(self):
        for child, parent in self.parent.items():
            if parent:
                yield parent, child

    def branch_length(self, child):
        if child not in self.parent:
            raise KeyError(child)
        length = self.length.get(child)
        if length is None:
            raise SystemExit(f"species_tree.tsv lacks a branch length for non-root node {child}")
        return float(length)
