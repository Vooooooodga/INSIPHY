"""Small rooted-tree and Sankoff routines for structural characters."""

import math
from collections import defaultdict

from .io import norm_state


class SpeciesTree:
    def __init__(self, rows):
        self.parent = {}
        self.children = defaultdict(list)
        self.label = {}
        for row in rows:
            node = row["node_id"]
            parent = row.get("parent_id", "")
            self.parent[node] = parent
            self.label[node] = row.get("label", node) or node
            if parent:
                self.children[parent].append(node)
        roots = [node for node, parent in self.parent.items() if not parent]
        if len(roots) != 1:
            raise SystemExit("species_tree.tsv must contain exactly one root")
        self.root = roots[0]
        self.leaves = [node for node in self.parent if not self.children.get(node)]
        self.leaf_by_label = {self.label[node]: node for node in self.leaves}

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


def transition_cost(layer, src, dst):
    if src == dst:
        return 0.0
    if layer in {"segment_presence", "adjacency_state"}:
        if src == "absent" and dst == "present":
            return 2.0
        if src == "present" and dst == "absent":
            return 1.5
        return 2.5
    if layer == "role_state":
        if "absent" in (src, dst):
            return 2.0
        return 1.0
    if layer == "source_mixture":
        return 2.0 if "multi_source" in (src, dst) else 1.0
    if layer == "copy_multiplicity":
        return 1.5
    return 1.0


def sankoff(tree, tips_by_label, states, layer):
    states = tuple(sorted(states))
    allowed = {}
    for label, value in tips_by_label.items():
        node = tree.leaf_by_label.get(label)
        if node is None:
            continue
        value = norm_state(value)
        allowed[node] = {value} if value in states else set(states)

    scores = {}
    for node in tree.postorder():
        scores[node] = {}
        if not tree.children.get(node):
            node_allowed = allowed.get(node, set(states))
            for state in states:
                scores[node][state] = 0.0 if state in node_allowed else math.inf
            continue
        for state in states:
            scores[node][state] = sum(
                min(scores[child][child_state] + transition_cost(layer, state, child_state) for child_state in states)
                for child in tree.children[node]
            )

    node_rows = []
    best_states = {}
    for node in tree.preorder():
        finite = {state: score for state, score in scores[node].items() if math.isfinite(score)}
        if not finite:
            finite = {state: 0.0 for state in states}
        best = min(finite.values())
        weights = {state: math.exp(-(score - best)) for state, score in finite.items()}
        total = sum(weights.values()) or 1.0
        best_set = {state for state, score in finite.items() if abs(score - best) < 1e-9}
        best_states[node] = best_set
        for state in states:
            node_rows.append(
                {
                    "node_id": node,
                    "node_label": tree.label[node],
                    "state": state,
                    "probability": f"{weights.get(state, 0.0) / total:.6g}",
                    "is_parsimony_best": int(state in best_set),
                }
            )

    edge_rows = []
    for parent, child in tree.edges():
        parent_states = best_states[parent]
        child_states = best_states[child]
        if parent_states & child_states:
            status = "unchanged_or_ambiguous"
            probability = 0.0
            change = "NA"
        elif len(parent_states) == 1 and len(child_states) == 1:
            status = "change_required"
            probability = 1.0
            change = f"{next(iter(parent_states))}->{next(iter(child_states))}"
        else:
            status = "change_possible"
            probability = 0.5
            change = f"{'|'.join(sorted(parent_states))}->{'|'.join(sorted(child_states))}"
        edge_rows.append(
            {
                "parent_node": parent,
                "child_node": child,
                "parent_label": tree.label[parent],
                "child_label": tree.label[child],
                "status": status,
                "change": change,
                "event_probability": f"{probability:.6g}",
            }
        )
    return min(scores[tree.root].values()), node_rows, edge_rows
