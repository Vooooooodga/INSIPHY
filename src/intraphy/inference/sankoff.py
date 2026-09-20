"""inference / sankoff: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

import math


STATES = (0, 1)


PARSIMONY_INFERENCE_METHOD = "equal_cost_maximum_parsimony"


def _fmt(value):
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):.8g}"


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
