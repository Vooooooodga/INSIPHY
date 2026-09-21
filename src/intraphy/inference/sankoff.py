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
    """Legacy two-state encoding routed through the generalized numerical owner."""
    import numpy as np
    from .configuration_dp import sankoff
    tips = {}
    for label in tree.leaf_by_label:
        value = observations.get(label, "unknown")
        tips[label] = np.ones(2) if value == "unknown" else np.eye(2)[int(value)]
    costs = {child: np.array([[0., 1.], [1., 0.]]) for _, child in tree.edges()}
    result = sankoff(tree, tips, costs)
    return result.cost, {n: set(v) for n, v in result.nodes.items()}, {
        branch: set(values) for branch, values in result.pairs.items()}
