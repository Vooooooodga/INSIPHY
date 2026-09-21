"""Directed edit distances, optimal edit-path support and compatible witnesses."""
from __future__ import annotations
from dataclasses import dataclass
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
import math
import numpy as np
from .space import StateSpace
from .types import ElementaryEdit
from .origins import permitted


@dataclass
class EditGraph:
    size: int
    edges: tuple[tuple[int, int, ElementaryEdit, float], ...]
    distance: np.ndarray
    outgoing: dict[int, list[tuple[int, ElementaryEdit, float]]]

    def witness_path(self, source: int, target: int) -> tuple[ElementaryEdit, ...]:
        if not math.isfinite(self.distance[source, target]):
            raise ValueError("Unreachable configuration pair")
        path = []
        cursor = source
        while cursor != target:
            choices = [(j, e, w) for j, e, w in self.outgoing.get(cursor, [])
                       if np.isclose(w+self.distance[j, target], self.distance[cursor, target], rtol=0, atol=1e-9)]
            if not choices:
                raise ArithmeticError("Shortest-path witness reconstruction failed")
            j, e, _ = min(choices, key=lambda x: (x[1].event_key, x[0]))
            path.append(e)
            cursor = j
        return tuple(path)

    def path_event_bounds(self, source: int, target: int) -> dict[str, tuple[int, int]]:
        """Minimum/maximum occurrences of each edit identity on all shortest paths."""
        best = self.distance[source, target]
        if not math.isfinite(best):
            return {}
        selected = [(i, j, e, w) for i, j, e, w in self.edges
                    if np.isclose(self.distance[source, i]+w+self.distance[j, target], best, rtol=0, atol=1e-9)]
        keys = sorted({e.event_key for _, _, e, _ in selected})
        if not keys:
            return {}
        key_index = {key: i for i, key in enumerate(keys)}
        lo, hi = {source: np.zeros(len(keys), dtype=int)}, {source: np.zeros(len(keys), dtype=int)}
        outgoing = {}
        for i, j, e, w in selected:
            outgoing.setdefault(i, []).append((j, e))
        order = sorted(set([source, target] + [i for i, _, _, _ in selected]), key=lambda i: self.distance[source, i])
        for i in order:
            if i not in lo:
                continue
            for j, edit in outgoing.get(i, []):
                increment = np.zeros(len(keys), dtype=int)
                increment[key_index[edit.event_key]] = 1
                lower, upper = lo[i]+increment, hi[i]+increment
                lo[j] = np.minimum(lo[j], lower) if j in lo else lower
                hi[j] = np.maximum(hi[j], upper) if j in hi else upper
        return {key: (int(lo[target][k]), int(hi[target][k])) for key, k in key_index.items()}


def make_graph(space: StateSpace, origins: dict[str, str], child: str,
               costs: dict[str, float] | None = None) -> EditGraph:
    if not space.complete:
        raise ValueError("state_space_incomplete: no formal edit graph")
    index = space.index
    costs = costs or {}
    if any(not np.isfinite(x) or x <= 0 for x in costs.values()):
        raise ValueError("Elementary edit costs must be positive and finite")
    edges = tuple((index[e.source], index[e.target], e, float(costs.get(e.kind, 1.)))
                  for e in space.edits if permitted(e, origins, child))
    outgoing = {}
    for i, j, e, w in edges:
        outgoing.setdefault(i, []).append((j, e, w))
    n = len(index)
    # Use sparse numerical shortest paths, but retain every parallel edit above
    # for all-optimal history identities. CSR duplicate entries must NOT be summed.
    minimum = {}
    for i, j, _, weight in edges:
        minimum[i, j] = min(weight, minimum.get((i, j), float("inf")))
    rows, cols, values = [], [], []
    for (i, j), value in sorted(minimum.items()):
        rows.append(i); cols.append(j); values.append(value)
    adjacency = csr_matrix((values, (rows, cols)), shape=(n, n))
    distances = dijkstra(adjacency, directed=True, return_predecessors=False)
    return EditGraph(n, edges, distances, outgoing)
