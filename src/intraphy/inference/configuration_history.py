"""All optimal origin/structure histories and non-additive elementary edits."""
from __future__ import annotations
from collections import defaultdict
import math
import numpy as np

from ..structure.space import StateSpace
from ..structure.origins import origin_scenarios
from ..structure.paths import make_graph
from .configuration_dp import sankoff
from ..structure.tree_context import canonical_tree


def reconstruct(space: StateSpace, tree, tips, *, max_origins=256, costs=None, normalize_tree=True):
    if normalize_tree:
        tree = canonical_tree(tree).tree
    scenarios, optimum = [], math.inf
    cache = {}
    for origins, root, _ in origin_scenarios(space, tree, max_origins, tips=tips):
        # All edges without an introduction share a graph.
        graphs = {}
        for _, child in tree.edges():
            signature = tuple(sorted(k for k, v in origins.items() if v == child))
            if signature not in cache:
                cache[signature] = make_graph(space, origins, child, costs)
            graphs[child] = cache[signature]
        result = sankoff(tree, tips, {child: graph.distance for child, graph in graphs.items()}, root)
        if result.cost < optimum-1e-9:
            optimum, scenarios = result.cost, []
        if math.isfinite(result.cost) and abs(result.cost-optimum) <= 1e-9:
            scenarios.append((origins, result, graphs))
    if not scenarios:
        return {"minimum_cost": None, "status": "no_compatible_history", "events": [], "witness": [], "nodes": {}}
    global_nodes = defaultdict(set)
    global_events = {}
    branch_pairs = defaultdict(set)
    for parent, child in tree.edges():
        alternatives = []
        details = {}
        path_cache = {}
        for origins, result, graphs in scenarios:
            for node, values in result.nodes.items():
                global_nodes[node].update(values)
            graph = graphs[child]
            for source, target in result.pairs[(parent, child)]:
                branch_pairs[(parent, child)].add((source, target))
                key = (id(graph), source, target)
                if key not in path_cache:
                    path_cache[key] = graph.path_event_bounds(source, target)
                alternatives.append(path_cache[key])
                best = graph.distance[source, target]
                for i, j, edit, weight in graph.edges:
                    if np.isclose(graph.distance[source, i]+weight+graph.distance[j, target], best, rtol=0, atol=1e-9):
                        details.setdefault(edit.event_key, []).append(edit)
        event_keys = sorted(set().union(*(a.keys() for a in alternatives)))
        for key in event_keys:
            lower = min(a.get(key, (0, 0))[0] for a in alternatives)
            upper = max(a.get(key, (0, 0))[1] for a in alternatives)
            edits = details[key]
            one = edits[0]
            global_events[(parent, child, key)] = {
                "parent": parent, "child": child, "edit_key": key, "operation": one.kind,
                "start": one.footprint[0], "end": one.footprint[1],
                "support": "required" if lower else "possible", "minimum_count": lower,
                "maximum_count": upper, "affected_spans": sorted({(e.start, e.end) for x in edits for e in x.affected}),
                "consequences": sorted({c for e in edits for c in e.consequences}),
                "molecular_mutation_count": "not_identified"}
    origins, representative, graphs = scenarios[0]
    witness = []
    for parent, child in tree.edges():
        for order, edit in enumerate(graphs[child].witness_path(representative.witness[parent], representative.witness[child]), 1):
            witness.append({"parent": parent, "child": child, "order": order, "edit_key": edit.event_key,
                "operation": edit.kind, "source": space.index[edit.source], "target": space.index[edit.target],
                "consequences": list(edit.consequences), "molecular_mutation_count": "not_identified"})
    unit_costs = not costs or all(float(v) == 1 for v in costs.values())
    if unit_costs and len(witness) != round(optimum):
        raise ArithmeticError("Representative edit history does not attain the global minimum")
    return {"status": "conditional_on_catalogue_and_observation", "minimum_cost": optimum,
            "minimum_structural_edits": int(round(optimum)) if unit_costs else None,
            "representative_edit_count": len(witness), "optimal_origin_scenarios": len(scenarios),
            "representative_origins": origins, "representative_states": representative.witness,
            "events": list(global_events.values()), "witness": witness,
            "nodes": {node: sorted(values) for node, values in global_nodes.items()},
            "pairs": [{"parent": p, "child": c, "pairs": sorted(values)} for (p, c), values in branch_pairs.items()]}
