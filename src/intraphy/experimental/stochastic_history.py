"""experimental / stochastic history: explicit implementation ownership."""
from __future__ import annotations

from collections import Counter
from intraphy.experimental.markov import discrete_likelihood_tables
from intraphy.experimental.markov import fit_invariant_test
from intraphy.experimental.markov import transition_cost
from intraphy.experimental.markov import transition_matrix
from intraphy.storage.values import norm_state
import math
import random


def _quantile(values, q):
    if not values:
        return 0.0
    values = sorted(values)
    pos = (len(values) - 1) * q
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return float(values[lo])
    frac = pos - lo
    return float(values[lo] * (1.0 - frac) + values[hi] * frac)


def _sample_weighted(items, rng):
    total = sum(max(0.0, weight) for _item, weight in items)
    if total <= 0:
        return items[0][0]
    draw = rng.random() * total
    acc = 0.0
    for item, weight in items:
        acc += max(0.0, weight)
        if draw <= acc:
            return item
    return items[-1][0]


def bootstrap_invariant_test(tree, tips_by_label, states, layer, replicates=100, seed=7, tip_error=1e-6):
    replicates = int(replicates or 0)
    if replicates <= 0:
        return None
    states = tuple(sorted(states))
    observed_test = fit_invariant_test(tree, tips_by_label, states, layer, tip_error=tip_error)
    observed_lrt = observed_test["lrt_statistic"]
    observed_states = [
        norm_state(value)
        for value in tips_by_label.values()
        if norm_state(value) in states and norm_state(value) != "unknown"
    ]
    root_candidates = sorted(set(observed_states)) or [state for state in states if state != "unknown"] or list(states)
    missing_labels = {label for label, value in tips_by_label.items() if norm_state(value) == "unknown"}
    rng = random.Random(seed)
    simulated_lrt = []
    for _rep in range(replicates):
        root_state = rng.choice(root_candidates)
        simulated = {}
        for leaf in tree.leaves:
            label = tree.label[leaf]
            if label in missing_labels:
                simulated[label] = "unknown"
                continue
            if tip_error > 0.0 and rng.random() < tip_error:
                alternatives = [state for state in states if state not in {root_state, "unknown"}]
                simulated[label] = rng.choice(alternatives or [root_state])
            else:
                simulated[label] = root_state
        simulated_lrt.append(fit_invariant_test(tree, simulated, states, layer, tip_error=tip_error)["lrt_statistic"])
    exceed = sum(1 for value in simulated_lrt if value >= observed_lrt - 1e-12)
    empirical_p = (exceed + 1.0) / (replicates + 1.0)
    mcse = math.sqrt(empirical_p * (1.0 - empirical_p) / (replicates + 1.0))
    return {
        "test_id": "ctmc_vs_invariant_parametric_bootstrap",
        "observed_lrt": observed_lrt,
        "bootstrap_replicates": replicates,
        "empirical_p_value": empirical_p,
        "monte_carlo_se": mcse,
        "null_lrt_mean": sum(simulated_lrt) / max(1, len(simulated_lrt)),
        "null_lrt_q025": _quantile(simulated_lrt, 0.025),
        "null_lrt_q500": _quantile(simulated_lrt, 0.5),
        "null_lrt_q975": _quantile(simulated_lrt, 0.975),
        "seed": seed,
        "tip_error": tip_error,
    }


def ctmc_posteriors(tree, tips_by_label, states, layer, rate=0.15, tip_error=1e-6):
    states = tuple(sorted(states))
    likelihoods, root_prior, total = discrete_likelihood_tables(tree, tips_by_label, states, layer, rate, tip_error)
    outside = {tree.root: {state: root_prior for state in states}}
    node_rows = []
    branch_rows = []

    for node in tree.preorder():
        denom = total or 1.0
        for state in states:
            posterior = outside[node][state] * likelihoods[node][state] / denom
            node_rows.append({"node_id": node, "state": state, "ctmc_probability": posterior})
        for child in tree.children.get(node, []):
            outside[child] = {}
            matrix = transition_matrix(layer, states, rate, tree.branch_length(child))
            sibling_terms = {}
            for state in states:
                prob = outside[node][state]
                for sibling in tree.children[node]:
                    if sibling == child:
                        continue
                    sibling_matrix = transition_matrix(layer, states, rate, tree.branch_length(sibling))
                    prob *= sum(sibling_matrix[(state, sibling_state)] * likelihoods[sibling][sibling_state] for sibling_state in states)
                sibling_terms[state] = prob
            for child_state in states:
                outside[child][child_state] = sum(sibling_terms[parent_state] * matrix[(parent_state, child_state)] for parent_state in states)

            joint = {}
            for parent_state in states:
                for child_state in states:
                    joint[(parent_state, child_state)] = sibling_terms[parent_state] * matrix[(parent_state, child_state)] * likelihoods[child][child_state] / denom
            change_probability = sum(prob for (parent_state, child_state), prob in joint.items() if parent_state != child_state)
            best_pair, best_prob = max(joint.items(), key=lambda item: item[1])
            branch_rows.append(
                {
                    "parent_node": node,
                    "child_node": child,
                    "ctmc_change_probability": change_probability,
                    "ctmc_most_likely_change": f"{best_pair[0]}->{best_pair[1]}",
                    "ctmc_most_likely_change_probability": best_prob,
                    "endpoint_joint": joint,
                }
            )
    return node_rows, branch_rows


def rate_matrix_array(layer, states, rate=0.15):
    import numpy as np  # type: ignore

    states = tuple(sorted(states))
    q = np.zeros((len(states), len(states)), dtype=float)
    for i, src in enumerate(states):
        weights = []
        for dst in states:
            weights.append(0.0 if dst == src else math.exp(-transition_cost(layer, src, dst)))
        denom = sum(weights) or 1.0
        for j, _dst in enumerate(states):
            if i == j:
                continue
            q[i, j] = rate * weights[j] / denom
        q[i, i] = -sum(q[i, j] for j in range(len(states)) if j != i)
    return q


def sample_ctmc_bridge(states, layer, rate, branch_length, start_state, end_state, rng):
    states = tuple(sorted(states))
    if start_state == end_state and rate <= 0:
        return [start_state]
    try:
        import numpy as np  # type: ignore

        q = rate_matrix_array(layer, states, rate)
        lam = max(float(-q[i, i]) for i in range(len(states)))
        if lam <= 0:
            return [start_state, end_state] if start_state != end_state else [start_state]
        r = np.eye(len(states)) + q / lam
        lamt = max(0.0, lam * branch_length)
        nmax = max(20, int(lamt + 10.0 * math.sqrt(lamt + 1.0) + 20.0))
        nmax = min(nmax, 1000)
        start = states.index(start_state)
        end = states.index(end_state)
        powers = [np.eye(len(states))]
        for _ in range(nmax):
            powers.append(powers[-1] @ r)
        poisson = math.exp(-lamt)
        weights = []
        for n in range(nmax + 1):
            if n > 0:
                poisson *= lamt / n
            weights.append(max(0.0, poisson * float(powers[n][start, end])))
        nsteps = _sample_weighted(list(enumerate(weights)), rng)
        if nsteps == 0:
            return [start_state]
        path = [start_state]
        current = start
        for step in range(1, nsteps):
            remaining = nsteps - step
            choices = []
            for state_idx, state in enumerate(states):
                weight = float(r[current, state_idx]) * float(powers[remaining][state_idx, end])
                choices.append((state_idx, max(0.0, weight)))
            current = _sample_weighted(choices, rng)
            path.append(states[current])
        path.append(end_state)
        return path
    except Exception:
        return [start_state, end_state] if start_state != end_state else [start_state]


def stochastic_map_summary(tree, tips_by_label, states, layer, rate=0.15, replicates=100, seed=7, tip_error=1e-6):
    replicates = int(replicates or 0)
    if replicates <= 0:
        return []
    rng = random.Random(seed)
    _nodes, edge_rows = ctmc_posteriors(tree, tips_by_label, states, layer, rate=rate, tip_error=tip_error)
    out = []
    for edge in edge_rows:
        joint = edge.get("endpoint_joint", {})
        endpoint_items = list(joint.items())
        counts = []
        any_changes = []
        transition_counter = Counter()
        for _rep in range(replicates):
            parent_state, child_state = _sample_weighted(endpoint_items, rng)
            path = sample_ctmc_bridge(states, layer, rate, tree.branch_length(edge["child_node"]), parent_state, child_state, rng)
            changes = [(a, b) for a, b in zip(path, path[1:]) if a != b]
            counts.append(len(changes))
            any_changes.append(1 if changes else 0)
            transition_counter.update(f"{a}->{b}" for a, b in changes)
        top_transition, top_count = ("NA", 0)
        if transition_counter:
            top_transition, top_count = transition_counter.most_common(1)[0]
        out.append(
            {
                "parent_node": edge["parent_node"],
                "child_node": edge["child_node"],
                "map_sample_count": replicates,
                "posterior_pr_any_change": sum(any_changes) / max(1, replicates),
                "posterior_expected_change_count": sum(counts) / max(1, replicates),
                "posterior_change_count_low": _quantile(counts, 0.025),
                "posterior_change_count_high": _quantile(counts, 0.975),
                "posterior_most_frequent_transition": top_transition,
                "posterior_transition_probability": top_count / max(1, sum(counts)),
            }
        )
    return out
