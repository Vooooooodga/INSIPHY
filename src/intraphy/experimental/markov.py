"""experimental / markov: explicit implementation ownership."""
from __future__ import annotations

from functools import lru_cache
from intraphy.storage.values import norm_state
import math


def transition_cost(layer, src, dst):
    if src == dst:
        return 0.0
    if layer in {"segment_presence", "adjacency_state", "element_presence", "element_adjacency_state"}:
        if src == "absent" and dst == "present":
            return 2.0
        if src == "present" and dst == "absent":
            return 1.5
        return 2.5
    if layer in {"role_state", "element_role_state"}:
        if "absent" in (src, dst):
            return 2.0
        return 1.0
    if layer == "source_mixture":
        return 2.0 if "multi_source" in (src, dst) else 1.0
    if layer == "copy_multiplicity":
        return 1.5
    return 1.0


def _transition_probability_approx(layer, src, dst, states, rate=0.15, branch_length=1.0):
    effective_rate = min(0.95, max(0.0, rate * branch_length))
    if src == dst:
        return max(1e-12, 1.0 - effective_rate)
    weights = []
    for state in states:
        if state != src:
            weights.append(math.exp(-transition_cost(layer, src, state)))
    denom = sum(weights) or 1.0
    return max(1e-12, effective_rate * math.exp(-transition_cost(layer, src, dst)) / denom)


@lru_cache(maxsize=4096)
def _transition_matrix_cached(layer, states, rate=0.15, branch_length=1.0):
    states = tuple(sorted(states))
    try:
        import numpy as np  # type: ignore
        from scipy.linalg import expm  # type: ignore

        q = np.zeros((len(states), len(states)))
        for i, src in enumerate(states):
            weights = []
            for dst in states:
                weights.append(0.0 if dst == src else math.exp(-transition_cost(layer, src, dst)))
            denom = sum(weights) or 1.0
            for j, dst in enumerate(states):
                if src == dst:
                    continue
                q[i, j] = rate * weights[j] / denom
            q[i, i] = -sum(q[i, j] for j in range(len(states)) if j != i)
        p = expm(q * branch_length)
        return {(src, dst): max(1e-12, float(p[i, j])) for i, src in enumerate(states) for j, dst in enumerate(states)}
    except Exception:
        return {(src, dst): _transition_probability_approx(layer, src, dst, states, rate, branch_length) for src in states for dst in states}


def transition_matrix(layer, states, rate=0.15, branch_length=1.0):
    return _transition_matrix_cached(layer, tuple(sorted(states)), rate, branch_length)


def transition_probability(layer, src, dst, states, rate=0.15, branch_length=1.0):
    return transition_matrix(layer, states, rate, branch_length)[(src, dst)]


def emission_probability(observed, state, states, tip_error=0.0):
    observed = norm_state(observed)
    if observed == "unknown" or observed not in states:
        return 1.0
    if observed == state:
        return max(1e-12, 1.0 - tip_error)
    if tip_error <= 0.0:
        return 0.0
    return max(1e-12, tip_error / max(1, len(states) - 1))


def _tip_observations(tree, tips_by_label):
    observations = {}
    for label, value in tips_by_label.items():
        node = tree.leaf_by_label.get(label)
        if node is not None:
            observations[node] = norm_state(value)
    return observations


def discrete_likelihood_tables(tree, tips_by_label, states, layer, rate=0.15, tip_error=0.0, rate_by_child=None):
    states = tuple(sorted(states))
    observations = _tip_observations(tree, tips_by_label)

    likelihoods = {}
    for node in tree.postorder():
        likelihoods[node] = {}
        if not tree.children.get(node):
            observed = observations.get(node, "unknown")
            for state in states:
                likelihoods[node][state] = emission_probability(observed, state, states, tip_error)
            continue
        for state in states:
            prob = 1.0
            for child in tree.children[node]:
                branch_length = tree.branch_length(child)
                child_sum = 0.0
                branch_rate = rate_by_child.get(child, rate) if rate_by_child else rate
                matrix = transition_matrix(layer, states, branch_rate, branch_length)
                for child_state in states:
                    child_sum += matrix[(state, child_state)] * likelihoods[child][child_state]
                prob *= max(child_sum, 1e-300)
            likelihoods[node][state] = prob
    root_prior = 1.0 / max(1, len(states))
    total = sum(root_prior * likelihoods[tree.root][state] for state in states)
    return likelihoods, root_prior, max(total, 1e-300)


def discrete_log_likelihood(tree, tips_by_label, states, layer, rate=0.15, tip_error=0.0, rate_by_child=None):
    _likelihoods, _root_prior, total = discrete_likelihood_tables(tree, tips_by_label, states, layer, rate, tip_error, rate_by_child=rate_by_child)
    return math.log(max(total, 1e-300))


def fit_discrete_ctmc(tree, tips_by_label, states, layer, tip_error=0.0):
    rates = [0.005, 0.01, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.4, 0.65, 1.0, 1.5]
    scored = [(discrete_log_likelihood(tree, tips_by_label, states, layer, rate, tip_error), rate) for rate in rates]
    log_likelihood, rate = max(scored, key=lambda item: item[0])
    observed = sum(1 for value in tips_by_label.values() if norm_state(value) != "unknown")
    k = 1
    aic = 2 * k - 2 * log_likelihood
    bic = math.log(max(1, observed)) * k - 2 * log_likelihood
    return {
        "model": "ctmc_mk_branch_length",
        "rate": rate,
        "log_likelihood": log_likelihood,
        "aic": aic,
        "bic": bic,
        "observed_tip_count": observed,
    }


def fit_foreground_rate_test(tree, tips_by_label, states, layer, foreground_edges, tip_error=0.0):
    foreground_edges = set(foreground_edges or [])
    if not foreground_edges:
        return None
    rates = [0.005, 0.01, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.4, 0.65, 1.0, 1.5]
    one = fit_discrete_ctmc(tree, tips_by_label, states, layer, tip_error=tip_error)
    best = None
    for background_rate in rates:
        for foreground_rate in rates:
            rate_by_child = {
                child: foreground_rate if (parent, child) in foreground_edges else background_rate
                for parent, child in tree.edges()
            }
            log_likelihood = discrete_log_likelihood(tree, tips_by_label, states, layer, rate=background_rate, tip_error=tip_error, rate_by_child=rate_by_child)
            candidate = (log_likelihood, background_rate, foreground_rate)
            if best is None or candidate[0] > best[0]:
                best = candidate
    if best is None:
        return None
    log_likelihood, background_rate, foreground_rate = best
    lrt = max(0.0, 2.0 * (log_likelihood - one["log_likelihood"]))
    observed = sum(1 for value in tips_by_label.values() if norm_state(value) != "unknown")
    return {
        "null_model": "one_rate_ctmc_mk",
        "alternative_model": "foreground_background_ctmc_mk",
        "null_log_likelihood": one["log_likelihood"],
        "alternative_log_likelihood": log_likelihood,
        "lrt_statistic": lrt,
        "df": 1,
        "p_value": chi_square_sf_df1(lrt),
        "p_value_method": "chi_square_df1_foreground_background_rate_test",
        "background_rate": background_rate,
        "foreground_rate": foreground_rate,
        "rate_ratio": foreground_rate / max(background_rate, 1e-12),
        "null_aic": one["aic"],
        "alternative_aic": 2 * 2 - 2 * log_likelihood,
        "null_bic": one["bic"],
        "alternative_bic": math.log(max(1, observed)) * 2 - 2 * log_likelihood,
        "observed_tip_count": observed,
    }


def chi_square_sf_df1(statistic):
    statistic = max(0.0, statistic)
    try:
        from scipy.stats import chi2  # type: ignore

        return float(chi2.sf(statistic, 1))
    except Exception:
        return math.erfc(math.sqrt(statistic / 2.0))


def fit_invariant_test(tree, tips_by_label, states, layer, tip_error=1e-6):
    observed = sum(1 for value in tips_by_label.values() if norm_state(value) != "unknown")
    alt = fit_discrete_ctmc(tree, tips_by_label, states, layer, tip_error=tip_error)
    null_log_likelihood = discrete_log_likelihood(tree, tips_by_label, states, layer, rate=0.0, tip_error=tip_error)
    lrt = max(0.0, 2.0 * (alt["log_likelihood"] - null_log_likelihood))
    boundary_mixture_p = 0.5 * chi_square_sf_df1(lrt)
    null_k = 0
    alt_k = 1
    return {
        "null_model": "invariant_no_structural_change",
        "alternative_model": alt["model"],
        "null_log_likelihood": null_log_likelihood,
        "alternative_log_likelihood": alt["log_likelihood"],
        "lrt_statistic": lrt,
        "df": 1,
        "p_value": boundary_mixture_p,
        "p_value_method": "0.5*chi_square_df1_boundary_rate_test",
        "null_aic": 2 * null_k - 2 * null_log_likelihood,
        "alternative_aic": 2 * alt_k - 2 * alt["log_likelihood"],
        "null_bic": math.log(max(1, observed)) * null_k - 2 * null_log_likelihood,
        "alternative_bic": math.log(max(1, observed)) * alt_k - 2 * alt["log_likelihood"],
        "fitted_rate": alt["rate"],
        "observed_tip_count": observed,
        "tip_error": tip_error,
    }


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
                min(scores[child][child_state] + tree.branch_length(child) * transition_cost(layer, state, child_state) for child_state in states)
                for child in tree.children[node]
            )

    def choose_state(candidates, preferred=None):
        if preferred in candidates:
            return preferred
        return sorted(candidates)[0]

    chosen_states = {}
    root_finite = {state: score for state, score in scores[tree.root].items() if math.isfinite(score)}
    if not root_finite:
        root_finite = {state: 0.0 for state in states}
    root_best = min(root_finite.values())
    chosen_states[tree.root] = choose_state({state for state, score in root_finite.items() if abs(score - root_best) < 1e-9})
    for parent in tree.preorder():
        parent_state = chosen_states[parent]
        for child in tree.children.get(parent, []):
            child_scores = {
                child_state: scores[child][child_state] + tree.branch_length(child) * transition_cost(layer, parent_state, child_state)
                for child_state in states
                if math.isfinite(scores[child][child_state])
            }
            if not child_scores:
                child_scores = {state: 0.0 for state in states}
            best = min(child_scores.values())
            candidates = {state for state, score in child_scores.items() if abs(score - best) < 1e-9}
            chosen_states[child] = choose_state(candidates, preferred=parent_state)

    node_rows = []
    for node in tree.preorder():
        finite = {state: score for state, score in scores[node].items() if math.isfinite(score)}
        if not finite:
            finite = {state: 0.0 for state in states}
        best = min(finite.values())
        weights = {state: math.exp(-(score - best)) for state, score in finite.items()}
        total = sum(weights.values()) or 1.0
        for state in states:
            node_rows.append(
                {
                    "node_id": node,
                    "node_label": tree.label[node],
                    "state": state,
                    "probability": f"{weights.get(state, 0.0) / total:.6g}",
                    "is_parsimony_best": int(state == chosen_states[node]),
                }
            )

    edge_rows = []
    for parent, child in tree.edges():
        parent_state = chosen_states[parent]
        child_state = chosen_states[child]
        if parent_state == child_state:
            status = "unchanged_or_ambiguous"
            probability = 0.0
            change = "NA"
        else:
            status = "change_required"
            probability = 1.0
            change = f"{parent_state}->{child_state}"
        edge_rows.append(
            {
                "parent_node": parent,
                "child_node": child,
                "parent_label": tree.label[parent],
                "child_label": tree.label[child],
                "branch_length": f"{tree.branch_length(child):.6g}",
                "status": status,
                "change": change,
                "event_probability": f"{probability:.6g}",
            }
        )
    return min(scores[tree.root].values()), node_rows, edge_rows
