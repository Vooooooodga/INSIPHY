"""Small rooted-tree and Sankoff routines for structural characters."""

import math
import random
from collections import Counter, defaultdict
from functools import lru_cache

from .io import norm_state


class SpeciesTree:
    def __init__(self, rows):
        self.parent = {}
        self.children = defaultdict(list)
        self.label = {}
        self.length = {}
        for row in rows:
            node = row["node_id"]
            parent = row.get("parent_id", "")
            self.parent[node] = parent
            self.label[node] = row.get("label", node) or node
            self.length[node] = float(row.get("branch_length") or row.get("length") or row.get("distance") or 1.0)
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

    def branch_length(self, child):
        return max(1e-9, float(self.length.get(child, 1.0)))


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
                "branch_length": f"{tree.branch_length(child):.6g}",
                "status": status,
                "change": change,
                "event_probability": f"{probability:.6g}",
            }
        )
    return min(scores[tree.root].values()), node_rows, edge_rows
