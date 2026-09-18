"""Likelihood analysis of single-copy intragenic structural sites."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq, minimize
from scipy.special import logsumexp
from scipy.stats import chi2

from .io import read_tsv, write_tsv
from .structural_sites import build_structural_site_matrix, _single_copy_families
from .tree import SpeciesTree


RATE_MIN = 1e-8
RATE_MAX = 100.0
MULTIPLIER_MIN = 1e-3
MULTIPLIER_MAX = 1e3
PROFILE_DROP_95 = 0.5 * float(chi2.ppf(0.95, 1))


def _fmt(value):
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):.8g}"


def _validated_tree_rows(path, branch_length_mode):
    if branch_length_mode not in {"supplied", "unit"}:
        raise SystemExit("branch_length_mode must be supplied or unit")
    rows = read_tsv(path, ["node_id", "parent_id", "label"])
    out = []
    for row in rows:
        item = dict(row)
        if not row.get("parent_id"):
            raw_root = row.get("branch_length") or row.get("length") or row.get("distance")
            item["branch_length"] = raw_root if raw_root not in {None, "", "NA"} else "NA"
        elif branch_length_mode == "unit":
            item["branch_length"] = 1.0
        else:
            raw = row.get("branch_length") or row.get("length") or row.get("distance")
            if raw in {None, "", "NA"}:
                raise SystemExit(
                    "species_tree.tsv has a branch without length; use --branch-length-mode unit to analyze branch counts"
                )
            length = float(raw)
            if not math.isfinite(length) or length < 0:
                raise SystemExit("all non-root branches must have finite non-negative lengths")
            item["branch_length"] = length
        out.append(item)
    return out


def _fit_valid_for_posterior(fit):
    return (
        bool(fit.get("converged"))
        and fit.get("fit_status") == "success"
        and fit.get("inference_status") != "no_observed_contrast"
        and bool(fit.get("identifiable"))
        and not bool(fit.get("boundary"))
    )


@lru_cache(maxsize=32768)
def _transition_matrix(gain, loss, branch_length, multiplier=1.0):
    gain_rate = float(gain) * float(multiplier)
    loss_rate = float(loss) * float(multiplier)
    total = gain_rate + loss_rate
    if total <= 0.0:
        return np.eye(2, dtype=float)
    changed = -math.expm1(-total * float(branch_length))
    stayed = 1.0 - changed
    return np.array(
        [
            [(loss_rate + gain_rate * stayed) / total, gain_rate * changed / total],
            [loss_rate * changed / total, (gain_rate + loss_rate * stayed) / total],
        ],
        dtype=float,
    )


def _root_prior(gain, loss, root_presence=None):
    if root_presence is not None:
        return np.array([1.0 - root_presence, root_presence], dtype=float)
    total = gain + loss
    if total <= 0:
        return np.array([0.5, 0.5], dtype=float)
    return np.array([loss / total, gain / total], dtype=float)


def _iter_weighted_patterns(patterns):
    for item in patterns:
        if isinstance(item, tuple) and len(item) == 2:
            observations, weight = item
            yield observations, int(weight)
        else:
            yield item, 1


def _compress_patterns(patterns, labels):
    counts = defaultdict(int)
    by_key = {}
    for observations in patterns:
        completed = {label: observations.get(label, "unknown") for label in labels}
        key = tuple((label, 2 if completed[label] == "unknown" else completed[label]) for label in labels)
        counts[key] += 1
        by_key[key] = completed
    return [(by_key[key], counts[key]) for key in sorted(by_key)]


def _site_key(row):
    return (row.get("family_id", "NA"), row.get("layer", "NA"), row.get("site_id", "NA"))


def _load_structural_site_universe(input_dir):
    path = Path(input_dir) / "structural_site_universe.tsv"
    rows = read_tsv(path, ["family_id", "layer", "site_id"], optional=True)
    if not rows:
        return None
    return {_site_key(row) for row in rows}


def _pattern_observed_values(pattern):
    return {value for value in pattern.values() if value in {0, 1}}


def _pattern_has_state_one(pattern):
    return any(value == 1 for value in pattern.values())


def _selected_for_ascertainment(pattern, ascertainment):
    observed = _pattern_observed_values(pattern)
    if len(observed) == 0:
        return False
    if ascertainment == "observed-at-least-one":
        return _pattern_has_state_one(pattern)
    if ascertainment == "variable-only":
        return len(observed) > 1
    if ascertainment == "complete-universe":
        return True
    raise SystemExit("unsupported ascertainment mode")


def _decode_parameters(theta, model, root_frequency="estimated", root_presence=0.5):
    values = np.exp(np.asarray(theta, dtype=float))
    if model == "ER":
        gain, loss, multiplier, offset = float(values[0]), float(values[0]), 1.0, 1
    elif model == "ARD":
        gain, loss, multiplier, offset = float(values[0]), float(values[1]), 1.0, 2
    elif model == "ARD_FOREGROUND":
        gain, loss, multiplier, offset = float(values[0]), float(values[1]), float(values[2]), 3
    else:
        raise ValueError(f"unknown model: {model}")
    if root_frequency == "stationary":
        rho = gain / max(gain + loss, 1e-300)
    elif root_frequency == "fixed":
        rho = float(root_presence)
    else:
        logit = float(theta[offset])
        rho = 1.0 / (1.0 + math.exp(-logit))
    return gain, loss, multiplier, rho


def _log_array(values):
    with np.errstate(divide="ignore"):
        return np.log(np.asarray(values, dtype=float))


def _inside_log_messages(tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence=None):
    inside = {}
    for node in tree.postorder():
        if not tree.children.get(node):
            observed = observations.get(tree.label[node], "unknown")
            if observed == 0:
                inside[node] = np.array([0.0, -math.inf], dtype=float)
            elif observed == 1:
                inside[node] = np.array([-math.inf, 0.0], dtype=float)
            else:
                inside[node] = np.array([0.0, 0.0], dtype=float)
            continue
        vector = np.zeros(2, dtype=float)
        for child in tree.children[node]:
            multiplier = foreground_multiplier if child in foreground_children else 1.0
            log_matrix = _log_array(_transition_matrix(gain, loss, tree.branch_length(child), multiplier))
            vector += np.logaddexp(
                log_matrix[:, 0] + inside[child][0],
                log_matrix[:, 1] + inside[child][1],
            )
        inside[node] = vector
    prior = _root_prior(gain, loss, root_presence)
    root_log = _log_array(prior) + inside[tree.root]
    return inside, float(np.logaddexp(root_log[0], root_log[1]))


def _inside_messages(tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence=None):
    log_inside, log_likelihood = _inside_log_messages(
        tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence
    )
    inside = {}
    for node, vector in log_inside.items():
        scale = float(np.logaddexp(vector[0], vector[1]))
        inside[node] = np.exp(vector - scale) if math.isfinite(scale) else np.zeros(2, dtype=float)
    return inside, log_likelihood


def _pattern_log_likelihood(tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence=None):
    _inside, log_likelihood = _inside_log_messages(
        tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence
    )
    return log_likelihood


def _ascertainment_log_probability(tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence, mode):
    observed_labels = {label for label, value in observations.items() if value in {0, 1}}
    if not observed_labels or (mode == "variable-only" and len(observed_labels) < 2):
        return -math.inf
    zero = {label: (0 if label in observed_labels else "unknown") for label in tree.leaf_by_label}
    log_excluded = _pattern_log_likelihood(
        tree, zero, gain, loss, foreground_multiplier, foreground_children, root_presence
    )
    if mode == "variable-only":
        one = {label: (1 if label in observed_labels else "unknown") for label in tree.leaf_by_label}
        log_one = _pattern_log_likelihood(
            tree, one, gain, loss, foreground_multiplier, foreground_children, root_presence
        )
        log_excluded = float(np.logaddexp(log_excluded, log_one))
    selected = -math.expm1(log_excluded)
    if selected <= 0.0:
        return -math.inf
    return math.log(selected)


def _dataset_log_likelihood(tree, patterns, model, theta, foreground_children, ascertainment, root_frequency, root_presence):
    gain, loss, multiplier, rho = _decode_parameters(theta, model, root_frequency, root_presence)
    total = 0.0
    for observations, weight in _iter_weighted_patterns(patterns):
        value = _pattern_log_likelihood(tree, observations, gain, loss, multiplier, foreground_children, rho)
        if ascertainment in {"observed-at-least-one", "variable-only"}:
            ascertainment_log_probability = _ascertainment_log_probability(
                tree, observations, gain, loss, multiplier, foreground_children, rho, ascertainment
            )
            if not math.isfinite(ascertainment_log_probability):
                return -math.inf
            value -= ascertainment_log_probability
        if not math.isfinite(value):
            return -math.inf
        total += int(weight) * value
    return total


def _model_bounds(model, root_frequency="estimated"):
    rate_bounds = [(math.log(RATE_MIN), math.log(RATE_MAX))]
    if model in {"ARD", "ARD_FOREGROUND"}:
        rate_bounds.append((math.log(RATE_MIN), math.log(RATE_MAX)))
    if model == "ARD_FOREGROUND":
        rate_bounds.append((math.log(MULTIPLIER_MIN), math.log(MULTIPLIER_MAX)))
    if root_frequency == "estimated":
        rate_bounds.append((-13.8155095579, 13.8155095579))
    return rate_bounds


def _model_starts(model, root_frequency="estimated", root_presence=0.5):
    if model == "ER":
        starts = [[math.log(value)] for value in (0.01, 0.1, 1.0)]
    elif model == "ARD":
        starts = [[math.log(a), math.log(b)] for a, b in ((0.01, 0.1), (0.1, 0.1), (0.1, 1.0), (1.0, 0.1))]
    else:
        starts = [
        [math.log(a), math.log(b), math.log(m)]
        for a, b, m in ((0.01, 0.1, 0.5), (0.1, 0.1, 1.0), (0.1, 1.0, 2.0), (1.0, 0.1, 5.0))
        ]
    if root_frequency == "estimated":
        rho = min(1.0 - 1e-6, max(1e-6, float(root_presence)))
        root_logit = math.log(rho / (1.0 - rho))
        starts = [start + [root_logit] for start in starts]
    return starts


def _profile_interval(objective, optimum, index, bounds, max_log_likelihood, transform=math.exp):
    target = max_log_likelihood - PROFILE_DROP_95
    optimum = np.asarray(optimum, dtype=float)

    class ProfileOptimizationFailed(Exception):
        pass

    def profile_at(fixed):
        free_indices = [idx for idx in range(len(optimum)) if idx != index]
        if not free_indices:
            trial = optimum.copy()
            trial[index] = fixed
            value = -float(objective(trial))
            if not math.isfinite(value):
                return None, "profile_optimization_failed"
            return value, "ok"

        def free_objective(free_values):
            trial = optimum.copy()
            trial[index] = fixed
            trial[free_indices] = free_values
            return objective(trial)

        start = optimum[free_indices]
        free_bounds = [bounds[idx] for idx in free_indices]
        result = minimize(free_objective, start, method="L-BFGS-B", bounds=free_bounds)
        if not result.success or not math.isfinite(float(result.fun)):
            return None, "profile_optimization_failed"
        return -float(result.fun), "ok"

    def profile_residual(fixed):
        profiled, status = profile_at(fixed)
        if status != "ok":
            raise ProfileOptimizationFailed
        return profiled - target

    def crossing(direction):
        edge = bounds[index][0] if direction < 0 else bounds[index][1]
        points = np.linspace(optimum[index], edge, 18)[1:]
        previous_x = optimum[index]
        previous_value = max_log_likelihood - target
        for point in points:
            profiled, status = profile_at(float(point))
            if status != "ok":
                return None, status
            value = profiled - target
            if value <= 0 <= previous_value:
                try:
                    root = brentq(profile_residual, float(point), float(previous_x))
                except ProfileOptimizationFailed:
                    return None, "profile_optimization_failed"
                except (ValueError, RuntimeError):
                    return None, "profile_root_failed"
                return root, "closed"
            previous_x = float(point)
            previous_value = value
        return edge, "range_limited"

    lower, lower_status = crossing(-1)
    upper, upper_status = crossing(1)
    status = "two_sided"
    failed_statuses = {lower_status, upper_status} & {"profile_optimization_failed", "profile_root_failed"}
    if failed_statuses:
        if lower_status == "range_limited":
            failed_statuses.add("lower_range_limited")
        if upper_status == "range_limited":
            failed_statuses.add("upper_range_limited")
        status = ";".join(sorted(failed_statuses))
    elif lower_status == "range_limited" and upper_status == "range_limited":
        status = "range_limited_both"
    elif lower_status == "range_limited":
        status = "lower_range_limited"
    elif upper_status == "range_limited":
        status = "upper_range_limited"
    return {
        "low": transform(lower) if lower is not None else None,
        "high": transform(upper) if upper is not None else None,
        "status": status,
        "raw_low": lower,
        "raw_high": upper,
    }


def _observed_information(objective, optimum):
    optimum = np.asarray(optimum, dtype=float)
    n = len(optimum)
    step = 1e-4
    hessian = np.zeros((n, n), dtype=float)
    center = float(objective(optimum))
    for i in range(n):
        unit_i = np.zeros(n)
        unit_i[i] = step
        hessian[i, i] = (objective(optimum + unit_i) - 2 * center + objective(optimum - unit_i)) / (step * step)
        for j in range(i + 1, n):
            unit_j = np.zeros(n)
            unit_j[j] = step
            value = (
                objective(optimum + unit_i + unit_j)
                - objective(optimum + unit_i - unit_j)
                - objective(optimum - unit_i + unit_j)
                + objective(optimum - unit_i - unit_j)
            ) / (4 * step * step)
            hessian[i, j] = hessian[j, i] = value
    eigenvalues = np.linalg.eigvalsh(hessian)
    identifiable = bool(np.all(np.isfinite(eigenvalues)) and np.min(eigenvalues) > 1e-7)
    condition = float(np.max(eigenvalues) / np.min(eigenvalues)) if identifiable else math.inf
    return identifiable and condition < 1e10, eigenvalues, condition


def fit_model(
    tree,
    patterns,
    model,
    foreground_children=frozenset(),
    ascertainment="observed-at-least-one",
    threads=1,
    root_frequency="estimated",
    root_presence=0.5,
    extra_starts=None,
):
    """Fit a CTMC with optional extra starts in log-rate/root-logit coordinates."""
    if root_frequency not in {"estimated", "stationary", "fixed"}:
        raise SystemExit("root_frequency must be estimated, stationary or fixed")
    if not 0 < float(root_presence) < 1:
        raise SystemExit("root_presence must lie strictly between 0 and 1")
    if ascertainment not in {"observed-at-least-one", "complete-universe", "variable-only"}:
        raise SystemExit("unsupported ascertainment mode")
    if ascertainment == "observed-at-least-one":
        for pattern, _weight in _iter_weighted_patterns(patterns):
            if not any(value == 1 for value in pattern.values()):
                raise SystemExit("observed-at-least-one ascertainment requires each structural site to be present in at least one observed tip")
    if ascertainment == "variable-only":
        for pattern, _weight in _iter_weighted_patterns(patterns):
            observed = {value for value in pattern.values() if value in {0, 1}}
            if len(observed) < 2:
                raise SystemExit("variable-only ascertainment requires every included site to vary among observed tips")
    bounds = _model_bounds(model, root_frequency)

    def objective(theta):
        return -_dataset_log_likelihood(
            tree, patterns, model, theta, frozenset(foreground_children), ascertainment,
            root_frequency, root_presence,
        )

    starts = _model_starts(model, root_frequency, root_presence)
    if extra_starts is not None:
        starts.extend(extra_starts)
    if threads > 1 and len(starts) > 1:
        with ThreadPoolExecutor(max_workers=min(threads, len(starts))) as executor:
            results = list(
                executor.map(
                    lambda start: minimize(objective, start, method="L-BFGS-B", bounds=bounds),
                    starts,
                )
            )
    else:
        results = [minimize(objective, start, method="L-BFGS-B", bounds=bounds) for start in starts]
    finite_results = [
        result for result in results
        if math.isfinite(float(result.fun)) and np.all(np.isfinite(result.x))
    ]
    converged_results = [result for result in finite_results if result.success]
    best = min(
        converged_results or finite_results or results,
        key=lambda result: float(result.fun) if math.isfinite(float(result.fun)) else math.inf,
    )
    converged = bool(converged_results)
    theta = np.asarray(best.x, dtype=float)
    log_likelihood = -float(best.fun)
    gain, loss, multiplier, rho = _decode_parameters(theta, model, root_frequency, root_presence)
    boundary = any(
        abs(theta[idx] - lower) < 1e-5 or abs(theta[idx] - upper) < 1e-5
        for idx, (lower, upper) in enumerate(bounds)
    )
    identifiable, information_eigenvalues, information_condition = (
        _observed_information(objective, theta) if converged else (False, np.array([], dtype=float), math.inf)
    )
    if not converged:
        intervals = [{"low": None, "high": None, "status": "optimizer_failed", "raw_low": None, "raw_high": None} for _idx in range(len(theta))]
        fit_status = "optimizer_failed"
    else:
        intervals = []
        for idx in range(len(theta)):
            transform = math.exp
            if root_frequency == "estimated" and idx == len(theta) - 1:
                transform = lambda value: 1.0 / (1.0 + math.exp(-value))
            intervals.append(_profile_interval(objective, theta, idx, bounds, log_likelihood, transform))
        if boundary:
            fit_status = "boundary_limited"
        elif not identifiable:
            fit_status = "nonidentifiable"
        else:
            fit_status = "success"
    inference_status = fit_status
    has_observed_contrast = any(
        weight > 0 and len(_pattern_observed_values(pattern)) > 1
        for pattern, weight in _iter_weighted_patterns(patterns)
    )
    if not has_observed_contrast:
        # Local curvature cannot establish an interior rate estimate without an observed contrast.
        fit_status = "not_estimable"
        inference_status = "no_observed_contrast"
        identifiable = False
    return {
        "model": model,
        "theta": theta,
        "gain_rate": gain,
        "loss_rate": loss,
        "foreground_multiplier": multiplier,
        "root_presence": rho,
        "log_likelihood": log_likelihood,
        "parameter_count": len(theta),
        "aic": 2 * len(theta) - 2 * log_likelihood,
        "converged": converged,
        "fit_status": fit_status,
        "inference_status": inference_status,
        "optimizer_message": str(best.message),
        "boundary": boundary,
        "intervals": intervals,
        "start_count": len(results),
        "identifiable": identifiable,
        "information_eigenvalues": information_eigenvalues,
        "information_condition": information_condition,
        "root_frequency": root_frequency,
        "profile_support_thetas": [],
    }


def _posterior_messages(tree, observations, gain, loss, multiplier, foreground_children, root_presence=None):
    inside, log_likelihood = _inside_log_messages(
        tree, observations, gain, loss, multiplier, foreground_children, root_presence
    )
    if not math.isfinite(log_likelihood):
        raise ValueError("posterior probabilities require a positive-probability observation pattern")
    prior = _root_prior(gain, loss, root_presence)
    outside = {tree.root: _log_array(prior)}
    node = {}
    edge = {}
    for current in tree.preorder():
        node_log = outside[current] + inside[current]
        node[current] = np.exp(node_log - np.logaddexp(node_log[0], node_log[1]))
        children = list(tree.children.get(current, []))
        if not children:
            continue
        child_messages = []
        child_log_matrices = {}
        for child in children:
            child_multiplier = multiplier if child in foreground_children else 1.0
            matrix = _transition_matrix(gain, loss, tree.branch_length(child), child_multiplier)
            log_matrix = _log_array(matrix)
            child_log_matrices[child] = log_matrix
            child_messages.append(
                np.logaddexp(
                    log_matrix[:, 0] + inside[child][0],
                    log_matrix[:, 1] + inside[child][1],
                )
            )
        prefix = [np.zeros(2, dtype=float)]
        for message in child_messages:
            prefix.append(prefix[-1] + message)
        suffix = [np.zeros(2, dtype=float) for _child in children]
        running = np.zeros(2, dtype=float)
        for idx in range(len(children) - 1, -1, -1):
            suffix[idx] = running
            running = running + child_messages[idx]
        for idx, child in enumerate(children):
            parent_context = outside[current] + prefix[idx] + suffix[idx]
            log_matrix = child_log_matrices[child]
            joint_log = parent_context[:, None] + log_matrix + inside[child][None, :]
            joint = np.exp(joint_log - logsumexp(joint_log))
            edge[(current, child)] = joint
            outside[child] = np.logaddexp(
                parent_context[0] + log_matrix[0, :],
                parent_context[1] + log_matrix[1, :],
            )
    return node, edge


@lru_cache(maxsize=32768)
def _conditional_transition_count(gain, loss, branch_length, multiplier, start, end, src, dst):
    matrix = _transition_matrix(gain, loss, branch_length, multiplier)
    denominator = float(matrix[start, end])
    if denominator <= 0.0:
        return 0.0
    gain_rate = float(gain) * float(multiplier)
    loss_rate = float(loss) * float(multiplier)
    q = np.array([[-gain_rate, gain_rate], [loss_rate, -loss_rate]], dtype=float)
    reward = np.zeros((2, 2), dtype=float)
    reward[src, dst] = gain_rate if (src, dst) == (0, 1) else loss_rate
    block = np.zeros((4, 4), dtype=float)
    block[:2, :2] = q
    block[:2, 2:] = reward
    block[2:, 2:] = q
    integral = float(expm(block * float(branch_length))[:2, 2:][start, end])
    return max(0.0, integral / denominator)


def _expected_transition_count(joint, gain, loss, branch_length, multiplier, src, dst):
    total = 0.0
    for start in (0, 1):
        for end in (0, 1):
            total += float(joint[start, end]) * _conditional_transition_count(
                gain, loss, branch_length, multiplier, start, end, src, dst
            )
    return total


def _bh_adjust(rows):
    by_test = defaultdict(list)
    for index, row in enumerate(rows):
        try:
            by_test[row["test_id"]].append((index, float(row["p_value"])))
        except (TypeError, ValueError):
            row["q_value"] = "NA"
            row["q_value_method"] = "not_available"
    for entries in by_test.values():
        if len(entries) == 1:
            rows[entries[0][0]]["q_value"] = _fmt(entries[0][1])
            rows[entries[0][0]]["q_value_method"] = "Benjamini-Hochberg_m_equals_1"
            continue
        ordered = sorted(entries, key=lambda item: item[1])
        adjusted = [0.0] * len(ordered)
        running = 1.0
        for rank in range(len(ordered), 0, -1):
            value = min(running, ordered[rank - 1][1] * len(ordered) / rank)
            adjusted[rank - 1] = value
            running = value
        for (index, _p), value in zip(ordered, adjusted):
            rows[index]["q_value"] = _fmt(value)
            rows[index]["q_value_method"] = "Benjamini-Hochberg_within_test_type"


def _read_foreground_children(path, tree):
    if not path:
        return frozenset()
    rows = read_tsv(path)
    label_to_node = {label: node for node, label in tree.label.items()}
    children = set()
    invalid = []
    edge_set = set(tree.edges())
    for row in rows:
        parent = row.get("parent_node") or row.get("parent_id") or ""
        child = row.get("child_node") or row.get("child_id") or ""
        if (parent, child) in edge_set:
            children.add(child)
            continue
        scope = row.get("branch_scope") or row.get("branch") or ""
        if "->" in scope:
            parent_label, child_label = (part.strip() for part in scope.split("->", 1))
            parent_node = label_to_node.get(parent_label)
            child_node = label_to_node.get(child_label)
            if (parent_node, child_node) in edge_set:
                children.add(child_node)
                continue
        invalid.append(scope or f"{parent}->{child}")
    if invalid:
        raise SystemExit("foreground branch file contains unmatched branches: " + ", ".join(invalid))
    if not children:
        raise SystemExit("foreground branch file did not match any branch in species_tree.tsv")
    if len(children) == len(edge_set):
        raise SystemExit("foreground model is unidentifiable when every branch is foreground")
    return frozenset(children)


def infer_single_copy_phylogeny(
    input_dir,
    output_dir,
    model="er-ard",
    foreground_branches=None,
    branch_length_mode="supplied",
    ascertainment="observed-at-least-one",
    threads=1,
    root_frequency="estimated",
    root_presence=0.5,
):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    site_rows, excluded_rows = build_structural_site_matrix(input_dir, output_dir)
    write_tsv(
        output_dir / "excluded_families.tsv",
        excluded_rows,
        ["family_id", "species", "copy_count", "gene_copy_ids", "reason"],
    )
    write_tsv(
        output_dir / "structural_site_matrix.tsv",
        site_rows,
        ["family_id", "layer", "site_id", "species", "state", "state_0", "state_1", "evidence"],
    )
    tree_rows = _validated_tree_rows(input_dir / "species_tree.tsv", branch_length_mode)
    tree = SpeciesTree(tree_rows)
    matrix_species = {row["species"] for row in site_rows}
    missing_tips = sorted(matrix_species - set(tree.leaf_by_label))
    if missing_tips:
        raise SystemExit(
            "species_tree.tsv lacks structural-matrix species: " + ", ".join(missing_tips)
        )
    foreground_children = _read_foreground_children(foreground_branches, tree)
    if model == "foreground" and not foreground_children:
        raise SystemExit("--model foreground requires --foreground-branches")
    site_universe = None
    if ascertainment == "complete-universe":
        site_universe = _load_structural_site_universe(input_dir)
        if site_universe is None:
            raise SystemExit(
                "complete-universe ascertainment requires input structural_site_universe.tsv with family_id, layer and site_id"
            )

    grouped = defaultdict(lambda: defaultdict(dict))
    state_labels = {}
    for row in site_rows:
        grouped[(row["family_id"], row["layer"])][row["site_id"]][row["species"]] = row["state"]
        state_labels[(row["family_id"], row["layer"])] = (row["state_0"], row["state_1"])

    fit_rows = []
    test_rows = []
    node_rows = []
    branch_rows = []
    change_rows = []
    for (family, layer), sites in sorted(grouped.items()):
        state_0, state_1 = state_labels[(family, layer)]
        raw_patterns = []
        encoded_sites = []
        for site_id, observations in sorted(sites.items()):
            if site_universe is not None and (family, layer, site_id) not in site_universe:
                continue
            encoded = {
                label: 0 if observations.get(label) == state_0 else 1 if observations.get(label) == state_1 else "unknown"
                for label in tree.leaf_by_label
            }
            observed_count = sum(value in {0, 1} for value in encoded.values())
            if observed_count >= 2 and _selected_for_ascertainment(encoded, ascertainment):
                raw_patterns.append(encoded)
                encoded_sites.append((site_id, encoded))
        patterns = _compress_patterns(raw_patterns, sorted(tree.leaf_by_label))
        if not patterns:
            continue
        site_count = sum(weight for _pattern, weight in patterns)
        compressed_pattern_count = len(patterns)
        informative = sum(
            weight
            for pattern, weight in patterns
            if len({value for value in pattern.values() if value in {0, 1}}) > 1
        )
        observed_taxa = len(
            {species for pattern, _weight in patterns for species, value in pattern.items() if value in {0, 1}}
        )
        if model == "foreground":
            null_fit = fit_model(tree, patterns, "ARD", ascertainment=ascertainment, threads=threads, root_frequency=root_frequency, root_presence=root_presence)
            extra_starts = []
            if null_fit["converged"]:
                null_theta = null_fit["theta"]
                extra_starts.append([null_theta[0], null_theta[1], 0.0, *null_theta[2:]])
            alternative_fit = fit_model(
                tree,
                patterns,
                "ARD_FOREGROUND",
                foreground_children=foreground_children,
                ascertainment=ascertainment,
                threads=threads,
                root_frequency=root_frequency,
                root_presence=root_presence,
                extra_starts=extra_starts,
            )
            test_id = "homogeneous_vs_foreground"
        else:
            null_fit = fit_model(tree, patterns, "ER", ascertainment=ascertainment, threads=threads, root_frequency=root_frequency, root_presence=root_presence)
            extra_starts = []
            if null_fit["converged"]:
                null_theta = null_fit["theta"]
                extra_starts.append([null_theta[0], null_theta[0], *null_theta[1:]])
            alternative_fit = fit_model(
                tree,
                patterns,
                "ARD",
                ascertainment=ascertainment,
                threads=threads,
                root_frequency=root_frequency,
                root_presence=root_presence,
                extra_starts=extra_starts,
            )
            test_id = "equal_rates_vs_gain_loss"

        for fit in (null_fit, alternative_fit):
            gain_ci = fit["intervals"][0]
            loss_ci = fit["intervals"][0] if fit["model"] == "ER" else fit["intervals"][1]
            multiplier_ci = fit["intervals"][2] if fit["model"] == "ARD_FOREGROUND" else {"low": 1.0, "high": 1.0, "status": "fixed"}
            root_ci = fit["intervals"][-1] if root_frequency == "estimated" else {"low": fit["root_presence"], "high": fit["root_presence"], "status": root_frequency}
            fit_rows.append(
                {
                    "family_id": family,
                    "layer": layer,
                    "model": fit["model"],
                    "n_taxa": observed_taxa,
                    "n_structural_sites": site_count,
                    "n_compressed_patterns": compressed_pattern_count,
                    "n_informative_patterns": informative,
                    "gain_rate": _fmt(fit["gain_rate"]),
                    "gain_rate_ci_low": _fmt(gain_ci["low"]),
                    "gain_rate_ci_high": _fmt(gain_ci["high"]),
                    "gain_rate_ci_status": gain_ci["status"],
                    "loss_rate": _fmt(fit["loss_rate"]),
                    "loss_rate_ci_low": _fmt(loss_ci["low"]),
                    "loss_rate_ci_high": _fmt(loss_ci["high"]),
                    "loss_rate_ci_status": loss_ci["status"],
                    "foreground_multiplier": _fmt(fit["foreground_multiplier"]),
                    "foreground_multiplier_ci_low": _fmt(multiplier_ci["low"]),
                    "foreground_multiplier_ci_high": _fmt(multiplier_ci["high"]),
                    "foreground_multiplier_ci_status": multiplier_ci["status"],
                    "root_presence": _fmt(fit["root_presence"]),
                    "root_presence_ci_low": _fmt(root_ci["low"]),
                    "root_presence_ci_high": _fmt(root_ci["high"]),
                    "root_presence_ci_status": root_ci["status"],
                    "root_frequency_mode": root_frequency,
                    "log_likelihood": _fmt(fit["log_likelihood"]),
                    "parameter_count": fit["parameter_count"],
                    "aic": _fmt(fit["aic"]),
                    "converged": str(fit["converged"]).lower(),
                    "fit_status": fit["fit_status"],
                    "inference_status": fit.get("inference_status", fit["fit_status"]),
                    "parameter_at_boundary": str(fit["boundary"]).lower(),
                    "identifiable": str(fit["identifiable"]).lower(),
                    "information_condition": _fmt(fit["information_condition"]),
                    "optimizer_starts": fit["start_count"],
                    "optimizer_message": fit["optimizer_message"],
                    "ascertainment": ascertainment,
                }
            )

        likelihood_difference = alternative_fit["log_likelihood"] - null_fit["log_likelihood"]
        lrt = 2.0 * likelihood_difference if likelihood_difference >= 0.0 else None
        estimable = (
            informative >= 1
            and null_fit["converged"]
            and alternative_fit["converged"]
            and not null_fit["boundary"]
            and not alternative_fit["boundary"]
            and null_fit["identifiable"]
            and alternative_fit["identifiable"]
            and null_fit["fit_status"] == "success"
            and alternative_fit["fit_status"] == "success"
            and lrt is not None
        )
        p_value = float(chi2.sf(lrt, 1)) if estimable else None
        if informative == 0:
            test_status = "parameters_not_estimable"
        elif likelihood_difference < -1e-7:
            test_status = "optimization_failure_alternative_below_null"
        elif estimable:
            test_status = "tested"
        else:
            test_status = "parameters_not_estimable"
        test_row = {
            "family_id": family,
            "layer": layer,
            "test_id": test_id,
            "null_model": null_fit["model"],
            "alternative_model": alternative_fit["model"],
            "null_log_likelihood": _fmt(null_fit["log_likelihood"]),
            "alternative_log_likelihood": _fmt(alternative_fit["log_likelihood"]),
            "lrt_statistic": _fmt(lrt),
            "df": 1,
            "p_value": _fmt(p_value),
            "q_value": "NA",
            "q_value_method": "not_available",
            "test_status": test_status,
            "inference_status": "no_observed_contrast" if informative == 0 else test_status,
            "reference_distribution": "chi_square_df1_asymptotic_regular_interior" if estimable else "not_available",
            "n_taxa": observed_taxa,
            "n_structural_sites": site_count,
            "n_compressed_patterns": compressed_pattern_count,
            "n_informative_patterns": informative,
        }
        test_rows.append(test_row)
        eligible_fits = [
            fit for fit in (null_fit, alternative_fit)
            if informative > 0 and _fit_valid_for_posterior(fit)
        ]
        selected_fit = min(eligible_fits, key=lambda fit: fit["aic"]) if eligible_fits else None
        if selected_fit is None:
            diagnostic_fit = alternative_fit if informative else null_fit
            status = diagnostic_fit.get("fit_status", "not_estimable")
            for site_id, observations in encoded_sites:
                known = sum(value in {0, 1} for value in observations.values())
                change_rows.append(
                    {
                        "family_id": family,
                        "layer": layer,
                        "site_id": site_id,
                        "parent_node": "NA",
                        "child_node": "NA",
                        "branch_scope": "NA",
                        "structural_change_type": "posterior_not_reported",
                        "structural_pattern": f"{state_0}<->{state_1}",
                        "endpoint_change_probability": "NA",
                        "gain_endpoint_probability": "NA",
                        "loss_endpoint_probability": "NA",
                        "direction_probability": "NA",
                        "expected_gain_count": "NA",
                        "expected_loss_count": "NA",
                        "model": diagnostic_fit.get("model", "NA"),
                        "rate_test_status": test_status,
                        "conditioning": (
                            f"posterior_not_reported;fit_status={status};"
                            f"inference_status={test_row['inference_status']};"
                            f"known_tip_count={known};requires=successful_identifiable_interior_fit"
                        ),
                    }
                )
            continue
        posterior_conditioning = (
            "conditional_MLE;"
            f"fit_status={selected_fit['fit_status']};"
            "uncertainty=sensitivity_not_estimated"
        )

        for site_id, observations in encoded_sites:
            site_log_likelihood = _pattern_log_likelihood(
                tree,
                observations,
                selected_fit["gain_rate"],
                selected_fit["loss_rate"],
                selected_fit["foreground_multiplier"],
                foreground_children,
                selected_fit["root_presence"],
            )
            if not math.isfinite(site_log_likelihood):
                known = sum(value in {0, 1} for value in observations.values())
                change_rows.append(
                    {
                        "family_id": family,
                        "layer": layer,
                        "site_id": site_id,
                        "parent_node": "NA",
                        "child_node": "NA",
                        "branch_scope": "NA",
                        "structural_change_type": "posterior_not_reported",
                        "structural_pattern": f"{state_0}<->{state_1}",
                        "endpoint_change_probability": "NA",
                        "gain_endpoint_probability": "NA",
                        "loss_endpoint_probability": "NA",
                        "direction_probability": "NA",
                        "expected_gain_count": "NA",
                        "expected_loss_count": "NA",
                        "model": selected_fit["model"],
                        "rate_test_status": test_status,
                        "conditioning": (
                            "posterior_not_reported;site_likelihood_zero;"
                            f"known_tip_count={known};requires=positive_probability_observation_pattern"
                        ),
                    }
                )
                continue
            node_posterior, edge_posterior = _posterior_messages(
                tree,
                observations,
                selected_fit["gain_rate"],
                selected_fit["loss_rate"],
                selected_fit["foreground_multiplier"],
                foreground_children,
                selected_fit["root_presence"],
            )
            for node_id, probabilities in node_posterior.items():
                for index, probability in enumerate(probabilities):
                    node_rows.append(
                        {
                            "family_id": family,
                            "layer": layer,
                            "site_id": site_id,
                            "node_id": node_id,
                            "node_label": tree.label[node_id],
                            "state": state_0 if index == 0 else state_1,
                            "posterior_probability": _fmt(probability),
                            "profile_probability_low": "NA",
                            "profile_probability_high": "NA",
                            "uncertainty_status": "sensitivity_not_estimated",
                            "model": selected_fit["model"],
                            "conditioning": posterior_conditioning,
                        }
                    )
            for (parent, child), joint in edge_posterior.items():
                branch_multiplier = (
                    selected_fit["foreground_multiplier"] if child in foreground_children else 1.0
                )
                expected_gain = _expected_transition_count(
                    joint,
                    selected_fit["gain_rate"],
                    selected_fit["loss_rate"],
                    tree.branch_length(child),
                    branch_multiplier,
                    0,
                    1,
                )
                expected_loss = _expected_transition_count(
                    joint,
                    selected_fit["gain_rate"],
                    selected_fit["loss_rate"],
                    tree.branch_length(child),
                    branch_multiplier,
                    1,
                    0,
                )
                gain_probability = float(joint[0, 1])
                loss_probability = float(joint[1, 0])
                change_probability = gain_probability + loss_probability
                for src, dst, probability in (
                    (state_0, state_1, gain_probability),
                    (state_1, state_0, loss_probability),
                ):
                    branch_rows.append(
                        {
                            "family_id": family,
                            "layer": layer,
                            "site_id": site_id,
                            "parent_node": parent,
                            "child_node": child,
                            "parent_label": tree.label[parent],
                            "child_label": tree.label[child],
                            "branch_length": _fmt(tree.branch_length(child)),
                            "from_state": src,
                            "to_state": dst,
                            "endpoint_transition_probability": _fmt(probability),
                            "profile_transition_probability_low": "NA",
                            "profile_transition_probability_high": "NA",
                            "total_endpoint_change_probability": _fmt(change_probability),
                            "profile_total_change_probability_low": "NA",
                            "profile_total_change_probability_high": "NA",
                            "expected_gain_count": _fmt(expected_gain),
                            "expected_loss_count": _fmt(expected_loss),
                            "uncertainty_status": "sensitivity_not_estimated",
                            "model": selected_fit["model"],
                            "conditioning": posterior_conditioning,
                        }
                    )
                change_rows.append(
                    {
                        "family_id": family,
                        "layer": layer,
                        "site_id": site_id,
                        "parent_node": parent,
                        "child_node": child,
                        "branch_scope": f"{tree.label[parent]}->{tree.label[child]}",
                        "structural_change_type": "bidirectional_transition_probabilities",
                        "structural_pattern": f"{state_0}<->{state_1}",
                        "endpoint_change_probability": _fmt(change_probability),
                        "gain_endpoint_probability": _fmt(gain_probability),
                        "loss_endpoint_probability": _fmt(loss_probability),
                        "direction_probability": "NA",
                        "expected_gain_count": _fmt(expected_gain),
                        "expected_loss_count": _fmt(expected_loss),
                        "model": selected_fit["model"],
                        "rate_test_status": test_status,
                        "conditioning": posterior_conditioning,
                    }
                )

    _bh_adjust(test_rows)
    write_tsv(
        output_dir / "model_fits.tsv",
        fit_rows,
        [
            "family_id", "layer", "model", "n_taxa", "n_structural_sites", "n_compressed_patterns",
            "n_informative_patterns",
            "gain_rate", "gain_rate_ci_low", "gain_rate_ci_high", "gain_rate_ci_status",
            "loss_rate", "loss_rate_ci_low", "loss_rate_ci_high", "loss_rate_ci_status",
            "foreground_multiplier", "foreground_multiplier_ci_low", "foreground_multiplier_ci_high",
            "foreground_multiplier_ci_status", "root_presence", "root_presence_ci_low",
            "root_presence_ci_high", "root_presence_ci_status", "root_frequency_mode", "log_likelihood",
            "parameter_count", "aic", "converged", "fit_status", "parameter_at_boundary", "identifiable",
            "information_condition", "optimizer_starts", "optimizer_message", "ascertainment", "inference_status",
        ],
    )
    write_tsv(
        output_dir / "model_tests.tsv",
        test_rows,
        [
            "family_id", "layer", "test_id", "null_model", "alternative_model", "null_log_likelihood",
            "alternative_log_likelihood", "lrt_statistic", "df", "p_value", "q_value", "q_value_method",
            "test_status", "reference_distribution", "n_taxa", "n_structural_sites", "n_compressed_patterns",
            "n_informative_patterns", "inference_status",
        ],
    )
    write_tsv(
        output_dir / "node_state_posteriors.tsv",
        node_rows,
        ["family_id", "layer", "site_id", "node_id", "node_label", "state", "posterior_probability",
         "profile_probability_low", "profile_probability_high", "uncertainty_status", "model", "conditioning"],
    )
    write_tsv(
        output_dir / "branch_transition_posteriors.tsv",
        branch_rows,
        [
            "family_id", "layer", "site_id", "parent_node", "child_node", "parent_label", "child_label",
            "branch_length", "from_state", "to_state", "endpoint_transition_probability",
            "profile_transition_probability_low", "profile_transition_probability_high",
            "total_endpoint_change_probability", "profile_total_change_probability_low",
            "profile_total_change_probability_high", "expected_gain_count", "expected_loss_count",
            "uncertainty_status", "model", "conditioning",
        ],
    )
    write_tsv(
        output_dir / "structural_changes.tsv",
        change_rows,
        [
            "family_id", "layer", "site_id", "parent_node", "child_node", "branch_scope",
            "structural_change_type", "structural_pattern", "endpoint_change_probability",
            "gain_endpoint_probability", "loss_endpoint_probability", "direction_probability",
            "expected_gain_count", "expected_loss_count", "model", "rate_test_status", "conditioning",
        ],
    )
    write_tsv(
        output_dir / "phylogeny_scope.tsv",
        [{
            "scope": "single_copy_structural_sites",
            "tree_file": "species_tree.tsv",
            "tree_scope": "species_tree",
            "layers": "exon_presence;exon_role;splice_junction",
            "note": "All homologous structural sites within each layer share fitted gain/loss parameters.",
        }],
        ["scope", "tree_file", "tree_scope", "layers", "note"],
    )
    parameters = {
        "analysis_scope": "single-copy",
        "model_test": model,
        "branch_length_mode": branch_length_mode,
        "ascertainment": ascertainment,
        "root_frequency": root_frequency,
        "root_presence_when_fixed": root_presence if root_frequency == "fixed" else None,
        "tree_file": str(input_dir / "species_tree.tsv"),
        "foreground_branches": str(foreground_branches) if foreground_branches else None,
        "structural_site_universe": str(input_dir / "structural_site_universe.tsv") if ascertainment == "complete-universe" else None,
        "threads": max(1, int(threads)),
        "fixed_inputs": ["species_tree", "ortholog_set", "structural_site_states"],
        "estimated_parameters": ["gain_rate", "loss_rate"] + (["root_presence"] if root_frequency == "estimated" else []) + (["foreground_multiplier"] if model == "foreground" else []),
        "excluded_family_count": len({row["family_id"] for row in excluded_rows}),
    }
    (output_dir / "run_parameters.json").write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")
    return site_rows, fit_rows, test_rows, change_rows
