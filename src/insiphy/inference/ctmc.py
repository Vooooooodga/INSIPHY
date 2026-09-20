"""inference / ctmc: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from insiphy.inference.parameters import _decode_parameters
import math
import numpy as np


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
