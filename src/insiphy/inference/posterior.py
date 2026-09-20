"""inference / posterior: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from functools import lru_cache
from insiphy.inference.ctmc import _inside_log_messages
from insiphy.inference.ctmc import _log_array
from insiphy.inference.ctmc import _root_prior
from insiphy.inference.ctmc import _transition_matrix
from insiphy.inference.diagnostics import _posterior_sensitivity_status
from insiphy.inference.parameters import _decode_parameters
from scipy.linalg import expm
from scipy.special import logsumexp
import math
import numpy as np


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


def _profile_posterior_envelope(
    tree,
    observations,
    fit,
    foreground_children,
    mle_node,
    mle_edge,
):
    sensitivity_status = fit.get(
        "posterior_sensitivity_status", _posterior_sensitivity_status(fit)
    )
    if sensitivity_status != "finite_profile_envelope":
        return None, None, sensitivity_status
    node_samples = {node: [np.asarray(values, dtype=float)] for node, values in mle_node.items()}
    edge_samples = {edge: [np.asarray(values, dtype=float)] for edge, values in mle_edge.items()}
    try:
        for theta in fit.get("profile_support_thetas", []):
            gain, loss, multiplier, rho = _decode_parameters(
                theta,
                fit["model"],
                fit.get("root_frequency", "estimated"),
                fit.get("root_presence", 0.5),
            )
            node, edge = _posterior_messages(
                tree, observations, gain, loss, multiplier, foreground_children, rho
            )
            for node_id, values in node.items():
                node_samples[node_id].append(np.asarray(values, dtype=float))
            for edge_id, values in edge.items():
                edge_samples[edge_id].append(np.asarray(values, dtype=float))
    except (ValueError, FloatingPointError):
        return None, None, "unavailable_profile_posterior_evaluation_failed"
    node_envelope = {
        node_id: (np.min(values, axis=0), np.max(values, axis=0))
        for node_id, values in node_samples.items()
    }
    edge_envelope = {}
    for edge_id, values in edge_samples.items():
        array = np.asarray(values, dtype=float)
        total_change = array[:, 0, 1] + array[:, 1, 0]
        edge_envelope[edge_id] = {
            "low": np.min(array, axis=0),
            "high": np.max(array, axis=0),
            "total_low": float(np.min(total_change)),
            "total_high": float(np.max(total_change)),
        }
    return node_envelope, edge_envelope, "finite_profile_envelope"
