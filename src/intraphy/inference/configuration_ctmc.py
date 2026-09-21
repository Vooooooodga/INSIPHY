"""Finite-state log-space pruning and inside/outside posteriors.

Supports non-reversible, branch-specific transition matrices. Expected edits use
marked transition integrals, not differences between endpoint states.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm
from scipy.special import logsumexp


@dataclass(frozen=True)
class ConfigurationLikelihood:
    log_likelihood: float
    nodes: dict[str, np.ndarray]
    endpoints: dict[tuple[str, str], np.ndarray]


def transition_matrix(q: np.ndarray, length: float) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    if q.ndim != 2 or q.shape[0] != q.shape[1] or not np.isfinite(q).all():
        raise ValueError("A finite square generator is required")
    off = q.copy()
    np.fill_diagonal(off, 0)
    if (off < 0).any() or not np.allclose(q.sum(axis=1), 0, atol=1e-10):
        raise ValueError("Invalid CTMC generator")
    if not np.isfinite(length) or length < 0:
        raise ValueError("Branch length must be finite and nonnegative")
    p = expm(q * length)
    if p.min() < -1e-10 or not np.allclose(p.sum(axis=1), 1, atol=1e-9):
        raise ArithmeticError("Transition calculation lost stochasticity")
    # Only roundoff is clipped; no truncation renormalization is performed.
    return np.maximum(p, 0)


def _log(values: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore"):
        return np.log(values)


def likelihood(tree, tips: dict[str, np.ndarray], transitions: dict[str, np.ndarray],
               root_prior: np.ndarray, posterior: bool = True) -> ConfigurationLikelihood:
    root_prior = np.asarray(root_prior, dtype=float)
    n = len(root_prior)
    if n == 0 or not np.isfinite(root_prior).all() or (root_prior < 0).any() or not np.isclose(root_prior.sum(), 1):
        raise ValueError("Root prior must be an explicit normalized nonnegative vector")
    if set(tips) != set(tree.leaf_by_label):
        raise ValueError("Exactly one likelihood vector for every tree tip is required")
    for label, values in tips.items():
        if len(values) != n or not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f"Invalid observation likelihood: {label}")
    logp = {}
    for _, child in tree.edges():
        p = np.asarray(transitions[child])
        if p.shape != (n, n) or not np.isfinite(p).all() or (p < 0).any() or not np.allclose(p.sum(axis=1), 1, atol=1e-9):
            raise ValueError("Invalid transition matrix")
        logp[child] = _log(p)
    inside, messages, siblings = {}, {}, {}
    for node in tree.postorder():
        children = tuple(tree.children.get(node, ()))
        if not children:
            inside[node] = _log(tips[tree.label[node]])
        else:
            terms = []
            for child in children:
                messages[child] = logsumexp(logp[child] + inside[child][None, :], axis=1)
                terms.append(messages[child])
            inside[node] = sum(terms, start=np.zeros(n))
            for j, child in enumerate(children):
                siblings[child] = sum((v for k, v in enumerate(terms) if k != j), start=np.zeros(n))
    ll = float(logsumexp(_log(root_prior) + inside[tree.root]))
    if not posterior or not np.isfinite(ll):
        return ConfigurationLikelihood(ll, {}, {})
    outside = {tree.root: _log(root_prior)}
    endpoints = {}
    for parent in tree.preorder():
        for child in tree.children.get(parent, ()):
            base = outside[parent] + siblings[child]
            outside[child] = logsumexp(base[:, None] + logp[child], axis=0)
            endpoints[(parent, child)] = np.exp(base[:, None] + logp[child] + inside[child][None, :] - ll)
    nodes = {node: np.exp(outside[node] + inside[node] - ll) for node in tree.preorder()}
    return ConfigurationLikelihood(ll, nodes, endpoints)


def marked_integral(q: np.ndarray, marked: np.ndarray, length: float) -> np.ndarray:
    """Integral exp(Qs) B exp(Q(t-s)) ds, with B containing selected edit rates."""
    n = len(q)
    if marked.shape != q.shape or np.any(marked < 0) or np.any(np.diag(marked) != 0):
        raise ValueError("A nonnegative off-diagonal marked-rate matrix is required")
    block = np.zeros((2*n, 2*n))
    block[:n, :n] = q
    block[n:, n:] = q
    block[:n, n:] = marked
    return expm(block*length)[:n, n:]


def expected_count(q, p, endpoints, length, marked=None) -> float:
    rates = q.copy() if marked is None else marked.copy()
    np.fill_diagonal(rates, 0)
    integral = marked_integral(q, rates, length)
    bridge = np.divide(integral, p, out=np.zeros_like(p), where=p > 0)
    value = float(np.sum(endpoints * bridge))
    if value < -1e-8 or not np.isfinite(value):
        raise ArithmeticError("Invalid expected edit count")
    return max(0., value)


def probability_any_change(q, p, endpoints, length) -> float:
    no_change = np.zeros_like(p)
    np.fill_diagonal(no_change, np.exp(np.diag(q)*length))
    bridge = np.divide(no_change, p, out=np.zeros_like(p), where=p > 0)
    value = 1. - float(np.sum(endpoints * bridge))
    if not -1e-8 <= value <= 1+1e-8:
        raise ArithmeticError("Invalid probability of any change")
    return min(1., max(0., value))
