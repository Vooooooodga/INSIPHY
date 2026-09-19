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

from .io import (
    STRUCTURAL_SITE_SCHEMA_VERSION,
    read_structural_site_matrix,
    read_tsv,
    structural_site_observed_state,
    validate_structural_site_tip_rows,
    write_structural_site_matrix,
    write_tsv,
)
from .structural_sites import (
    _single_copy_families,
    build_structural_site_matrix,
    structural_matrix_annotation_view,
)
from .tree import SpeciesTree


RATE_MIN = 1e-8
RATE_MAX = 100.0
MULTIPLIER_MIN = 1e-3
MULTIPLIER_MAX = 1e3
PROFILE_DROP_95 = 0.5 * float(chi2.ppf(0.95, 1))
COMPLETE_UNIVERSE_RULES = {"independent_catalogue", "curated_complete_universe"}
GENERATED_ALL_ZERO_MASKS = {"generated_all_zero", "explicit_all_zero"}


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


def _posterior_unavailable_reason(fit):
    if fit.get("inference_status") == "no_observed_states":
        return "no_observed_states"
    if fit.get("inference_status") == "no_observed_contrast":
        return "no_observed_contrast"
    if not bool(fit.get("converged")):
        return "optimizer_failed"
    if bool(fit.get("boundary")):
        return "parameter_at_boundary"
    if not bool(fit.get("identifiable")):
        return "parameters_not_identifiable"
    if fit.get("fit_status") != "success":
        return str(fit.get("fit_status") or "fit_not_successful")
    return "NA"


def _fit_valid_for_posterior(fit):
    return _posterior_unavailable_reason(fit) == "NA"


def _posterior_sensitivity_status(fit):
    statuses = {interval.get("status", "") for interval in fit.get("intervals", [])}
    if any("range_limited" in status for status in statuses):
        return "unavailable_open_profile_interval"
    if any(status in {"profile_optimization_failed", "profile_root_failed", "optimizer_failed"}
           or status.startswith("profile_optimization_failed;")
           or status.startswith("profile_root_failed;")
           for status in statuses):
        return "unavailable_profile_failure"
    endpoints_complete = bool(fit.get("intervals")) and all(
        interval.get("theta_low") is not None and interval.get("theta_high") is not None
        for interval in fit.get("intervals", [])
    )
    if statuses and statuses == {"two_sided"} and endpoints_complete:
        return "finite_profile_envelope"
    if statuses and statuses == {"two_sided"}:
        return "unavailable_profile_endpoints_missing"
    return "sensitivity_not_estimated"


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


def _transition_event_type(layer, source_index, target_index):
    if source_index == target_index:
        return "no_change"
    if layer == "splice_junction":
        return "intron_gain" if (source_index, target_index) == (0, 1) else "intron_loss"
    if layer == "exon_presence":
        return "sequence_gain" if (source_index, target_index) == (0, 1) else "sequence_loss"
    if layer == "exon_role":
        return "exon_role_gain" if (source_index, target_index) == (0, 1) else "exon_role_loss"
    return "state_gain" if (source_index, target_index) == (0, 1) else "state_loss"


def _transition_structural_relation(layer, site_kind, source_index, target_index):
    if layer != "splice_junction" or site_kind not in {
        "within_element_junction",
        "within_exon_boundary",
    }:
        return "NA"
    if (source_index, target_index) == (0, 1):
        return "exon_split"
    if (source_index, target_index) == (1, 0):
        return "exon_fusion"
    return "NA"


def _load_site_matrix(
    input_dir,
    output_dir,
    structural_site_matrix_path,
    annotation_view,
):
    if structural_site_matrix_path is not None:
        source = Path(structural_site_matrix_path)
        site_rows = read_structural_site_matrix(source)
        excluded_rows = read_tsv(source.parent / "excluded_families.tsv", optional=True)
        output_matrix = Path(output_dir) / "structural_site_matrix.tsv"
        if source.resolve() != output_matrix.resolve():
            write_structural_site_matrix(output_matrix, site_rows)
        return site_rows, excluded_rows, source, "frozen_matrix"
    site_rows, excluded_rows = build_structural_site_matrix(
        input_dir,
        output_dir,
        annotation_view=annotation_view,
    )
    source = Path(output_dir) / "structural_site_matrix.tsv"
    write_structural_site_matrix(source, site_rows)
    return read_structural_site_matrix(source), excluded_rows, source, "prepared_in_analysis"


def _masked_state(row):
    return structural_site_observed_state(row)


def _validate_complete_universe(site_rows, tree):
    rows_by_site = defaultdict(dict)
    for row in site_rows:
        rows_by_site[_site_key(row)][row.get("species", "NA")] = row
    missing = []
    for key, by_species in sorted(rows_by_site.items()):
        for species in sorted(tree.leaf_by_label):
            row = by_species.get(species)
            if row is None:
                missing.append("/".join((*key, species)))
                continue
            rule = str(row.get("discovery_rule", "")).strip()
            if rule not in COMPLETE_UNIVERSE_RULES:
                missing.append("/".join((*key, species)) + ":nonindependent_discovery_rule")
                continue
            state = row.get("state", "unknown")
            state_0 = row.get("state_0", "NA")
            state_1 = row.get("state_1", "NA")
            mask = str(row.get("observation_mask", "observed")).strip().lower()
            explicit_observation = mask == "observed" and state in {state_0, state_1}
            generated_all_zero = mask in GENERATED_ALL_ZERO_MASKS and state == state_0
            if not (explicit_observation or generated_all_zero):
                missing.append("/".join((*key, species)) + ":state_not_explicitly_enumerated")
    if missing:
        preview = ", ".join(missing[:8])
        suffix = " ..." if len(missing) > 8 else ""
        raise SystemExit(
            "complete-universe ascertainment requires independent_catalogue or "
            "curated_complete_universe rows with an observed state, or an explicitly generated "
            "all-zero state, for every family/layer/site/species combination; unavailable: "
            + preview + suffix
        )


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


def _comparison_models(model):
    if model == "foreground":
        return "ARD", "ARD_FOREGROUND", "homogeneous_vs_foreground"
    return "ER", "ARD", "equal_rates_vs_gain_loss"


def _unavailable_fit(model, root_frequency, root_presence, reason):
    parameter_count = len(_model_bounds(model, root_frequency))
    intervals = [
        {
            "low": None, "high": None, "status": reason, "raw_low": None,
            "raw_high": None, "theta_low": None, "theta_high": None,
        }
        for _index in range(parameter_count)
    ]
    fit = {
        "model": model,
        "theta": np.array([], dtype=float),
        "gain_rate": None,
        "loss_rate": None,
        "foreground_multiplier": None,
        "root_presence": None,
        "log_likelihood": None,
        "parameter_count": parameter_count,
        "aic": None,
        "converged": False,
        "fit_status": "not_fitted",
        "inference_status": reason,
        "optimizer_message": reason,
        "boundary": False,
        "intervals": intervals,
        "start_count": 0,
        "identifiable": False,
        "information_eigenvalues": np.array([], dtype=float),
        "information_condition": math.inf,
        "root_frequency": root_frequency,
        "profile_support_thetas": [],
        "contrast_profile": {
            "name": "NA", "low": None, "high": None, "status": reason,
            "raw_low": None, "raw_high": None, "theta_low": None, "theta_high": None,
        },
        "posterior_available": False,
        "posterior_unavailable_reason": reason,
        "posterior_sensitivity_status": "unavailable_no_posterior",
    }
    return fit


def _fit_output_row(
    family,
    layer,
    fit,
    observed_taxa,
    site_count,
    compressed_pattern_count,
    informative,
    root_frequency,
    ascertainment,
    analysis_summary=None,
):
    analysis_summary = analysis_summary or {}
    gain_ci = fit["intervals"][0]
    loss_ci = fit["intervals"][0] if fit["model"] == "ER" else fit["intervals"][1]
    multiplier_ci = (
        fit["intervals"][2]
        if fit["model"] == "ARD_FOREGROUND"
        else {"low": 1.0, "high": 1.0, "status": "fixed"}
    )
    root_ci = (
        fit["intervals"][-1]
        if root_frequency == "estimated"
        else {"low": fit["root_presence"], "high": fit["root_presence"], "status": root_frequency}
    )
    posterior_available = _fit_valid_for_posterior(fit)
    contrast = fit.get("contrast_profile", {})
    return {
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
        "posterior_available": str(posterior_available).lower(),
        "posterior_unavailable_reason": _posterior_unavailable_reason(fit),
        "posterior_sensitivity_status": fit.get(
            "posterior_sensitivity_status", _posterior_sensitivity_status(fit)
        ),
        "parameter_at_boundary": str(fit["boundary"]).lower(),
        "identifiable": str(fit["identifiable"]).lower(),
        "information_condition": _fmt(fit["information_condition"]),
        "optimizer_starts": fit["start_count"],
        "optimizer_message": fit["optimizer_message"],
        "ascertainment": ascertainment,
        "total_structural_sites": analysis_summary.get("total_structural_sites", site_count),
        "included_structural_sites": site_count,
        "variable_structural_sites": informative,
        "known_tip_count_min": analysis_summary.get("known_tip_count_min", "NA"),
        "known_tip_count_median": analysis_summary.get("known_tip_count_median", "NA"),
        "known_tip_count_max": analysis_summary.get("known_tip_count_max", "NA"),
        "known_tip_count_distribution": analysis_summary.get("known_tip_count_distribution", "NA"),
        "linked_group_count": analysis_summary.get("linked_group_count", 0),
        "correlated_linked_group_count": analysis_summary.get("correlated_linked_group_count", 0),
        "discovery_rules": analysis_summary.get("discovery_rules", "NA"),
        "observation_masks": analysis_summary.get("observation_masks", "NA"),
        "branch_length_mode": analysis_summary.get("branch_length_mode", "NA"),
        "matrix_schema_version": analysis_summary.get("matrix_schema_version", "NA"),
        "matrix_source": analysis_summary.get("matrix_source", "NA"),
        "contrast_name": contrast.get("name", "NA"),
        "contrast_ci_low": _fmt(contrast.get("low")),
        "contrast_ci_high": _fmt(contrast.get("high")),
        "contrast_profile_status": contrast.get("status", "not_applicable"),
    }


def _no_patterns_reason(encoded_sites):
    observed_by_site = [
        {value for value in observations.values() if value in {0, 1}}
        for _site_id, observations in encoded_sites
    ]
    if not encoded_sites:
        return "no_sites_in_analysis_universe"
    if not any(states for states in observed_by_site):
        return "no_observed_states"
    if not any(sum(value in {0, 1} for value in observations.values()) >= 2
               for _site_id, observations in encoded_sites):
        return "insufficient_observed_tips"
    if not any(len(states) > 1 for states in observed_by_site):
        return "no_observed_contrast"
    return "no_sites_selected_by_ascertainment"


def _lrt_unavailable_reason(
    informative,
    null_fit,
    alternative_fit,
    likelihood_difference,
    correlated_linked_sites=False,
):
    if informative == 0:
        return "no_observed_contrast"
    if correlated_linked_sites:
        return "correlated_linked_sites_not_modelled"
    if likelihood_difference < -1e-7:
        return "alternative_log_likelihood_below_null"
    null_reason = _posterior_unavailable_reason(null_fit)
    if null_reason != "NA":
        return f"null_model_{null_reason}"
    alternative_reason = _posterior_unavailable_reason(alternative_fit)
    if alternative_reason != "NA":
        return f"alternative_model_{alternative_reason}"
    contrast = alternative_fit.get("contrast_profile") or {}
    status = contrast.get("status", "missing")
    if status != "two_sided":
        return f"tested_contrast_profile_{status}"
    low = contrast.get("low")
    high = contrast.get("high")
    if alternative_fit.get("model") == "ARD":
        estimate = alternative_fit["gain_rate"] / alternative_fit["loss_rate"]
    else:
        estimate = alternative_fit.get("foreground_multiplier")
    if (
        low is None
        or high is None
        or estimate is None
        or not all(math.isfinite(float(value)) for value in (low, high, estimate))
        or not float(low) < float(estimate) < float(high)
        or contrast.get("theta_low") is None
        or contrast.get("theta_high") is None
    ):
        return "tested_contrast_profile_not_finite_internal"
    return "NA"


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
                return None, None, "profile_optimization_failed"
            return value, trial, "ok"

        def free_objective(free_values):
            trial = optimum.copy()
            trial[index] = fixed
            trial[free_indices] = free_values
            return objective(trial)

        start = optimum[free_indices]
        free_bounds = [bounds[idx] for idx in free_indices]
        result = minimize(free_objective, start, method="L-BFGS-B", bounds=free_bounds)
        if not result.success or not math.isfinite(float(result.fun)):
            return None, None, "profile_optimization_failed"
        trial = optimum.copy()
        trial[index] = fixed
        trial[free_indices] = result.x
        return -float(result.fun), trial, "ok"

    def profile_residual(fixed):
        profiled, _theta, status = profile_at(fixed)
        if status != "ok":
            raise ProfileOptimizationFailed
        return profiled - target

    def crossing(direction):
        edge = bounds[index][0] if direction < 0 else bounds[index][1]
        points = np.linspace(optimum[index], edge, 18)[1:]
        previous_x = optimum[index]
        previous_value = max_log_likelihood - target
        for point in points:
            profiled, endpoint_theta, status = profile_at(float(point))
            if status != "ok":
                return None, status, None
            value = profiled - target
            if value <= 0 <= previous_value:
                try:
                    root = brentq(profile_residual, float(point), float(previous_x))
                except ProfileOptimizationFailed:
                    return None, "profile_optimization_failed", None
                except (ValueError, RuntimeError):
                    return None, "profile_root_failed", None
                _value, root_theta, root_status = profile_at(root)
                if root_status != "ok":
                    return None, root_status, None
                return root, "closed", root_theta
            previous_x = float(point)
            previous_value = value
        return edge, "range_limited", endpoint_theta

    lower, lower_status, theta_low = crossing(-1)
    upper, upper_status, theta_high = crossing(1)
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
        "theta_low": theta_low,
        "theta_high": theta_high,
    }


def _profile_contrast_interval(objective, optimum, bounds, max_log_likelihood, weights, name):
    """Profile a declared linear contrast in optimization coordinates."""

    optimum = np.asarray(optimum, dtype=float)
    weights = np.asarray(weights, dtype=float)
    target = max_log_likelihood - PROFILE_DROP_95
    contrast_optimum = float(weights @ optimum)
    corners_low = np.array([bound[0] if weight >= 0 else bound[1] for weight, bound in zip(weights, bounds)])
    corners_high = np.array([bound[1] if weight >= 0 else bound[0] for weight, bound in zip(weights, bounds)])
    contrast_bounds = (float(weights @ corners_low), float(weights @ corners_high))

    class ProfileOptimizationFailed(Exception):
        pass

    def profile_at(fixed):
        constraint = {
            "type": "eq",
            "fun": lambda theta: float(weights @ np.asarray(theta, dtype=float) - fixed),
        }
        displacement = (fixed - contrast_optimum) * weights / float(weights @ weights)
        start = np.clip(optimum + displacement, [item[0] for item in bounds], [item[1] for item in bounds])
        result = minimize(objective, start, method="SLSQP", bounds=bounds, constraints=[constraint])
        if (
            not result.success
            or not math.isfinite(float(result.fun))
            or abs(float(weights @ np.asarray(result.x, dtype=float) - fixed)) > 1e-6
        ):
            return None, None, "profile_optimization_failed"
        return -float(result.fun), np.asarray(result.x, dtype=float), "ok"

    def residual(fixed):
        profiled, _theta, status = profile_at(fixed)
        if status != "ok":
            raise ProfileOptimizationFailed
        return profiled - target

    def crossing(direction):
        edge = contrast_bounds[0] if direction < 0 else contrast_bounds[1]
        points = np.linspace(contrast_optimum, edge, 18)[1:]
        previous_x = contrast_optimum
        previous_value = max_log_likelihood - target
        endpoint_theta = None
        for point in points:
            profiled, endpoint_theta, status = profile_at(float(point))
            if status != "ok":
                return None, status, None
            value = profiled - target
            if value <= 0 <= previous_value:
                try:
                    root = brentq(residual, float(point), float(previous_x))
                except ProfileOptimizationFailed:
                    return None, "profile_optimization_failed", None
                except (ValueError, RuntimeError):
                    return None, "profile_root_failed", None
                _value, root_theta, root_status = profile_at(root)
                if root_status != "ok":
                    return None, root_status, None
                return root, "closed", root_theta
            previous_x = float(point)
            previous_value = value
        return edge, "range_limited", endpoint_theta

    low, low_status, theta_low = crossing(-1)
    high, high_status, theta_high = crossing(1)
    failures = {low_status, high_status} & {"profile_optimization_failed", "profile_root_failed"}
    if failures:
        status = ";".join(sorted(failures))
    elif low_status == "range_limited" and high_status == "range_limited":
        status = "range_limited_both"
    elif low_status == "range_limited":
        status = "lower_range_limited"
    elif high_status == "range_limited":
        status = "upper_range_limited"
    else:
        status = "two_sided"
    return {
        "name": name,
        "low": math.exp(low) if low is not None else None,
        "high": math.exp(high) if high is not None else None,
        "status": status,
        "raw_low": low,
        "raw_high": high,
        "theta_low": theta_low,
        "theta_high": theta_high,
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
    has_observed_contrast = any(
        weight > 0 and len(_pattern_observed_values(pattern)) > 1
        for pattern, weight in _iter_weighted_patterns(patterns)
    )
    if not has_observed_contrast:
        return _unavailable_fit(model, root_frequency, root_presence, "no_observed_contrast")
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
        intervals = [
            {
                "low": None, "high": None, "status": "optimizer_failed", "raw_low": None,
                "raw_high": None, "theta_low": None, "theta_high": None,
            }
            for _idx in range(len(theta))
        ]
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
    contrast_profile = {
        "name": "NA", "low": None, "high": None, "status": "not_applicable",
        "raw_low": None, "raw_high": None, "theta_low": None, "theta_high": None,
    }
    if converged and model == "ARD":
        weights = [1.0, -1.0] + [0.0] * (len(theta) - 2)
        contrast_profile = _profile_contrast_interval(
            objective, theta, bounds, log_likelihood, weights, "gain_loss_rate_ratio"
        )
    elif converged and model == "ARD_FOREGROUND":
        weights = [0.0, 0.0, 1.0] + [0.0] * (len(theta) - 3)
        contrast_profile = _profile_contrast_interval(
            objective, theta, bounds, log_likelihood, weights, "foreground_multiplier"
        )
    profile_support_thetas = []
    for interval in intervals:
        if interval.get("status") == "two_sided":
            profile_support_thetas.extend(
                endpoint
                for endpoint in (interval.get("theta_low"), interval.get("theta_high"))
                if endpoint is not None and np.all(np.isfinite(endpoint))
            )
    if contrast_profile.get("status") == "two_sided":
        profile_support_thetas.extend(
            endpoint
            for endpoint in (
                contrast_profile.get("theta_low"), contrast_profile.get("theta_high")
            )
            if endpoint is not None and np.all(np.isfinite(endpoint))
        )
    fit = {
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
        "profile_support_thetas": profile_support_thetas,
        "contrast_profile": contrast_profile,
    }
    fit["posterior_available"] = _fit_valid_for_posterior(fit)
    fit["posterior_unavailable_reason"] = _posterior_unavailable_reason(fit)
    fit["posterior_sensitivity_status"] = _posterior_sensitivity_status(fit)
    return fit


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


def _analysis_summary(
    site_rows,
    all_encoded_sites,
    included_site_ids,
    branch_length_mode,
    matrix_schema_version,
    matrix_source,
):
    known_counts = sorted(
        sum(value in {0, 1} for value in observations.values())
        for _site_id, observations in all_encoded_sites
    )
    grouped_links = defaultdict(set)
    discovery_rules = set()
    observation_masks = set()
    for row in site_rows:
        discovery_rules.add(str(row.get("discovery_rule", "NA")))
        observation_masks.add(str(row.get("observation_mask", "NA")))
        linked_ids = {
            token
            for token in str(row.get("linked_group_id", "NA")).split(";")
            if token not in {"", "NA", "None", "none"}
        }
        for linked in linked_ids:
            grouped_links[linked].add(row.get("site_id", "NA"))
    included = set(included_site_ids)
    correlated = sum(1 for sites in grouped_links.values() if len(sites & included) > 1)
    if known_counts:
        midpoint = len(known_counts) // 2
        median = (
            known_counts[midpoint]
            if len(known_counts) % 2
            else 0.5 * (known_counts[midpoint - 1] + known_counts[midpoint])
        )
        distribution = ";".join(
            f"{count}:{known_counts.count(count)}" for count in sorted(set(known_counts))
        )
    else:
        median = None
        distribution = "NA"
    return {
        "total_structural_sites": len(all_encoded_sites),
        "known_tip_count_min": min(known_counts) if known_counts else "NA",
        "known_tip_count_median": _fmt(median),
        "known_tip_count_max": max(known_counts) if known_counts else "NA",
        "known_tip_count_distribution": distribution,
        "linked_group_count": len(grouped_links),
        "correlated_linked_group_count": correlated,
        "discovery_rules": ";".join(sorted(discovery_rules)) or "NA",
        "observation_masks": ";".join(sorted(observation_masks)) or "NA",
        "branch_length_mode": branch_length_mode,
        "matrix_schema_version": matrix_schema_version,
        "matrix_source": str(matrix_source),
    }


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
    structural_site_matrix_path=None,
    annotation_view="repertoire",
):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    site_rows, excluded_rows, matrix_source, matrix_mode = _load_site_matrix(
        input_dir,
        output_dir,
        structural_site_matrix_path,
        annotation_view,
    )
    matrix_annotation_view = structural_matrix_annotation_view(site_rows)
    effective_annotation_view = (
        annotation_view if matrix_annotation_view == "view_independent" else matrix_annotation_view
    )
    if matrix_annotation_view in {"canonical", "repertoire"} and matrix_annotation_view != annotation_view:
        raise SystemExit(
            "requested annotation view does not match the frozen structural-site matrix: "
            f"requested={annotation_view}, matrix={matrix_annotation_view}"
        )
    write_tsv(
        output_dir / "excluded_families.tsv",
        excluded_rows,
        ["family_id", "species", "copy_count", "gene_copy_ids", "reason"],
    )
    tree_rows = _validated_tree_rows(input_dir / "species_tree.tsv", branch_length_mode)
    tree = SpeciesTree(tree_rows)
    validate_structural_site_tip_rows(site_rows, tree.leaf_by_label)
    foreground_children = _read_foreground_children(foreground_branches, tree)
    if model == "foreground" and not foreground_children:
        raise SystemExit("--model foreground requires --foreground-branches")
    matrix_versions = sorted({row.get("schema_version", "NA") for row in site_rows})
    matrix_schema_version = ";".join(matrix_versions) if matrix_versions else STRUCTURAL_SITE_SCHEMA_VERSION
    if ascertainment == "complete-universe":
        _validate_complete_universe(site_rows, tree)

    grouped = defaultdict(lambda: defaultdict(dict))
    state_labels = {}
    site_kinds = {}
    for row in site_rows:
        grouped[(row["family_id"], row["layer"])][row["site_id"]][row["species"]] = row
        site_key = _site_key(row)
        state_labels[site_key] = (row["state_0"], row["state_1"])
        site_kinds[site_key] = row.get("site_kind", "NA")

    fit_rows = []
    test_rows = []
    node_rows = []
    branch_rows = []
    change_rows = []
    for (family, layer), sites in sorted(grouped.items()):
        raw_patterns = []
        all_encoded_sites = []
        encoded_sites = []
        included_site_ids = []
        family_layer_rows = [
            row for row in site_rows
            if row.get("family_id") == family and row.get("layer") == layer
        ]
        for site_id, observation_rows in sorted(sites.items()):
            state_0, state_1 = state_labels[(family, layer, site_id)]
            encoded = {
                label: (
                    0 if _masked_state(observation_rows[label]) == state_0
                    else 1 if _masked_state(observation_rows[label]) == state_1
                    else "unknown"
                )
                for label in tree.leaf_by_label
            }
            all_encoded_sites.append((site_id, encoded))
            observed_count = sum(value in {0, 1} for value in encoded.values())
            if observed_count >= 2 and _selected_for_ascertainment(encoded, ascertainment):
                raw_patterns.append(encoded)
                encoded_sites.append((site_id, encoded))
                included_site_ids.append(site_id)
        patterns = _compress_patterns(raw_patterns, sorted(tree.leaf_by_label))
        analysis_summary = _analysis_summary(
            family_layer_rows,
            all_encoded_sites,
            included_site_ids,
            branch_length_mode,
            matrix_schema_version,
            matrix_source,
        )
        null_model, alternative_model, test_id = _comparison_models(model)
        if patterns:
            site_count = sum(weight for _pattern, weight in patterns)
            compressed_pattern_count = len(patterns)
            informative = sum(
                weight
                for pattern, weight in patterns
                if len({value for value in pattern.values() if value in {0, 1}}) > 1
            )
            observed_taxa = len(
                {
                    species
                    for pattern, _weight in patterns
                    for species, value in pattern.items()
                    if value in {0, 1}
                }
            )
            null_fit = fit_model(
                tree,
                patterns,
                null_model,
                ascertainment=ascertainment,
                threads=threads,
                root_frequency=root_frequency,
                root_presence=root_presence,
            )
            extra_starts = []
            if null_fit["converged"]:
                null_theta = null_fit["theta"]
                if model == "foreground":
                    extra_starts.append([null_theta[0], null_theta[1], 0.0, *null_theta[2:]])
                else:
                    extra_starts.append([null_theta[0], null_theta[0], *null_theta[1:]])
            alternative_fit = fit_model(
                tree,
                patterns,
                alternative_model,
                foreground_children=foreground_children,
                ascertainment=ascertainment,
                threads=threads,
                root_frequency=root_frequency,
                root_presence=root_presence,
                extra_starts=extra_starts,
            )
            analysis_unavailable_reason = "NA" if informative > 0 else "no_observed_contrast"
        else:
            analysis_unavailable_reason = _no_patterns_reason(all_encoded_sites)
            null_fit = _unavailable_fit(
                null_model, root_frequency, root_presence, analysis_unavailable_reason
            )
            alternative_fit = _unavailable_fit(
                alternative_model, root_frequency, root_presence, analysis_unavailable_reason
            )
            site_count = 0
            compressed_pattern_count = 0
            informative = 0
            observed_taxa = len(
                {
                    species
                    for _site_id, observations in all_encoded_sites
                    for species, value in observations.items()
                    if value in {0, 1}
                }
            )
            encoded_sites = all_encoded_sites

        for fit in (null_fit, alternative_fit):
            fit_rows.append(
                _fit_output_row(
                    family,
                    layer,
                    fit,
                    observed_taxa,
                    site_count,
                    compressed_pattern_count,
                    informative,
                    root_frequency,
                    ascertainment,
                    analysis_summary,
                )
            )

        if (
            patterns
            and null_fit.get("log_likelihood") is not None
            and alternative_fit.get("log_likelihood") is not None
        ):
            likelihood_difference = alternative_fit["log_likelihood"] - null_fit["log_likelihood"]
            lrt_unavailable_reason = _lrt_unavailable_reason(
                informative,
                null_fit,
                alternative_fit,
                likelihood_difference,
                correlated_linked_sites=(
                    analysis_summary["correlated_linked_group_count"] > 0
                ),
            )
            lrt = 2.0 * max(0.0, likelihood_difference) if likelihood_difference >= -1e-7 else None
        elif patterns:
            likelihood_difference = None
            lrt = None
            lrt_unavailable_reason = (
                "no_observed_contrast"
                if informative == 0
                else f"alternative_model_{_posterior_unavailable_reason(alternative_fit)}"
            )
        else:
            likelihood_difference = None
            lrt = None
            lrt_unavailable_reason = analysis_unavailable_reason
        lrt_df = (
            int(alternative_fit["parameter_count"])
            - int(null_fit["parameter_count"])
        )
        if lrt_unavailable_reason == "NA" and lrt_df <= 0:
            lrt_unavailable_reason = "nonpositive_model_dimension_difference"
        lrt_available = lrt_unavailable_reason == "NA" and lrt is not None
        p_value = float(chi2.sf(lrt, lrt_df)) if lrt_available else None
        if lrt_available:
            test_status = "tested"
        elif lrt_unavailable_reason == "correlated_linked_sites_not_modelled":
            test_status = "not_tested_correlated_sites"
        elif lrt_unavailable_reason in {"no_observed_states", "no_observed_contrast"}:
            test_status = "parameters_not_estimable"
        elif lrt_unavailable_reason == "alternative_log_likelihood_below_null":
            test_status = "optimization_failure_alternative_below_null"
        else:
            test_status = "parameters_not_estimable"
        eligible_fits = [
            fit for fit in (null_fit, alternative_fit)
            if informative > 0 and _fit_valid_for_posterior(fit)
        ]
        selected_fit = min(eligible_fits, key=lambda fit: fit["aic"]) if eligible_fits else None
        posterior_available = selected_fit is not None
        if posterior_available:
            posterior_unavailable_reason = "NA"
        elif analysis_unavailable_reason != "NA":
            posterior_unavailable_reason = analysis_unavailable_reason
        else:
            posterior_unavailable_reason = (
                f"no_valid_fit;null={_posterior_unavailable_reason(null_fit)};"
                f"alternative={_posterior_unavailable_reason(alternative_fit)}"
            )
        test_row = {
            "family_id": family,
            "layer": layer,
            "test_id": test_id,
            "null_model": null_fit["model"],
            "alternative_model": alternative_fit["model"],
            "null_log_likelihood": _fmt(null_fit["log_likelihood"]),
            "alternative_log_likelihood": _fmt(alternative_fit["log_likelihood"]),
            "lrt_statistic": _fmt(lrt),
            "df": lrt_df,
            "p_value": _fmt(p_value),
            "q_value": "NA",
            "q_value_method": "not_available",
            "test_status": test_status,
            "inference_status": (
                lrt_unavailable_reason
                if lrt_unavailable_reason in {"no_observed_states", "no_observed_contrast"}
                else test_status
            ),
            "lrt_available": str(lrt_available).lower(),
            "lrt_unavailable_reason": lrt_unavailable_reason,
            "reference_distribution": "asymptotic_chi_square",
            "small_sample_accuracy": "unassessed",
            "posterior_available": str(posterior_available).lower(),
            "posterior_model": selected_fit["model"] if selected_fit is not None else "NA",
            "posterior_unavailable_reason": posterior_unavailable_reason,
            "n_taxa": observed_taxa,
            "n_structural_sites": site_count,
            "n_compressed_patterns": compressed_pattern_count,
            "n_informative_patterns": informative,
            "tested_contrast": alternative_fit.get("contrast_profile", {}).get("name", "NA"),
            "tested_contrast_estimate": _fmt(
                alternative_fit["gain_rate"] / alternative_fit["loss_rate"]
                if alternative_fit.get("model") == "ARD"
                and alternative_fit.get("gain_rate") is not None
                and alternative_fit.get("loss_rate") not in {None, 0}
                else alternative_fit.get("foreground_multiplier")
            ),
            "tested_contrast_ci_low": _fmt(
                alternative_fit.get("contrast_profile", {}).get("low")
            ),
            "tested_contrast_ci_high": _fmt(
                alternative_fit.get("contrast_profile", {}).get("high")
            ),
            "tested_contrast_profile_status": alternative_fit.get(
                "contrast_profile", {}
            ).get("status", "unavailable"),
            **analysis_summary,
            "included_structural_sites": site_count,
            "variable_structural_sites": informative,
        }
        test_rows.append(test_row)
        if selected_fit is None:
            diagnostic_fit = alternative_fit if informative else null_fit
            status = diagnostic_fit.get("fit_status", "not_estimable")
            for site_id, observations in encoded_sites:
                state_0, state_1 = state_labels[(family, layer, site_id)]
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
                        "structural_relation": "NA",
                        "structural_pattern": f"{state_0}<->{state_1}",
                        "endpoint_change_probability": "NA",
                        "gain_endpoint_probability": "NA",
                        "loss_endpoint_probability": "NA",
                        "direction_probability": "NA",
                        "expected_gain_count": "NA",
                        "expected_loss_count": "NA",
                        "model": diagnostic_fit.get("model", "NA"),
                        "rate_test_status": test_status,
                        "posterior_available": "false",
                        "posterior_unavailable_reason": posterior_unavailable_reason,
                        "conditioning": (
                            f"posterior_not_reported;fit_status={status};"
                            f"inference_status={test_row['inference_status']};"
                            f"known_tip_count={known};requires=successful_identifiable_interior_fit"
                        ),
                    }
                )
            continue
        sensitivity_status = selected_fit.get(
            "posterior_sensitivity_status", _posterior_sensitivity_status(selected_fit)
        )
        posterior_conditioning = (
            "conditional_MLE;"
            f"fit_status={selected_fit['fit_status']};"
            f"uncertainty={sensitivity_status}"
        )

        for site_id, observations in encoded_sites:
            state_0, state_1 = state_labels[(family, layer, site_id)]
            site_kind = site_kinds[(family, layer, site_id)]
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
                        "structural_relation": "NA",
                        "structural_pattern": f"{state_0}<->{state_1}",
                        "endpoint_change_probability": "NA",
                        "gain_endpoint_probability": "NA",
                        "loss_endpoint_probability": "NA",
                        "direction_probability": "NA",
                        "expected_gain_count": "NA",
                        "expected_loss_count": "NA",
                        "model": selected_fit["model"],
                        "rate_test_status": test_status,
                        "posterior_available": "false",
                        "posterior_unavailable_reason": "site_likelihood_zero",
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
            node_envelope, edge_envelope, site_sensitivity_status = _profile_posterior_envelope(
                tree,
                observations,
                selected_fit,
                foreground_children,
                node_posterior,
                edge_posterior,
            )
            for node_id, probabilities in node_posterior.items():
                for index, probability in enumerate(probabilities):
                    profile_low = (
                        node_envelope[node_id][0][index] if node_envelope is not None else None
                    )
                    profile_high = (
                        node_envelope[node_id][1][index] if node_envelope is not None else None
                    )
                    node_rows.append(
                        {
                            "family_id": family,
                            "layer": layer,
                            "site_id": site_id,
                            "node_id": node_id,
                            "node_label": tree.label[node_id],
                            "state": state_0 if index == 0 else state_1,
                            "posterior_probability": _fmt(probability),
                            "profile_probability_low": _fmt(profile_low),
                            "profile_probability_high": _fmt(profile_high),
                            "uncertainty_status": site_sensitivity_status,
                            "posterior_available": "true",
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
                envelope = edge_envelope.get((parent, child)) if edge_envelope is not None else None
                for src_index, dst_index, src, dst, probability in (
                    (0, 1, state_0, state_1, gain_probability),
                    (1, 0, state_1, state_0, loss_probability),
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
                            "event_type": _transition_event_type(layer, src_index, dst_index),
                            "structural_relation": _transition_structural_relation(
                                layer,
                                site_kind,
                                src_index,
                                dst_index,
                            ),
                            "endpoint_transition_probability": _fmt(probability),
                            "profile_transition_probability_low": _fmt(
                                envelope["low"][src_index, dst_index] if envelope is not None else None
                            ),
                            "profile_transition_probability_high": _fmt(
                                envelope["high"][src_index, dst_index] if envelope is not None else None
                            ),
                            "total_endpoint_change_probability": _fmt(change_probability),
                            "profile_total_change_probability_low": _fmt(
                                envelope["total_low"] if envelope is not None else None
                            ),
                            "profile_total_change_probability_high": _fmt(
                                envelope["total_high"] if envelope is not None else None
                            ),
                            "expected_gain_count": _fmt(expected_gain),
                            "expected_loss_count": _fmt(expected_loss),
                            "uncertainty_status": site_sensitivity_status,
                            "posterior_available": "true",
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
                        "structural_relation": (
                            "exon_split_or_exon_fusion"
                            if layer == "splice_junction" and site_kind in {
                                "within_element_junction",
                                "within_exon_boundary",
                            }
                            else "NA"
                        ),
                        "structural_pattern": f"{state_0}<->{state_1}",
                        "endpoint_change_probability": _fmt(change_probability),
                        "gain_endpoint_probability": _fmt(gain_probability),
                        "loss_endpoint_probability": _fmt(loss_probability),
                        "direction_probability": "NA",
                        "expected_gain_count": _fmt(expected_gain),
                        "expected_loss_count": _fmt(expected_loss),
                        "model": selected_fit["model"],
                        "rate_test_status": test_status,
                        "posterior_available": "true",
                        "posterior_unavailable_reason": "NA",
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
            "parameter_count", "aic", "converged", "fit_status", "posterior_available",
            "posterior_unavailable_reason", "posterior_sensitivity_status", "parameter_at_boundary",
            "identifiable", "information_condition", "optimizer_starts", "optimizer_message",
            "ascertainment", "inference_status", "total_structural_sites",
            "included_structural_sites", "variable_structural_sites", "known_tip_count_min",
            "known_tip_count_median", "known_tip_count_max", "known_tip_count_distribution",
            "linked_group_count", "correlated_linked_group_count", "discovery_rules",
            "observation_masks", "branch_length_mode", "matrix_schema_version", "matrix_source",
            "contrast_name", "contrast_ci_low", "contrast_ci_high", "contrast_profile_status",
        ],
    )
    write_tsv(
        output_dir / "model_tests.tsv",
        test_rows,
        [
            "family_id", "layer", "test_id", "null_model", "alternative_model", "null_log_likelihood",
            "alternative_log_likelihood", "lrt_statistic", "df", "p_value", "q_value", "q_value_method",
            "test_status", "lrt_available", "lrt_unavailable_reason", "reference_distribution",
            "small_sample_accuracy", "posterior_available", "posterior_model",
            "posterior_unavailable_reason", "n_taxa", "n_structural_sites", "n_compressed_patterns",
            "n_informative_patterns", "inference_status", "tested_contrast",
            "tested_contrast_estimate", "tested_contrast_ci_low", "tested_contrast_ci_high",
            "tested_contrast_profile_status", "total_structural_sites",
            "included_structural_sites", "variable_structural_sites", "known_tip_count_min",
            "known_tip_count_median", "known_tip_count_max", "known_tip_count_distribution",
            "linked_group_count", "correlated_linked_group_count", "discovery_rules",
            "observation_masks", "branch_length_mode", "matrix_schema_version", "matrix_source",
        ],
    )
    write_tsv(
        output_dir / "node_state_posteriors.tsv",
        node_rows,
        ["family_id", "layer", "site_id", "node_id", "node_label", "state", "posterior_probability",
         "profile_probability_low", "profile_probability_high", "uncertainty_status",
         "posterior_available", "model", "conditioning"],
    )
    write_tsv(
        output_dir / "branch_transition_posteriors.tsv",
        branch_rows,
        [
            "family_id", "layer", "site_id", "parent_node", "child_node", "parent_label", "child_label",
            "branch_length", "from_state", "to_state", "event_type", "structural_relation",
            "endpoint_transition_probability",
            "profile_transition_probability_low", "profile_transition_probability_high",
            "total_endpoint_change_probability", "profile_total_change_probability_low",
            "profile_total_change_probability_high", "expected_gain_count", "expected_loss_count",
            "uncertainty_status", "posterior_available", "model", "conditioning",
        ],
    )
    write_tsv(
        output_dir / "structural_changes.tsv",
        change_rows,
        [
            "family_id", "layer", "site_id", "parent_node", "child_node", "branch_scope",
            "structural_change_type", "structural_relation", "structural_pattern",
            "endpoint_change_probability",
            "gain_endpoint_probability", "loss_endpoint_probability", "direction_probability",
            "expected_gain_count", "expected_loss_count", "model", "rate_test_status",
            "posterior_available", "posterior_unavailable_reason", "conditioning",
        ],
    )
    write_tsv(
        output_dir / "phylogeny_scope.tsv",
        [{
            "scope": "single_copy_structural_sites",
            "tree_file": "species_tree.tsv",
            "tree_scope": "species_tree",
            "layers": "exon_presence;exon_role;splice_junction",
            "matrix_source": str(matrix_source),
            "matrix_mode": matrix_mode,
            "matrix_schema_version": matrix_schema_version,
            "annotation_view": effective_annotation_view,
            "note": "All homologous structural sites within each layer share fitted gain/loss parameters.",
        }],
        [
            "scope", "tree_file", "tree_scope", "layers", "matrix_source", "matrix_mode",
            "matrix_schema_version", "annotation_view", "note",
        ],
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
        "structural_site_matrix": str(matrix_source),
        "structural_site_matrix_mode": matrix_mode,
        "structural_site_matrix_schema_version": matrix_schema_version,
        "annotation_view": effective_annotation_view,
        "threads": max(1, int(threads)),
        "fixed_inputs": ["species_tree", "ortholog_set", "structural_site_states"],
        "estimated_parameters": ["gain_rate", "loss_rate"] + (["root_presence"] if root_frequency == "estimated" else []) + (["foreground_multiplier"] if model == "foreground" else []),
        "excluded_family_count": len({row["family_id"] for row in excluded_rows}),
    }
    (output_dir / "run_parameters.json").write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")
    return site_rows, fit_rows, test_rows, change_rows
