"""inference / fitting: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from intraphy.inference.ctmc import _dataset_log_likelihood
from intraphy.inference.ctmc import _iter_weighted_patterns
from intraphy.inference.ctmc import _pattern_observed_values
from intraphy.inference.diagnostics import _fit_valid_for_posterior
from intraphy.inference.diagnostics import _posterior_sensitivity_status
from intraphy.inference.diagnostics import _posterior_unavailable_reason
from intraphy.inference.diagnostics import _unavailable_fit
from intraphy.inference.parameters import PROFILE_DROP_95
from intraphy.inference.parameters import _decode_parameters
from intraphy.inference.parameters import _model_bounds
from intraphy.inference.parameters import _model_starts
from scipy.optimize import brentq
from scipy.optimize import minimize
import math
import numpy as np


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
