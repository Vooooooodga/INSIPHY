"""inference / diagnostics: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.inference.formatting import _fmt
from intraphy.inference.parameters import _model_bounds
import math
import numpy as np


def _posterior_unavailable_reason(fit):
    if fit.get("inference_status") == "correlated_linked_sites_not_modelled":
        return "correlated_linked_sites_not_modelled"
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
