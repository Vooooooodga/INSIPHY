"""inference / model comparison: explicit implementation ownership."""
from __future__ import annotations

from intraphy.inference.diagnostics import _fit_valid_for_posterior
from intraphy.inference.diagnostics import _lrt_unavailable_reason
from intraphy.inference.diagnostics import _posterior_unavailable_reason
from intraphy.inference.formatting import _fmt
from scipy.stats import chi2


def _compare_structural_models(dependent, patterns, null_fit, alternative_fit, informative, analysis_summary, analysis_unavailable_reason, family, layer, test_id, observed_taxa, site_count, compressed_pattern_count):
    if dependent and patterns:
        likelihood_difference = lrt = None
        lrt_unavailable_reason = "correlated_linked_sites_not_modelled"
    elif (
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
    return test_row, selected_fit, test_status, posterior_unavailable_reason
