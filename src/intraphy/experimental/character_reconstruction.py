"""experimental / character reconstruction: explicit implementation ownership."""
from __future__ import annotations

from intraphy.experimental.markov import discrete_log_likelihood
from intraphy.experimental.markov import fit_discrete_ctmc
from intraphy.experimental.markov import fit_foreground_rate_test
from intraphy.experimental.markov import fit_invariant_test
from intraphy.experimental.markov import sankoff
from intraphy.experimental.stochastic_history import bootstrap_invariant_test
from intraphy.experimental.stochastic_history import ctmc_posteriors
from intraphy.experimental.stochastic_history import stochastic_map_summary
from intraphy.io import norm_state


def add_character(
    tree,
    layer,
    object_id,
    tips,
    states,
    state_rows,
    branch_rows,
    model_score_rows,
    model_fit_rows,
    hypothesis_rows,
    bootstrap_rows,
    stochastic_rows,
    foreground_rows,
    bootstrap_replicates=0,
    stochastic_maps=0,
    seed=7,
    foreground_edges=None,
):
    score, probs, edges = sankoff(tree, tips, states, layer)
    observed_tip_count = sum(1 for value in tips.values() if norm_state(value) != "unknown")
    if observed_tip_count < 2:
        model_score_rows.append(
            {
                "layer": layer,
                "object_id": object_id,
                "model": "insufficient_observed_tips",
                "parsimony_score": f"{score:.6g}",
                "log_likelihood": "NA",
                "state_count": len(states),
                "observed_tip_count": observed_tip_count,
                "fitted_rate": "NA",
                "aic": "NA",
                "bic": "NA",
            }
        )
        hypothesis_rows.append(
            {
                "layer": layer,
                "object_id": object_id,
                "test_id": "ctmc_vs_invariant",
                "null_model": "invariant_no_structural_change",
                "alternative_model": "ctmc_mk_branch_length",
                "null_log_likelihood": "NA",
                "alternative_log_likelihood": "NA",
                "lrt_statistic": "NA",
                "df": 1,
                "p_value": "NA",
                "p_value_method": "insufficient_observed_tips",
                "fitted_rate": "NA",
                "null_aic": "NA",
                "alternative_aic": "NA",
                "null_bic": "NA",
                "alternative_bic": "NA",
                "observed_tip_count": observed_tip_count,
                "tip_error": "NA",
            }
        )
        model_fit_rows.append(
            {
                "layer": layer,
                "object_id": object_id,
                "model": "insufficient_observed_tips",
                "fitted_rate": "NA",
                "log_likelihood": "NA",
                "aic": "NA",
                "bic": "NA",
                "observed_tip_count": observed_tip_count,
            }
        )
        for row in probs:
            state_rows.append({"layer": layer, "object_id": object_id, "score": f"{score:.6g}", **row})
        return
    log_likelihood = discrete_log_likelihood(tree, tips, states, layer)
    fit = fit_discrete_ctmc(tree, tips, states, layer)
    test = fit_invariant_test(tree, tips, states, layer)
    _ctmc_nodes, ctmc_edges = ctmc_posteriors(tree, tips, states, layer, rate=fit["rate"])
    ctmc_by_edge = {(row["parent_node"], row["child_node"]): row for row in ctmc_edges}
    bootstrap = bootstrap_invariant_test(tree, tips, states, layer, replicates=bootstrap_replicates, seed=seed, tip_error=1e-6)
    if bootstrap:
        bootstrap_rows.append(
            {
                "layer": layer,
                "object_id": object_id,
                "test_id": bootstrap["test_id"],
                "observed_lrt": f"{bootstrap['observed_lrt']:.6g}",
                "bootstrap_replicates": bootstrap["bootstrap_replicates"],
                "empirical_p_value": f"{bootstrap['empirical_p_value']:.6g}",
                "monte_carlo_se": f"{bootstrap['monte_carlo_se']:.6g}",
                "null_lrt_mean": f"{bootstrap['null_lrt_mean']:.6g}",
                "null_lrt_q025": f"{bootstrap['null_lrt_q025']:.6g}",
                "null_lrt_q500": f"{bootstrap['null_lrt_q500']:.6g}",
                "null_lrt_q975": f"{bootstrap['null_lrt_q975']:.6g}",
                "seed": bootstrap["seed"],
                "tip_error": f"{bootstrap['tip_error']:.6g}",
            }
        )
    for row in stochastic_map_summary(tree, tips, states, layer, rate=fit["rate"], replicates=stochastic_maps, seed=seed, tip_error=1e-6):
        stochastic_rows.append(
            {
                "layer": layer,
                "object_id": object_id,
                "parent_node": row["parent_node"],
                "child_node": row["child_node"],
                "parent_label": tree.label[row["parent_node"]],
                "child_label": tree.label[row["child_node"]],
                "map_sample_count": row["map_sample_count"],
                "posterior_pr_any_change": f"{row['posterior_pr_any_change']:.6g}",
                "posterior_expected_change_count": f"{row['posterior_expected_change_count']:.6g}",
                "posterior_change_count_low": f"{row['posterior_change_count_low']:.6g}",
                "posterior_change_count_high": f"{row['posterior_change_count_high']:.6g}",
                "posterior_most_frequent_transition": row["posterior_most_frequent_transition"],
                "posterior_transition_probability": f"{row['posterior_transition_probability']:.6g}",
            }
        )
    foreground = fit_foreground_rate_test(tree, tips, states, layer, foreground_edges, tip_error=1e-6) if foreground_edges else None
    if foreground:
        foreground_rows.append(
            {
                "layer": layer,
                "object_id": object_id,
                "test_id": "foreground_background_rate",
                "null_model": foreground["null_model"],
                "alternative_model": foreground["alternative_model"],
                "null_log_likelihood": f"{foreground['null_log_likelihood']:.6g}",
                "alternative_log_likelihood": f"{foreground['alternative_log_likelihood']:.6g}",
                "lrt_statistic": f"{foreground['lrt_statistic']:.6g}",
                "df": foreground["df"],
                "p_value": f"{foreground['p_value']:.6g}",
                "p_value_method": foreground["p_value_method"],
                "background_rate": f"{foreground['background_rate']:.6g}",
                "foreground_rate": f"{foreground['foreground_rate']:.6g}",
                "rate_ratio": f"{foreground['rate_ratio']:.6g}",
                "null_aic": f"{foreground['null_aic']:.6g}",
                "alternative_aic": f"{foreground['alternative_aic']:.6g}",
                "null_bic": f"{foreground['null_bic']:.6g}",
                "alternative_bic": f"{foreground['alternative_bic']:.6g}",
                "observed_tip_count": foreground["observed_tip_count"],
            }
        )
    model_score_rows.append(
        {
            "layer": layer,
            "object_id": object_id,
            "model": "likelihood_like_discrete_character",
            "parsimony_score": f"{score:.6g}",
            "log_likelihood": f"{log_likelihood:.6g}",
            "state_count": len(states),
            "observed_tip_count": sum(1 for value in tips.values() if norm_state(value) != "unknown"),
            "fitted_rate": f"{fit['rate']:.6g}",
            "aic": f"{fit['aic']:.6g}",
            "bic": f"{fit['bic']:.6g}",
        }
    )
    hypothesis_rows.append(
        {
            "layer": layer,
            "object_id": object_id,
            "test_id": "ctmc_vs_invariant",
            "null_model": test["null_model"],
            "alternative_model": test["alternative_model"],
            "null_log_likelihood": f"{test['null_log_likelihood']:.6g}",
            "alternative_log_likelihood": f"{test['alternative_log_likelihood']:.6g}",
            "lrt_statistic": f"{test['lrt_statistic']:.6g}",
            "df": test["df"],
            "p_value": f"{test['p_value']:.6g}",
            "p_value_method": test["p_value_method"],
            "fitted_rate": f"{test['fitted_rate']:.6g}",
            "null_aic": f"{test['null_aic']:.6g}",
            "alternative_aic": f"{test['alternative_aic']:.6g}",
            "null_bic": f"{test['null_bic']:.6g}",
            "alternative_bic": f"{test['alternative_bic']:.6g}",
            "observed_tip_count": test["observed_tip_count"],
            "tip_error": f"{test['tip_error']:.6g}",
        }
    )
    model_fit_rows.append(
        {
            "layer": layer,
            "object_id": object_id,
            "model": fit["model"],
            "fitted_rate": f"{fit['rate']:.6g}",
            "log_likelihood": f"{fit['log_likelihood']:.6g}",
            "aic": f"{fit['aic']:.6g}",
            "bic": f"{fit['bic']:.6g}",
            "observed_tip_count": fit["observed_tip_count"],
        }
    )
    for row in probs:
        state_rows.append({"layer": layer, "object_id": object_id, "score": f"{score:.6g}", **row})
    for row in edges:
        ctmc = ctmc_by_edge.get((row["parent_node"], row["child_node"]), {})
        event_type = f"{layer}_change" if row["status"] == "change_required" else "state_change"
        branch_rows.append(
            {
                "layer": layer,
                "object_id": object_id,
                "event_type": event_type,
                **row,
                "ctmc_change_probability": f"{ctmc.get('ctmc_change_probability', 0.0):.6g}",
                "ctmc_most_likely_change": ctmc.get("ctmc_most_likely_change", "NA"),
                "ctmc_most_likely_change_probability": f"{ctmc.get('ctmc_most_likely_change_probability', 0.0):.6g}",
            }
        )
