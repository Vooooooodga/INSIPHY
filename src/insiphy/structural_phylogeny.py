"""structural_phylogeny: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.inference.ctmc import _compress_patterns
from insiphy.inference.ctmc import _pattern_log_likelihood
from insiphy.inference.ctmc import _selected_for_ascertainment
from insiphy.inference.diagnostics import _fit_output_row
from insiphy.inference.diagnostics import _fit_valid_for_posterior
from insiphy.inference.diagnostics import _lrt_unavailable_reason
from insiphy.inference.diagnostics import _posterior_sensitivity_status
from insiphy.inference.diagnostics import _posterior_unavailable_reason
from insiphy.inference.diagnostics import _unavailable_fit
from insiphy.inference.events import _transition_event_type
from insiphy.inference.events import _transition_structural_relation
from insiphy.inference.fitting import _no_patterns_reason
from insiphy.inference.fitting import fit_model
from insiphy.inference.formatting import _bh_adjust
from insiphy.inference.formatting import _fmt
from insiphy.inference.inputs import _analysis_summary
from insiphy.inference.inputs import _load_site_matrix
from insiphy.inference.inputs import _masked_state
from insiphy.inference.inputs import _read_foreground_children
from insiphy.inference.inputs import _site_key
from insiphy.inference.inputs import _validate_complete_universe
from insiphy.inference.inputs import _validated_tree_rows
from insiphy.inference.parameters import _comparison_models
from insiphy.inference.posterior import _expected_transition_count
from insiphy.inference.posterior import _posterior_messages
from insiphy.inference.posterior import _profile_posterior_envelope
from insiphy.observations.schema import STRUCTURAL_SITE_SCHEMA_VERSION
from insiphy.observations.schema import validate_structural_site_tip_rows
from insiphy.run_result import record_run_result
from insiphy.storage.tabular import write_tsv
from insiphy.structural_sites import structural_matrix_annotation_view
from insiphy.topology import SpeciesTree
from pathlib import Path
from scipy.stats import chi2
import json
import math


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
    *, observation_matrix=None,
):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    if observation_matrix is None:
        site_rows, excluded_rows, matrix_source, matrix_mode = _load_site_matrix(
            input_dir, output_dir, structural_site_matrix_path, annotation_view,
        )
    else:
        site_rows, excluded_rows, matrix_source, matrix_mode = observation_matrix.as_legacy_tuple()
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
    record_run_result(output_dir, input_dir, model, "single-copy")
    return site_rows, fit_rows, test_rows, change_rows


# Backward-compatible symbol exports; no alternate implementations.
from insiphy.inference.ctmc import (
    _transition_matrix,
    _root_prior,
    _iter_weighted_patterns,
    _compress_patterns,
    _pattern_observed_values,
    _pattern_has_state_one,
    _selected_for_ascertainment,
    _log_array,
    _inside_log_messages,
    _inside_messages,
    _pattern_log_likelihood,
    _ascertainment_log_probability,
    _dataset_log_likelihood,
)
from insiphy.inference.diagnostics import (
    _posterior_unavailable_reason,
    _fit_valid_for_posterior,
    _posterior_sensitivity_status,
    _unavailable_fit,
    _fit_output_row,
    _lrt_unavailable_reason,
)
from insiphy.inference.events import (
    _transition_event_type,
    _transition_structural_relation,
)
from insiphy.inference.fitting import (
    _no_patterns_reason,
    _profile_interval,
    _profile_contrast_interval,
    _observed_information,
    fit_model,
)
from insiphy.inference.formatting import (
    _fmt,
    _bh_adjust,
)
from insiphy.inference.inputs import (
    _validated_tree_rows,
    _site_key,
    _load_site_matrix,
    _masked_state,
    _validate_complete_universe,
    _read_foreground_children,
    _analysis_summary,
)
from insiphy.inference.parameters import (
    RATE_MIN,
    RATE_MAX,
    MULTIPLIER_MIN,
    MULTIPLIER_MAX,
    PROFILE_DROP_95,
    COMPLETE_UNIVERSE_RULES,
    GENERATED_ALL_ZERO_MASKS,
    _decode_parameters,
    _model_bounds,
    _comparison_models,
    _model_starts,
)
from insiphy.inference.posterior import (
    _posterior_messages,
    _conditional_transition_count,
    _expected_transition_count,
    _profile_posterior_envelope,
)
import json
import math
import numpy as np

from insiphy.structural_sites import build_structural_site_matrix
