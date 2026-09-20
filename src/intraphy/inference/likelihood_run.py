"""inference / likelihood run: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.inference.diagnostics import _fit_output_row
from intraphy.inference.formatting import _bh_adjust
from intraphy.inference.inputs import _load_site_matrix
from intraphy.inference.inputs import _read_foreground_children
from intraphy.inference.inputs import _site_key
from intraphy.inference.inputs import _validate_complete_universe
from intraphy.inference.inputs import _validated_tree_rows
from intraphy.inference.layer_fitting import _fit_structural_layer
from intraphy.inference.model_comparison import _compare_structural_models
from intraphy.inference.posterior_export import _append_structural_posteriors
from intraphy.observations.schema import STRUCTURAL_SITE_SCHEMA_VERSION
from intraphy.observations.schema import validate_structural_site_tip_rows
from intraphy.run_result import record_run_result
from intraphy.storage.tabular import write_tsv
from intraphy.structural_sites import structural_matrix_annotation_view
from intraphy.topology import SpeciesTree
from pathlib import Path
import json


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
    *, observation_matrix=None, analysis_range="all", min_callable_fraction=0.70,
):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    if observation_matrix is None:
        site_rows, excluded_rows, matrix_source, matrix_mode = _load_site_matrix(
            input_dir, output_dir, structural_site_matrix_path, annotation_view,
            analysis_range=analysis_range, min_callable_fraction=min_callable_fraction,
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
        null_fit, alternative_fit, observed_taxa, site_count, compressed_pattern_count, informative, analysis_summary, dependent, patterns, analysis_unavailable_reason, test_id, encoded_sites = _fit_structural_layer(site_rows, family, layer, sites, state_labels, tree, ascertainment, branch_length_mode, matrix_schema_version, matrix_source, model, threads, root_frequency, root_presence, foreground_children)

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

        test_row, selected_fit, test_status, posterior_unavailable_reason = _compare_structural_models(dependent, patterns, null_fit, alternative_fit, informative, analysis_summary, analysis_unavailable_reason, family, layer, test_id, observed_taxa, site_count, compressed_pattern_count)
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
        _append_structural_posteriors(selected_fit, encoded_sites, state_labels, family, layer, site_kinds, tree, foreground_children, change_rows, test_status, node_rows, branch_rows)

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
