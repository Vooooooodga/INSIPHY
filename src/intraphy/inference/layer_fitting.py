"""inference / layer fitting: explicit implementation ownership."""
from __future__ import annotations

from intraphy.inference.ctmc import _compress_patterns
from intraphy.inference.ctmc import _selected_for_ascertainment
from intraphy.inference.diagnostics import _unavailable_fit
from intraphy.inference.fitting import _no_patterns_reason
from intraphy.inference.fitting import fit_model
from intraphy.inference.inputs import _analysis_summary
from intraphy.inference.inputs import _masked_state
from intraphy.inference.parameters import _comparison_models


def _fit_structural_layer(site_rows, family, layer, sites, state_labels, tree, ascertainment, branch_length_mode, matrix_schema_version, matrix_source, model, threads, root_frequency, root_presence, foreground_children):
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
    dependent = analysis_summary["correlated_linked_group_count"] > 0
    if patterns and not dependent:
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
    elif patterns:
        analysis_unavailable_reason = "correlated_linked_sites_not_modelled"
        null_fit = _unavailable_fit(null_model, root_frequency, root_presence, analysis_unavailable_reason)
        alternative_fit = _unavailable_fit(alternative_model, root_frequency, root_presence, analysis_unavailable_reason)
        site_count = len(raw_patterns)
        compressed_pattern_count = len(patterns)
        informative = sum(len({v for v in p.values() if v in {0, 1}}) > 1 for p in raw_patterns)
        observed_taxa = len({s for p in raw_patterns for s, v in p.items() if v in {0, 1}})
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
    return null_fit, alternative_fit, observed_taxa, site_count, compressed_pattern_count, informative, analysis_summary, dependent, patterns, analysis_unavailable_reason, test_id, encoded_sites
