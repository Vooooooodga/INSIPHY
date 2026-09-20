"""inference / posterior export: explicit implementation ownership."""
from __future__ import annotations

from intraphy.inference.ctmc import _pattern_log_likelihood
from intraphy.inference.diagnostics import _posterior_sensitivity_status
from intraphy.inference.events import _transition_event_type
from intraphy.inference.events import _transition_structural_relation
from intraphy.inference.formatting import _fmt
from intraphy.inference.posterior import _expected_transition_count
from intraphy.inference.posterior import _posterior_messages
from intraphy.inference.posterior import _profile_posterior_envelope
import math


def _append_structural_posteriors(selected_fit, encoded_sites, state_labels, family, layer, site_kinds, tree, foreground_children, change_rows, test_status, node_rows, branch_rows):
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
