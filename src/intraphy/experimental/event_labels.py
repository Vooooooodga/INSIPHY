"""experimental / event labels: explicit implementation ownership."""
from __future__ import annotations

from intraphy.elements import element_role_from_occurrences
from intraphy.io import norm_state
from intraphy.io import to_float


EXONIC = {"CDS", "UTR", "noncoding_exon", "exon"}


NONCODING = {"intron", "regulatory", "intergenic", "noncoding"}


COPY_CONTEXT_PATTERNS = {
    "copy_duplication_or_expansion",
    "copy_duplication_or_relocation",
    "copy_loss_or_collapse",
    "copy_multiplicity_shift",
    "retrocopy_or_dispersed_duplication_candidate",
}


ANNOTATION_EVIDENCE_PATTERNS = {
    "sequence_supported_annotation_gap",
    "annotation_or_alignment_evidence",
}


AMBIGUOUS_EVIDENCE_PATTERNS = {
    "ambiguous_paralogous_similarity",
}


UNINFORMATIVE_SOURCE_LABELS = {"", "NA", "unknown", "unknown_source", "ancestral_source", "intronic_source", "te_source"}


def call_scope_for_pattern(pattern):
    if pattern in AMBIGUOUS_EVIDENCE_PATTERNS:
        return "ambiguous_evidence"
    if pattern in COPY_CONTEXT_PATTERNS:
        return "copy_context"
    if pattern in ANNOTATION_EVIDENCE_PATTERNS:
        return "annotation_evidence"
    return "core_structural_event"


def informative_source_labels(labels):
    out = set()
    for label in labels:
        text = str(label or "").strip()
        low = text.lower()
        if low in UNINFORMATIVE_SOURCE_LABELS:
            continue
        if "intronic" in low or low.startswith("te_"):
            continue
        out.add(text)
    return out


def interpretation_hint_for_pattern(pattern, inferred_event=""):
    if inferred_event in {"te_associated_exonization", "transposable_element_exonization"}:
        return "transposable_element_associated_exonization"
    if inferred_event in {"introner_insertion", "te_intron_gain"}:
        return "introner_or_transposable_element_intron_gain"
    return {
        "segment_gain": "segment_insertion_or_annotation_recovery",
        "segment_loss": "segment_deletion_or_unresolved_annotation",
        "exonization_candidate": "exonization_or_role_gain",
        "coding_or_exonic_role_loss": "loss_of_exonic_or_coding_role",
        "segment_fusion_or_new_adjacency": "segment_fusion_rearrangement_or_new_adjacency",
        "segment_split_or_adjacency_loss": "segment_split_rearrangement_or_adjacency_loss",
        "chimeric_origin_or_source_mixing": "source_mixture_in_derived_copy",
        "chimeric_source_join_candidate": "source_joining_or_chimeric_gene_structure",
        "copy_duplication_or_expansion": "copy_number_expansion_context",
        "copy_duplication_or_relocation": "copy_number_or_location_context",
        "retrocopy_or_dispersed_duplication_candidate": "retrocopy_like_or_dispersed_copy_context",
        "copy_loss_or_collapse": "copy_loss_or_assembly_collapse_context",
        "sequence_supported_annotation_gap": "annotation_dropout_or_hidden_segment",
        "splice_boundary_shift": "splice_boundary_shift_or_intron_sliding",
        "te_associated_exonization": "transposable_element_associated_exonization",
        "introner_or_te_intron_gain": "introner_or_transposable_element_intron_gain",
        "pseudogenization_or_frame_disruption": "frame_disruption_or_pseudogenization_candidate",
        "ambiguous_paralogous_similarity": "gene_conversion_or_recent_duplication_or_unresolved_paralogy",
    }.get(pattern, "unresolved_structural_mechanism")


def interpretation_caveat_for_scope(call_scope):
    return {
        "core_structural_event": "structural evidence only; mechanism requires external biological evidence",
        "copy_context": "copy-number context; copy origin requires gene-tree and locus-level evidence",
        "annotation_evidence": "sequence-supported annotation evidence; expression or transcript evidence is not used",
        "ambiguous_evidence": "ambiguous correspondence evidence; do not assign a mechanism from this table alone",
    }.get(call_scope, "interpret as structural evidence with unresolved biological mechanism")


def make_event(
    family_id,
    event_type,
    structural_pattern,
    evidence_layer,
    object_id,
    branch_scope,
    event_probability,
    ctmc_change_probability,
    change,
    alternative_explanation,
    call_scope=None,
):
    call_scope = call_scope or call_scope_for_pattern(structural_pattern)
    return {
        "family_id": family_id,
        "event_type": event_type,
        "event_class": structural_pattern,
        "structural_change_type": structural_pattern,
        "structural_pattern": structural_pattern,
        "call_scope": call_scope,
        "evidence_layer": evidence_layer,
        "object_id": object_id,
        "branch_scope": branch_scope,
        "event_probability": event_probability,
        "ctmc_change_probability": ctmc_change_probability,
        "change": change,
        "alternative_explanation": alternative_explanation,
    }


def interpretation_hints(event_rows):
    rows = []
    for event in event_rows:
        pattern = event.get("structural_change_type") or event.get("structural_pattern") or event.get("event_class", "NA")
        call_scope = event.get("call_scope", "core_structural_event")
        rows.append(
            {
                "family_id": event.get("family_id", "NA"),
                "object_id": event.get("object_id", "NA"),
                "branch_scope": event.get("branch_scope", "NA"),
                "evidence_layer": event.get("evidence_layer", "NA"),
                "structural_change_type": pattern,
                "possible_interpretation": interpretation_hint_for_pattern(pattern, event.get("event_type", "")),
                "interpretation_caveat": interpretation_caveat_for_scope(call_scope),
                "call_scope": call_scope,
            }
        )
    return rows


def state_from_occurrences(rows):
    states = [norm_state(row.get("presence_status")) for row in rows]
    if "present" in states:
        return "present"
    if states and all(state == "absent" for state in states):
        return "absent"
    return "unknown"


def role_from_occurrences(rows):
    return element_role_from_occurrences(rows)


def classify_branch_event(layer, change):
    if "->" not in change:
        return "ambiguous_structural_change"
    src, dst = change.split("->", 1)
    if layer in {"segment_presence", "element_presence"}:
        if src == "absent" and dst == "present":
            return "segment_gain"
        if src == "present" and dst == "absent":
            return "segment_loss"
        return "element_presence_shift"
    if layer in {"role_state", "element_role_state"}:
        if src in {"non_exonic_source", "intron_or_noncoding", "absent"} and dst in {"CDS", "exon_or_UTR"}:
            return "exonization_candidate"
        if src in {"CDS", "exon_or_UTR"} and dst in {"non_exonic_source", "intron_or_noncoding"}:
            return "coding_or_exonic_role_loss"
        return "element_role_shift"
    if layer in {"adjacency_state", "element_adjacency_state"}:
        if src == "absent" and dst in {"present", "copy_variable"}:
            return "segment_fusion_or_new_adjacency"
        if src == "present" and dst in {"absent", "copy_variable"}:
            return "segment_split_or_adjacency_loss"
        if src == "copy_variable" and dst == "absent":
            return "segment_split_or_adjacency_loss"
        if src == "copy_variable" and dst == "present":
            return "segment_fusion_or_new_adjacency"
        return "adjacency_shift"
    if layer == "source_mixture":
        if dst == "multi_source":
            return "chimeric_origin_or_source_mixing"
        return "source_mixture_shift"
    if layer == "copy_multiplicity":
        if dst == "tandem_multi_copy":
            return "copy_duplication_or_expansion"
        if dst in {"same_contig_multi_copy", "dispersed_multi_copy"}:
            return "copy_duplication_or_relocation"
        if src == "tandem_multi_copy" and dst == "single_copy":
            return "copy_loss_or_collapse"
        return "copy_multiplicity_shift"
    return "structural_state_change"


def completion_event_class(call, inferred_event=""):
    if inferred_event in {"te_associated_exonization", "transposable_element_exonization"}:
        return "te_associated_exonization"
    if inferred_event in {"introner_insertion", "te_intron_gain"}:
        return "introner_or_te_intron_gain"
    return {
        "hidden_segment_candidate": "sequence_supported_annotation_gap",
        "shifted_splice_site_candidate": "splice_boundary_shift",
        "joined_exon_candidate": "segment_fusion_or_new_adjacency",
        "hidden_segment_with_frame_disruption": "pseudogenization_or_frame_disruption",
    }.get(call, "annotation_or_alignment_evidence")


def add_q_values(rows):
    numeric = []
    for idx, row in enumerate(rows):
        p_value = to_float(row.get("p_value"), None)
        if p_value is not None:
            numeric.append((idx, p_value))
    for row in rows:
        row["q_value"] = "NA"
        row["q_value_method"] = "NA"
    if not numeric:
        return rows
    numeric.sort(key=lambda item: item[1])
    m = len(numeric)
    adjusted = {}
    running = 1.0
    for rank, (idx, p_value) in reversed(list(enumerate(numeric, start=1))):
        running = min(running, p_value * m / rank)
        adjusted[idx] = min(1.0, running)
    for idx, q_value in adjusted.items():
        rows[idx]["q_value"] = f"{q_value:.6g}"
        rows[idx]["q_value_method"] = "benjamini_hochberg"
    return rows


def event_support_summary(event_rows, hypothesis_rows, bootstrap_rows, stochastic_rows):
    tests = {(row.get("layer"), row.get("object_id")): row for row in hypothesis_rows}
    boot = {(row.get("layer"), row.get("object_id")): row for row in bootstrap_rows}
    histories = {}
    for row in stochastic_rows:
        scope = f"{row.get('parent_label')}->{row.get('child_label')}"
        histories[(row.get("layer"), row.get("object_id"), scope)] = row
    out = []
    for event in event_rows:
        key = (event.get("evidence_layer"), event.get("object_id"))
        test = tests.get(key, {})
        bootstrap = boot.get(key, {})
        history = histories.get((event.get("evidence_layer"), event.get("object_id"), event.get("branch_scope")), {})
        p_value = to_float(test.get("p_value"), None)
        q_value = to_float(test.get("q_value"), None)
        empirical = to_float(bootstrap.get("empirical_p_value"), None)
        ctmc = to_float(event.get("ctmc_change_probability"), None)
        stochastic = to_float(history.get("posterior_pr_any_change"), None)
        event_probability = to_float(event.get("event_probability"), None)
        tier = "qualitative"
        if any(value is not None and value >= 0.8 for value in [ctmc, stochastic]) or (empirical is not None and empirical <= 0.05) or (q_value is not None and q_value <= 0.10):
            tier = "high"
        elif any(value is not None and value >= 0.5 for value in [ctmc, stochastic, event_probability]) or (p_value is not None and p_value <= 0.10):
            tier = "moderate"
        evidence_count = 1 + sum(value is not None for value in [p_value, q_value, empirical, ctmc, stochastic, event_probability])
        out.append(
            {
                "family_id": event.get("family_id", "NA"),
                "event_class": event.get("event_class", "NA"),
                "structural_change_type": event.get("structural_change_type", event.get("structural_pattern", event.get("event_class", "NA"))),
                "structural_pattern": event.get("structural_pattern", event.get("event_class", "NA")),
                "call_scope": event.get("call_scope", "core_structural_event"),
                "object_id": event.get("object_id", "NA"),
                "branch_scope": event.get("branch_scope", "NA"),
                "evidence_layer": event.get("evidence_layer", "NA"),
                "support_tier": tier,
                "lrt_p_value": test.get("p_value", "NA"),
                "lrt_q_value": test.get("q_value", "NA"),
                "empirical_p_value": bootstrap.get("empirical_p_value", "NA"),
                "fitted_rate": test.get("fitted_rate", "NA"),
                "ctmc_change_probability": event.get("ctmc_change_probability", "NA"),
                "stochastic_pr_any_change": history.get("posterior_pr_any_change", "NA"),
                "evidence_count": evidence_count,
                "alternative_explanation": event.get("alternative_explanation", "NA"),
            }
        )
    return out


def best_branch_support(branch_rows):
    out = {}
    for row in branch_rows:
        if row.get("status") != "change_required":
            continue
        key = (row.get("layer"), row.get("object_id"))
        score = to_float(row.get("ctmc_change_probability"), None)
        if score is None:
            score = to_float(row.get("event_probability"), 0.0)
        current = out.get(key)
        if current is None or score > current["score"]:
            out[key] = {
                "branch_scope": f"{row.get('parent_label')}->{row.get('child_label')}",
                "event_probability": row.get("event_probability", "NA"),
                "ctmc_change_probability": row.get("ctmc_change_probability", "NA"),
                "change": row.get("change", "NA"),
                "score": score,
            }
    return out
