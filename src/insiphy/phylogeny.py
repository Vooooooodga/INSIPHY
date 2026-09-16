"""Fixed-tree structural inference for intragenic structural evolution."""

from collections import defaultdict

from .elements import (
    collect_element_profiles,
    element_id_from_homology,
    element_role_from_occurrences,
)
from .io import norm_state, read_tsv, to_float, write_tsv
from .tree import (
    SpeciesTree,
    bootstrap_invariant_test,
    ctmc_posteriors,
    discrete_log_likelihood,
    fit_discrete_ctmc,
    fit_foreground_rate_test,
    fit_invariant_test,
    sankoff,
    stochastic_map_summary,
)


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


def read_foreground_edges(path, tree):
    if not path:
        return set()
    rows = read_tsv(path, optional=True)
    label_to_node = {label: node for node, label in tree.label.items()}
    out = set()
    for row in rows:
        parent = row.get("parent_node") or row.get("parent_id") or ""
        child = row.get("child_node") or row.get("child_id") or ""
        scope = row.get("branch_scope") or row.get("branch") or ""
        if (parent, child) in set(tree.edges()):
            out.add((parent, child))
            continue
        if scope and "->" in scope:
            parent_label, child_label = scope.split("->", 1)
            parent_node = label_to_node.get(parent_label.strip())
            child_node = label_to_node.get(child_label.strip())
            if parent_node and child_node:
                out.add((parent_node, child_node))
    return out


def copy_tip_label(occ):
    return f"{occ.get('species', 'NA')}:{occ.get('gene_copy_id', 'NA')}"


def copy_tip_label_from_adjacency(row):
    return f"{row.get('species', 'NA')}:{row.get('gene_copy_id', 'NA')}"


def tree_tip_label_for_occ(occ, tree, tree_scope):
    if tree_scope == "species_tree":
        return occ.get("species", "NA")
    preferred = copy_tip_label(occ)
    if preferred in tree.leaf_by_label:
        return preferred
    copy_id = occ.get("gene_copy_id", "NA")
    if copy_id in tree.leaf_by_label:
        return copy_id
    return preferred


def tree_tip_label_for_adjacency(row, tree, tree_scope):
    if tree_scope == "species_tree":
        return row.get("species", "NA")
    preferred = copy_tip_label_from_adjacency(row)
    if preferred in tree.leaf_by_label:
        return preferred
    copy_id = row.get("gene_copy_id", "NA")
    if copy_id in tree.leaf_by_label:
        return copy_id
    return preferred


def read_structural_tree(input_dir, species_tree):
    copy_rows = read_tsv(f"{input_dir}/copy_tree.tsv", ["node_id", "parent_id", "label"], optional=True)
    if copy_rows:
        return SpeciesTree(copy_rows), "copy_tree", "copy_tree.tsv"
    gene_rows = read_tsv(f"{input_dir}/gene_tree.tsv", ["node_id", "parent_id", "label"], optional=True)
    if gene_rows:
        return SpeciesTree(gene_rows), "gene_tree", "gene_tree.tsv"
    return species_tree, "species_tree", "species_tree.tsv"


def has_multicopy_species(occurrences):
    copies_by_family_species = defaultdict(set)
    for row in occurrences:
        copies_by_family_species[(row.get("family_id", "NA"), row.get("species", "NA"))].add(row.get("gene_copy_id", "NA"))
    return any(len(copies) > 1 for copies in copies_by_family_species.values())


def tree_depths(tree):
    depths = {tree.root: 0}
    for node in tree.preorder():
        for child in tree.children.get(node, []):
            depths[child] = depths[node] + 1
    return depths


def ancestor_set(tree, node):
    ancestors = []
    cur = node
    while cur:
        ancestors.append(cur)
        cur = tree.parent.get(cur, "")
    return ancestors


def mrca_node(tree, nodes):
    if not nodes:
        return "NA"
    depths = tree_depths(tree)
    common = set(ancestor_set(tree, nodes[0]))
    for node in nodes[1:]:
        common &= set(ancestor_set(tree, node))
    if not common:
        return "NA"
    return max(common, key=lambda node: depths.get(node, 0))


def phylogenetic_coverage(tree, occ_by_object, occ_by_id, object_field):
    rows = []
    total_tips = len(tree.leaves)
    for object_id, occ_ids in sorted(occ_by_object.items()):
        species = set()
        copies = set()
        families = set()
        for occ_id in occ_ids:
            occ = occ_by_id.get(occ_id)
            if not occ or norm_state(occ.get("presence_status")) != "present":
                continue
            species.add(occ.get("species", "NA"))
            copies.add(f"{occ.get('species', 'NA')}:{occ.get('gene_copy_id', 'NA')}")
            families.add(occ.get("family_id", "NA"))
        nodes = [tree.leaf_by_label[sp] for sp in species if sp in tree.leaf_by_label]
        mrca = mrca_node(tree, nodes)
        count = len({sp for sp in species if sp in tree.leaf_by_label})
        ratio = count / max(1, total_tips)
        if count == 0:
            coverage_class = "absent_or_unplaced"
        elif count == total_tips:
            coverage_class = "tree_spanning"
        elif count == 1:
            coverage_class = "tip_specific"
        else:
            coverage_class = "partial"
        rows.append(
            {
                "family_id": ";".join(sorted(families)) if families else "NA",
                object_field: object_id,
                "present_species_count": count,
                "tree_tip_count": total_tips,
                "coverage_ratio": f"{ratio:.6g}",
                "mrca_node": mrca,
                "mrca_label": tree.label.get(mrca, "NA"),
                "coverage_class": coverage_class,
                "present_species": ";".join(sorted(species)) if species else "NA",
                "present_copy_count": len(copies),
                "present_copies": ";".join(sorted(copies)) if copies else "NA",
            }
        )
    return rows


def internal_homology_phylogenetic_coverage(tree, occ_by_homology, occ_by_id):
    return phylogenetic_coverage(tree, occ_by_homology, occ_by_id, "homology_id")


def element_phylogenetic_coverage(tree, occ_by_element, occ_by_id):
    return phylogenetic_coverage(tree, occ_by_element, occ_by_id, "element_id")


def fallback_element_rows(homology, occurrences, input_dir, output_dir):
    rows = read_tsv(f"{output_dir}/element_correspondence.tsv", optional=True)
    if rows:
        return rows
    evidence = read_tsv(f"{input_dir}/sequence_synteny_evidence.tsv", optional=True)
    annotation = read_tsv(f"{output_dir}/annotation_completion_candidates.tsv", optional=True)
    profiles, element_by_homology = collect_element_profiles(homology, occurrences, evidence + annotation)
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    out = []
    for row in homology:
        element_id = element_by_homology.get(row["homology_id"])
        if not element_id:
            continue
        occ = occ_by_id.get(row["occurrence_id"], {})
        out.append(
            {
                "element_id": element_id,
                "family_id": occ.get("family_id", "NA"),
                "homology_id": row["homology_id"],
                "occurrence_id": row["occurrence_id"],
                "species": occ.get("species", "NA"),
                "gene_copy_id": occ.get("gene_copy_id", "NA"),
                "element_class": "exon_like" if occ.get("role") in EXONIC else "candidate_source",
                "display_role": role_from_occurrences([occ]),
                "source_label": row.get("source_label", "NA"),
                "support_type": row.get("support_type", "NA"),
                "confidence": row.get("confidence", "NA"),
                "membership_score": row.get("confidence", "NA"),
                "membership_call": "core_member",
            }
        )
    return out


def infer_phylogeny(
    input_dir,
    output_dir,
    bootstrap_replicates=0,
    stochastic_maps=0,
    seed=7,
    foreground_branches=None,
    analysis_scope="single-copy",
    model="parsimony",
    branch_length_mode="supplied",
    ascertainment="observed-at-least-one",
    threads=1,
    root_frequency="estimated",
    root_presence=0.5,
):
    if analysis_scope == "single-copy":
        if bootstrap_replicates or stochastic_maps:
            raise SystemExit(
                "single-copy analysis uses parsimony or analytic likelihood; "
                "--bootstrap-replicates and --stochastic-maps must be 0"
            )
        if model != "foreground" and foreground_branches:
            raise SystemExit("--foreground-branches requires --model foreground")
        if int(threads) < 1:
            raise SystemExit("--threads must be at least 1")
        if model == "parsimony":
            from .parsimony import infer_single_copy_parsimony

            return infer_single_copy_parsimony(input_dir, output_dir, threads=threads)
        from .structural_phylogeny import infer_single_copy_phylogeny

        return infer_single_copy_phylogeny(
            input_dir,
            output_dir,
            model=model,
            foreground_branches=foreground_branches,
            branch_length_mode=branch_length_mode,
            ascertainment=ascertainment,
            threads=threads,
            root_frequency=root_frequency,
            root_presence=root_presence,
        )
    occurrences = read_tsv(f"{input_dir}/segment_occurrences.tsv", ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"])
    homology = read_tsv(f"{input_dir}/segment_homology.tsv", ["homology_id", "occurrence_id", "support_type", "confidence"])
    adjacencies = read_tsv(f"{input_dir}/physical_adjacencies.tsv", ["adjacency_id", "family_id", "species", "gene_copy_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status"])
    species_tree = SpeciesTree(read_tsv(f"{input_dir}/species_tree.tsv", ["node_id", "parent_id", "label"]))
    structural_tree, structural_tree_scope, structural_tree_file = read_structural_tree(input_dir, species_tree)
    phylogeny_scope_rows = [
        {
            "scope": "structural_characters",
            "tree_file": structural_tree_file,
            "tree_scope": structural_tree_scope,
            "layers": "element_presence;element_role_state;element_adjacency_state;source_mixture",
            "note": "EG structural characters are evaluated on copy_tree.tsv or gene_tree.tsv when supplied; species_tree.tsv is the fallback for single-copy cases.",
        },
        {
            "scope": "copy_multiplicity",
            "tree_file": "species_tree.tsv",
            "tree_scope": "species_tree",
            "layers": "copy_multiplicity",
            "note": "Copy multiplicity is a species-level character and remains evaluated on species_tree.tsv.",
        },
    ]
    if structural_tree_scope == "species_tree" and has_multicopy_species(occurrences):
        phylogeny_scope_rows.append(
            {
                "scope": "warning",
                "tree_file": "species_tree.tsv",
                "tree_scope": "species_tree_fallback",
                "layers": "element_presence;element_role_state;element_adjacency_state;source_mixture",
                "note": "Multiple gene copies occur in at least one species, but no copy_tree.tsv or gene_tree.tsv was supplied; copy-specific structural histories are collapsed to species-level tips.",
            }
        )
    foreground_edges = read_foreground_edges(foreground_branches, structural_tree)
    species_foreground_edges = read_foreground_edges(foreground_branches, species_tree)
    copy_context = read_tsv(f"{input_dir}/copy_context.tsv", ["family_id", "species", "gene_copy_id", "copy_class"], optional=True)
    copy_relationships = read_tsv(f"{input_dir}/copy_relationships.tsv", ["family_id", "species", "query_copy_id", "subject_copy_id", "relationship_class"], optional=True)
    segment_matches = read_tsv(f"{input_dir}/segment_matches.tsv", optional=True)
    annotation_candidates = read_tsv(f"{output_dir}/annotation_completion_candidates.tsv", optional=True)

    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    homology_by_occ = defaultdict(list)
    occ_by_homology = defaultdict(list)
    source_by_occ = defaultdict(set)
    for row in homology:
        homology_by_occ[row["occurrence_id"]].append(row["homology_id"])
        occ_by_homology[row["homology_id"]].append(row["occurrence_id"])
        if row.get("source_label"):
            source_by_occ[row["occurrence_id"]].add(row["source_label"])
    element_rows = fallback_element_rows(homology, occurrences, input_dir, output_dir)
    element_by_occ = defaultdict(list)
    occ_by_element = defaultdict(list)
    for row in element_rows:
        if row.get("element_class") == "context":
            continue
        element_id = row["element_id"]
        occ_id = row["occurrence_id"]
        element_by_occ[occ_id].append(element_id)
        occ_by_element[element_id].append(occ_id)

    state_rows = []
    branch_rows = []
    model_score_rows = []
    model_fit_rows = []
    hypothesis_rows = []
    bootstrap_rows = []
    stochastic_rows = []
    foreground_rows = []
    object_family = {}
    for element_id, occ_ids in sorted(occ_by_element.items()):
        by_tip = defaultdict(list)
        families = set()
        for occ_id in occ_ids:
            occ = occ_by_id.get(occ_id)
            if occ:
                by_tip[tree_tip_label_for_occ(occ, structural_tree, structural_tree_scope)].append(occ)
                families.add(occ["family_id"])
        object_family[element_id] = ";".join(sorted(families)) if families else element_id
        add_character(structural_tree, "element_presence", element_id, {tip: state_from_occurrences(rows) for tip, rows in by_tip.items()}, {"present", "absent"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)
        add_character(structural_tree, "element_role_state", element_id, {tip: role_from_occurrences(rows) for tip, rows in by_tip.items()}, {"absent", "CDS", "exon_or_UTR", "non_exonic_source", "unknown"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)

    graph_edges = []
    adj_by_pair = defaultdict(lambda: defaultdict(list))
    pair_family = defaultdict(set)
    for row in adjacencies:
        left_elements = sorted(element_by_occ.get(row["left_occurrence_id"], []))
        right_elements = sorted(element_by_occ.get(row["right_occurrence_id"], []))
        if not left_elements or not right_elements:
            continue
        pair_id = "|".join(left_elements) + "__" + "|".join(right_elements)
        adj_by_pair[pair_id][tree_tip_label_for_adjacency(row, structural_tree, structural_tree_scope)].append(row["adjacency_status"])
        pair_family[pair_id].add(row["family_id"])
        left_sources = informative_source_labels(source_by_occ.get(row["left_occurrence_id"], set()))
        right_sources = informative_source_labels(source_by_occ.get(row["right_occurrence_id"], set()))
        graph_edges.append(
            {
                "family_id": row["family_id"],
                "species": row["species"],
                "gene_copy_id": row["gene_copy_id"],
                "edge_id": row["adjacency_id"],
                "left_element_id": "|".join(left_elements),
                "right_element_id": "|".join(right_elements),
                "left_occurrence_id": row["left_occurrence_id"],
                "right_occurrence_id": row["right_occurrence_id"],
                "adjacency_status": norm_state(row["adjacency_status"]),
                "left_source_labels": ";".join(sorted(left_sources)) or "NA",
                "right_source_labels": ";".join(sorted(right_sources)) or "NA",
            }
        )
    for pair_id, by_tip in sorted(adj_by_pair.items()):
        object_family[pair_id] = ";".join(sorted(pair_family[pair_id]))
        tips = {}
        for tip_label, vals in by_tip.items():
            vals = [norm_state(value) for value in vals]
            if "present" in vals and "absent" in vals:
                tips[tip_label] = "copy_variable"
            elif "present" in vals:
                tips[tip_label] = "present"
            elif vals and all(value == "absent" for value in vals):
                tips[tip_label] = "absent"
            else:
                tips[tip_label] = "unknown"
        add_character(structural_tree, "element_adjacency_state", pair_id, tips, {"present", "absent", "copy_variable"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)

    source_states_by_family_tip = defaultdict(list)
    for key in sorted({(row["family_id"], row["species"], row["gene_copy_id"]) for row in occurrences}):
        family, species, copy = key
        rows = [row for row in occurrences if (row["family_id"], row["species"], row["gene_copy_id"]) == key and norm_state(row["presence_status"]) == "present"]
        sources = informative_source_labels(source for row in rows for source in source_by_occ.get(row["occurrence_id"], set()))
        tip_label = tree_tip_label_for_occ({"species": species, "gene_copy_id": copy}, structural_tree, structural_tree_scope)
        source_states_by_family_tip[(family, tip_label)].append("multi_source" if len(sources) > 1 else "single_source" if len(sources) == 1 else "unknown")
    source_tips = defaultdict(dict)
    for (family, tip_label), states in source_states_by_family_tip.items():
        source_tips[family][tip_label] = "multi_source" if "multi_source" in states else "single_source" if "single_source" in states else "unknown"
    for family, tips in sorted(source_tips.items()):
        object_family[family] = family
        add_character(structural_tree, "source_mixture", family, tips, {"single_source", "multi_source", "unknown"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)

    copy_states_by_family_species = defaultdict(list)
    for row in copy_context:
        copy_states_by_family_species[(row["family_id"], row["species"])].append(row.get("copy_class", "unknown"))
    copy_tips = defaultdict(dict)
    for (family, species), states in copy_states_by_family_species.items():
        if "tandem_multi_copy" in states:
            copy_tips[family][species] = "tandem_multi_copy"
        elif "same_contig_multi_copy" in states:
            copy_tips[family][species] = "same_contig_multi_copy"
        elif "dispersed_multi_copy" in states:
            copy_tips[family][species] = "dispersed_multi_copy"
        elif "unresolved" in states:
            copy_tips[family][species] = "unresolved"
        elif "single_copy" in states:
            copy_tips[family][species] = "single_copy"
        else:
            copy_tips[family][species] = "unknown"
    for family, tips in sorted(copy_tips.items()):
        object_family[family] = family
        add_character(species_tree, "copy_multiplicity", family, tips, {"single_copy", "tandem_multi_copy", "same_contig_multi_copy", "dispersed_multi_copy", "unresolved", "unknown"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, species_foreground_edges)

    hidden_calls = {"hidden_segment_candidate", "shifted_splice_site_candidate", "joined_exon_candidate", "hidden_segment_with_frame_disruption"}
    hidden_support = sum(to_float(row.get("support_score")) for row in annotation_candidates if row.get("completion_call") in hidden_calls)
    hidden_count = sum(1 for row in annotation_candidates if row.get("completion_call") in hidden_calls)
    multi_source_count = sum(1 for tips in source_tips.values() for state in tips.values() if state == "multi_source")
    source_join_count = sum(1 for edge in graph_edges if edge["adjacency_status"] == "present" and edge["left_source_labels"] != "NA" and edge["right_source_labels"] != "NA" and edge["left_source_labels"] != edge["right_source_labels"])

    model_rows = [
        {"comparison_id": "annotation_error", "model": "strict_annotation", "score": f"{hidden_count:.6g}", "delta_vs_best": "NA", "interpretation": "penalizes hidden segments as missing_or_unknown"},
        {"comparison_id": "annotation_error", "model": "annotation_error", "score": f"{max(0.0, hidden_count - hidden_support):.6g}", "delta_vs_best": "NA", "interpretation": "allows sequence_synteny_supported_hidden_segments"},
        {"comparison_id": "compound_event", "model": "independent_character", "score": f"{(2.0 * source_join_count + multi_source_count):.6g}", "delta_vs_best": "NA", "interpretation": "counts source mixture and source joining separately"},
        {"comparison_id": "compound_event", "model": "compound_chimeric_or_copy_event", "score": f"{(2.5 if source_join_count or multi_source_count else 0.0):.6g}", "delta_vs_best": "NA", "interpretation": "one copy-level event explains coordinated source mixture and adjacency"},
    ]
    for comparison in sorted({row["comparison_id"] for row in model_rows}):
        rows = [row for row in model_rows if row["comparison_id"] == comparison]
        best = min(to_float(row["score"]) for row in rows)
        for row in rows:
            row["delta_vs_best"] = f"{to_float(row['score']) - best:.6g}"

    event_rows = []
    branch_support = best_branch_support(branch_rows)
    for row in branch_rows:
        if row["status"] == "change_required":
            pattern = classify_branch_event(row["layer"], row["change"])
            event_rows.append(
                make_event(
                    object_family.get(row["object_id"], row["object_id"].split("__")[0]),
                    row["event_type"],
                    pattern,
                    row["layer"],
                    row["object_id"],
                    f"{row['parent_label']}->{row['child_label']}",
                    row["event_probability"],
                    row.get("ctmc_change_probability", "NA"),
                    row["change"],
                    "annotation_gap_or_mapping_ambiguity_if_sequence_support_low",
                )
            )
    for edge in graph_edges:
        if edge["adjacency_status"] == "present" and edge["left_source_labels"] != "NA" and edge["right_source_labels"] != "NA" and edge["left_source_labels"] != edge["right_source_labels"]:
            pair_id = f"{edge['left_element_id']}__{edge['right_element_id']}"
            support = branch_support.get(("element_adjacency_state", pair_id), {})
            event_rows.append(
                make_event(
                    edge["family_id"],
                    "source_join_candidate",
                    "chimeric_source_join_candidate",
                    "element_adjacency_state",
                    pair_id,
                    support.get("branch_scope", "estimated_from_adjacency_state_history"),
                    support.get("event_probability", "NA"),
                    support.get("ctmc_change_probability", "NA"),
                    f"{edge['left_source_labels']}->{edge['right_source_labels']}",
                    "paralogy_or_homology_assignment_error_if_low_support",
                    call_scope="core_structural_event",
                )
            )
    for row in copy_relationships:
        rel = row.get("relationship_class", "")
        if rel in {"tandem_duplication_candidate", "same_contig_duplication_candidate", "dispersed_or_retrocopy_candidate"}:
            event_class = "retrocopy_or_dispersed_duplication_candidate" if rel == "dispersed_or_retrocopy_candidate" else "copy_duplication_or_expansion"
            support = branch_support.get(("copy_multiplicity", row["family_id"]), {})
            event_rows.append(
                make_event(
                    row["family_id"],
                    rel,
                    event_class,
                    "copy_multiplicity",
                    row["family_id"],
                    support.get("branch_scope", "tip_copy_relationship_requires_phylogenetic_placement"),
                    row.get("synteny_score", "NA"),
                    support.get("ctmc_change_probability", "NA"),
                    f"{row['species']}:{row['query_copy_id']}--{row['subject_copy_id']}:{rel}",
                    "assembly_fragmentation_or_unresolved_paralogy_if_low_support",
                    call_scope="copy_context",
                )
            )

    for row in annotation_candidates:
        call = row.get("completion_call", "")
        if call in {"hidden_segment_candidate", "shifted_splice_site_candidate", "joined_exon_candidate", "hidden_segment_with_frame_disruption"}:
            pattern = completion_event_class(call, row.get("inferred_event", ""))
            object_id = element_id_from_homology(row.get("homology_id", "")) if row.get("homology_id") else row.get("evidence_id", "NA")
            event_rows.append(
                make_event(
                    row.get("family_id", "NA"),
                    call,
                    pattern,
                    "annotation_completion",
                    object_id,
                    "tip_or_alignment_evidence_requires_phylogenetic_context",
                    row.get("support_score", "NA"),
                    "NA",
                    row.get("inferred_event", call),
                    "annotation_dropout_or_shifted_boundary_if_sequence_support_is_partial",
                )
            )

    for row in segment_matches:
        if row.get("match_status") != "mapped":
            continue
        score = to_float(row.get("total_score"), 0.0)
        if score < 0.90:
            continue
        left = occ_by_id.get(row.get("query_occurrence_id"))
        right = occ_by_id.get(row.get("subject_occurrence_id"))
        if not left or not right:
            continue
        if left["species"] == right["species"] and left["gene_copy_id"] != right["gene_copy_id"]:
            event_rows.append(
                make_event(
                    left.get("family_id", "NA"),
                    "high_identity_paralogous_segment_match",
                    "ambiguous_paralogous_similarity",
                    "segment_correspondence_graph",
                    row.get("match_id", "NA"),
                    "tip_paralogous_similarity_requires_phylogenetic_context",
                    row.get("total_score", "NA"),
                    "NA",
                    f"{left['gene_copy_id']}<->{right['gene_copy_id']}",
                    "recent_duplication_or_unresolved_paralogy_if_context_support_low",
                    call_scope="ambiguous_evidence",
                )
            )

    case_summary = []
    for family in sorted({row["family_id"] for row in occurrences}):
        family_events = [row for row in event_rows if row["family_id"] == family or row["object_id"].startswith(family)]
        case_summary.append(
            {
                "family_id": family,
                "hidden_segment_candidates": hidden_count,
                "source_join_candidates": source_join_count,
                "multi_source_tip_count": multi_source_count,
                "branch_event_candidates": sum(1 for row in family_events if row.get("call_scope") == "core_structural_event"),
                "copy_context_candidates": sum(1 for row in family_events if row.get("call_scope") == "copy_context"),
                "ambiguous_evidence_candidates": sum(1 for row in family_events if row.get("call_scope") == "ambiguous_evidence"),
                "best_annotation_model": min([row for row in model_rows if row["comparison_id"] == "annotation_error"], key=lambda row: to_float(row["score"]))["model"],
                "best_compound_model": min([row for row in model_rows if row["comparison_id"] == "compound_event"], key=lambda row: to_float(row["score"]))["model"],
            }
        )

    write_tsv(f"{output_dir}/ancestral_state_probabilities.tsv", state_rows, ["layer", "object_id", "score", "node_id", "node_label", "state", "probability", "is_parsimony_best"])
    write_tsv(f"{output_dir}/branch_event_probabilities.tsv", branch_rows, ["layer", "object_id", "event_type", "parent_node", "child_node", "parent_label", "child_label", "branch_length", "status", "change", "event_probability", "ctmc_change_probability", "ctmc_most_likely_change", "ctmc_most_likely_change_probability"])
    add_q_values(hypothesis_rows)
    add_q_values(foreground_rows)
    support_rows = event_support_summary(event_rows, hypothesis_rows, bootstrap_rows, stochastic_rows)
    hint_rows = interpretation_hints(event_rows)
    homology_coverage_rows = internal_homology_phylogenetic_coverage(species_tree, occ_by_homology, occ_by_id)
    element_coverage_rows = element_phylogenetic_coverage(species_tree, occ_by_element, occ_by_id)
    write_tsv(f"{output_dir}/character_model_scores.tsv", model_score_rows, ["layer", "object_id", "model", "parsimony_score", "log_likelihood", "state_count", "observed_tip_count", "fitted_rate", "aic", "bic"])
    write_tsv(f"{output_dir}/model_fit.tsv", model_fit_rows, ["layer", "object_id", "model", "fitted_rate", "log_likelihood", "aic", "bic", "observed_tip_count"])
    write_tsv(f"{output_dir}/hypothesis_tests.tsv", hypothesis_rows, ["layer", "object_id", "test_id", "null_model", "alternative_model", "null_log_likelihood", "alternative_log_likelihood", "lrt_statistic", "df", "p_value", "p_value_method", "q_value", "q_value_method", "fitted_rate", "null_aic", "alternative_aic", "null_bic", "alternative_bic", "observed_tip_count", "tip_error"])
    write_tsv(f"{output_dir}/hypothesis_bootstrap.tsv", bootstrap_rows, ["layer", "object_id", "test_id", "observed_lrt", "bootstrap_replicates", "empirical_p_value", "monte_carlo_se", "null_lrt_mean", "null_lrt_q025", "null_lrt_q500", "null_lrt_q975", "seed", "tip_error"])
    write_tsv(f"{output_dir}/branch_history_posteriors.tsv", stochastic_rows, ["layer", "object_id", "parent_node", "child_node", "parent_label", "child_label", "map_sample_count", "posterior_pr_any_change", "posterior_expected_change_count", "posterior_change_count_low", "posterior_change_count_high", "posterior_most_frequent_transition", "posterior_transition_probability"])
    write_tsv(f"{output_dir}/foreground_tests.tsv", foreground_rows, ["layer", "object_id", "test_id", "null_model", "alternative_model", "null_log_likelihood", "alternative_log_likelihood", "lrt_statistic", "df", "p_value", "p_value_method", "q_value", "q_value_method", "background_rate", "foreground_rate", "rate_ratio", "null_aic", "alternative_aic", "null_bic", "alternative_bic", "observed_tip_count"])
    write_tsv(f"{output_dir}/candidate_structural_events.tsv", event_rows, ["family_id", "event_type", "event_class", "structural_change_type", "structural_pattern", "call_scope", "evidence_layer", "object_id", "branch_scope", "event_probability", "ctmc_change_probability", "change", "alternative_explanation"])
    write_tsv(f"{output_dir}/event_support_summary.tsv", support_rows, ["family_id", "event_class", "structural_change_type", "structural_pattern", "call_scope", "object_id", "branch_scope", "evidence_layer", "support_tier", "lrt_p_value", "lrt_q_value", "empirical_p_value", "fitted_rate", "ctmc_change_probability", "stochastic_pr_any_change", "evidence_count", "alternative_explanation"])
    write_tsv(f"{output_dir}/interpretation_hints.tsv", hint_rows, ["family_id", "object_id", "branch_scope", "evidence_layer", "structural_change_type", "possible_interpretation", "interpretation_caveat", "call_scope"])
    write_tsv(f"{output_dir}/element_phylogenetic_coverage.tsv", element_coverage_rows, ["family_id", "element_id", "present_species_count", "tree_tip_count", "coverage_ratio", "mrca_node", "mrca_label", "coverage_class", "present_species", "present_copy_count", "present_copies"])
    write_tsv(f"{output_dir}/internal_homology_phylogenetic_coverage.tsv", homology_coverage_rows, ["family_id", "homology_id", "present_species_count", "tree_tip_count", "coverage_ratio", "mrca_node", "mrca_label", "coverage_class", "present_species", "present_copy_count", "present_copies"])
    write_tsv(f"{output_dir}/phylogeny_scope.tsv", phylogeny_scope_rows, ["scope", "tree_file", "tree_scope", "layers", "note"])
    write_tsv(f"{output_dir}/model_comparison.tsv", model_rows, ["comparison_id", "model", "score", "delta_vs_best", "interpretation"])
    write_tsv(f"{output_dir}/intragenic_graph_edges.tsv", graph_edges, ["family_id", "species", "gene_copy_id", "edge_id", "left_element_id", "right_element_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status", "left_source_labels", "right_source_labels"])
    write_tsv(f"{output_dir}/case_summary.tsv", case_summary, ["family_id", "hidden_segment_candidates", "source_join_candidates", "multi_source_tip_count", "branch_event_candidates", "copy_context_candidates", "ambiguous_evidence_candidates", "best_annotation_model", "best_compound_model"])
    return state_rows, branch_rows, event_rows, model_rows
