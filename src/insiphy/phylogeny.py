"""Fixed-tree structural inference for intragenic structural evolution."""

from collections import defaultdict

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


def state_from_occurrences(rows):
    states = [norm_state(row.get("presence_status")) for row in rows]
    if "present" in states:
        return "present"
    if states and all(state == "absent" for state in states):
        return "absent"
    return "unknown"


def role_from_occurrences(rows):
    roles = {row.get("role", "unknown") for row in rows if norm_state(row.get("presence_status")) == "present"}
    if not roles:
        return "absent"
    if "CDS" in roles:
        return "CDS"
    if roles & {"UTR", "noncoding_exon", "exon"}:
        return "exon_or_UTR"
    if roles & NONCODING:
        return "intron_or_noncoding"
    return "unknown"


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
    if layer == "segment_presence":
        if src == "absent" and dst == "present":
            return "segment_gain"
        if src == "present" and dst == "absent":
            return "segment_loss"
        return "segment_presence_shift"
    if layer == "role_state":
        if src in {"intron_or_noncoding", "absent"} and dst in {"CDS", "exon_or_UTR"}:
            return "exonization_candidate"
        if src in {"CDS", "exon_or_UTR"} and dst == "intron_or_noncoding":
            return "coding_or_exonic_role_loss"
        return "segment_role_shift"
    if layer == "adjacency_state":
        if src == "absent" and dst == "present":
            return "segment_fusion_or_new_adjacency"
        if src == "present" and dst == "absent":
            return "segment_split_or_adjacency_loss"
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


def infer_phylogeny(input_dir, output_dir, bootstrap_replicates=0, stochastic_maps=0, seed=7, foreground_branches=None):
    occurrences = read_tsv(f"{input_dir}/segment_occurrences.tsv", ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"])
    homology = read_tsv(f"{input_dir}/segment_homology.tsv", ["homology_id", "occurrence_id", "support_type", "confidence"])
    adjacencies = read_tsv(f"{input_dir}/physical_adjacencies.tsv", ["adjacency_id", "family_id", "species", "gene_copy_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status"])
    tree = SpeciesTree(read_tsv(f"{input_dir}/species_tree.tsv", ["node_id", "parent_id", "label"]))
    foreground_edges = read_foreground_edges(foreground_branches, tree)
    copy_context = read_tsv(f"{input_dir}/copy_context.tsv", ["family_id", "species", "gene_copy_id", "copy_class"], optional=True)
    copy_relationships = read_tsv(f"{input_dir}/copy_relationships.tsv", ["family_id", "species", "query_copy_id", "subject_copy_id", "relationship_class"], optional=True)
    segment_matches = read_tsv(f"{input_dir}/segment_matches.tsv", optional=True)
    annotation_candidates = read_tsv(f"{output_dir}/annotation_completion_candidates.tsv", optional=True)

    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    hsg_by_occ = defaultdict(list)
    occ_by_hsg = defaultdict(list)
    source_by_occ = defaultdict(set)
    for row in homology:
        hsg_by_occ[row["occurrence_id"]].append(row["homology_id"])
        occ_by_hsg[row["homology_id"]].append(row["occurrence_id"])
        if row.get("source_label"):
            source_by_occ[row["occurrence_id"]].add(row["source_label"])

    state_rows = []
    branch_rows = []
    model_score_rows = []
    model_fit_rows = []
    hypothesis_rows = []
    bootstrap_rows = []
    stochastic_rows = []
    foreground_rows = []
    object_family = {}
    for hsg, occ_ids in sorted(occ_by_hsg.items()):
        by_species = defaultdict(list)
        families = set()
        for occ_id in occ_ids:
            occ = occ_by_id.get(occ_id)
            if occ:
                by_species[occ["species"]].append(occ)
                families.add(occ["family_id"])
        object_family[hsg] = ";".join(sorted(families)) if families else hsg
        add_character(tree, "segment_presence", hsg, {sp: state_from_occurrences(rows) for sp, rows in by_species.items()}, {"present", "absent"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)
        add_character(tree, "role_state", hsg, {sp: role_from_occurrences(rows) for sp, rows in by_species.items()}, {"absent", "CDS", "exon_or_UTR", "intron_or_noncoding", "unknown"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)

    graph_edges = []
    adj_by_pair = defaultdict(lambda: defaultdict(list))
    pair_family = defaultdict(set)
    for row in adjacencies:
        left_h = sorted(hsg_by_occ.get(row["left_occurrence_id"], []))
        right_h = sorted(hsg_by_occ.get(row["right_occurrence_id"], []))
        if not left_h or not right_h:
            continue
        pair_id = "|".join(left_h) + "__" + "|".join(right_h)
        adj_by_pair[pair_id][row["species"]].append(row["adjacency_status"])
        pair_family[pair_id].add(row["family_id"])
        graph_edges.append(
            {
                "family_id": row["family_id"],
                "species": row["species"],
                "gene_copy_id": row["gene_copy_id"],
                "edge_id": row["adjacency_id"],
                "left_hsg": "|".join(left_h),
                "right_hsg": "|".join(right_h),
                "left_occurrence_id": row["left_occurrence_id"],
                "right_occurrence_id": row["right_occurrence_id"],
                "adjacency_status": norm_state(row["adjacency_status"]),
                "left_source_labels": ";".join(sorted(source_by_occ.get(row["left_occurrence_id"], set()))) or "NA",
                "right_source_labels": ";".join(sorted(source_by_occ.get(row["right_occurrence_id"], set()))) or "NA",
            }
        )
    for pair_id, by_species in sorted(adj_by_pair.items()):
        object_family[pair_id] = ";".join(sorted(pair_family[pair_id]))
        tips = {}
        for species, vals in by_species.items():
            vals = [norm_state(value) for value in vals]
            tips[species] = "present" if "present" in vals else "absent" if vals and all(value == "absent" for value in vals) else "unknown"
        add_character(tree, "adjacency_state", pair_id, tips, {"present", "absent"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)

    source_states_by_family_species = defaultdict(list)
    for key in sorted({(row["family_id"], row["species"], row["gene_copy_id"]) for row in occurrences}):
        family, species, copy = key
        rows = [row for row in occurrences if (row["family_id"], row["species"], row["gene_copy_id"]) == key and norm_state(row["presence_status"]) == "present"]
        sources = {source for row in rows for source in source_by_occ.get(row["occurrence_id"], set()) if source not in {"unknown_source", "NA"}}
        source_states_by_family_species[(family, species)].append("multi_source" if len(sources) > 1 else "single_source" if len(sources) == 1 else "unknown")
    source_tips = defaultdict(dict)
    for (family, species), states in source_states_by_family_species.items():
        source_tips[family][species] = "multi_source" if "multi_source" in states else "single_source" if "single_source" in states else "unknown"
    for family, tips in sorted(source_tips.items()):
        object_family[family] = family
        add_character(tree, "source_mixture", family, tips, {"single_source", "multi_source", "unknown"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)

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
        add_character(tree, "copy_multiplicity", family, tips, {"single_copy", "tandem_multi_copy", "same_contig_multi_copy", "dispersed_multi_copy", "unresolved", "unknown"}, state_rows, branch_rows, model_score_rows, model_fit_rows, hypothesis_rows, bootstrap_rows, stochastic_rows, foreground_rows, bootstrap_replicates, stochastic_maps, seed, foreground_edges)

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
            event_rows.append(
                {
                    "family_id": object_family.get(row["object_id"], row["object_id"].split("__")[0]),
                    "event_type": row["event_type"],
                    "event_class": classify_branch_event(row["layer"], row["change"]),
                    "evidence_layer": row["layer"],
                    "object_id": row["object_id"],
                    "branch_scope": f"{row['parent_label']}->{row['child_label']}",
                    "event_probability": row["event_probability"],
                    "ctmc_change_probability": row.get("ctmc_change_probability", "NA"),
                    "change": row["change"],
                    "alternative_explanation": "annotation_gap_or_mapping_ambiguity_if_sequence_support_low",
                }
            )
    for edge in graph_edges:
        if edge["adjacency_status"] == "present" and edge["left_source_labels"] != "NA" and edge["right_source_labels"] != "NA" and edge["left_source_labels"] != edge["right_source_labels"]:
            pair_id = f"{edge['left_hsg']}__{edge['right_hsg']}"
            support = branch_support.get(("adjacency_state", pair_id), {})
            event_rows.append(
                {
                    "family_id": edge["family_id"],
                    "event_type": "source_join_candidate",
                    "event_class": "chimeric_source_join_candidate",
                    "evidence_layer": "adjacency_state",
                    "object_id": pair_id,
                    "branch_scope": support.get("branch_scope", "estimated_from_adjacency_state_history"),
                    "event_probability": support.get("event_probability", "NA"),
                    "ctmc_change_probability": support.get("ctmc_change_probability", "NA"),
                    "change": f"{edge['left_source_labels']}->{edge['right_source_labels']}",
                    "alternative_explanation": "paralogy_or_homology_assignment_error_if_low_support",
                }
            )
    for row in copy_relationships:
        rel = row.get("relationship_class", "")
        if rel in {"tandem_duplication_candidate", "same_contig_duplication_candidate", "dispersed_or_retrocopy_candidate"}:
            event_class = "retrocopy_or_dispersed_duplication_candidate" if rel == "dispersed_or_retrocopy_candidate" else "copy_duplication_or_expansion"
            support = branch_support.get(("copy_multiplicity", row["family_id"]), {})
            event_rows.append(
                {
                    "family_id": row["family_id"],
                    "event_type": rel,
                    "event_class": event_class,
                    "evidence_layer": "copy_multiplicity",
                    "object_id": row["family_id"],
                    "branch_scope": support.get("branch_scope", "tip_copy_relationship_requires_phylogenetic_placement"),
                    "event_probability": row.get("synteny_score", "NA"),
                    "ctmc_change_probability": support.get("ctmc_change_probability", "NA"),
                    "change": f"{row['species']}:{row['query_copy_id']}--{row['subject_copy_id']}:{rel}",
                    "alternative_explanation": "assembly_fragmentation_or_unresolved_paralogy_if_low_support",
                }
            )

    for row in annotation_candidates:
        call = row.get("completion_call", "")
        if call in {"hidden_segment_candidate", "shifted_splice_site_candidate", "joined_exon_candidate", "hidden_segment_with_frame_disruption"}:
            event_rows.append(
                {
                    "family_id": row.get("family_id", "NA"),
                    "event_type": call,
                    "event_class": completion_event_class(call, row.get("inferred_event", "")),
                    "evidence_layer": "annotation_completion",
                    "object_id": row.get("evidence_id", "NA"),
                    "branch_scope": "tip_or_alignment_evidence_requires_phylogenetic_context",
                    "event_probability": row.get("support_score", "NA"),
                    "ctmc_change_probability": "NA",
                    "change": row.get("inferred_event", call),
                    "alternative_explanation": "annotation_dropout_or_shifted_boundary_if_sequence_support_is_partial",
                }
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
                {
                    "family_id": left.get("family_id", "NA"),
                    "event_type": "high_identity_paralogous_segment_match",
                    "event_class": "gene_conversion_candidate",
                    "evidence_layer": "segment_correspondence_graph",
                    "object_id": row.get("match_id", "NA"),
                    "branch_scope": "tip_paralogous_similarity_requires_phylogenetic_context",
                    "event_probability": row.get("total_score", "NA"),
                    "ctmc_change_probability": "NA",
                    "change": f"{left['gene_copy_id']}<->{right['gene_copy_id']}",
                    "alternative_explanation": "recent_duplication_or_unresolved_paralogy_if_context_support_low",
                }
            )

    case_summary = []
    for family in sorted({row["family_id"] for row in occurrences}):
        case_summary.append(
            {
                "family_id": family,
                "hidden_segment_candidates": hidden_count,
                "source_join_candidates": source_join_count,
                "multi_source_tip_count": multi_source_count,
                "branch_event_candidates": sum(1 for row in event_rows if row["family_id"] == family or row["object_id"].startswith(family)),
                "best_annotation_model": min([row for row in model_rows if row["comparison_id"] == "annotation_error"], key=lambda row: to_float(row["score"]))["model"],
                "best_compound_model": min([row for row in model_rows if row["comparison_id"] == "compound_event"], key=lambda row: to_float(row["score"]))["model"],
            }
        )

    write_tsv(f"{output_dir}/ancestral_state_probabilities.tsv", state_rows, ["layer", "object_id", "score", "node_id", "node_label", "state", "probability", "is_parsimony_best"])
    write_tsv(f"{output_dir}/branch_event_probabilities.tsv", branch_rows, ["layer", "object_id", "event_type", "parent_node", "child_node", "parent_label", "child_label", "branch_length", "status", "change", "event_probability", "ctmc_change_probability", "ctmc_most_likely_change", "ctmc_most_likely_change_probability"])
    add_q_values(hypothesis_rows)
    add_q_values(foreground_rows)
    support_rows = event_support_summary(event_rows, hypothesis_rows, bootstrap_rows, stochastic_rows)
    write_tsv(f"{output_dir}/character_model_scores.tsv", model_score_rows, ["layer", "object_id", "model", "parsimony_score", "log_likelihood", "state_count", "observed_tip_count", "fitted_rate", "aic", "bic"])
    write_tsv(f"{output_dir}/model_fit.tsv", model_fit_rows, ["layer", "object_id", "model", "fitted_rate", "log_likelihood", "aic", "bic", "observed_tip_count"])
    write_tsv(f"{output_dir}/hypothesis_tests.tsv", hypothesis_rows, ["layer", "object_id", "test_id", "null_model", "alternative_model", "null_log_likelihood", "alternative_log_likelihood", "lrt_statistic", "df", "p_value", "p_value_method", "q_value", "q_value_method", "fitted_rate", "null_aic", "alternative_aic", "null_bic", "alternative_bic", "observed_tip_count", "tip_error"])
    write_tsv(f"{output_dir}/hypothesis_bootstrap.tsv", bootstrap_rows, ["layer", "object_id", "test_id", "observed_lrt", "bootstrap_replicates", "empirical_p_value", "monte_carlo_se", "null_lrt_mean", "null_lrt_q025", "null_lrt_q500", "null_lrt_q975", "seed", "tip_error"])
    write_tsv(f"{output_dir}/branch_history_posteriors.tsv", stochastic_rows, ["layer", "object_id", "parent_node", "child_node", "parent_label", "child_label", "map_sample_count", "posterior_pr_any_change", "posterior_expected_change_count", "posterior_change_count_low", "posterior_change_count_high", "posterior_most_frequent_transition", "posterior_transition_probability"])
    write_tsv(f"{output_dir}/foreground_tests.tsv", foreground_rows, ["layer", "object_id", "test_id", "null_model", "alternative_model", "null_log_likelihood", "alternative_log_likelihood", "lrt_statistic", "df", "p_value", "p_value_method", "q_value", "q_value_method", "background_rate", "foreground_rate", "rate_ratio", "null_aic", "alternative_aic", "null_bic", "alternative_bic", "observed_tip_count"])
    write_tsv(f"{output_dir}/candidate_structural_events.tsv", event_rows, ["family_id", "event_type", "event_class", "evidence_layer", "object_id", "branch_scope", "event_probability", "ctmc_change_probability", "change", "alternative_explanation"])
    write_tsv(f"{output_dir}/event_support_summary.tsv", support_rows, ["family_id", "event_class", "object_id", "branch_scope", "evidence_layer", "support_tier", "lrt_p_value", "lrt_q_value", "empirical_p_value", "fitted_rate", "ctmc_change_probability", "stochastic_pr_any_change", "evidence_count", "alternative_explanation"])
    write_tsv(f"{output_dir}/model_comparison.tsv", model_rows, ["comparison_id", "model", "score", "delta_vs_best", "interpretation"])
    write_tsv(f"{output_dir}/intragenic_graph_edges.tsv", graph_edges, ["family_id", "species", "gene_copy_id", "edge_id", "left_hsg", "right_hsg", "left_occurrence_id", "right_occurrence_id", "adjacency_status", "left_source_labels", "right_source_labels"])
    write_tsv(f"{output_dir}/case_summary.tsv", case_summary, ["family_id", "hidden_segment_candidates", "source_join_candidates", "multi_source_tip_count", "branch_event_candidates", "best_annotation_model", "best_compound_model"])
    return state_rows, branch_rows, event_rows, model_rows
