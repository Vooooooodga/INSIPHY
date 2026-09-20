"""experimental / reconstruction: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.elements import element_id_from_homology
from intraphy.experimental.character_reconstruction import add_character
from intraphy.experimental.event_labels import add_q_values
from intraphy.experimental.event_labels import best_branch_support
from intraphy.experimental.event_labels import classify_branch_event
from intraphy.experimental.event_labels import completion_event_class
from intraphy.experimental.event_labels import event_support_summary
from intraphy.experimental.event_labels import informative_source_labels
from intraphy.experimental.event_labels import interpretation_hints
from intraphy.experimental.event_labels import make_event
from intraphy.experimental.event_labels import role_from_occurrences
from intraphy.experimental.event_labels import state_from_occurrences
from intraphy.experimental.tree_coverage import element_phylogenetic_coverage
from intraphy.experimental.tree_coverage import fallback_element_rows
from intraphy.experimental.tree_coverage import has_multicopy_species
from intraphy.experimental.tree_coverage import internal_homology_phylogenetic_coverage
from intraphy.experimental.tree_coverage import read_foreground_edges
from intraphy.experimental.tree_coverage import read_structural_tree
from intraphy.experimental.tree_coverage import tree_tip_label_for_adjacency
from intraphy.experimental.tree_coverage import tree_tip_label_for_occ
from intraphy.experimental.tree_model import SpeciesTree
from intraphy.io import norm_state
from intraphy.io import read_tsv
from intraphy.io import to_float
from intraphy.io import write_tsv


def infer_experimental_phylogeny(
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
    structural_site_matrix_path=None,
    annotation_view="repertoire",
):
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
