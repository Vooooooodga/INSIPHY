"""Fixed-tree structural inference for duplicated/chimeric gene demos."""

from collections import defaultdict

from .io import norm_state, read_tsv, to_float, write_tsv
from .tree import SpeciesTree, sankoff


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


def add_character(tree, layer, object_id, tips, states, state_rows, branch_rows):
    score, probs, edges = sankoff(tree, tips, states, layer)
    for row in probs:
        state_rows.append({"layer": layer, "object_id": object_id, "score": f"{score:.6g}", **row})
    for row in edges:
        event_type = f"{layer}_change" if row["status"] == "change_required" else "state_change"
        branch_rows.append({"layer": layer, "object_id": object_id, "event_type": event_type, **row})


def infer_phylogeny(input_dir, output_dir):
    occurrences = read_tsv(f"{input_dir}/segment_occurrences.tsv", ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"])
    homology = read_tsv(f"{input_dir}/segment_homology.tsv", ["homology_id", "occurrence_id", "support_type", "confidence"])
    adjacencies = read_tsv(f"{input_dir}/physical_adjacencies.tsv", ["adjacency_id", "family_id", "species", "gene_copy_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status"])
    tree = SpeciesTree(read_tsv(f"{input_dir}/species_tree.tsv", ["node_id", "parent_id", "label"]))
    copy_context = read_tsv(f"{input_dir}/copy_context.tsv", ["family_id", "species", "gene_copy_id", "copy_class"], optional=True)
    annotation_candidates = read_tsv(f"{output_dir}/annotation_completion_candidates.tsv", ["evidence_id", "completion_call"], optional=True)

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
        add_character(tree, "segment_presence", hsg, {sp: state_from_occurrences(rows) for sp, rows in by_species.items()}, {"present", "absent"}, state_rows, branch_rows)
        add_character(tree, "role_state", hsg, {sp: role_from_occurrences(rows) for sp, rows in by_species.items()}, {"absent", "CDS", "exon_or_UTR", "intron_or_noncoding", "unknown"}, state_rows, branch_rows)

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
        add_character(tree, "adjacency_state", pair_id, tips, {"present", "absent"}, state_rows, branch_rows)

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
        add_character(tree, "source_mixture", family, tips, {"single_source", "multi_source", "unknown"}, state_rows, branch_rows)

    copy_states_by_family_species = defaultdict(list)
    for row in copy_context:
        copy_states_by_family_species[(row["family_id"], row["species"])].append(row.get("copy_class", "unknown"))
    copy_tips = defaultdict(dict)
    for (family, species), states in copy_states_by_family_species.items():
        if "tandem_multi_copy" in states:
            copy_tips[family][species] = "tandem_multi_copy"
        elif "unresolved" in states:
            copy_tips[family][species] = "unresolved"
        elif "single_copy" in states:
            copy_tips[family][species] = "single_copy"
        else:
            copy_tips[family][species] = "unknown"
    for family, tips in sorted(copy_tips.items()):
        object_family[family] = family
        add_character(tree, "copy_multiplicity", family, tips, {"single_copy", "tandem_multi_copy", "unresolved", "unknown"}, state_rows, branch_rows)

    hidden_support = sum(to_float(row.get("support_score")) for row in annotation_candidates if row.get("completion_call") == "hidden_segment_candidate")
    hidden_count = sum(1 for row in annotation_candidates if row.get("completion_call") == "hidden_segment_candidate")
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
    for row in branch_rows:
        if row["status"] == "change_required":
            event_rows.append(
                {
                    "family_id": object_family.get(row["object_id"], row["object_id"].split("__")[0]),
                    "event_type": row["event_type"],
                    "object_id": row["object_id"],
                    "branch_scope": f"{row['parent_label']}->{row['child_label']}",
                    "event_probability": row["event_probability"],
                    "change": row["change"],
                    "alternative_explanation": "annotation_gap_or_mapping_ambiguity_if_sequence_support_low",
                }
            )
    for edge in graph_edges:
        if edge["adjacency_status"] == "present" and edge["left_source_labels"] != "NA" and edge["right_source_labels"] != "NA" and edge["left_source_labels"] != edge["right_source_labels"]:
            event_rows.append(
                {
                    "family_id": edge["family_id"],
                    "event_type": "source_join_candidate",
                    "object_id": edge["edge_id"],
                    "branch_scope": "estimated_from_adjacency_state_history",
                    "event_probability": "NA",
                    "change": f"{edge['left_source_labels']}->{edge['right_source_labels']}",
                    "alternative_explanation": "paralogy_or_homology_assignment_error_if_low_support",
                }
            )

    demo_summary = []
    for family in sorted({row["family_id"] for row in occurrences}):
        demo_summary.append(
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
    write_tsv(f"{output_dir}/branch_event_probabilities.tsv", branch_rows, ["layer", "object_id", "event_type", "parent_node", "child_node", "parent_label", "child_label", "status", "change", "event_probability"])
    write_tsv(f"{output_dir}/candidate_structural_events.tsv", event_rows, ["family_id", "event_type", "object_id", "branch_scope", "event_probability", "change", "alternative_explanation"])
    write_tsv(f"{output_dir}/model_comparison.tsv", model_rows, ["comparison_id", "model", "score", "delta_vs_best", "interpretation"])
    write_tsv(f"{output_dir}/intragenic_graph_edges.tsv", graph_edges, ["family_id", "species", "gene_copy_id", "edge_id", "left_hsg", "right_hsg", "left_occurrence_id", "right_occurrence_id", "adjacency_status", "left_source_labels", "right_source_labels"])
    write_tsv(f"{output_dir}/demo_summary.tsv", demo_summary, ["family_id", "hidden_segment_candidates", "source_join_candidates", "multi_source_tip_count", "branch_event_candidates", "best_annotation_model", "best_compound_model"])
    return state_rows, branch_rows, event_rows, model_rows
