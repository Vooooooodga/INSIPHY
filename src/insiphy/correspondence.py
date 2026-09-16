"""Homologous segment group assignment and segment correspondence scoring."""

from collections import Counter, defaultdict

from .elements import (
    collect_element_profiles,
    element_class_for_occurrence,
    element_role_from_occurrences,
    membership_call,
)
from .alignment import global_alignment_stats
from .io import norm_state, parse_fasta, read_tsv, to_float, uniq, write_tsv
from .tree import SpeciesTree


def simple_identity(seq_a, seq_b):
    return global_alignment_stats(seq_a, seq_b).identity


def match_total(row):
    if row.get("total_score") not in ("", "NA", None):
        return to_float(row.get("total_score"))
    boundary = to_float(row.get("boundary_score"), None)
    if boundary is None:
        boundary = 0.5 * (to_float(row.get("left_boundary_score"), 0.5) + to_float(row.get("right_boundary_score"), 0.5))
    return (
        0.40 * to_float(row.get("alignment_score"))
        + 0.15 * to_float(row.get("coverage_score"), to_float(row.get("size_ratio"), 1.0))
        + 0.10 * to_float(row.get("left_context_score"), 0.5)
        + 0.10 * to_float(row.get("right_context_score"), 0.5)
        + 0.10 * boundary
        + 0.10 * to_float(row.get("phase_score"), 0.5)
        + 0.05 * to_float(row.get("order_score"), 0.5)
    )


def normalized_components(row):
    boundary = to_float(row.get("boundary_score"), None)
    if boundary is None:
        boundary = 0.5 * (to_float(row.get("left_boundary_score"), 0.5) + to_float(row.get("right_boundary_score"), 0.5))
    return {
        "alignment_score": to_float(row.get("alignment_score")),
        "coverage_score": to_float(row.get("coverage_score"), to_float(row.get("size_ratio"), 1.0)),
        "left_context_score": to_float(row.get("left_context_score"), 0.5),
        "right_context_score": to_float(row.get("right_context_score"), 0.5),
        "boundary_score": boundary,
        "phase_score": to_float(row.get("phase_score"), 0.5),
        "order_score": to_float(row.get("order_score"), 0.5),
        "strand_score": to_float(row.get("strand_score"), 0.5),
        "splice_score": to_float(row.get("splice_score"), 0.5),
    }


def tree_distances(input_dir):
    rows = read_tsv(f"{input_dir}/species_tree.tsv", optional=True)
    if not rows:
        return {}
    tree = SpeciesTree(rows)
    depth = {tree.root: 0.0}
    for node in tree.preorder():
        for child in tree.children.get(node, []):
            depth[child] = depth[node] + tree.branch_length(child)

    ancestors = {}
    for node in tree.parent:
        vals = []
        cur = node
        while cur:
            vals.append(cur)
            cur = tree.parent.get(cur)
        ancestors[node] = vals

    distances = {}
    labels = sorted(tree.leaf_by_label)
    for left_label in labels:
        for right_label in labels:
            left = tree.leaf_by_label[left_label]
            right = tree.leaf_by_label[right_label]
            right_ancestors = set(ancestors[right])
            lca = next(node for node in ancestors[left] if node in right_ancestors)
            distances[(left_label, right_label)] = depth[left] + depth[right] - 2 * depth[lca]
    return distances


def progressive_tier(distance):
    if distance == "NA":
        return "tree_not_available"
    value = to_float(distance, 0.0)
    if value <= 2.0:
        return "nearest_species_support"
    if value <= 4.0:
        return "within_clade_support"
    return "deep_tree_support"


def reciprocal_status(scored):
    best_subject_by_query = defaultdict(list)
    best_query_by_subject = defaultdict(list)
    for row in scored:
        score = to_float(row["total_score"])
        q = row["query_occurrence_id"]
        s = row["subject_occurrence_id"]
        best_subject_by_query[q].append((score, s, row["match_id"]))
        best_query_by_subject[s].append((score, q, row["match_id"]))
        best_subject_by_query[s].append((score, q, row["match_id"]))
        best_query_by_subject[q].append((score, s, row["match_id"]))

    query_best = {}
    for query, vals in best_subject_by_query.items():
        best = max(score for score, _node, _mid in vals)
        query_best[query] = {_node for score, _node, _mid in vals if abs(score - best) < 1e-12}
    subject_best = {}
    for subject, vals in best_query_by_subject.items():
        best = max(score for score, _node, _mid in vals)
        subject_best[subject] = {_node for score, _node, _mid in vals if abs(score - best) < 1e-12}

    out = {}
    for row in scored:
        q = row["query_occurrence_id"]
        s = row["subject_occurrence_id"]
        if s in query_best.get(q, set()) and q in subject_best.get(s, set()):
            out[row["match_id"]] = "reciprocal_best"
        elif s in query_best.get(q, set()):
            out[row["match_id"]] = "query_best_only"
        else:
            out[row["match_id"]] = "not_reciprocal_best"
    return out


def _tree_depths(tree):
    depths = {tree.root: 0}
    for node in tree.preorder():
        for child in tree.children.get(node, []):
            depths[child] = depths[node] + 1
    return depths


def _ancestor_chain(tree, node):
    chain = []
    cur = node
    while cur:
        chain.append(cur)
        cur = tree.parent.get(cur, "")
    return chain


def _mrca_label(tree, labels):
    nodes = [tree.leaf_by_label[label] for label in labels if label in tree.leaf_by_label]
    if not nodes:
        return "NA"
    depths = _tree_depths(tree)
    common = set(_ancestor_chain(tree, nodes[0]))
    for node in nodes[1:]:
        common &= set(_ancestor_chain(tree, node))
    if not common:
        return "NA"
    node = max(common, key=lambda item: depths.get(item, 0))
    return tree.label.get(node, node)


def _coverage_class(present_count, total_count):
    if present_count == 0:
        return "absent_or_unplaced"
    if present_count == total_count:
        return "tree_spanning"
    if present_count == 1:
        return "tip_specific"
    return "partial"


def summarize_progressive_elements(input_dir, occurrences, element_rows, scored):
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    element_by_occ = {row["occurrence_id"]: row["element_id"] for row in element_rows}
    class_by_occ = {row["occurrence_id"]: row.get("element_class", "context") for row in element_rows}
    distances = tree_distances(input_dir)
    rows_by_element = defaultdict(list)
    homology_by_element = defaultdict(set)
    families_by_element = defaultdict(set)
    for row in element_rows:
        element_id = row["element_id"]
        rows_by_element[element_id].append(row)
        homology_by_element[element_id].add(row.get("homology_id", "NA"))
        if row.get("family_id"):
            families_by_element[element_id].add(row["family_id"])

    pair_scores = defaultdict(list)
    support_by_tier = defaultdict(Counter)
    for row in scored:
        left_element = element_by_occ.get(row["query_occurrence_id"])
        right_element = element_by_occ.get(row["subject_occurrence_id"])
        if not left_element or left_element != right_element:
            continue
        if row.get("correspondence_call") in {"split", "fusion", "deleted", "unmapped", "nonintersecting"}:
            continue
        score = to_float(row.get("total_score"), 0.0)
        left = occ_by_id.get(row["query_occurrence_id"], {})
        right = occ_by_id.get(row["subject_occurrence_id"], {})
        distance = distances.get((left.get("species"), right.get("species")), "NA")
        tier = progressive_tier(f"{distance:.6g}" if distance != "NA" else "NA")
        pair_scores[left_element].append(score)
        support_by_tier[left_element][tier] += 1

    out = []
    for element_id, rows in sorted(rows_by_element.items()):
        species = {row.get("species", "NA") for row in rows if row.get("species")}
        copies = {f"{row.get('species', 'NA')}:{row.get('gene_copy_id', 'NA')}" for row in rows if row.get("gene_copy_id")}
        scores = pair_scores.get(element_id, [])
        classes = Counter(class_by_occ.get(row["occurrence_id"], "context") for row in rows)
        high_near = support_by_tier[element_id]["nearest_species_support"]
        within = support_by_tier[element_id]["within_clade_support"]
        deep = support_by_tier[element_id]["deep_tree_support"]
        if high_near or within:
            call = "progressive_supported"
        elif deep:
            call = "deep_only_support"
        elif len(species) >= 2:
            call = "membership_supported_without_pair_score"
        else:
            call = "single_tip_or_low_support"
        out.append(
            {
                "element_id": element_id,
                "family_id": ";".join(sorted(families_by_element[element_id])) or "NA",
                "homology_ids": ";".join(sorted(homology_by_element[element_id])) or "NA",
                "member_count": len(rows),
                "species_count": len(species),
                "copy_count": len(copies),
                "exon_like_members": classes.get("exon_like", 0),
                "candidate_source_members": classes.get("candidate_source", 0),
                "nearest_species_support": high_near,
                "within_clade_support": within,
                "deep_tree_support": deep,
                "best_pair_score": f"{max(scores):.6g}" if scores else "NA",
                "mean_pair_score": f"{sum(scores) / len(scores):.6g}" if scores else "NA",
                "progressive_call": call,
            }
        )
    return out


def summarize_ancestral_element_graph(input_dir, element_rows):
    tree_rows = read_tsv(f"{input_dir}/species_tree.tsv", optional=True)
    if not tree_rows:
        return []
    tree = SpeciesTree(tree_rows)
    total_tips = len(tree.leaves)
    rows_by_element = defaultdict(list)
    families_by_element = defaultdict(set)
    for row in element_rows:
        rows_by_element[row["element_id"]].append(row)
        if row.get("family_id"):
            families_by_element[row["element_id"]].add(row["family_id"])
    out = []
    for element_id, rows in sorted(rows_by_element.items()):
        species = sorted({row.get("species", "NA") for row in rows if row.get("species") and row.get("element_class") != "absent"})
        copies = sorted({f"{row.get('species', 'NA')}:{row.get('gene_copy_id', 'NA')}" for row in rows if row.get("gene_copy_id")})
        present_count = len([sp for sp in species if sp in tree.leaf_by_label])
        out.append(
            {
                "family_id": ";".join(sorted(families_by_element[element_id])) or "NA",
                "element_id": element_id,
                "mrca_label": _mrca_label(tree, species),
                "present_species_count": present_count,
                "tree_tip_count": total_tips,
                "coverage_class": _coverage_class(present_count, total_tips),
                "present_species": ";".join(species) or "NA",
                "present_copies": ";".join(copies) or "NA",
            }
        )
    return out


def summarize_intragenic_paths(occurrences, element_rows):
    element_by_occ = {row["occurrence_id"]: row["element_id"] for row in element_rows if row.get("element_class") not in {"context", "absent"}}
    grouped = defaultdict(list)
    for row in occurrences:
        grouped[(row.get("family_id", "NA"), row.get("species", "NA"), row.get("gene_copy_id", "NA"))].append(row)
    out = []
    path_counts = defaultdict(Counter)
    for (family, species, copy), rows in sorted(grouped.items()):
        rows = sorted(rows, key=lambda row: (row.get("contig", ""), int(row.get("start", "0")), int(row.get("end", "0"))))
        path = []
        context_count = 0
        for row in rows:
            if norm_state(row.get("presence_status")) != "present":
                continue
            element_id = element_by_occ.get(row["occurrence_id"])
            if element_id:
                path.append(element_id)
            else:
                context_count += 1
        element_path = ">".join(path) if path else "NA"
        path_counts[family][element_path] += 1
        out.append(
            {
                "family_id": family,
                "path_scope": "extant_copy",
                "node_label": f"{species}:{copy}",
                "species": species,
                "gene_copy_id": copy,
                "path_type": "observed",
                "element_path": element_path,
                "context_count": context_count,
                "path_support": "1",
            }
        )
    for family, counts in sorted(path_counts.items()):
        if not counts:
            continue
        path, count = counts.most_common(1)[0]
        out.append(
            {
                "family_id": family,
                "path_scope": "family_consensus",
                "node_label": family,
                "species": "NA",
                "gene_copy_id": "NA",
                "path_type": "majority_observed_path",
                "element_path": path,
                "context_count": "NA",
                "path_support": f"{count}/{sum(counts.values())}",
            }
        )
    return out


def infer_correspondence(input_dir, output_dir):
    homology = read_tsv(f"{input_dir}/segment_homology.tsv", ["homology_id", "occurrence_id", "support_type", "confidence"])
    occurrences = read_tsv(f"{input_dir}/segment_occurrences.tsv", ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"])
    matches = read_tsv(f"{input_dir}/segment_matches.tsv", ["match_id", "query_occurrence_id", "subject_occurrence_id", "match_status"], optional=True)
    evidence_rows = read_tsv(f"{input_dir}/sequence_synteny_evidence.tsv", optional=True)
    annotation_rows = read_tsv(f"{output_dir}/annotation_completion_candidates.tsv", optional=True)
    seqs = parse_fasta(f"{input_dir}/segment_sequences.fasta")
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    distances = tree_distances(input_dir)
    profiles, element_by_homology = collect_element_profiles(homology, occurrences, evidence_rows + annotation_rows)

    component_rows = []
    element_rows = []
    by_component = defaultdict(list)
    for row in homology:
        occ = occ_by_id.get(row["occurrence_id"], {})
        by_component[row["homology_id"]].append(row["occurrence_id"])
        score, call = membership_call(row.get("confidence"), 0)
        component_rows.append(
            {
                "homology_id": row["homology_id"],
                "occurrence_id": row["occurrence_id"],
                "family_id": occ.get("family_id", "NA"),
                "species": occ.get("species", "NA"),
                "gene_copy_id": occ.get("gene_copy_id", "NA"),
                "source_label": row.get("source_label", "NA"),
                "support_type": row.get("support_type", "NA"),
                "confidence": row.get("confidence", "NA"),
                "membership_score": f"{score:.6g}",
                "membership_call": call,
            }
        )
        element_id = element_by_homology.get(row["homology_id"])
        if element_id:
            element_class = element_class_for_occurrence(occ, row["homology_id"], profiles, element_by_homology)
            element_rows.append(
                {
                    "element_id": element_id,
                    "family_id": occ.get("family_id", "NA"),
                    "homology_id": row["homology_id"],
                    "occurrence_id": row["occurrence_id"],
                    "species": occ.get("species", "NA"),
                    "gene_copy_id": occ.get("gene_copy_id", "NA"),
                    "element_class": element_class,
                    "display_role": element_role_from_occurrences([occ]),
                    "source_label": row.get("source_label", "NA"),
                    "support_type": row.get("support_type", "NA"),
                    "confidence": row.get("confidence", "NA"),
                    "membership_score": f"{score:.6g}",
                    "membership_call": call,
                }
            )

    conservation = []
    for component_id, occ_ids in sorted(by_component.items()):
        ids_with_seq = [occ_id for occ_id in occ_ids if occ_id in seqs]
        identities = []
        for i, left in enumerate(ids_with_seq):
            for right in ids_with_seq[i + 1 :]:
                identities.append(simple_identity(seqs[left], seqs[right]))
        mean_identity = sum(identities) / len(identities) if identities else 0.0
        roles = Counter(occ_by_id.get(occ_id, {}).get("role", "unknown") for occ_id in occ_ids)
        conservation.append(
            {
                "homology_id": component_id,
                "occurrence_count": len(occ_ids),
                "species_count": uniq(occ_by_id.get(occ_id, {}).get("species", "") for occ_id in occ_ids),
                "mean_pairwise_identity": f"{mean_identity:.6g}",
                "role_spectrum": ";".join(f"{key}:{roles[key]}" for key in sorted(roles)),
                "conservation_call": "conserved" if mean_identity >= 0.75 or len(occ_ids) >= 2 else "single_or_low_sequence_support",
            }
        )

    scored = []
    grouped = defaultdict(list)
    for row in matches:
        score = match_total(row)
        components = normalized_components(row)
        grouped[row["query_occurrence_id"]].append((score, row["match_id"]))
        scored.append(
            {
                "match_id": row["match_id"],
                "query_occurrence_id": row["query_occurrence_id"],
                "subject_occurrence_id": row["subject_occurrence_id"],
                "alignment_score": f"{components['alignment_score']:.6g}",
                "coverage_score": f"{components['coverage_score']:.6g}",
                "left_context_score": f"{components['left_context_score']:.6g}",
                "right_context_score": f"{components['right_context_score']:.6g}",
                "boundary_score": f"{components['boundary_score']:.6g}",
                "phase_score": f"{components['phase_score']:.6g}",
                "order_score": f"{components['order_score']:.6g}",
                "strand_score": f"{components['strand_score']:.6g}",
                "splice_score": f"{components['splice_score']:.6g}",
                "total_score": f"{score:.6g}",
                "match_status": row.get("match_status", "ambiguous"),
                "distance_class": row.get("distance_class", "NA"),
                "threshold": row.get("threshold", "NA"),
                "alignment_cigar": row.get("alignment_cigar", "NA"),
                "reciprocal_status": "unclassified",
                "correspondence_call": "unclassified",
            }
        )
    best_status = {}
    for vals in grouped.values():
        best = max(score for score, _ in vals)
        best_ids = [match_id for score, match_id in vals if abs(score - best) < 1e-12]
        for match_id in best_ids:
            best_status[match_id] = "best_tie" if len(best_ids) > 1 else "best_unique"
    reciprocal = reciprocal_status(scored)
    graph_edges = []
    degree = Counter()
    for row in scored:
        row["reciprocal_status"] = reciprocal.get(row["match_id"], "not_reciprocal_best")
        if row["match_status"] == "mapped":
            degree[row["query_occurrence_id"]] += 1
            degree[row["subject_occurrence_id"]] += 1
            graph_edges.append(
                {
                    "edge_id": row["match_id"],
                    "query_occurrence_id": row["query_occurrence_id"],
                    "subject_occurrence_id": row["subject_occurrence_id"],
                    "total_score": row["total_score"],
                    "reciprocal_status": row["reciprocal_status"],
                    "edge_call": "high_confidence_correspondence" if row["reciprocal_status"] == "reciprocal_best" and to_float(row["total_score"]) >= 0.7 else "supporting_correspondence",
                }
            )
        if row["match_status"] != "mapped":
            row["correspondence_call"] = row["match_status"]
        elif row["reciprocal_status"] == "reciprocal_best":
            row["correspondence_call"] = "reciprocal_best"
        else:
            row["correspondence_call"] = best_status.get(row["match_id"], "not_best")

    progressive_rows = []
    for row in scored:
        left = occ_by_id.get(row["query_occurrence_id"], {})
        right = occ_by_id.get(row["subject_occurrence_id"], {})
        pair_distance = distances.get((left.get("species"), right.get("species")), "NA")
        pair_distance_text = f"{pair_distance:.6g}" if pair_distance != "NA" else "NA"
        progressive_rows.append(
            {
                "match_id": row["match_id"],
                "query_species": left.get("species", "NA"),
                "subject_species": right.get("species", "NA"),
                "query_gene_copy_id": left.get("gene_copy_id", "NA"),
                "subject_gene_copy_id": right.get("gene_copy_id", "NA"),
                "tree_distance": pair_distance_text,
                "progressive_tier": progressive_tier(pair_distance_text),
                "total_score": row["total_score"],
                "correspondence_call": row["correspondence_call"],
            }
        )
    progressive_element_rows = summarize_progressive_elements(input_dir, occurrences, element_rows, scored)
    ancestral_graph_rows = summarize_ancestral_element_graph(input_dir, element_rows)
    path_rows = summarize_intragenic_paths(occurrences, element_rows)

    for row in component_rows:
        deg = degree.get(row["occurrence_id"], 0)
        adjusted, call = membership_call(row.get("confidence"), deg)
        row["membership_score"] = f"{adjusted:.6g}"
        row["membership_call"] = call
    element_by_occ = {row["occurrence_id"]: row for row in element_rows}
    for occ_id, row in element_by_occ.items():
        adjusted, call = membership_call(row.get("confidence"), degree.get(occ_id, 0))
        row["membership_score"] = f"{adjusted:.6g}"
        row["membership_call"] = call

    write_tsv(f"{output_dir}/internal_homology_assignments.tsv", component_rows, ["homology_id", "occurrence_id", "family_id", "species", "gene_copy_id", "source_label", "support_type", "confidence", "membership_score", "membership_call"])
    write_tsv(f"{output_dir}/element_correspondence.tsv", element_rows, ["element_id", "family_id", "homology_id", "occurrence_id", "species", "gene_copy_id", "element_class", "display_role", "source_label", "support_type", "confidence", "membership_score", "membership_call"])
    write_tsv(f"{output_dir}/segment_conservation.tsv", conservation, ["homology_id", "occurrence_count", "species_count", "mean_pairwise_identity", "role_spectrum", "conservation_call"])
    write_tsv(
        f"{output_dir}/segment_correspondence.tsv",
        scored,
        [
            "match_id",
            "query_occurrence_id",
            "subject_occurrence_id",
            "alignment_score",
            "coverage_score",
            "left_context_score",
            "right_context_score",
            "boundary_score",
            "phase_score",
            "order_score",
            "strand_score",
            "splice_score",
            "total_score",
            "distance_class",
            "threshold",
            "alignment_cigar",
            "match_status",
            "reciprocal_status",
            "correspondence_call",
        ],
    )
    write_tsv(f"{output_dir}/internal_homology_graph_edges.tsv", graph_edges, ["edge_id", "query_occurrence_id", "subject_occurrence_id", "total_score", "reciprocal_status", "edge_call"])
    write_tsv(f"{output_dir}/progressive_correspondence.tsv", progressive_rows, ["match_id", "query_species", "subject_species", "query_gene_copy_id", "subject_gene_copy_id", "tree_distance", "progressive_tier", "total_score", "correspondence_call"])
    write_tsv(f"{output_dir}/progressive_element_correspondence.tsv", progressive_element_rows, ["element_id", "family_id", "homology_ids", "member_count", "species_count", "copy_count", "exon_like_members", "candidate_source_members", "nearest_species_support", "within_clade_support", "deep_tree_support", "best_pair_score", "mean_pair_score", "progressive_call"])
    write_tsv(f"{output_dir}/ancestral_element_graph.tsv", ancestral_graph_rows, ["family_id", "element_id", "mrca_label", "present_species_count", "tree_tip_count", "coverage_class", "present_species", "present_copies"])
    write_tsv(f"{output_dir}/ancestral_intragenic_paths.tsv", path_rows, ["family_id", "path_scope", "node_label", "species", "gene_copy_id", "path_type", "element_path", "context_count", "path_support"])
    return component_rows, conservation, scored
