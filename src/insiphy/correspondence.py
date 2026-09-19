"""Homologous segment group assignment and segment correspondence scoring."""

import json
from collections import Counter, defaultdict

from .coordinates import ClosedInterval1, Interval0, parse_legacy_blocks
from .elements import (
    collect_element_profiles,
    element_class_for_occurrence,
    element_role_from_occurrences,
    membership_call,
)
from .alignment import global_alignment_stats
from .io import norm_state, parse_fasta, read_tsv, to_float, uniq, write_tsv
from .tree import SpeciesTree


CORRESPONDENCE_EVIDENCE_FIELDS = (
    "matched_blocks",
    "query_genomic_matched_blocks",
    "subject_genomic_matched_blocks",
    "query_parent_feature_ids",
    "subject_parent_feature_ids",
    "query_transcript_ids",
    "subject_transcript_ids",
    "candidate_id",
    "alternative_candidate_ids",
    "aligned_blocks",
    "gap_blocks",
    "sequence_kind",
    "backend",
    "backend_version",
    "raw_score",
    "nt_identity",
    "aa_identity",
    "known_aligned_pairs",
    "unknown_aligned_pairs",
    "query_covered_bases",
    "target_covered_bases",
    "query_length",
    "target_length",
    "relative_strand",
    "is_secondary",
    "left_anchor_id",
    "right_anchor_id",
    "search_interval",
    "search_interval_side",
    "short_sequence_coverage",
    "alignment_input_transposed",
    "mapq",
    "mapping_quality",
    "hit_count",
    "ambiguous_hit_count",
    "alternative_hits",
    "candidate_assessments",
    "dna_candidate_assessments",
    "candidate_ids",
    "score_scheme",
    "raw_alignment_score",
    "enumeration_complete",
    "candidate_enumeration_status",
    "incomplete_reason",
    "short_context_route",
    "local_boundary_range",
    "candidate_resolution",
    "retained_candidate_ids",
    "best_path_candidate_ids",
    "chain_best_score",
    "chain_score_delta",
    "chain_configuration",
    "chain_delta_rule",
    "chain_local_mode",
    "chain_status",
    "chain_ambiguity",
    "chain_start_anchor_ids",
    "chain_end_anchor_ids",
    "chain_retained_edges",
    "chain_best_path_count_capped",
    "chain_near_optimal_path_count_capped",
    "flanking_anchor_status",
    "membership_edge_eligible",
    "membership_edge_reason",
    "position_edge_eligible",
    "position_edge_reason",
    "true_absence_eligible",
    "true_absence_evidence_status",
    "true_absence_reason",
    "protein_mapping_status",
    "protein_candidate_evidence_available",
    "protein_membership_eligible",
    "protein_position_eligible",
    "protein_hard_observation_eligible",
    "protein_known_aa_pairs",
    "protein_blosum62_score",
    "protein_gap_fraction",
    "protein_left_anchor_pairs",
    "protein_right_anchor_pairs",
    "protein_left_anchor_score",
    "protein_right_anchor_score",
    "protein_left_anchor_supported",
    "protein_right_anchor_supported",
    "protein_terminal_side",
    "protein_msa_column_start0",
    "protein_msa_column_end0",
    "protein_msa_column_interval",
    "protein_msa_mode",
    "protein_candidate_mapping_count",
    "protein_candidate_coordinate_consensus",
    "protein_candidate_details",
    "protein_competing_occurrences",
    "protein_query_source_features",
    "protein_target_source_features",
)


def simple_identity(seq_a, seq_b):
    return global_alignment_stats(seq_a, seq_b, backend="auto").identity


def match_total(row):
    if row.get("correspondence_score") not in ("", "NA", None):
        return to_float(row["correspondence_score"])
    if row.get("total_score") not in ("", "NA", None):
        return to_float(row.get("total_score"))
    boundary = to_float(row.get("boundary_score"), None)
    if boundary is None:
        boundary = 0.5 * (to_float(row.get("left_boundary_score"), 0.5) + to_float(row.get("right_boundary_score"), 0.5))
    return (
        0.34 * to_float(row.get("alignment_score"))
        + 0.14 * to_float(row.get("coverage_score"), to_float(row.get("size_ratio"), 1.0))
        + 0.10 * to_float(row.get("left_context_score"), 0.5)
        + 0.10 * to_float(row.get("right_context_score"), 0.5)
        + 0.10 * boundary
        + 0.08 * to_float(row.get("phase_score"), 0.5)
        + 0.06 * to_float(row.get("order_score"), 0.5)
        + 0.04 * to_float(row.get("strand_score"), 0.5)
        + 0.04 * to_float(row.get("splice_score"), 0.5)
    )


def normalized_components(row):
    boundary = to_float(row.get("boundary_score"), None)
    if boundary is None:
        boundary = 0.5 * (to_float(row.get("left_boundary_score"), 0.5) + to_float(row.get("right_boundary_score"), 0.5))
    return {
        "alignment_score": to_float(row.get("alignment_score")),
        "coverage_score": to_float(row.get("coverage_score"), to_float(row.get("size_ratio"), 1.0)),
        "sequence_score": to_float(row.get("sequence_score"), None),
        "structural_context_score": to_float(row.get("structural_context_score"), None),
        "left_context_score": to_float(row.get("left_context_score"), 0.5),
        "right_context_score": to_float(row.get("right_context_score"), 0.5),
        "boundary_score": boundary,
        "phase_score": to_float(row.get("phase_score"), 0.5),
        "order_score": to_float(row.get("order_score"), 0.5),
        "strand_score": to_float(row.get("strand_score"), 0.5),
        "splice_score": to_float(row.get("splice_score"), 0.5),
    }


def tree_distances(input_dir):
    """Prioritize correspondence by supplied lengths, or topology if any are missing.

    Unit edges apply to this ordering only; CTMC branch lengths are unchanged.
    """
    rows = read_tsv(f"{input_dir}/species_tree.tsv", optional=True)
    if not rows:
        return {}
    tree = SpeciesTree(rows)
    use_topology = any(tree.length[node] is None for node in tree.parent if node != tree.root)
    depth = {tree.root: 0.0}
    for node in tree.preorder():
        for child in tree.children.get(node, []):
            depth[child] = depth[node] + (1.0 if use_topology else tree.branch_length(child))

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
    elements_by_occ = defaultdict(set)
    for row in element_rows:
        elements_by_occ[row["occurrence_id"]].add(row["element_id"])
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
        common_elements = (
            elements_by_occ.get(row["query_occurrence_id"], set())
            & elements_by_occ.get(row["subject_occurrence_id"], set())
        )
        if not common_elements:
            continue
        if row.get("correspondence_call") in {"split", "fusion", "deleted", "unmapped", "nonintersecting"}:
            continue
        score = to_float(row.get("total_score"), 0.0)
        left = occ_by_id.get(row["query_occurrence_id"], {})
        right = occ_by_id.get(row["subject_occurrence_id"], {})
        distance = distances.get((left.get("species"), right.get("species")), "NA")
        tier = progressive_tier(f"{distance:.6g}" if distance != "NA" else "NA")
        for element_id in common_elements:
            pair_scores[element_id].append(score)
            support_by_tier[element_id][tier] += 1

    out = []
    for element_id, rows in sorted(rows_by_element.items()):
        species = {row.get("species", "NA") for row in rows if row.get("species")}
        copies = {f"{row.get('species', 'NA')}:{row.get('gene_copy_id', 'NA')}" for row in rows if row.get("gene_copy_id")}
        scores = pair_scores.get(element_id, [])
        classes = Counter(row.get("element_class", "context") for row in rows)
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


def summarize_observed_element_tree_coverage(input_dir, element_rows):
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
    memberships_by_occ = defaultdict(list)
    for row in element_rows:
        if row.get("element_class") not in {"context", "absent"}:
            memberships_by_occ[row["occurrence_id"]].append(row)
    grouped = defaultdict(list)
    for row in occurrences:
        grouped[(row.get("family_id", "NA"), row.get("species", "NA"), row.get("gene_copy_id", "NA"))].append(row)
    out = []
    path_counts = defaultdict(Counter)
    for (family, species, copy), rows in sorted(grouped.items()):
        rows = sorted(rows, key=lambda row: (
            row.get("contig", ""),
            -int(row.get("end", "0"))
            if row.get("strand") == "-"
            else int(row.get("start", "0")),
            row.get("occurrence_id", ""),
        ))
        path = []
        context_count = 0
        for row in rows:
            if norm_state(row.get("presence_status")) != "present":
                continue
            memberships = memberships_by_occ.get(row["occurrence_id"], ())
            if memberships:
                def membership_coordinate(membership):
                    blocks = _parse_genomic_blocks(membership.get("genomic_matched_blocks"))
                    starts = [interval.start0 for _contig, _strand, interval in blocks]
                    if not starts:
                        return int(row.get("start", 0)) - 1
                    return min(starts) if row.get("strand") != "-" else -max(starts)

                path.extend(
                    membership["element_id"]
                    for membership in sorted(
                        memberships,
                        key=lambda membership: (
                            membership_coordinate(membership), membership["element_id"],
                        ),
                    )
                )
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


def _tokens(value):
    return {
        token
        for token in str(value or "").replace(",", ";").split(";")
        if token and token != "NA"
    }


def _edge_eligible(match, field, legacy_default=False):
    value = match.get(field)
    if value in {1, "1", True}:
        return True
    if value in {0, "0", False}:
        return False
    return bool(legacy_default)


def _hard_membership_eligible(match):
    legacy_resolved = (
        match.get("membership_edge_eligible") in {None, "", "NA"}
        and match.get("match_status") == "mapped"
        and match.get("candidate_resolution", "resolved") == "resolved"
    )
    if not _edge_eligible(match, "membership_edge_eligible", legacy_resolved):
        return False
    if "annotated_CDS_protein" not in str(match.get("correspondence_basis", "")):
        return True
    return _edge_eligible(match, "protein_hard_observation_eligible", False)


def _hard_position_eligible(match):
    legacy_resolved = (
        match.get("position_edge_eligible") in {None, "", "NA"}
        and match.get("match_status") == "mapped"
        and match.get("candidate_resolution", "resolved") == "resolved"
    )
    return (
        _hard_membership_eligible(match)
        and _edge_eligible(match, "position_edge_eligible", legacy_resolved)
    )


def _genomic_block_records(blocks):
    records = []
    for index, token in enumerate(sorted(blocks), start=1):
        try:
            contig, bounds, strand = token.rsplit(":", 2)
            start, end = (int(value) for value in bounds.split("-", 1))
        except (TypeError, ValueError):
            continue
        records.append(
            {
                "block_id": f"membership_block_{index:03d}",
                "target_contig": contig,
                "target_start": min(start, end),
                "target_end": max(start, end),
                "target_strand": strand,
            }
        )
    return records


def _parse_genomic_blocks(value):
    records = []
    for token in str(value or "").split(";"):
        if not token or token == "NA":
            continue
        try:
            contig, bounds, strand = token.rsplit(":", 2)
            start, end = (int(item) for item in bounds.split("-", 1))
            interval = ClosedInterval1(
                min(start, end), max(start, end),
            ).to_interval0()
        except (TypeError, ValueError):
            continue
        records.append((contig, strand, interval))
    return tuple(records)


def _blocks_overlap(left, right):
    return any(
        left_contig == right_contig
        and left_strand == right_strand
        and left_interval.overlaps(right_interval)
        for left_contig, left_strand, left_interval in left
        for right_contig, right_strand, right_interval in right
    )


def _block_union_length(blocks):
    grouped = defaultdict(list)
    for contig, strand, interval in blocks:
        grouped[(contig, strand)].append(interval)
    total = 0
    for intervals in grouped.values():
        intervals = sorted(intervals)
        start0, end0 = intervals[0].start0, intervals[0].end0
        for interval in intervals[1:]:
            if interval.start0 <= end0:
                end0 = max(end0, interval.end0)
            else:
                total += end0 - start0
                start0, end0 = interval.start0, interval.end0
        total += end0 - start0
    return total


def _tagged_membership_blocks(nodes, node_info):
    """Preserve each alignment block and the occurrence side it supports."""

    blocks_by_match = defaultdict(set)
    for node in nodes:
        match = node_info[node]["match"]
        side = node_info[node]["side"]
        try:
            blocks = parse_legacy_blocks(match.get("matched_blocks"))
        except (TypeError, ValueError):
            continue
        for block in blocks:
            query = ClosedInterval1.from_interval0(block.query)
            target = ClosedInterval1.from_interval0(block.target)
            blocks_by_match[(match["match_id"], side)].add(
                (query.start, query.end, target.start, target.end)
            )
    tokens = []
    for (match_id, side), blocks in sorted(blocks_by_match.items()):
        for index, (query_start, query_end, target_start, target_end) in enumerate(
            sorted(blocks)
        ):
            prefix = f"{match_id}:{side}:" if index == 0 else ""
            tokens.append(
                f"{prefix}{query_start}-{query_end}:"
                f"{target_start}-{target_end}"
            )
    return ";".join(tokens) or "NA"


def _candidate_block_evidence(match, side):
    field = (
        "query_genomic_matched_blocks"
        if side == "query"
        else "subject_genomic_matched_blocks"
    )
    evidence = []
    retained_ids = _tokens(match.get("retained_candidate_ids"))
    for encoded_field in ("candidate_assessments", "dna_candidate_assessments"):
        encoded = match.get(encoded_field)
        if encoded is None or encoded == "" or encoded == "NA":
            continue
        try:
            records = json.loads(encoded) if isinstance(encoded, str) else encoded
        except (TypeError, ValueError):
            continue
        if not isinstance(records, list):
            continue
        for record in records:
            if not isinstance(record, dict):
                continue
            if record.get("accepted", 1) not in {1, "1", True}:
                continue
            record_ids = _tokens(record.get("candidate_id"))
            if retained_ids and record_ids.isdisjoint(retained_ids):
                continue
            blocks = _parse_genomic_blocks(record.get(field))
            if not blocks:
                continue
            evidence.append((blocks, record_ids))
    if not evidence:
        top_blocks = _parse_genomic_blocks(match.get(field))
        if top_blocks:
            evidence.append((top_blocks, _tokens(match.get("candidate_ids"))))
    distinct = {}
    for blocks, candidate_ids in evidence:
        signature = tuple(
            (contig, strand, interval.start0, interval.end0)
            for contig, strand, interval in blocks
        )
        distinct.setdefault(signature, [blocks, set()])[1].update(candidate_ids)
    return tuple(
        (blocks, candidate_ids)
        for blocks, candidate_ids in distinct.values()
    )


def _membership_match_details(element_rows, occurrences, scored):
    """Replace whole-feature membership rows with resolved local memberships."""

    occurrence_by_id = {row["occurrence_id"]: row for row in occurrences}
    templates_by_occurrence = defaultdict(dict)
    for row in element_rows:
        templates_by_occurrence[row["occurrence_id"]].setdefault(
            row["element_id"], dict(row),
        )
    base_elements_by_occurrence = {
        occurrence_id: set(templates)
        for occurrence_id, templates in templates_by_occurrence.items()
    }
    parent = {}

    def find(item):
        parent.setdefault(item, item)
        if parent[item] != item:
            parent[item] = find(parent[item])
        return parent[item]

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    node_info = {}
    nodes_by_occurrence = defaultdict(list)
    matches_by_occurrence = defaultdict(list)
    for match in scored:
        query_id = match.get("query_occurrence_id")
        subject_id = match.get("subject_occurrence_id")
        matches_by_occurrence[query_id].append((match, "query"))
        matches_by_occurrence[subject_id].append((match, "subject"))
        if not _hard_position_eligible(match):
            continue
        shared_base_elements = (
            base_elements_by_occurrence.get(query_id, set())
            & base_elements_by_occurrence.get(subject_id, set())
        )
        if not shared_base_elements:
            continue
        query_blocks = _parse_genomic_blocks(match.get("query_genomic_matched_blocks"))
        subject_blocks = _parse_genomic_blocks(match.get("subject_genomic_matched_blocks"))
        if not query_blocks or not subject_blocks:
            continue
        for base_element in sorted(shared_base_elements):
            query_node = (match["match_id"], "query", base_element)
            subject_node = (match["match_id"], "subject", base_element)
            node_info[query_node] = {
                "occurrence_id": query_id,
                "base_element": base_element,
                "blocks": query_blocks,
                "match": match,
                "side": "query",
            }
            node_info[subject_node] = {
                "occurrence_id": subject_id,
                "base_element": base_element,
                "blocks": subject_blocks,
                "match": match,
                "side": "subject",
            }
            nodes_by_occurrence[(query_id, base_element)].append(query_node)
            nodes_by_occurrence[(subject_id, base_element)].append(subject_node)
            union(query_node, subject_node)

    for nodes in nodes_by_occurrence.values():
        for index, left in enumerate(nodes):
            for right in nodes[index + 1:]:
                if _blocks_overlap(node_info[left]["blocks"], node_info[right]["blocks"]):
                    union(left, right)

    components = defaultdict(list)
    for node in node_info:
        components[find(node)].append(node)
    components_by_base = defaultdict(list)
    for nodes in components.values():
        bases = {
            node_info[node]["base_element"]
            for node in nodes
        } - {None}
        if len(bases) == 1:
            components_by_base[next(iter(bases))].append(nodes)

    def component_key(nodes):
        coordinates = [
            (contig, interval.start0, interval.end0, node_info[node]["occurrence_id"])
            for node in nodes
            for contig, _strand, interval in node_info[node]["blocks"]
        ]
        return min(coordinates) if coordinates else ("", 0, 0, "")

    localized_rows = []
    rows_by_occurrence = defaultdict(list)
    for base_element, base_components in sorted(components_by_base.items()):
        ordered_components = sorted(base_components, key=component_key)
        for component_index, nodes in enumerate(ordered_components, start=1):
            element_id = (
                base_element if len(ordered_components) == 1
                else f"{base_element}.local_{component_index:03d}"
            )
            by_occurrence = defaultdict(list)
            for node in nodes:
                by_occurrence[node_info[node]["occurrence_id"]].append(node)
            for occurrence_index, (occurrence_id, occurrence_nodes) in enumerate(
                sorted(by_occurrence.items()), start=1,
            ):
                template = dict(
                    templates_by_occurrence[occurrence_id][base_element]
                )
                occurrence = occurrence_by_id.get(occurrence_id, {})
                blocks = sorted({
                    (contig, strand, interval)
                    for node in occurrence_nodes
                    for contig, strand, interval in node_info[node]["blocks"]
                }, key=lambda item: (item[0], item[2].start0, item[2].end0, item[1]))
                block_tokens = []
                actual_blocks = []
                for block_index, (contig, strand, interval) in enumerate(blocks, start=1):
                    public = ClosedInterval1.from_interval0(interval)
                    block_tokens.append(
                        f"{contig}:{public.start}-{public.end}:{strand}"
                    )
                    actual_blocks.append({
                        "block_id": f"{element_id}.{occurrence_id}.block_{block_index:03d}",
                        "target_contig": contig,
                        "target_start": public.start,
                        "target_end": public.end,
                        "target_strand": strand,
                    })
                matches = [node_info[node]["match"] for node in occurrence_nodes]
                candidate_ids = {
                    candidate_id
                    for match in matches
                    for candidate_id in (
                        _tokens(match.get("candidate_ids"))
                        | _tokens(match.get("retained_candidate_ids"))
                    )
                }
                template.update({
                    "element_id": element_id,
                    "parent_element_id": base_element,
                    "member_interval_id": f"{element_id}.{occurrence_id}.member_{occurrence_index:03d}",
                    "matched_blocks": _tagged_membership_blocks(
                        occurrence_nodes, node_info,
                    ),
                    "genomic_matched_blocks": ";".join(block_tokens) or "NA",
                    "actual_matched_blocks": json.dumps(
                        actual_blocks, sort_keys=True, separators=(",", ":"),
                    ),
                    "parent_feature_ids": occurrence.get("source_feature_id", "NA"),
                    "transcript_ids": occurrence.get("transcript_id", "NA"),
                    "path_roles": occurrence.get("path_roles", "unknown"),
                    "coding_roles": occurrence.get("coding_roles", "unknown"),
                    "position_roles": occurrence.get("position_roles", "unknown"),
                    "annotation_source": occurrence.get("annotation_source", "NA"),
                    "original_attributes": occurrence.get("original_attributes", "NA"),
                    "partial_start": occurrence.get("partial_start", "NA"),
                    "partial_end": occurrence.get("partial_end", "NA"),
                    "path_role_records": occurrence.get("path_role_records", "[]"),
                    "repeat_instance_id": "NA",
                    "correspondence_status": "resolved",
                    "candidate_ids": ";".join(sorted(candidate_ids)) or "NA",
                    "membership_edge_eligible": 1,
                    "membership_edge_reason": "resolved_local_matched_subinterval",
                    "position_edge_eligible": 1,
                    "position_edge_reason": "unique_resolved_actual_coordinates",
                    "reference_occurrence_id": "NA",
                    "reference_length": "NA",
                    "reference_length_source": "localized_matched_block_union",
                    "reference_coverage_relation": "localized_actual_blocks",
                    "reference_covered_bases": _block_union_length(blocks),
                    "reference_overlap_bases": 0,
                    "reference_uncovered_bases": "NA",
                    "_interval_records": tuple(blocks),
                })
                localized_rows.append(template)
                rows_by_occurrence[occurrence_id].append(template)

    localized_keys = {
        (row["occurrence_id"], row["parent_element_id"])
        for row in localized_rows
    }
    for occurrence_id, templates in sorted(templates_by_occurrence.items()):
        for base_element, template_source in sorted(templates.items()):
            if (occurrence_id, base_element) in localized_keys:
                continue
            occurrence = occurrence_by_id.get(occurrence_id, {})
            relevant = matches_by_occurrence.get(occurrence_id, ())
            candidate_ids = {
                candidate_id
                for match, _side in relevant
                for candidate_id in (
                    _tokens(match.get("candidate_ids"))
                    | _tokens(match.get("retained_candidate_ids"))
                )
            }
            ambiguous = any(
                match.get("candidate_resolution") == "ambiguous"
                or match.get("match_status") == "candidate_ambiguous"
                or match.get("protein_mapping_status") == "ambiguous_mapping"
                for match, _side in relevant
            )
            template = dict(template_source)
            membership_supported = any(
                _hard_membership_eligible(match) for match, _side in relevant
            )
            template.update({
                "parent_element_id": base_element,
                "member_interval_id": "NA",
                "matched_blocks": "NA",
                "genomic_matched_blocks": "NA",
                "actual_matched_blocks": "[]",
                "parent_feature_ids": occurrence.get("source_feature_id", "NA"),
                "transcript_ids": occurrence.get("transcript_id", "NA"),
                "path_roles": occurrence.get("path_roles", "unknown"),
                "coding_roles": occurrence.get("coding_roles", "unknown"),
                "position_roles": occurrence.get("position_roles", "unknown"),
                "annotation_source": occurrence.get("annotation_source", "NA"),
                "original_attributes": occurrence.get("original_attributes", "NA"),
                "partial_start": occurrence.get("partial_start", "NA"),
                "partial_end": occurrence.get("partial_end", "NA"),
                "path_role_records": occurrence.get("path_role_records", "[]"),
                "repeat_instance_id": "NA",
                "correspondence_status": "ambiguous" if ambiguous else "candidate",
                "candidate_ids": ";".join(sorted(candidate_ids)) or "NA",
                "membership_edge_eligible": int(membership_supported),
                "membership_edge_reason": (
                    "retained_membership_with_unresolved_local_position"
                    if membership_supported else "local_matched_subinterval_unresolved"
                ),
                "position_edge_eligible": 0,
                "position_edge_reason": "actual_local_coordinates_unavailable",
                "reference_occurrence_id": "NA",
                "reference_length": "NA",
                "reference_length_source": "localized_matched_block_union",
                "reference_coverage_relation": "uncovered",
                "reference_covered_bases": 0,
                "reference_overlap_bases": 0,
                "reference_uncovered_bases": "NA",
                "_interval_records": tuple(),
            })
            localized_rows.append(template)
            rows_by_occurrence[occurrence_id].append(template)

    for match in scored:
        if _hard_position_eligible(match):
            continue
        if match.get("match_status") in {
            "candidate_chain_excluded", "candidate_low_similarity", "low_similarity",
        }:
            continue
        if not (
            str(match.get("match_status", "")).startswith("candidate_")
            or match.get("candidate_resolution") in {"candidate", "ambiguous"}
            or _edge_eligible(match, "protein_candidate_evidence_available", False)
        ):
            continue
        for occurrence_id, side in (
            (match.get("query_occurrence_id"), "query"),
            (match.get("subject_occurrence_id"), "subject"),
        ):
            for candidate_blocks, local_candidate_ids in _candidate_block_evidence(
                match, side,
            ):
                for row in rows_by_occurrence.get(occurrence_id, ()):
                    resolved_blocks = row.get("_interval_records", ())
                    if resolved_blocks and not _blocks_overlap(
                        resolved_blocks, candidate_blocks,
                    ):
                        continue
                    row["correspondence_status"] = "ambiguous"
                    if not _hard_membership_eligible(match):
                        row["membership_edge_eligible"] = 0
                    row["position_edge_eligible"] = 0
                    if not _hard_membership_eligible(match):
                        row["membership_edge_reason"] = "overlapping_competing_candidate"
                    row["position_edge_reason"] = "overlapping_competing_candidate"
                    row["candidate_ids"] = ";".join(sorted(
                        _tokens(row.get("candidate_ids"))
                        | local_candidate_ids
                        | _tokens(match.get("retained_candidate_ids"))
                    )) or "NA"

    by_element_copy = defaultdict(list)
    by_element = defaultdict(list)
    for row in localized_rows:
        by_element[row["element_id"]].append(row)
        by_element_copy[
            (row["element_id"], row.get("species", "NA"), row.get("gene_copy_id", "NA"))
        ].append(row)
    for element_id, rows in by_element.items():
        occurrence_ids = {row["occurrence_id"] for row in rows}
        reference_id = max(
            sorted(occurrence_ids),
            key=lambda item: int(occurrence_by_id.get(item, {}).get("end", 0))
            - int(occurrence_by_id.get(item, {}).get("start", 1)) + 1,
        )
        reference_length = (
            int(occurrence_by_id.get(reference_id, {}).get("end", 0))
            - int(occurrence_by_id.get(reference_id, {}).get("start", 1)) + 1
        )
        intervals_by_occurrence = defaultdict(list)
        for match in scored:
            query_id = match.get("query_occurrence_id")
            subject_id = match.get("subject_occurrence_id")
            if reference_id not in {query_id, subject_id}:
                continue
            try:
                blocks = parse_legacy_blocks(match.get("matched_blocks"))
            except (TypeError, ValueError):
                continue
            member_id = subject_id if query_id == reference_id else query_id
            if member_id not in occurrence_ids:
                continue
            intervals_by_occurrence[member_id].extend(
                block.query if query_id == reference_id else block.target
                for block in blocks
            )
        copy_groups = defaultdict(list)
        for row in rows:
            copy_groups[(row.get("species", "NA"), row.get("gene_copy_id", "NA"))].append(row)
        for copy_rows in copy_groups.values():
            intervals = sorted(
                interval
                for row in copy_rows
                for interval in intervals_by_occurrence.get(row["occurrence_id"], ())
            )
            if not intervals:
                continue
            covered = _block_union_length([
                ("reference", "+", interval) for interval in intervals
            ])
            overlap = sum(
                max(0, left.end0 - right.start0)
                for left, right in zip(intervals, intervals[1:])
            )
            relation = (
                "repeated_overlap" if overlap > 0
                else "complementary_complete" if covered >= reference_length and len(copy_rows) > 1
                else "single_partial"
            )
            for row in copy_rows:
                row["reference_occurrence_id"] = reference_id
                row["reference_length"] = reference_length
                row["reference_length_source"] = "alignment_sequence_length"
                row["reference_coverage_relation"] = relation
                row["reference_covered_bases"] = covered
                row["reference_overlap_bases"] = overlap
                row["reference_uncovered_bases"] = max(0, reference_length - covered)
            if relation == "repeated_overlap":
                for index, row in enumerate(
                    sorted(copy_rows, key=lambda item: item["occurrence_id"]), start=1,
                ):
                    row["repeat_instance_id"] = f"{element_id}.repeat_{index:03d}"
    for (element_id, _species, _copy), rows in by_element_copy.items():
        resolved = [row for row in rows if row.get("_interval_records")]
        if len({row["occurrence_id"] for row in resolved}) < 2:
            continue
        locus_groups = []
        for row in sorted(resolved, key=lambda item: item["occurrence_id"]):
            overlapping_groups = [
                group
                for group in locus_groups
                if any(
                    _blocks_overlap(
                        row["_interval_records"], member["_interval_records"],
                    )
                    for member in group
                )
            ]
            if not overlapping_groups:
                locus_groups.append([row])
                continue
            primary = overlapping_groups[0]
            primary.append(row)
            for group in overlapping_groups[1:]:
                primary.extend(group)
                locus_groups.remove(group)
        if len(locus_groups) < 2:
            continue
        locus_groups.sort(key=lambda group: min(
            (
                contig, interval.start0, interval.end0, row["occurrence_id"],
            )
            for row in group
            for contig, _strand, interval in row["_interval_records"]
        ))
        for index, group in enumerate(locus_groups, start=1):
            repeat_id = f"{element_id}.repeat_{index:03d}"
            for row in group:
                row["repeat_instance_id"] = repeat_id
                row["reference_coverage_relation"] = "repeated_overlap"
    for element_id, rows in by_element.items():
        if any(row.get("reference_length_source") == "alignment_sequence_length" for row in rows):
            reference_ids = {
                row.get("reference_occurrence_id") for row in rows
                if row.get("reference_occurrence_id") not in {None, "", "NA"}
            }
            if len(reference_ids) == 1:
                reference_id = next(iter(reference_ids))
                for row in rows:
                    if (
                        row["occurrence_id"] == reference_id
                        and row.get("reference_coverage_relation") != "repeated_overlap"
                    ):
                        row["reference_coverage_relation"] = "reference_axis"
                continue
        reference = max(
            rows,
            key=lambda row: (int(row.get("reference_covered_bases", 0)), row["occurrence_id"]),
        )
        reference_id = reference["occurrence_id"]
        reference_length = reference.get("reference_covered_bases", 0) or "NA"
        for row in rows:
            row["reference_occurrence_id"] = reference_id
            row["reference_length"] = reference_length
            if (
                row["occurrence_id"] == reference_id
                and row.get("reference_coverage_relation") != "repeated_overlap"
            ):
                row["reference_coverage_relation"] = "reference_axis"

    for row in localized_rows:
        row.pop("_interval_records", None)
    element_rows[:] = sorted(
        localized_rows,
        key=lambda row: (
            row.get("family_id", ""), row.get("element_id", ""),
            row.get("species", ""), row.get("gene_copy_id", ""),
            row.get("occurrence_id", ""), row.get("member_interval_id", ""),
        ),
    )


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
        observed_species = {
            occ_by_id[occ_id]["species"] for occ_id in occ_ids
            if norm_state(occ_by_id.get(occ_id, {}).get("presence_status")) == "present"
            and norm_state(occ_by_id.get(occ_id, {}).get("species")) != "unknown"
        }
        observation_call = (
            "observed_in_multiple_species" if len(observed_species) > 1
            else "observed_single_species" if observed_species
            else "no_species_presence_observed"
        )
        conservation.append(
            {
                "homology_id": component_id,
                "occurrence_count": len(occ_ids),
                "species_count": len(observed_species),
                "mean_pairwise_identity": f"{mean_identity:.6g}",
                "role_spectrum": ";".join(f"{key}:{roles[key]}" for key in sorted(roles)),
                "conservation_call": observation_call,
            }
        )

    scored = []
    grouped = defaultdict(list)
    for row in matches:
        score = match_total(row)
        components = normalized_components(row)
        if components["sequence_score"] is None:
            components["sequence_score"] = 0.70 * components["alignment_score"] + 0.30 * components["coverage_score"]
        if components["structural_context_score"] is None:
            components["structural_context_score"] = 0.5 * components["left_context_score"] + 0.5 * components["right_context_score"]
        grouped[row["query_occurrence_id"]].append((score, row["match_id"]))
        scored.append(
            {
                "match_id": row["match_id"],
                "query_occurrence_id": row["query_occurrence_id"],
                "subject_occurrence_id": row["subject_occurrence_id"],
                "alignment_score": f"{components['alignment_score']:.6g}",
                "coverage_score": f"{components['coverage_score']:.6g}",
                "sequence_score": f"{components['sequence_score']:.6g}",
                "structural_context_score": f"{components['structural_context_score']:.6g}",
                "left_context_score": f"{components['left_context_score']:.6g}",
                "right_context_score": f"{components['right_context_score']:.6g}",
                "boundary_score": f"{components['boundary_score']:.6g}",
                "phase_score": f"{components['phase_score']:.6g}",
                "order_score": f"{components['order_score']:.6g}",
                "strand_score": f"{components['strand_score']:.6g}",
                "splice_score": f"{components['splice_score']:.6g}",
                "total_score": f"{score:.6g}",
                "dna_total_score": row.get("total_score", "NA"),
                "correspondence_basis": row.get("correspondence_basis", "DNA"),
                "protein_status": row.get("protein_status", "not_requested"),
                "protein_aa_identity": row.get("protein_aa_identity", "NA"),
                "protein_query_cds_coverage": row.get("protein_query_cds_coverage", "NA"),
                "protein_target_cds_coverage": row.get("protein_target_cds_coverage", "NA"),
                "protein_metrics_scope": row.get("protein_metrics_scope", "NA"),
                "protein_best_query_transcript": row.get("protein_best_query_transcript", "NA"),
                "protein_best_target_transcript": row.get("protein_best_target_transcript", "NA"),
                "protein_supporting_transcripts": row.get("protein_supporting_transcripts", "NA"),
                "match_status": row.get("match_status", "ambiguous"),
                "distance_class": row.get("distance_class", "NA"),
                "threshold": row.get("threshold", "NA"),
                "alignment_cigar": row.get("alignment_cigar", "NA"),
                "reciprocal_status": "unclassified",
                "correspondence_call": "unclassified",
                **{
                    field: row.get(field, "NA")
                    for field in CORRESPONDENCE_EVIDENCE_FIELDS
                },
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
        if row["match_status"] == "mapped" and _hard_membership_eligible(row):
            degree[row["query_occurrence_id"]] += 1
            degree[row["subject_occurrence_id"]] += 1
            graph_edges.append(
                {
                    "edge_id": row["match_id"],
                    "query_occurrence_id": row["query_occurrence_id"],
                    "subject_occurrence_id": row["subject_occurrence_id"],
                    "total_score": row["total_score"],
                    "reciprocal_status": row["reciprocal_status"],
                    "edge_call": (
                        "high_confidence_correspondence"
                        if row.get("candidate_resolution") == "resolved"
                        and row["reciprocal_status"] == "reciprocal_best"
                        and to_float(row["total_score"]) >= 0.7
                        else "supporting_correspondence"
                    ),
                }
            )
        if row["match_status"] != "mapped":
            row["correspondence_call"] = row["match_status"]
        elif row.get("candidate_resolution") == "ambiguous":
            row["correspondence_call"] = "ambiguous_candidate_mapping"
        elif row.get("candidate_resolution") == "excluded":
            row["correspondence_call"] = "outside_near_optimal_chain"
        elif row["reciprocal_status"] == "reciprocal_best":
            row["correspondence_call"] = "reciprocal_best"
        else:
            row["correspondence_call"] = best_status.get(row["match_id"], "not_best")

    _membership_match_details(element_rows, occurrences, scored)

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
    observed_coverage_rows = summarize_observed_element_tree_coverage(input_dir, element_rows)
    path_rows = summarize_intragenic_paths(occurrences, element_rows)

    for row in component_rows:
        deg = degree.get(row["occurrence_id"], 0)
        adjusted, call = membership_call(row.get("confidence"), deg)
        row["membership_score"] = f"{adjusted:.6g}"
        row["membership_call"] = call
    for row in element_rows:
        occ_id = row["occurrence_id"]
        adjusted, call = membership_call(row.get("confidence"), degree.get(occ_id, 0))
        row["membership_score"] = f"{adjusted:.6g}"
        row["membership_call"] = call

    write_tsv(f"{output_dir}/internal_homology_assignments.tsv", component_rows, ["homology_id", "occurrence_id", "family_id", "species", "gene_copy_id", "source_label", "support_type", "confidence", "membership_score", "membership_call"])
    write_tsv(
        f"{output_dir}/element_correspondence.tsv",
        element_rows,
        [
            "element_id",
            "parent_element_id",
            "member_interval_id",
            "family_id",
            "homology_id",
            "occurrence_id",
            "species",
            "gene_copy_id",
            "element_class",
            "display_role",
            "source_label",
            "support_type",
            "confidence",
            "membership_score",
            "membership_call",
            "matched_blocks",
            "genomic_matched_blocks",
            "actual_matched_blocks",
            "parent_feature_ids",
            "transcript_ids",
            "path_roles",
            "coding_roles",
            "position_roles",
            "annotation_source",
            "original_attributes",
            "partial_start",
            "partial_end",
            "path_role_records",
            "repeat_instance_id",
            "correspondence_status",
            "candidate_ids",
            "membership_edge_eligible",
            "membership_edge_reason",
            "position_edge_eligible",
            "position_edge_reason",
            "reference_occurrence_id",
            "reference_length",
            "reference_length_source",
            "reference_coverage_relation",
            "reference_covered_bases",
            "reference_overlap_bases",
            "reference_uncovered_bases",
        ],
    )
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
            "sequence_score",
            "structural_context_score",
            "left_context_score",
            "right_context_score",
            "boundary_score",
            "phase_score",
            "order_score",
            "strand_score",
            "splice_score",
            "total_score",
            "dna_total_score",
            "correspondence_basis",
            "protein_status",
            "protein_aa_identity",
            "protein_query_cds_coverage",
            "protein_target_cds_coverage",
            "protein_metrics_scope",
            "protein_best_query_transcript",
            "protein_best_target_transcript",
            "protein_supporting_transcripts",
            "distance_class",
            "threshold",
            "alignment_cigar",
            "match_status",
            "reciprocal_status",
            "correspondence_call",
            *CORRESPONDENCE_EVIDENCE_FIELDS,
        ],
    )
    write_tsv(f"{output_dir}/internal_homology_graph_edges.tsv", graph_edges, ["edge_id", "query_occurrence_id", "subject_occurrence_id", "total_score", "reciprocal_status", "edge_call"])
    write_tsv(f"{output_dir}/progressive_correspondence.tsv", progressive_rows, ["match_id", "query_species", "subject_species", "query_gene_copy_id", "subject_gene_copy_id", "tree_distance", "progressive_tier", "total_score", "correspondence_call"])
    write_tsv(f"{output_dir}/progressive_element_correspondence.tsv", progressive_element_rows, ["element_id", "family_id", "homology_ids", "member_count", "species_count", "copy_count", "exon_like_members", "candidate_source_members", "nearest_species_support", "within_clade_support", "deep_tree_support", "best_pair_score", "mean_pair_score", "progressive_call"])
    write_tsv(f"{output_dir}/observed_element_tree_coverage.tsv", observed_coverage_rows, ["family_id", "element_id", "mrca_label", "present_species_count", "tree_tip_count", "coverage_class", "present_species", "present_copies"])
    write_tsv(f"{output_dir}/observed_intragenic_paths.tsv", path_rows, ["family_id", "path_scope", "node_label", "species", "gene_copy_id", "path_type", "element_path", "context_count", "path_support"])
    return component_rows, conservation, scored
