"""Descriptive extant structure summaries; not inferred ancestral paths."""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from intraphy.coordinates import ClosedInterval1, Interval0, parse_legacy_blocks
from intraphy.elements import collect_element_profiles, element_class_for_occurrence, element_role_from_occurrences, membership_call
from intraphy.alignment import global_alignment_stats
from intraphy.io import norm_state, parse_fasta, read_tsv, to_float, uniq, write_tsv
from intraphy.topology import SpeciesTree
from intraphy.mapping.membership import _parse_genomic_blocks
from intraphy.mapping.profiles import progressive_tier
from intraphy.mapping.profiles import tree_distances


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
