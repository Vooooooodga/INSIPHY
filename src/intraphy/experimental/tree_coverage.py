"""experimental / tree coverage: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.elements import collect_element_profiles
from intraphy.experimental.event_labels import EXONIC
from intraphy.experimental.event_labels import role_from_occurrences
from intraphy.experimental.tree_model import SpeciesTree
from intraphy.io import norm_state
from intraphy.io import read_tsv


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
