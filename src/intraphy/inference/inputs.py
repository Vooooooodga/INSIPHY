"""inference / inputs: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from intraphy.inference.formatting import _fmt
from intraphy.inference.parameters import COMPLETE_UNIVERSE_RULES
from intraphy.inference.parameters import GENERATED_ALL_ZERO_MASKS
from intraphy.observations.matrix import load_site_matrix
from intraphy.observations.schema import structural_site_observed_state
from intraphy.storage.tabular import read_tsv
import math


def _validated_tree_rows(path, branch_length_mode):
    if branch_length_mode not in {"supplied", "unit"}:
        raise SystemExit("branch_length_mode must be supplied or unit")
    rows = read_tsv(path, ["node_id", "parent_id", "label"])
    out = []
    for row in rows:
        item = dict(row)
        if not row.get("parent_id"):
            raw_root = row.get("branch_length") or row.get("length") or row.get("distance")
            item["branch_length"] = raw_root if raw_root not in {None, "", "NA"} else "NA"
        elif branch_length_mode == "unit":
            item["branch_length"] = 1.0
        else:
            raw = row.get("branch_length") or row.get("length") or row.get("distance")
            if raw in {None, "", "NA"}:
                raise SystemExit(
                    "species_tree.tsv has a branch without length; use --branch-length-mode unit to analyze branch counts"
                )
            length = float(raw)
            if not math.isfinite(length) or length < 0:
                raise SystemExit("all non-root branches must have finite non-negative lengths")
            item["branch_length"] = length
        out.append(item)
    return out


def _site_key(row):
    return (row.get("family_id", "NA"), row.get("layer", "NA"), row.get("site_id", "NA"))


def _load_site_matrix(input_dir, output_dir, structural_site_matrix_path, annotation_view,
                      *, analysis_range="all", min_callable_fraction=0.70):
    return load_site_matrix(input_dir, output_dir, structural_site_matrix_path, annotation_view,
                           analysis_range=analysis_range, min_callable_fraction=min_callable_fraction)


def _masked_state(row):
    return structural_site_observed_state(row)


def _validate_complete_universe(site_rows, tree):
    rows_by_site = defaultdict(dict)
    for row in site_rows:
        rows_by_site[_site_key(row)][row.get("species", "NA")] = row
    missing = []
    for key, by_species in sorted(rows_by_site.items()):
        for species in sorted(tree.leaf_by_label):
            row = by_species.get(species)
            if row is None:
                missing.append("/".join((*key, species)))
                continue
            rule = str(row.get("discovery_rule", "")).strip()
            if rule not in COMPLETE_UNIVERSE_RULES:
                missing.append("/".join((*key, species)) + ":nonindependent_discovery_rule")
                continue
            state = row.get("state", "unknown")
            state_0 = row.get("state_0", "NA")
            state_1 = row.get("state_1", "NA")
            mask = str(row.get("observation_mask", "observed")).strip().lower()
            explicit_observation = mask == "observed" and state in {state_0, state_1}
            generated_all_zero = mask in GENERATED_ALL_ZERO_MASKS and state == state_0
            missing_observation = (mask == "missing" and state == "unknown"
                and str(row.get("observation_reason", "")).strip() not in {"", "NA"})
            if not (explicit_observation or generated_all_zero or missing_observation):
                missing.append("/".join((*key, species)) + ":state_not_explicitly_enumerated")
    if missing:
        preview = ", ".join(missing[:8])
        suffix = " ..." if len(missing) > 8 else ""
        raise SystemExit(
            "complete-universe ascertainment requires independent_catalogue or "
            "curated_complete_universe rows with an observed state, or an explicitly generated "
            "all-zero state, or an explicitly missing observation with a reason, for every family/layer/site/species combination; unavailable: "
            + preview + suffix
        )


def _read_foreground_children(path, tree):
    if not path:
        return frozenset()
    rows = read_tsv(path)
    label_to_node = {label: node for node, label in tree.label.items()}
    children = set()
    invalid = []
    edge_set = set(tree.edges())
    for row in rows:
        parent = row.get("parent_node") or row.get("parent_id") or ""
        child = row.get("child_node") or row.get("child_id") or ""
        if (parent, child) in edge_set:
            children.add(child)
            continue
        scope = row.get("branch_scope") or row.get("branch") or ""
        if "->" in scope:
            parent_label, child_label = (part.strip() for part in scope.split("->", 1))
            parent_node = label_to_node.get(parent_label)
            child_node = label_to_node.get(child_label)
            if (parent_node, child_node) in edge_set:
                children.add(child_node)
                continue
        invalid.append(scope or f"{parent}->{child}")
    if invalid:
        raise SystemExit("foreground branch file contains unmatched branches: " + ", ".join(invalid))
    if not children:
        raise SystemExit("foreground branch file did not match any branch in species_tree.tsv")
    if len(children) == len(edge_set):
        raise SystemExit("foreground model is unidentifiable when every branch is foreground")
    return frozenset(children)


def _analysis_summary(
    site_rows,
    all_encoded_sites,
    included_site_ids,
    branch_length_mode,
    matrix_schema_version,
    matrix_source,
):
    known_counts = sorted(
        sum(value in {0, 1} for value in observations.values())
        for _site_id, observations in all_encoded_sites
    )
    grouped_links = defaultdict(set)
    discovery_rules = set()
    observation_masks = set()
    for row in site_rows:
        discovery_rules.add(str(row.get("discovery_rule", "NA")))
        observation_masks.add(str(row.get("observation_mask", "NA")))
        linked_ids = {
            token
            for token in str(row.get("linked_group_id", "NA")).split(";")
            if token not in {"", "NA", "None", "none"}
        }
        for linked in linked_ids:
            grouped_links[linked].add(row.get("site_id", "NA"))
    included = set(included_site_ids)
    correlated = sum(1 for sites in grouped_links.values() if len(sites & included) > 1)
    if known_counts:
        midpoint = len(known_counts) // 2
        median = (
            known_counts[midpoint]
            if len(known_counts) % 2
            else 0.5 * (known_counts[midpoint - 1] + known_counts[midpoint])
        )
        distribution = ";".join(
            f"{count}:{known_counts.count(count)}" for count in sorted(set(known_counts))
        )
    else:
        median = None
        distribution = "NA"
    return {
        "total_structural_sites": len(all_encoded_sites),
        "known_tip_count_min": min(known_counts) if known_counts else "NA",
        "known_tip_count_median": _fmt(median),
        "known_tip_count_max": max(known_counts) if known_counts else "NA",
        "known_tip_count_distribution": distribution,
        "linked_group_count": len(grouped_links),
        "correlated_linked_group_count": correlated,
        "discovery_rules": ";".join(sorted(discovery_rules)) or "NA",
        "observation_masks": ";".join(sorted(observation_masks)) or "NA",
        "branch_length_mode": branch_length_mode,
        "matrix_schema_version": matrix_schema_version,
        "matrix_source": str(matrix_source),
    }
