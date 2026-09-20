"""One owner for callable scope. No coordinates, paths, or discoveries are trimmed.

A known absence is a call. An unknown is not an absence. Selection uses the
fraction of the supplied tree panel with a known state, not state-1 frequency.
This conditions inference on a measured mask; it does not model state-dependent
annotation failure or replace the model's ascertainment correction.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import json
from math import isfinite
from pathlib import Path
from typing import Iterable, Mapping, Tuple

from intraphy.observations.schema import normalize_observation_mask
from intraphy.storage.tabular import write_tsv


@dataclass(frozen=True)
class ScopePolicy:
    analysis_range: str = "all"
    min_callable_fraction: float = 0.70
    report_thresholds: Tuple[float, ...] = (0.50, 0.70, 0.90)

    def __post_init__(self):
        if self.analysis_range not in {"all", "high-coverage"}:
            raise ValueError("analysis_range must be all or high-coverage")
        for value in (self.min_callable_fraction, *self.report_thresholds):
            if not isfinite(float(value)) or not 0 <= value <= 1:
                raise ValueError("callable fractions must be finite and within [0,1]")

    def selected(self, known: int, panel: int) -> bool:
        return self.analysis_range == "all" or (
            panel > 0 and known / panel >= self.min_callable_fraction)


def site_key(row):
    return tuple(str(row.get(field, "NA")) for field in ("family_id", "layer", "site_id"))


def _tree_scope(tree, called):
    if tree is None or not called:
        return {"called_mrca": "NA", "represented_root_children": 0,
                "root_children_total": len(tree.children.get(tree.root, ())) if tree else 0}
    def ancestors(node):
        values = []
        while node:
            values.append(node)
            node = tree.parent.get(node)
        return values
    chains = [ancestors(tree.leaf_by_label[label]) for label in sorted(called)]
    common = set(chains[0]).intersection(*(set(chain) for chain in chains[1:]))
    mrca = next(node for node in chains[0] if node in common)
    root_children = set(tree.children.get(tree.root, ()))
    represented = {node for chain in chains for node in chain if node in root_children}
    return {"called_mrca": mrca, "represented_root_children": len(represented),
            "root_children_total": len(root_children)}


def assess_scope(rows: Iterable[Mapping], panel_species=None, policy=None, tree=None):
    """Return effective full rows, selected rows, cell audit, character summaries.

    Does not generate new rows or links. Explicitly masked known values are
    preserved in the audit but become unknown for inference. Existing evidence
    adjudication owns scientific validity; this layer never invents a new score.
    """
    policy = policy or ScopePolicy()
    rows = [dict(row) for row in rows]
    panel = tuple(sorted(panel_species if panel_species is not None else
                         {row["species"] for row in rows}))
    if len(set(panel)) != len(panel):
        raise ValueError("scope panel contains duplicate species")
    by_site = defaultdict(list)
    audit = []
    for row in rows:
        if row["species"] not in panel:
            raise ValueError(f"observation outside scope panel: {row['species']}")
        key = site_key(row)
        raw_state = row["state"]
        mask = normalize_observation_mask(row["observation_mask"])
        applicable = row["applicability"]
        called = (applicable == "applicable" and mask != "missing"
                  and raw_state in (row["state_0"], row["state_1"]))
        if not called:
            row["state"] = "unknown"
            row["observation_mask"] = "missing"
        disposition = "called" if called else (
            "inapplicable" if applicable == "inapplicable" else "unknown")
        if not str(row.get("observation_reason", "")).strip() or row.get("observation_reason") == "NA":
            row["observation_reason"] = "declared_observation" if called else "scope_observation_unresolved"
        audit.append(dict(zip(("family_id", "layer", "site_id"), key),
                          species=row["species"], input_state=raw_state,
                          effective_state=row["state"], disposition=disposition,
                          applicability=applicable, observation_mask=row["observation_mask"],
                          reason=row.get("observation_reason", "NA"),
                          annotation_view=row.get("annotation_view", "NA")))
        by_site[key].append(row)
    summaries = []
    selected_keys = set()
    for key, site_rows in sorted(by_site.items()):
        labels = [row["species"] for row in site_rows]
        if len(labels) != len(set(labels)) or set(labels) != set(panel):
            raise ValueError("scope requires exactly one row per panel species for " + "/".join(key))
        called = [row for row in site_rows if row["state"] != "unknown"]
        n0 = sum(row["state"] == row["state_0"] for row in called)
        n1 = len(called) - n0
        na = sum(row["applicability"] == "inapplicable" for row in site_rows)
        applicable = sum(row["applicability"] == "applicable" for row in site_rows)
        undetermined = len(panel) - applicable - na
        keep = policy.selected(len(called), len(panel))
        if keep:
            selected_keys.add(key)
        reasons = Counter(row["observation_reason"] for row in site_rows if row["state"] == "unknown")
        summary = dict(zip(("family_id", "layer", "site_id"), key),
                       panel_n=len(panel), applicable_n=applicable,
                       inapplicable_n=na, applicability_undetermined_n=undetermined,
                       state0_n=n0, state1_n=n1, callable_n=len(called),
                       unknown_n=len(panel)-na-len(called),
                       callable_fraction_panel=len(called)/len(panel) if panel else 0.0,
                       callable_fraction_applicable=len(called)/applicable if applicable else "NA",
                       state1_fraction_called=n1/len(called) if called else "NA",
                       selected_for_analysis=int(keep), analysis_range=policy.analysis_range,
                       min_callable_fraction=policy.min_callable_fraction,
                       called_species=";".join(sorted(row["species"] for row in called)) or "NA",
                       unknown_reasons=json.dumps(reasons, ensure_ascii=False, sort_keys=True),
                       linked_group_id=site_rows[0].get("linked_group_id", "NA"))
        summary.update(_tree_scope(tree, {row["species"] for row in called}))
        summaries.append(summary)
    selected = [row for row in rows if site_key(row) in selected_keys]
    return rows, selected, audit, summaries


def write_scope_reports(output_dir, audit, summaries, policy):
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "structural_observation_eligibility.tsv", audit,
              ["family_id","layer","site_id","species","input_state","effective_state",
               "disposition","applicability","observation_mask","reason","annotation_view"])
    fields = ["family_id","layer","site_id","panel_n","applicable_n","inapplicable_n",
              "applicability_undetermined_n","state0_n","state1_n","callable_n","unknown_n",
              "callable_fraction_panel","callable_fraction_applicable","state1_fraction_called",
              "selected_for_analysis","analysis_range","min_callable_fraction","called_species",
              "called_mrca","represented_root_children","root_children_total","unknown_reasons",
              "linked_group_id"]
    write_tsv(output / "structural_character_eligibility.tsv", summaries, fields)
    coverage = []
    layers = sorted({row["layer"] for row in summaries})
    for layer in layers:
        subset = [row for row in summaries if row["layer"] == layer]
        for threshold in sorted(set((*policy.report_thresholds, policy.min_callable_fraction))):
            coverage.append({"layer": layer, "threshold": threshold,
                "characters_total": len(subset),
                "characters_retained": sum(row["callable_fraction_panel"] >= threshold for row in subset),
                "purpose": "coverage_summary_not_additional_model_fit"})
    write_tsv(output / "scope_coverage_views.tsv", coverage,
              ["layer","threshold","characters_total","characters_retained","purpose"])
    metadata = {
        "analysis_range": policy.analysis_range,
        "selection_rule": "all evidence-admissible observations" if policy.analysis_range == "all"
             else f"known state 0 or 1 / full supplied tree panel >= {policy.min_callable_fraction}",
        "min_callable_fraction": policy.min_callable_fraction,
        "total_characters": len(summaries),
        "selected_characters": sum(row["selected_for_analysis"] for row in summaries),
        "full_matrix": "structural_site_matrix.tsv",
        "analysis_matrix": "analysis_structural_site_matrix.tsv",
        "coordinates_and_adjacencies_changed": False,
        "ascertainment": "unchanged; handled separately by selected inference model",
        "limitation": "Conditional on observed masks; no state-dependent detection model or threshold calibration.",
    }
    (output / "analysis_scope.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata
