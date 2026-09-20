"""observations / catalogue: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.observations.schema import structural_site_annotation_view
from insiphy.observations.support import MISSING
from insiphy.observations.support import _finalize_site_metadata
from insiphy.observations.support import _observation_row
from insiphy.observations.support import _tokens
from insiphy.storage.tabular import read_tsv
from insiphy.storage.values import norm_state
from pathlib import Path


def _structural_site_universe(input_dir, output_dir):
    for base in (Path(input_dir), Path(output_dir)):
        rows = read_tsv(base / "structural_site_universe.tsv", optional=True)
        if rows:
            return rows
    return []


def _state_labels_for_layer(layer):
    return ("not_exonic", "exonic") if layer == "exon_role" else ("absent", "present")


def _merge_catalogue_observations(
    rows, catalogue_rows, valid_families, species_by_family, annotation_view
):
    if annotation_view not in {"canonical", "repertoire"}:
        raise ValueError("annotation_view must be 'repertoire' or 'canonical'")
    merged = {(row["family_id"], row["layer"], row["site_id"], row["species"]): dict(row)
              for row in rows}
    for catalogue in catalogue_rows:
        family, layer = catalogue.get("family_id"), catalogue.get("layer")
        site_id, species = catalogue.get("site_id"), catalogue.get("species")
        state = norm_state(catalogue.get("state"))
        if not all(value and value != "NA" for value in (family, layer, site_id)):
            continue
        catalogue_view = structural_site_annotation_view(catalogue)
        if layer == "exon_presence":
            if catalogue_view != "view_independent":
                raise ValueError("structural-site presence catalogue entries are view-independent")
        elif catalogue_view not in {"canonical", "repertoire"}:
            raise ValueError(
                "structural-site role and junction catalogue entries require an explicit "
                "or inferable canonical/repertoire view"
            )
        elif catalogue_view != annotation_view:
            raise ValueError(
                "structural-site catalogue annotation view conflicts with this analysis: "
                f"catalogue={catalogue_view}, analysis={annotation_view}"
            )
        if not species or species == "NA":
            continue
        if state == "unknown" and catalogue.get("state") in MISSING:
            continue
        if family not in valid_families or species not in species_by_family.get(family, set()):
            continue
        state_0, state_1 = catalogue.get("state_0"), catalogue.get("state_1")
        if state_0 in MISSING or state_1 in MISSING:
            state_0, state_1 = _state_labels_for_layer(layer)
        discovery_rule = catalogue.get("discovery_rule")
        if discovery_rule in MISSING:
            discovery_rule = "independent_catalogue"
        row = _observation_row(
            family=family, layer=layer, site_id=site_id, species=species, state=state,
            state_0=state_0, state_1=state_1,
            evidence={catalogue.get("evidence", "structural_site_universe_explicit_observation")},
            applicability=catalogue.get("applicability", "applicable" if state != "unknown" else "undetermined"),
            reason=catalogue.get("observation_reason", "explicit_catalogue_observation"),
            transcript_scope=catalogue.get("transcript_scope", "catalogue_defined"),
            parent_ids=_tokens(catalogue.get("parent_feature_ids")),
            interval_ids=_tokens(catalogue.get("member_interval_ids")),
            evidence_ids=_tokens(catalogue.get("evidence_ids")),
            discovery_rule=discovery_rule,
            linked_group=catalogue.get("linked_group_id", "NA"),
            site_kind=catalogue.get("site_kind", "catalogue_explicit_observation"),
            confidence=catalogue.get("confidence_flag", "curated"),
            conclusion=catalogue.get("conclusion_flag", "explicit_catalogue_observation"),
            annotation_completeness=catalogue.get("annotation_completeness", "catalogue_reported"))
        row["observation_source"] = catalogue.get("observation_source", "independent_structural_site_catalogue")
        row["discovery_species"] = catalogue.get("discovery_species", "NA")
        row["annotation_view"] = catalogue_view
        catalogue_mask = catalogue.get("observation_mask")
        if catalogue_mask not in {None, ""}:
            row["observation_mask"] = catalogue_mask
        merged[(family, layer, site_id, species)] = row
    return _finalize_site_metadata([merged[key] for key in sorted(merged)])


def _complete_tree_tip_observations(rows, tip_labels):
    """Materialize unknown observations instead of letting inference invent missing tips."""

    by_site = defaultdict(list)
    for row in rows:
        by_site[(row["family_id"], row["layer"], row["site_id"])].append(row)

    completed = list(rows)
    for site_key, site_rows in sorted(by_site.items()):
        template = site_rows[0]
        observed_species = {row["species"] for row in site_rows}
        for species in sorted(set(tip_labels) - observed_species):
            missing = dict(template)
            missing.update({
                "species": species,
                "state": "unknown",
                "evidence": "no_structural_observation_for_tree_tip",
                "applicability": "undetermined",
                "observation_reason": "tree_tip_missing_from_structural_evidence",
                "parent_feature_ids": "NA",
                "member_interval_ids": "NA",
                "evidence_ids": "NA",
                "observation_mask": "missing",
                "observation_source": "no_observation_for_tree_tip",
                "annotation_completeness": "tree_tip_not_observed",
                "confidence_flag": "low",
                "conclusion_flag": "unknown",
            })
            completed.append(missing)
    return completed
