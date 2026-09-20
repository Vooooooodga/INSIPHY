"""Count structural character changes without claiming independent mutations."""
from collections import defaultdict
from pathlib import Path

from ..observations.characters import character_id, site_key
from ..storage.tabular import write_tsv

HISTORY_FIELDS = ["family_id", "layer", "site_id", "character_id", "node_id",
                  "node_label", "state_index", "history_id", "selection_rule",
                  "probability", "minimum_changes"]


def annotate_changes(branch_rows):
    for row in branch_rows:
        key = site_key(row)
        row["character_id"] = row["count_unit_id"] = character_id(key)
        row["event_id"] = "NA" if row["event_type"] == "no_change" else character_id(
            (*key, row["parent_id"], row["child_id"], row["parent_state_index"], row["child_state_index"]))
        row["count_interpretation"] = "structural_character_state_transition"
    return branch_rows


def write_change_summary(output_dir, site_rows, branch_rows, histories, node_rows):
    """Summaries are stratified by layer. Possible placements are never added."""
    families = defaultdict(list)
    unique_sites = {site_key(row): row for row in site_rows}
    for key, row in unique_sites.items():
        families[key[:2]].append(row)
    summaries = []
    for (family, layer), sites in sorted(families.items()):
        informative = [r for r in sites if r["pattern_class"] != "all_missing"]
        summaries.append({"family_id": family, "layer": layer,
            "structural_characters": len(sites), "observed_characters": len(informative),
            "variable_characters": sum(r["pattern_class"] == "observed_contrast" for r in sites),
            "minimum_character_changes": sum(int(float(r["min_changes"])) for r in informative) if informative else "NA",
            "count_interpretation": "sum_of_minimum_changes_within_layer",
            "mutation_event_count": "not_estimated",
            "possible_placements_summed": "false"})
    directory = Path(output_dir)
    write_tsv(directory / "gene_change_summary.tsv", summaries,
        ["family_id", "layer", "structural_characters", "observed_characters", "variable_characters",
         "minimum_character_changes", "count_interpretation", "mutation_event_count", "possible_placements_summed"])
    write_tsv(directory / "minimum_change_history.tsv", histories, HISTORY_FIELDS)
    write_ancestral_consistency(directory, node_rows)


def write_ancestral_consistency(output_dir, node_rows):
    """Flag incompatible singleton estimates; do not assemble an ancestral transcript."""
    presence = {(r["family_id"], r["site_id"], r["node_id"]): r["state_indices"]
                for r in node_rows if r["layer"] == "exon_presence"}
    rows = []
    for row in node_rows:
        if row["layer"] != "exon_role":
            continue
        dna = presence.get((row["family_id"], row["site_id"], row["node_id"]))
        if dna is None:
            continue
        role = row["state_indices"]
        status = "incompatible_singleton_states" if dna == "0" and role in {"0", "1"} else (
            "dna_presence_unresolved" if dna != "1" else "no_applicability_conflict")
        rows.append({"family_id": row["family_id"], "site_id": row["site_id"],
                     "node_id": row["node_id"], "sequence_state_indices": dna,
                     "exon_role_state_indices": role, "status": status,
                     "ancestral_transcript_reconstructed": "false"})
    write_tsv(Path(output_dir) / "ancestral_state_consistency.tsv", rows,
              ["family_id", "site_id", "node_id", "sequence_state_indices", "exon_role_state_indices",
               "status", "ancestral_transcript_reconstructed"])
