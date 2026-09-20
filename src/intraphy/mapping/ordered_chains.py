"""mapping / ordered chains: explicit implementation ownership."""
from __future__ import annotations

from intraphy.mapping.chain_candidates import _collect_chain_candidates
from intraphy.mapping.chain_qualification import _qualify_chain_edges
from intraphy.mapping.path_evaluation import _evaluate_candidate_paths


def _apply_ordered_candidate_chains(rows, occurrence_by_id, occurrences, transcript_paths=None):
    for row in rows:
        row["membership_edge_eligible"] = 0
        row["position_edge_eligible"] = 0
        row["_membership_edge_eligible"] = False
        row["_position_edge_eligible"] = False
        row["retained_candidate_ids"] = "NA"
        row["best_path_candidate_ids"] = "NA"
        row["chain_best_path_count_capped"] = 0
        row["chain_near_optimal_path_count_capped"] = 0
        row["chain_status"] = "unassessed"
        row["chain_ambiguity"] = "unassessed"
        row["chain_start_anchor_ids"] = "NA"
        row["chain_end_anchor_ids"] = "NA"
        row["chain_retained_edges"] = "NA"
        row["left_anchor_id"] = "NA"
        row["right_anchor_id"] = "NA"
    original_intervals, record_by_candidate, owner_by_candidate, groups, candidate_by_id = _collect_chain_candidates(occurrences, rows, occurrence_by_id, transcript_paths)

    physical_retained, physical_metadata, retained_ids, physical_best, best_ids, chain_metadata, left_flank_ids_by_candidate, right_flank_ids_by_candidate, anchor_status_by_candidate = _evaluate_candidate_paths(original_intervals, record_by_candidate, owner_by_candidate, groups, transcript_paths, candidate_by_id)

    for candidate_id in physical_retained:
        candidate = candidate_by_id[candidate_id]
        if candidate.path_memberships:
            continue
        dna, delta = physical_metadata[candidate_id]
        retained_ids.add(candidate_id)
        if candidate_id in physical_best:
            best_ids.add(candidate_id)
        chain_metadata[candidate_id] = {
            "status": "genomic_DNA_only", "best_score": f"{dna.best_score:.6g}",
            "score_delta": f"{delta:.6g}", "local_mode": 1,
            "configuration": dna.configuration_name, "ambiguity": dna.ambiguity_status,
            "candidate_ambiguous": candidate_id in dna.ambiguous_ids,
            "start_anchor_ids": "NA", "end_anchor_ids": "NA",
            "retained_edges": ";".join(f"{a}>{b}" for a,b in sorted(dna.retained_edges)) or "NA",
            "best_path_count_capped": dna.best_path_count_capped,
            "near_optimal_path_count_capped": dna.near_optimal_path_count_capped,
        }

    _qualify_chain_edges(rows, retained_ids, best_ids, physical_retained, candidate_by_id, chain_metadata, left_flank_ids_by_candidate, right_flank_ids_by_candidate, anchor_status_by_candidate, occurrence_by_id, original_intervals, owner_by_candidate, record_by_candidate)
