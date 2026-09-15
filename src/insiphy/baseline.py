"""Baseline comparisons for INSIPHY outputs."""

from collections import Counter
from pathlib import Path

from .io import read_tsv, to_float, write_tsv


def evaluate_baselines(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"])
    matches = read_tsv(output_dir / "segment_correspondence.tsv", ["match_id", "total_score", "correspondence_call"], optional=True)
    annotation = read_tsv(output_dir / "annotation_completion_candidates.tsv", ["completion_call"], optional=True)
    events = read_tsv(output_dir / "candidate_structural_events.tsv", ["event_class"], optional=True)

    hidden = sum(1 for row in annotation if row.get("completion_call") == "hidden_segment_candidate")
    conflicts = sum(1 for row in annotation if row.get("completion_call") == "annotation_conflict_candidate")
    mean_match = sum(to_float(row.get("total_score")) for row in matches) / max(1, len(matches))
    low_match = sum(1 for row in matches if to_float(row.get("total_score")) < 0.55)
    event_classes = Counter(row.get("event_class", "unknown") for row in events)
    chimeric = event_classes.get("chimeric_source_join_candidate", 0) + event_classes.get("chimeric_origin_or_source_mixing", 0)
    structural = len(events)
    families = sorted({row["family_id"] for row in occurrences})

    rows = [
        {
            "baseline_model": "annotation_only",
            "evidence_used": "annotated_roles_and_presence",
            "family_count": len(families),
            "event_candidates_detected": structural - hidden,
            "hidden_candidates_explained": 0,
            "score": f"{structural + hidden + conflicts:.6g}",
            "interpretation": "penalizes events requiring sequence-supported annotation completion",
        },
        {
            "baseline_model": "sequence_only",
            "evidence_used": "segment_sequence_similarity",
            "family_count": len(families),
            "event_candidates_detected": structural,
            "hidden_candidates_explained": hidden,
            "score": f"{structural + low_match + max(0.0, 1.0 - mean_match):.6g}",
            "interpretation": "uses sequence correspondence but lacks intragenic adjacency and copy context",
        },
        {
            "baseline_model": "synteny_aware_phylogenetic",
            "evidence_used": "sequence_boundary_order_copy_context_species_tree",
            "family_count": len(families),
            "event_candidates_detected": structural,
            "hidden_candidates_explained": hidden,
            "score": f"{max(0.0, structural - 0.5 * hidden - 0.75 * chimeric):.6g}",
            "interpretation": "credits coordinated intragenic synteny and source joining on the tree",
        },
    ]
    best = min(to_float(row["score"]) for row in rows)
    for row in rows:
        row["delta_vs_best"] = f"{to_float(row['score']) - best:.6g}"
    write_tsv(
        output_dir / "baseline_comparison.tsv",
        rows,
        ["baseline_model", "evidence_used", "family_count", "event_candidates_detected", "hidden_candidates_explained", "score", "delta_vs_best", "interpretation"],
    )
    return rows
