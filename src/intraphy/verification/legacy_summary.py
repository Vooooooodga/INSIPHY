"""verification / legacy_summary: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from pathlib import Path


def evaluate_baselines(input_dir, output_dir):
    """Deprecated output summary, explicitly not a comparison of methods.

    Old versions rewarded the same output with hand-written method scores.
    Those numbers are intentionally not reproduced as scientific evidence.
    """
    output_dir = Path(output_dir)
    events = read_tsv(output_dir / "candidate_structural_events.tsv", optional=True)
    occurrences = read_tsv(Path(input_dir) / "segment_occurrences.tsv", optional=True)
    rows = [{"baseline_model": "not_evaluated", "status": "not_a_method_comparison",
             "family_count": len({row.get("family_id") for row in occurrences}),
             "event_candidates_detected": len(events), "score": "NA", "delta_vs_best": "NA",
             "interpretation": "Separate predictions and independent truth are required; heuristic method rewards were removed in 0.16."}]
    write_tsv(output_dir / "baseline_comparison.tsv", rows, list(rows[0]))
    return rows
