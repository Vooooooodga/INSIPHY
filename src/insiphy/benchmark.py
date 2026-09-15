"""Benchmark INSIPHY event calls against simulated truth."""

from pathlib import Path

from .io import read_tsv, write_tsv


def event_key(row):
    return (row.get("family_id", ""), row.get("event_class", ""))


def benchmark_events(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    truth = read_tsv(input_dir / "truth_events.tsv", ["family_id", "event_class"], optional=True)
    calls = read_tsv(output_dir / "candidate_structural_events.tsv", ["family_id", "event_class"], optional=True)
    truth_set = {event_key(row) for row in truth}
    call_set = {event_key(row) for row in calls}
    tp = len(truth_set & call_set)
    fp = len(call_set - truth_set)
    fn = len(truth_set - call_set)
    precision = tp / max(1, tp + fp)
    recall = tp / max(1, tp + fn)
    f1 = 2 * precision * recall / max(1e-12, precision + recall)
    rows = [
        {
            "truth_events": len(truth_set),
            "called_events": len(call_set),
            "true_positive": tp,
            "false_positive": fp,
            "false_negative": fn,
            "precision": f"{precision:.6g}",
            "recall": f"{recall:.6g}",
            "f1": f"{f1:.6g}",
        }
    ]
    write_tsv(output_dir / "benchmark_summary.tsv", rows, ["truth_events", "called_events", "true_positive", "false_positive", "false_negative", "precision", "recall", "f1"])
    return rows
