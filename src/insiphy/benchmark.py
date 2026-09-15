"""Benchmark INSIPHY event calls against simulated truth."""

from pathlib import Path

from .io import read_tsv, to_float, write_tsv


def event_key(row, include_branch=False):
    base = (row.get("family_id", ""), row.get("event_class", ""))
    if include_branch:
        return base + (row.get("branch_scope", ""),)
    return base


def benchmark_events(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    truth = read_tsv(input_dir / "truth_events.tsv", ["family_id", "event_class"], optional=True)
    calls = read_tsv(output_dir / "candidate_structural_events.tsv", ["family_id", "event_class"], optional=True)
    bootstraps = read_tsv(output_dir / "hypothesis_bootstrap.tsv", ["empirical_p_value"], optional=True)
    truth_set = {event_key(row) for row in truth}
    call_set = {event_key(row) for row in calls}
    truth_branch_set = {event_key(row, include_branch=True) for row in truth if row.get("branch_scope")}
    call_branch_set = {event_key(row, include_branch=True) for row in calls if row.get("branch_scope")}
    tp = len(truth_set & call_set)
    fp = len(call_set - truth_set)
    fn = len(truth_set - call_set)
    precision = tp / max(1, tp + fp)
    recall = tp / max(1, tp + fn)
    f1 = 2 * precision * recall / max(1e-12, precision + recall)
    branch_tp = len(truth_branch_set & call_branch_set)
    branch_accuracy = branch_tp / max(1, len(truth_branch_set))
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
            "branch_true_positive": branch_tp,
            "branch_accuracy": f"{branch_accuracy:.6g}",
        }
    ]
    details = []
    for key in sorted(truth_set | call_set):
        details.append(
            {
                "family_id": key[0],
                "event_class": key[1],
                "truth_status": "truth_present" if key in truth_set else "truth_absent",
                "call_status": "called" if key in call_set else "not_called",
                "benchmark_call": "true_positive" if key in truth_set and key in call_set else "false_positive" if key in call_set else "false_negative",
            }
        )
    write_tsv(output_dir / "benchmark_summary.tsv", rows, ["truth_events", "called_events", "true_positive", "false_positive", "false_negative", "precision", "recall", "f1", "branch_true_positive", "branch_accuracy"])
    write_tsv(output_dir / "benchmark_detailed.tsv", details, ["family_id", "event_class", "truth_status", "call_status", "benchmark_call"])
    empirical = [to_float(row.get("empirical_p_value"), None) for row in bootstraps]
    empirical = [value for value in empirical if value is not None]
    mcse = [to_float(row.get("monte_carlo_se"), None) for row in bootstraps]
    mcse = [value for value in mcse if value is not None]
    calibration = [
        {
            "bootstrap_tests": len(empirical),
            "empirical_p_le_0_05": sum(1 for value in empirical if value <= 0.05),
            "empirical_p_le_0_10": sum(1 for value in empirical if value <= 0.10),
            "mean_empirical_p": f"{(sum(empirical) / len(empirical)):.6g}" if empirical else "NA",
            "mean_monte_carlo_se": f"{(sum(mcse) / len(mcse)):.6g}" if mcse else "NA",
        }
    ]
    write_tsv(output_dir / "benchmark_calibration.tsv", calibration, ["bootstrap_tests", "empirical_p_le_0_05", "empirical_p_le_0_10", "mean_empirical_p", "mean_monte_carlo_se"])
    return rows
