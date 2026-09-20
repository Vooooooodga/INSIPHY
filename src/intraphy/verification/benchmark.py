"""verification / benchmark: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from intraphy.storage.values import to_float
from pathlib import Path


NON_BIOLOGICAL_BENCHMARK_CLASSES = {
    "sequence_supported_annotation_gap",
    "annotation_or_alignment_evidence",
}


COPY_CONTEXT_CLASSES = {
    "copy_duplication_or_expansion",
    "copy_duplication_or_relocation",
    "copy_loss_or_collapse",
    "copy_multiplicity_shift",
    "retrocopy_or_dispersed_duplication_candidate",
}


AMBIGUOUS_CLASSES = {
    "gene_conversion_candidate",
    "ambiguous_paralogous_similarity",
}


def event_key(row, include_branch=False):
    pattern = row.get("structural_pattern") or row.get("event_class", "")
    object_id = next((str(row.get(key)) for key in ("site_id", "object_id", "element_id")
                      if row.get(key) not in {None, "", "NA"}), "unresolved_object")
    base = (row.get("family_id", ""), pattern, object_id)
    return base + (row.get("branch_scope", ""),) if include_branch else base


def call_scope(row):
    scope = row.get("call_scope", "")
    if scope:
        return scope
    event_class = row.get("structural_pattern") or row.get("event_class", "")
    if event_class in AMBIGUOUS_CLASSES:
        return "ambiguous_evidence"
    if event_class in COPY_CONTEXT_CLASSES:
        return "copy_context"
    if event_class in NON_BIOLOGICAL_BENCHMARK_CLASSES:
        return "annotation_evidence"
    return "core_structural_event"


def benchmark_events(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    if not (input_dir / "truth_events.tsv").exists():
        rows = [{"status": "unavailable", "reason": "truth_file_not_provided"}]
        write_tsv(output_dir / "benchmark_summary.tsv", rows, ["status", "reason"])
        return rows
    truth = read_tsv(input_dir / "truth_events.tsv", ["family_id", "event_class"], optional=True)
    calls = read_tsv(output_dir / "candidate_structural_events.tsv", ["family_id", "event_class"], optional=True)
    core_truth = [row for row in truth if call_scope(row) == "core_structural_event"]
    core_calls = [row for row in calls if call_scope(row) == "core_structural_event"]
    bootstraps = read_tsv(output_dir / "hypothesis_bootstrap.tsv", ["empirical_p_value"], optional=True)
    truth_set = {event_key(row) for row in core_truth}
    call_set = {event_key(row) for row in core_calls}
    truth_branch_set = {event_key(row, include_branch=True) for row in core_truth if "->" in row.get("branch_scope", "")}
    call_branch_set = {event_key(row, include_branch=True) for row in core_calls if "->" in row.get("branch_scope", "")}
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
            "copy_context_truth": sum(1 for row in truth if call_scope(row) == "copy_context"),
            "copy_context_called": sum(1 for row in calls if call_scope(row) == "copy_context"),
            "ambiguous_evidence_truth": sum(1 for row in truth if call_scope(row) == "ambiguous_evidence"),
            "ambiguous_evidence_called": sum(1 for row in calls if call_scope(row) == "ambiguous_evidence"),
            "unplaced_core_calls": sum(1 for row in core_calls if "->" not in row.get("branch_scope", "")),
        }
    ]
    details = []
    for key in sorted(truth_set | call_set):
        details.append(
            {
                "family_id": key[0],
                "event_class": key[1],
                "object_id": key[2],
                "structural_pattern": key[1],
                "call_scope": "core_structural_event",
                "truth_status": "truth_present" if key in truth_set else "truth_absent",
                "call_status": "called" if key in call_set else "not_called",
                "benchmark_call": "true_positive" if key in truth_set and key in call_set else "false_positive" if key in call_set else "false_negative",
            }
        )
    write_tsv(output_dir / "benchmark_summary.tsv", rows, ["truth_events", "called_events", "true_positive", "false_positive", "false_negative", "precision", "recall", "f1", "branch_true_positive", "branch_accuracy", "copy_context_truth", "copy_context_called", "ambiguous_evidence_truth", "ambiguous_evidence_called", "unplaced_core_calls"])
    write_tsv(output_dir / "benchmark_detailed.tsv", details, ["family_id", "event_class", "object_id", "structural_pattern", "call_scope", "truth_status", "call_status", "benchmark_call"])
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
