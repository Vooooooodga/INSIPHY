"""verification / calibration: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.correspondence import infer_correspondence
from intraphy.evidence.completion import complete_annotation
from intraphy.phylogeny import infer_phylogeny
from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from intraphy.storage.values import to_float
from intraphy.verification.benchmark import benchmark_events
from intraphy.verification.fixtures import simulate_dataset
from intraphy.verification.legacy_summary import evaluate_baselines
from pathlib import Path


DEFAULT_SCENARIOS = [
    "exonization",
    "source_join",
    "segment_split",
    "segment_fusion",
    "tandem_duplication",
    "processed_copy_or_intron_loss",
    "annotation_dropout",
    "negative_control",
]


def _run_case(input_dir, output_dir, bootstrap_replicates=0, stochastic_maps=0, seed=7):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    complete_annotation(input_dir, output_dir)
    infer_correspondence(input_dir, output_dir, annotation_completion_path=Path(output_dir) / "annotation_completion_candidates.tsv")
    infer_phylogeny(
        input_dir,
        output_dir,
        bootstrap_replicates=bootstrap_replicates,
        stochastic_maps=stochastic_maps,
        seed=seed,
        analysis_scope="experimental-multicopy",
    )
    evaluate_baselines(input_dir, output_dir)
    benchmark_events(input_dir, output_dir)


def _mean(values):
    values = [value for value in values if value is not None]
    return sum(values) / len(values) if values else None


def _fmt(value):
    return f"{value:.6g}" if value is not None else "NA"


def calibrate_simulations(
    output_dir,
    scenarios=None,
    replicates=10,
    bootstrap_replicates=0,
    stochastic_maps=0,
    seed=101,
):
    """Run fixed simulated scenarios and summarize statistical behavior.

    The calibration output describes method behavior under known simulated
    histories. It is not used to alter real-data event calls.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    scenarios = scenarios or DEFAULT_SCENARIOS
    replicate_rows = []

    for scenario in scenarios:
        for idx in range(1, int(replicates) + 1):
            run_seed = int(seed) + idx - 1
            case_dir = output_dir / scenario / f"rep_{idx:03d}"
            input_dir = case_dir / "input"
            result_dir = case_dir / "result"
            simulate_dataset(input_dir, seed=run_seed, scenario=scenario)
            _run_case(input_dir, result_dir, bootstrap_replicates, stochastic_maps, run_seed)
            summary = read_tsv(result_dir / "benchmark_summary.tsv")[0]
            calibration = read_tsv(result_dir / "benchmark_calibration.tsv")[0]
            replicate_rows.append(
                {
                    "scenario": scenario,
                    "replicate": idx,
                    "seed": run_seed,
                    "truth_events": summary["truth_events"],
                    "called_events": summary["called_events"],
                    "true_positive": summary["true_positive"],
                    "false_positive": summary["false_positive"],
                    "false_negative": summary["false_negative"],
                    "precision": summary["precision"],
                    "recall": summary["recall"],
                    "f1": summary["f1"],
                    "branch_accuracy": summary["branch_accuracy"],
                    "bootstrap_tests": calibration["bootstrap_tests"],
                    "empirical_p_le_0_05": calibration["empirical_p_le_0_05"],
                    "mean_empirical_p": calibration["mean_empirical_p"],
                }
            )

    summary_rows = []
    for scenario in scenarios:
        rows = [row for row in replicate_rows if row["scenario"] == scenario]
        positive_rows = [row for row in rows if int(row["truth_events"]) > 0]
        null_rows = [row for row in rows if int(row["truth_events"]) == 0]
        precision = _mean([to_float(row["precision"], None) for row in rows])
        recall = _mean([to_float(row["recall"], None) for row in rows])
        f1 = _mean([to_float(row["f1"], None) for row in rows])
        branch_accuracy = _mean([to_float(row["branch_accuracy"], None) for row in positive_rows])
        false_positive_rate = _mean([1.0 if int(row["called_events"]) > 0 else 0.0 for row in null_rows])
        power = _mean([1.0 if int(row["true_positive"]) > 0 else 0.0 for row in positive_rows])
        empirical_p_rate = _mean(
            [
                to_float(row["empirical_p_le_0_05"], 0.0) / max(1.0, to_float(row["bootstrap_tests"], 0.0))
                for row in rows
                if int(row["bootstrap_tests"]) > 0
            ]
        )
        summary_rows.append(
            {
                "scenario": scenario,
                "replicates": len(rows),
                "mean_precision": _fmt(precision),
                "mean_recall": _fmt(recall),
                "mean_f1": _fmt(f1),
                "power_any_true_positive": _fmt(power),
                "false_positive_rate_on_null": _fmt(false_positive_rate),
                "mean_branch_accuracy": _fmt(branch_accuracy),
                "empirical_p_le_0_05_rate": _fmt(empirical_p_rate),
            }
        )

    write_tsv(
        output_dir / "calibration_replicates.tsv",
        replicate_rows,
        [
            "scenario",
            "replicate",
            "seed",
            "truth_events",
            "called_events",
            "true_positive",
            "false_positive",
            "false_negative",
            "precision",
            "recall",
            "f1",
            "branch_accuracy",
            "bootstrap_tests",
            "empirical_p_le_0_05",
            "mean_empirical_p",
        ],
    )
    write_tsv(
        output_dir / "calibration_operating_characteristics.tsv",
        summary_rows,
        [
            "scenario",
            "replicates",
            "mean_precision",
            "mean_recall",
            "mean_f1",
            "power_any_true_positive",
            "false_positive_rate_on_null",
            "mean_branch_accuracy",
            "empirical_p_le_0_05_rate",
        ],
    )
    return summary_rows
