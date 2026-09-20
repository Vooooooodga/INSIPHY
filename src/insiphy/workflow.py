"""Run-level orchestration. Stages consume explicit inputs, not stale outputs."""
import json
from pathlib import Path
from .annotation import complete_annotation, generate_sequence_evidence
from .correspondence import infer_correspondence
from .phylogeny import infer_phylogeny
from .baseline import evaluate_baselines


def run_all(
    input_dir,
    output_dir,
    bootstrap_replicates=0,
    stochastic_maps=0,
    seed=7,
    foreground_branches=None,
    analysis_scope="single-copy",
    model="parsimony",
    branch_length_mode="supplied",
    ascertainment="observed-at-least-one",
    threads=1,
    root_frequency="estimated",
    root_presence=0.5,
    evidence_aligner="minimap2",
    short_context_max_length=300,
    annotation_view="repertoire",
    structural_site_matrix_path=None,
    analysis_range="all",
    min_callable_fraction=0.70,
):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    from .run_result import RunResult
    RunResult(model, analysis_scope, "species_tree.tsv", (), status="running").write(output_dir)
    infer_correspondence(input_dir, output_dir)
    generate_sequence_evidence(
        input_dir,
        output_dir,
        threads=threads,
        aligner=evidence_aligner,
        short_context_max_length=short_context_max_length,
    )
    complete_annotation(input_dir, output_dir)
    infer_correspondence(
        input_dir, output_dir,
        annotation_completion_path=Path(output_dir) / "annotation_completion_candidates.tsv",
    )
    infer_phylogeny(
        input_dir,
        output_dir,
        bootstrap_replicates=bootstrap_replicates,
        stochastic_maps=stochastic_maps,
        seed=seed,
        foreground_branches=foreground_branches,
        analysis_scope=analysis_scope,
        model=model,
        branch_length_mode=branch_length_mode,
        ascertainment=ascertainment,
        threads=threads,
        root_frequency=root_frequency,
        root_presence=root_presence,
        annotation_view=annotation_view,
        structural_site_matrix_path=structural_site_matrix_path,
        analysis_range=analysis_range,
        min_callable_fraction=min_callable_fraction,
    )
    _record_evidence_aligner(
        output_dir,
        evidence_aligner,
        short_context_max_length=short_context_max_length,
    )
    if analysis_scope == "experimental-multicopy":
        evaluate_baselines(input_dir, output_dir)

def _record_evidence_aligner(
    output_dir,
    evidence_aligner,
    short_context_max_length=None,
):
    parameter_path = Path(output_dir) / "run_parameters.json"
    if not parameter_path.exists():
        return
    parameters = json.loads(parameter_path.read_text())
    if not isinstance(parameters, dict):
        raise ValueError(f"{parameter_path} does not contain a JSON object")
    parameters["evidence_aligner"] = evidence_aligner
    if short_context_max_length is not None:
        parameters["short_context_max_length"] = int(short_context_max_length)
    parameter_path.write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")
