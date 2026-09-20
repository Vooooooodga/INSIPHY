"""The IntraPhy command-line interface."""

from .commands.parser import build_parser
import json
from pathlib import Path

from . import __version__
from .alignment import available_alignment_backends
from .annotation import complete_annotation, generate_sequence_evidence
from .baseline import evaluate_baselines
from .benchmark import benchmark_events
from .calibration import DEFAULT_SCENARIOS, calibrate_simulations
from .case import build_case, inspect_annotation, scan_hidden_segments
from .correspondence import infer_correspondence
from .orthofinder import import_orthofinder
from .phylogeny import infer_phylogeny
from .preprocess import assess_short_candidate_thresholds, derive_tables, extract_gene
from .simulate import simulate_dataset
from .visualize import visualize_results


from .workflow import run_all, _record_evidence_aligner



def _dispatch(args):
    if args.command == "example":
        from .verification.native_cases import build_native_example
        manifest = build_native_example(args.output_dir, args.scenario, args.seed)
        print(f"Synthetic inputs written: {manifest}")
    elif args.command == "check":
        from .commands.environment import environment_report
        print(json.dumps({"status": "valid", "environment": environment_report()}, indent=2))
    elif args.command == "extract-gene":
        extract_gene(args.genome, args.annotation, args.gene_id, args.family_id, args.species, args.gene_copy_id, args.output_dir, args.append, args.transcript_policy, args.canonical_rule, args.source_label, args.copy_role, flank=args.flank, max_extension=args.max_extension)
    elif args.command == "derive-tables":
        derive_tables(
            args.input_dir,
            args.output_dir,
            args.identity_threshold,
            args.distance_table,
            args.aligner,
            args.threads,
            args.min_size_ratio,
            context_aligner=args.context_aligner,
            coding_msa_mode=args.coding_msa_mode,
            short_context_max_length=args.short_context_max_length,
        )
    elif args.command == "simulate":
        simulate_dataset(args.output_dir, args.seed, args.scenario)
    elif args.command == "benchmark":
        benchmark_events(args.input_dir, args.output_dir)
    elif args.command == "calibrate":
        calibrate_simulations(args.output_dir, args.scenario, args.replicates, args.bootstrap_replicates, args.stochastic_maps, args.seed)
    elif args.command == "visualize":
        encoding = args.correspondence_encoding or "color"
        visualize_results(args.input_dir, args.result_dir, args.output_dir, encoding,
                          layout=args.layout, targets=args.target,
                          target_manifest=args.target_manifest)
    elif args.command == "inspect-aligners":
        from .commands.environment import inspect_tools
        print("tool\tavailable\tpath\tversion")
        for row in inspect_tools():
            print(f"{row['tool']}\t{row['available']}\t{row['path'] or 'NA'}\t{row['version'] or 'NA'}")
    elif args.command == "inspect-annotation":
        inspect_annotation(args.annotation, args.output_dir, args.query, args.alias_file, args.species, args.case_id)
    elif args.command == "build-case":
        build_case(
            args.manifest,
            args.output_dir,
            args.identity_threshold,
            args.species_tree,
            args.transcript_policy,
            args.canonical_rule,
            args.aligner,
            args.threads,
            args.min_size_ratio,
            args.copy_tree,
            args.gene_tree,
            flank=args.flank,
            max_extension=args.max_extension,
            context_aligner=args.context_aligner,
            coding_msa_mode=args.coding_msa_mode,
            short_context_max_length=args.short_context_max_length,
        )
    elif args.command == "import-orthofinder":
        import_orthofinder(args.orthofinder_dir, args.orthogroup, args.genome_manifest, args.output_dir, args.species_tree)
    elif args.command == "scan-hidden-segments":
        scan_hidden_segments(args.source_fasta, args.target_fasta, args.output_dir, args.family_id, args.species, args.gene_copy_id, args.min_identity, args.min_coverage, args.aligner, args.threads)
    elif args.command == "candidate-sensitivity":
        assess_short_candidate_thresholds(
            args.segment_matches,
            args.output,
            args.identity,
            args.coverage,
        )
    elif args.command == "complete-annotation":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        complete_annotation(args.input_dir, args.output_dir)
    elif args.command == "segment-correspondence":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_correspondence(args.input_dir, args.output_dir)
    elif args.command == "infer-phylogeny":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_phylogeny(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            bootstrap_replicates=args.bootstrap_replicates,
            stochastic_maps=args.stochastic_maps,
            seed=args.seed,
            foreground_branches=args.foreground_branches,
            analysis_scope=args.analysis_scope,
            model=args.model,
            branch_length_mode=args.branch_length_mode,
            ascertainment=args.ascertainment,
            threads=args.threads,
            root_frequency=args.root_frequency,
            root_presence=args.root_presence,
            structural_site_matrix_path=args.structural_site_matrix,
            analysis_range=args.analysis_range,
            min_callable_fraction=args.min_callable_fraction,
            annotation_view=args.annotation_view,
        )
    elif args.command == "compare-baselines":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        evaluate_baselines(args.input_dir, args.output_dir)
    elif args.command == "run":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        run_all(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            bootstrap_replicates=args.bootstrap_replicates,
            stochastic_maps=args.stochastic_maps,
            seed=args.seed,
            foreground_branches=args.foreground_branches,
            analysis_scope=args.analysis_scope,
            model=args.model,
            branch_length_mode=args.branch_length_mode,
            ascertainment=args.ascertainment,
            threads=args.threads,
            root_frequency=args.root_frequency,
            root_presence=args.root_presence,
            evidence_aligner=args.evidence_aligner,
            short_context_max_length=args.short_context_max_length,
            annotation_view=args.annotation_view,
            structural_site_matrix_path=args.structural_site_matrix,
            analysis_range=args.analysis_range,
            min_callable_fraction=args.min_callable_fraction,
        )


from .commands.session import command_session


def main(argv=None):
    import sys
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        with command_session(args):
            _dispatch(args)
        return 0
    except KeyboardInterrupt:
        print("IntraPhy interrupted; see execution.json for the failed run.", file=sys.stderr)
        return 130
    except (OSError, ValueError, RuntimeError, SystemExit) as exc:
        if getattr(args, "debug", False):
            raise
        if isinstance(exc, SystemExit) and exc.code in {None, 0}:
            return 0
        print(f"IntraPhy error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
