"""Command line interface for INSIPHY."""

import argparse
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
from .preprocess import derive_tables, extract_gene
from .simulate import simulate_dataset
from .visualize import visualize_results


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
):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    infer_correspondence(input_dir, output_dir)
    generate_sequence_evidence(input_dir, output_dir, threads=threads, aligner=evidence_aligner)
    complete_annotation(input_dir, output_dir)
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
    )
    _record_evidence_aligner(output_dir, evidence_aligner)
    if analysis_scope == "experimental-multicopy":
        evaluate_baselines(input_dir, output_dir)


def _record_evidence_aligner(output_dir, evidence_aligner):
    parameter_path = Path(output_dir) / "run_parameters.json"
    if not parameter_path.exists():
        return
    parameters = json.loads(parameter_path.read_text())
    if not isinstance(parameters, dict):
        raise ValueError(f"{parameter_path} does not contain a JSON object")
    parameters["evidence_aligner"] = evidence_aligner
    parameter_path.write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="insiphy")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    extract = sub.add_parser("extract-gene")
    extract.add_argument("--genome", required=True)
    extract.add_argument("--annotation", required=True)
    extract.add_argument("--gene-id", required=True)
    extract.add_argument("--family-id", required=True)
    extract.add_argument("--species", required=True)
    extract.add_argument("--gene-copy-id", required=True)
    extract.add_argument("--output-dir", required=True)
    extract.add_argument("--append", action="store_true")
    extract.add_argument(
        "--transcript-policy",
        choices=["canonical", "all"],
        default="all",
        help="Transcript structures to extract; default keeps all annotated transcript paths.",
    )
    extract.add_argument("--canonical-rule", choices=["longest_cds", "longest_span"], default="longest_cds")
    extract.add_argument("--source-label", default="unknown_source")
    extract.add_argument("--copy-role", choices=["source", "background", "derived", "candidate"], default="candidate")
    extract.add_argument("--flank", type=int, default=1000)
    extract.add_argument("--max-extension", type=int, default=10000)

    derive = sub.add_parser("derive-tables")
    derive.add_argument("--input-dir", required=True)
    derive.add_argument("--output-dir")
    derive.add_argument("--identity-threshold", type=float, default=0.7)
    derive.add_argument("--distance-table")
    derive.add_argument("--aligner", choices=["auto", "internal", "mafft", "minimap2", "lastz"], default="mafft", help="Exon-pair backend: mafft/auto uses overlap projection; others use local alignment.")
    derive.add_argument("--context-aligner", choices=["internal", "minimap2", "lastz"], default="minimap2", help="Local backend for pairs involving non-exon sequence.")
    derive.add_argument("--threads", type=int, default=1)
    derive.add_argument("--min-size-ratio", type=float, default=0.25)

    sim = sub.add_parser(
        "simulate",
        description="Legacy experimental-multicopy simulator; it does not validate formal single-copy statistics.",
    )
    sim.add_argument("--output-dir", required=True)
    sim.add_argument("--seed", type=int, default=7)
    sim.add_argument(
        "--scenario",
        choices=[
            "compound",
            "exonization",
            "source_join",
            "tandem_duplication",
            "processed_copy_or_intron_loss",
            "segment_split",
            "segment_fusion",
            "segment_split_fusion",
            "splice_boundary_shift",
            "te_exonization",
            "gene_conversion",
            "negative_control",
            "annotation_dropout",
        ],
        default="compound",
    )

    bench = sub.add_parser("benchmark")
    bench.add_argument("--input-dir", required=True)
    bench.add_argument("--output-dir", required=True)

    cal = sub.add_parser(
        "calibrate",
        description="Legacy experimental-multicopy calibration wrapper; it cannot validate formal single-copy CTMC results.",
    )
    cal.add_argument("--output-dir", required=True)
    cal.add_argument("--scenario", action="append", choices=DEFAULT_SCENARIOS)
    cal.add_argument("--replicates", type=int, default=10)
    cal.add_argument("--bootstrap-replicates", type=int, default=0)
    cal.add_argument("--stochastic-maps", type=int, default=0)
    cal.add_argument("--seed", type=int, default=101)

    viz = sub.add_parser("visualize")
    viz.add_argument("--input-dir", required=True)
    viz.add_argument("--result-dir", required=True)
    viz.add_argument("--output-dir", required=True)
    viz.add_argument("--correspondence-encoding", choices=["pattern", "color"], default=None)

    sub.add_parser("inspect-aligners")

    inspect = sub.add_parser("inspect-annotation")
    inspect.add_argument("--annotation", required=True)
    inspect.add_argument("--output-dir", required=True)
    inspect.add_argument("--query", action="append")
    inspect.add_argument("--alias-file")
    inspect.add_argument("--species", default="NA")
    inspect.add_argument("--case-id", default="case")

    case = sub.add_parser("build-case")
    case.add_argument("--manifest", required=True)
    case.add_argument("--output-dir", required=True)
    case.add_argument("--identity-threshold", type=float, default=0.7)
    case.add_argument("--species-tree")
    case.add_argument("--copy-tree")
    case.add_argument("--gene-tree")
    case.add_argument(
        "--transcript-policy",
        choices=["canonical", "all"],
        default="all",
        help="Transcript structures to extract; default keeps all annotated transcript paths.",
    )
    case.add_argument("--canonical-rule", choices=["longest_cds", "longest_span"], default="longest_cds")
    case.add_argument("--aligner", choices=["auto", "internal", "mafft", "minimap2", "lastz"], default="mafft", help="Exon-pair backend: mafft/auto uses overlap projection; others use local alignment.")
    case.add_argument("--context-aligner", choices=["internal", "minimap2", "lastz"], default="minimap2", help="Local backend for pairs involving non-exon sequence.")
    case.add_argument("--threads", type=int, default=1)
    case.add_argument("--min-size-ratio", type=float, default=0.25)
    case.add_argument("--flank", type=int, default=1000)
    case.add_argument("--max-extension", type=int, default=10000)

    orthofinder = sub.add_parser(
        "import-orthofinder",
        description=(
            "Import one upstream orthogroup: every selected species must map to exactly "
            "one annotated gene locus. Multiple isoform members of that locus are accepted "
            "and their source IDs retained. Ambiguous or unresolved members are excluded."
        ),
    )
    orthofinder.add_argument(
        "--orthofinder-dir", required=True,
        help="Completed results directory containing Orthogroups.tsv or Orthogroups.txt; uses WorkingDirectory/SequenceIDs.txt when available.",
    )
    orthofinder.add_argument("--orthogroup", required=True)
    orthofinder.add_argument(
        "--genome-manifest", required=True,
        help="TSV selecting species and their genome_fasta/annotation_file resources; gene IDs are resolved from the orthogroup.",
    )
    orthofinder.add_argument("--species-tree", help="Supplied TSV or Newick tree, retained without pruning.")
    orthofinder.add_argument("--output-dir", required=True)

    hidden = sub.add_parser("scan-hidden-segments")
    hidden.add_argument("--source-fasta", required=True)
    hidden.add_argument("--target-fasta", required=True)
    hidden.add_argument("--output-dir", required=True)
    hidden.add_argument("--family-id", default="NA")
    hidden.add_argument("--species", default="NA")
    hidden.add_argument("--gene-copy-id", default="NA")
    hidden.add_argument("--min-identity", type=float, default=0.75)
    hidden.add_argument("--min-coverage", type=float, default=0.5)
    hidden.add_argument("--aligner", choices=["internal", "minimap2", "miniprot"], default="internal")
    hidden.add_argument("--threads", type=int, default=1)

    for name in ["complete-annotation", "segment-correspondence", "compare-baselines"]:
        cmd = sub.add_parser(name)
        cmd.add_argument("--input-dir", required=True)
        cmd.add_argument("--output-dir", required=True)
    for name in ["infer-phylogeny", "run"]:
        cmd = sub.add_parser(name)
        cmd.add_argument("--input-dir", required=True)
        cmd.add_argument("--output-dir", required=True)
        cmd.add_argument(
            "--analysis-scope",
            choices=["single-copy", "experimental-multicopy"],
            default="single-copy",
        )
        cmd.add_argument("--model", choices=["parsimony", "er-ard", "foreground"], default="parsimony")
        cmd.add_argument("--branch-length-mode", choices=["supplied", "unit"], default="supplied")
        cmd.add_argument(
            "--ascertainment",
            choices=["observed-at-least-one", "complete-universe", "variable-only"],
            default="observed-at-least-one",
        )
        cmd.add_argument(
            "--root-frequency",
            choices=["estimated", "stationary", "fixed"],
            default="estimated",
        )
        cmd.add_argument("--root-presence", type=float, default=0.5)
        cmd.add_argument("--threads", type=int, default=1)
        cmd.add_argument("--bootstrap-replicates", type=int, default=0)
        cmd.add_argument("--stochastic-maps", type=int, default=0)
        cmd.add_argument("--seed", type=int, default=7)
        cmd.add_argument("--foreground-branches")
        if name == "run":
            cmd.add_argument(
                "--evidence-aligner",
                choices=["internal", "mafft", "minimap2", "lastz", "miniprot"],
                default="minimap2",
                help="Aligner for sequence/protein evidence projection.",
            )

    args = parser.parse_args(argv)
    if args.command == "extract-gene":
        extract_gene(args.genome, args.annotation, args.gene_id, args.family_id, args.species, args.gene_copy_id, args.output_dir, args.append, args.transcript_policy, args.canonical_rule, args.source_label, args.copy_role, flank=args.flank, max_extension=args.max_extension)
    elif args.command == "derive-tables":
        derive_tables(args.input_dir, args.output_dir, args.identity_threshold, args.distance_table, args.aligner, args.threads, args.min_size_ratio, context_aligner=args.context_aligner)
    elif args.command == "simulate":
        simulate_dataset(args.output_dir, args.seed, args.scenario)
    elif args.command == "benchmark":
        benchmark_events(args.input_dir, args.output_dir)
    elif args.command == "calibrate":
        calibrate_simulations(args.output_dir, args.scenario, args.replicates, args.bootstrap_replicates, args.stochastic_maps, args.seed)
    elif args.command == "visualize":
        encoding = args.correspondence_encoding or "color"
        visualize_results(args.input_dir, args.result_dir, args.output_dir, encoding)
    elif args.command == "inspect-aligners":
        for row in available_alignment_backends():
            print(f"{row['aligner']}\t{row['available']}\t{row['notes']}")
    elif args.command == "inspect-annotation":
        inspect_annotation(args.annotation, args.output_dir, args.query, args.alias_file, args.species, args.case_id)
    elif args.command == "build-case":
        build_case(args.manifest, args.output_dir, args.identity_threshold, args.species_tree, args.transcript_policy, args.canonical_rule, args.aligner, args.threads, args.min_size_ratio, args.copy_tree, args.gene_tree, flank=args.flank, max_extension=args.max_extension, context_aligner=args.context_aligner)
    elif args.command == "import-orthofinder":
        import_orthofinder(args.orthofinder_dir, args.orthogroup, args.genome_manifest, args.output_dir, args.species_tree)
    elif args.command == "scan-hidden-segments":
        scan_hidden_segments(args.source_fasta, args.target_fasta, args.output_dir, args.family_id, args.species, args.gene_copy_id, args.min_identity, args.min_coverage, args.aligner, args.threads)
    elif args.command == "complete-annotation":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        complete_annotation(args.input_dir, args.output_dir)
    elif args.command == "segment-correspondence":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_correspondence(args.input_dir, args.output_dir)
    elif args.command == "infer-phylogeny":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_phylogeny(
            args.input_dir,
            args.output_dir,
            args.bootstrap_replicates,
            args.stochastic_maps,
            args.seed,
            args.foreground_branches,
            args.analysis_scope,
            args.model,
            args.branch_length_mode,
            args.ascertainment,
            args.threads,
            args.root_frequency,
            args.root_presence,
        )
    elif args.command == "compare-baselines":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        evaluate_baselines(args.input_dir, args.output_dir)
    elif args.command == "run":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        run_all(
            args.input_dir,
            args.output_dir,
            args.bootstrap_replicates,
            args.stochastic_maps,
            args.seed,
            args.foreground_branches,
            args.analysis_scope,
            args.model,
            args.branch_length_mode,
            args.ascertainment,
            args.threads,
            args.root_frequency,
            args.root_presence,
            evidence_aligner=args.evidence_aligner,
        )


if __name__ == "__main__":
    main()
