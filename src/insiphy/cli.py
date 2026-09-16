"""Command line interface for INSIPHY."""

import argparse
from pathlib import Path

from .alignment import available_alignment_backends
from .annotation import complete_annotation
from .baseline import evaluate_baselines
from .benchmark import benchmark_events
from .case import build_case, inspect_annotation, scan_hidden_segments
from .correspondence import infer_correspondence
from .phylogeny import infer_phylogeny
from .preprocess import derive_tables, extract_gene
from .simulate import simulate_dataset
from .visualize import visualize_results


def run_all(input_dir, output_dir, bootstrap_replicates=0, stochastic_maps=0, seed=7, foreground_branches=None):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    complete_annotation(input_dir, output_dir)
    infer_correspondence(input_dir, output_dir)
    infer_phylogeny(input_dir, output_dir, bootstrap_replicates=bootstrap_replicates, stochastic_maps=stochastic_maps, seed=seed, foreground_branches=foreground_branches)
    evaluate_baselines(input_dir, output_dir)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="insiphy")
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
    extract.add_argument("--transcript-policy", choices=["canonical", "all"], default="canonical")
    extract.add_argument("--canonical-rule", choices=["longest_cds", "longest_span"], default="longest_cds")
    extract.add_argument("--source-label", default="unknown_source")
    extract.add_argument("--copy-role", choices=["source", "background", "derived", "candidate"], default="candidate")

    derive = sub.add_parser("derive-tables")
    derive.add_argument("--input-dir", required=True)
    derive.add_argument("--output-dir")
    derive.add_argument("--identity-threshold", type=float, default=0.7)
    derive.add_argument("--distance-table")
    derive.add_argument("--aligner", choices=["internal", "minimap2", "miniprot"], default="internal")
    derive.add_argument("--threads", type=int, default=1)
    derive.add_argument("--min-size-ratio", type=float, default=0.25)

    sim = sub.add_parser("simulate")
    sim.add_argument("--output-dir", required=True)
    sim.add_argument("--seed", type=int, default=7)
    sim.add_argument("--scenario", choices=["compound", "exonization", "source_join", "tandem_duplication", "segment_split_fusion", "splice_boundary_shift", "te_exonization", "gene_conversion", "negative_control", "annotation_dropout"], default="compound")

    bench = sub.add_parser("benchmark")
    bench.add_argument("--input-dir", required=True)
    bench.add_argument("--output-dir", required=True)

    viz = sub.add_parser("visualize")
    viz.add_argument("--input-dir", required=True)
    viz.add_argument("--result-dir", required=True)
    viz.add_argument("--output-dir", required=True)

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
    case.add_argument("--transcript-policy", choices=["canonical", "all"], default="canonical")
    case.add_argument("--canonical-rule", choices=["longest_cds", "longest_span"], default="longest_cds")
    case.add_argument("--aligner", choices=["internal", "minimap2", "miniprot"], default="internal")
    case.add_argument("--threads", type=int, default=1)
    case.add_argument("--min-size-ratio", type=float, default=0.25)

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
        cmd.add_argument("--bootstrap-replicates", type=int, default=0)
        cmd.add_argument("--stochastic-maps", type=int, default=0)
        cmd.add_argument("--seed", type=int, default=7)
        cmd.add_argument("--foreground-branches")

    args = parser.parse_args(argv)
    if args.command == "extract-gene":
        extract_gene(args.genome, args.annotation, args.gene_id, args.family_id, args.species, args.gene_copy_id, args.output_dir, args.append, args.transcript_policy, args.canonical_rule, args.source_label, args.copy_role)
    elif args.command == "derive-tables":
        derive_tables(args.input_dir, args.output_dir, args.identity_threshold, args.distance_table, args.aligner, args.threads, args.min_size_ratio)
    elif args.command == "simulate":
        simulate_dataset(args.output_dir, args.seed, args.scenario)
    elif args.command == "benchmark":
        benchmark_events(args.input_dir, args.output_dir)
    elif args.command == "visualize":
        visualize_results(args.input_dir, args.result_dir, args.output_dir)
    elif args.command == "inspect-aligners":
        for row in available_alignment_backends():
            print(f"{row['aligner']}\t{row['available']}\t{row['notes']}")
    elif args.command == "inspect-annotation":
        inspect_annotation(args.annotation, args.output_dir, args.query, args.alias_file, args.species, args.case_id)
    elif args.command == "build-case":
        build_case(args.manifest, args.output_dir, args.identity_threshold, args.species_tree, args.transcript_policy, args.canonical_rule, args.aligner, args.threads, args.min_size_ratio)
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
        infer_phylogeny(args.input_dir, args.output_dir, args.bootstrap_replicates, args.stochastic_maps, args.seed, args.foreground_branches)
    elif args.command == "compare-baselines":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        evaluate_baselines(args.input_dir, args.output_dir)
    elif args.command == "run":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        run_all(args.input_dir, args.output_dir, args.bootstrap_replicates, args.stochastic_maps, args.seed, args.foreground_branches)


if __name__ == "__main__":
    main()
