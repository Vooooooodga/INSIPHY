"""Command line interface for INSIPHY."""

import argparse
from pathlib import Path

from .annotation import complete_annotation
from .baseline import evaluate_baselines
from .benchmark import benchmark_events
from .case import build_case, inspect_annotation, scan_hidden_segments
from .correspondence import infer_correspondence
from .phylogeny import infer_phylogeny
from .preprocess import derive_tables, extract_gene
from .simulate import simulate_dataset


def run_all(input_dir, output_dir):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    complete_annotation(input_dir, output_dir)
    infer_correspondence(input_dir, output_dir)
    infer_phylogeny(input_dir, output_dir)
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

    derive = sub.add_parser("derive-tables")
    derive.add_argument("--input-dir", required=True)
    derive.add_argument("--output-dir")
    derive.add_argument("--identity-threshold", type=float, default=0.7)

    sim = sub.add_parser("simulate")
    sim.add_argument("--output-dir", required=True)
    sim.add_argument("--seed", type=int, default=7)

    bench = sub.add_parser("benchmark")
    bench.add_argument("--input-dir", required=True)
    bench.add_argument("--output-dir", required=True)

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

    hidden = sub.add_parser("scan-hidden-segments")
    hidden.add_argument("--source-fasta", required=True)
    hidden.add_argument("--target-fasta", required=True)
    hidden.add_argument("--output-dir", required=True)
    hidden.add_argument("--family-id", default="NA")
    hidden.add_argument("--species", default="NA")
    hidden.add_argument("--gene-copy-id", default="NA")
    hidden.add_argument("--min-identity", type=float, default=0.75)
    hidden.add_argument("--min-coverage", type=float, default=0.5)

    for name in ["complete-annotation", "segment-correspondence", "infer-phylogeny", "compare-baselines", "run", "run-demo"]:
        cmd = sub.add_parser(name)
        cmd.add_argument("--input-dir", required=True)
        cmd.add_argument("--output-dir", required=True)

    args = parser.parse_args(argv)
    if args.command == "extract-gene":
        extract_gene(args.genome, args.annotation, args.gene_id, args.family_id, args.species, args.gene_copy_id, args.output_dir, args.append)
    elif args.command == "derive-tables":
        derive_tables(args.input_dir, args.output_dir, args.identity_threshold)
    elif args.command == "simulate":
        simulate_dataset(args.output_dir, args.seed)
    elif args.command == "benchmark":
        benchmark_events(args.input_dir, args.output_dir)
    elif args.command == "inspect-annotation":
        inspect_annotation(args.annotation, args.output_dir, args.query, args.alias_file, args.species, args.case_id)
    elif args.command == "build-case":
        build_case(args.manifest, args.output_dir, args.identity_threshold, args.species_tree)
    elif args.command == "scan-hidden-segments":
        scan_hidden_segments(args.source_fasta, args.target_fasta, args.output_dir, args.family_id, args.species, args.gene_copy_id, args.min_identity, args.min_coverage)
    elif args.command == "complete-annotation":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        complete_annotation(args.input_dir, args.output_dir)
    elif args.command == "segment-correspondence":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_correspondence(args.input_dir, args.output_dir)
    elif args.command == "infer-phylogeny":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_phylogeny(args.input_dir, args.output_dir)
    elif args.command == "compare-baselines":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        evaluate_baselines(args.input_dir, args.output_dir)
    elif args.command in {"run", "run-demo"}:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        run_all(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
