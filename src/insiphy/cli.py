"""Command line interface for INSIPHY."""

import argparse
from pathlib import Path

from .annotation import complete_annotation
from .correspondence import infer_correspondence
from .phylogeny import infer_phylogeny
from .preprocess import derive_tables, extract_gene


def run_all(input_dir, output_dir):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    complete_annotation(input_dir, output_dir)
    infer_correspondence(input_dir, output_dir)
    infer_phylogeny(input_dir, output_dir)


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

    for name in ["complete-annotation", "segment-correspondence", "infer-phylogeny", "run", "run-demo"]:
        cmd = sub.add_parser(name)
        cmd.add_argument("--input-dir", required=True)
        cmd.add_argument("--output-dir", required=True)

    args = parser.parse_args(argv)
    if args.command == "extract-gene":
        extract_gene(args.genome, args.annotation, args.gene_id, args.family_id, args.species, args.gene_copy_id, args.output_dir, args.append)
    elif args.command == "derive-tables":
        derive_tables(args.input_dir, args.output_dir, args.identity_threshold)
    elif args.command == "complete-annotation":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        complete_annotation(args.input_dir, args.output_dir)
    elif args.command == "segment-correspondence":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_correspondence(args.input_dir, args.output_dir)
    elif args.command == "infer-phylogeny":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_phylogeny(args.input_dir, args.output_dir)
    elif args.command in {"run", "run-demo"}:
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        run_all(args.input_dir, args.output_dir)


if __name__ == "__main__":
    main()
