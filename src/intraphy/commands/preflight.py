"""Fail before expensive analysis when arguments or declared inputs are invalid."""
import math
from pathlib import Path
import shutil

from ..storage.tabular import read_tsv
from ..topology import SpeciesTree


def required_tools(args):
    names = set()
    if args.command in {"build-case", "derive-tables"}:
        names.add("mafft")  # Family protein alignment and exon-pair alignment.
        names.add(getattr(args, "context_aligner", "minimap2"))
        names.add(getattr(args, "aligner", "mafft"))
    if args.command == "run":
        names.update({"mafft", getattr(args, "evidence_aligner", "minimap2")})
        if getattr(args, "evidence_aligner", "minimap2") == "miniprot":
            names.add("minimap2")  # Nucleotide evidence is a separate channel.
    if args.command == "normalize-annotation":
        names.add("agat_convert_sp_gxf2gxf.pl")
    return sorted(names - {"internal", "auto"})


def validate_arguments(args):
    if getattr(args, "threads", 1) < 1:
        raise ValueError("--threads must be at least 1")
    for field in ("min_callable_fraction", "root_presence", "identity_threshold", "min_identity", "min_coverage"):
        value = getattr(args, field, None)
        if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError(f"--{field.replace('_', '-')} must be finite and in [0, 1]")
    for field in ("flank", "max_extension", "short_context_max_length"):
        if getattr(args, field, 0) < 0:
            raise ValueError(f"--{field.replace('_', '-')} must be nonnegative")
    if args.command in {"run", "infer-phylogeny"} and args.analysis_scope == "single-copy":
        if args.bootstrap_replicates or args.stochastic_maps:
            raise ValueError("Bootstrap and stochastic-map flags belong to the experimental model; formal values must be 0")
        if (args.model == "foreground") != bool(args.foreground_branches):
            raise ValueError("--model foreground requires --foreground-branches; other models must omit it")


def validate_input_paths(args):
    for field in ("manifest", "genome", "annotation", "species_tree", "structural_site_matrix", "foreground_branches"):
        value = getattr(args, field, None)
        if value and not Path(value).is_file():
            raise FileNotFoundError(f"--{field.replace('_', '-')}: file does not exist: {value}")
    if args.command == "normalize-annotation":
        from ..inputs.resources import expand_files, GFF_SUFFIXES
        expand_files(args.gff, GFF_SUFFIXES)
        if args.config and not Path(args.config).is_file():
            raise FileNotFoundError(f"AGAT configuration does not exist: {args.config}")
        if args.timeout < 1:
            raise ValueError("--timeout must be positive")
    if args.command in {"build-case", "check", "extract-loci"}:
        from ..inputs.selection import resolve_inputs
        args._input_selection = resolve_inputs(args)
    if args.command in {"run", "infer-phylogeny"}:
        directory = Path(args.input_dir)
        tree_path = directory / "species_tree.tsv"
        if not tree_path.is_file():
            raise FileNotFoundError(f"Prepared input lacks {tree_path.name}; run build-case with --species-tree")
        SpeciesTree(read_tsv(tree_path))
        if args.command == "run" or not args.structural_site_matrix:
            for name in ("segment_occurrences.tsv", "segment_homology.tsv", "segment_sequences.fasta"):
                if not (directory / name).is_file():
                    raise FileNotFoundError(f"Prepared input lacks {name}")


def preflight(args):
    validate_arguments(args)
    validate_input_paths(args)
    missing = [name for name in required_tools(args) if not shutil.which(name)]
    if missing:
        raise RuntimeError("Required external tools are unavailable on PATH: " + ", ".join(missing)
                           + ". See docs/installation.md and run intraphy inspect-aligners.")
