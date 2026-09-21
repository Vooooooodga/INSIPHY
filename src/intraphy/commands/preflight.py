"""Fail before expensive analysis when arguments or declared inputs are invalid."""
import math
from pathlib import Path
import shutil

from ..storage.tabular import read_tsv
from ..topology import SpeciesTree


def required_tools(args):
    names = set()
    if args.command in {"build-case", "derive-tables", "analyze"}:
        names.add("mafft")  # Family protein alignment and exon-pair alignment.
        names.add(getattr(args, "context_aligner", "minimap2"))
        names.add(getattr(args, "aligner", "mafft"))
    new_model = str(getattr(args, "model", "")).startswith("exon-")
    if (args.command == "run" and not (new_model and getattr(args, "exon_configurations", None))) or (args.command == "infer-phylogeny" and new_model and not getattr(args, "exon_configurations", None)):
        names.update({"mafft", getattr(args, "evidence_aligner", "minimap2")})
        if getattr(args, "evidence_aligner", "minimap2") == "miniprot":
            names.add("minimap2")  # Nucleotide evidence is a separate channel.
    if args.command == "realign-exons":
        names.add(args.cesar)
    if args.command == "normalize-annotation":
        names.add("agat_convert_sp_gxf2gxf.pl")
    return sorted(names - {"internal", "auto"})


def validate_arguments(args):
    if args.command == "realign-exons" and args.timeout < 1:
        raise ValueError("--timeout must be positive")
    if args.command == "fit-exon-rates":
        if args.max_states < 1 or args.gene_bootstrap < 0 or args.parametric_bootstrap < 0:
            raise ValueError("State budget must be positive and replicate counts nonnegative")
    if getattr(args, "threads", 1) < 1:
        raise ValueError("--threads must be at least 1")
    for field in ("min_callable_fraction", "root_presence", "identity_threshold", "min_identity", "min_coverage"):
        value = getattr(args, field, None)
        if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError(f"--{field.replace('_', '-')} must be finite and in [0, 1]")
    for field in ("flank", "max_extension", "short_context_max_length"):
        if getattr(args, field, 0) < 0:
            raise ValueError(f"--{field.replace('_', '-')} must be nonnegative")
    if str(getattr(args, "model", "")).startswith("exon-"):
        for name in ("max_states", "max_origin_scenarios", "max_observation_scenarios", "max_locus_bases", "alignment_timeout", "anchor_bases"):
            if getattr(args, name) < 1:
                raise ValueError(f"--{name.replace('_', '-')} must be positive")
        for name in ("exon_identity", "anchor_identity"):
            value = getattr(args, name)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"--{name.replace('_', '-')} must lie in [0,1]")
        if getattr(args, "analysis_scope", "single-copy") != "single-copy":
            raise ValueError("V19 exon configurations do not implement multicopy genealogy")
        if getattr(args, "structural_site_matrix", None):
            raise ValueError("Legacy structural-site matrices are not V19 exon configurations")
        if getattr(args, "analysis_range", "all") != "all":
            raise ValueError("V19 does not use a percentage trim rule")
        if getattr(args, "annotation_view", "repertoire") != "repertoire":
            raise ValueError("V19 uses explicit local annotation scenarios, not the legacy --annotation-view selector")
        if getattr(args, "root_frequency", "estimated") != "estimated" or getattr(args, "root_presence", .5) != .5:
            raise ValueError("Legacy binary root-frequency flags do not apply to V19 configurations")
        if getattr(args, "ascertainment", "observed-at-least-one") != "observed-at-least-one":
            raise ValueError("V19 discovery is declared in its catalogue, not by a legacy ascertainment flag")
        if getattr(args, "evidence_aligner", "minimap2") != "minimap2":
            raise ValueError("V19 uses the MAFFT/minimap2 evidence route; legacy evidence backends are not substituted")
        if getattr(args, "expected_edits", False) and args.model != "exon-ctmc":
            raise ValueError("--expected-edits requires exon-ctmc")
        if args.command == "analyze" and args.exon_configurations:
            raise ValueError("Use infer-phylogeny to read an existing V19 catalogue without rebuilding raw inputs")
        if (args.model == "exon-ctmc") != bool(args.exon_rates):
            raise ValueError("exon-ctmc requires --exon-rates; parsimony must omit it")
    if args.command in {"run", "infer-phylogeny"} and args.analysis_scope == "single-copy":
        if args.bootstrap_replicates or args.stochastic_maps:
            raise ValueError("Bootstrap and stochastic-map flags belong to the experimental model; formal values must be 0")
        if (args.model == "foreground") != bool(args.foreground_branches):
            raise ValueError("--model foreground requires --foreground-branches; other models must omit it")


def validate_input_paths(args):
    for field in ("manifest", "genome", "annotation", "species_tree", "structural_site_matrix", "foreground_branches", "exon_configurations", "exon_rates"):
        value = getattr(args, field, None)
        values = value if isinstance(value, (tuple, list)) else [value]
        for item in values:
            if item and not Path(item).is_file():
                raise FileNotFoundError(f"--{field.replace('_', '-')}: file does not exist: {item}")
    if args.command == "normalize-annotation":
        from ..inputs.resources import expand_files, GFF_SUFFIXES
        expand_files(args.gff, GFF_SUFFIXES)
        if args.config and not Path(args.config).is_file():
            raise FileNotFoundError(f"AGAT configuration does not exist: {args.config}")
        if args.timeout < 1:
            raise ValueError("--timeout must be positive")
    if args.command in {"build-case", "check", "extract-loci", "analyze"}:
        from ..inputs.selection import resolve_inputs
        args._input_selection = resolve_inputs(args)
    if args.command in {"run", "infer-phylogeny"}:
        directory = Path(args.input_dir)
        tree_path = directory / "species_tree.tsv"
        if not tree_path.is_file():
            raise FileNotFoundError(f"Prepared input lacks {tree_path.name}; run build-case with --species-tree")
        SpeciesTree(read_tsv(tree_path))
        if not getattr(args, "exon_configurations", None) and (args.command == "run" or not args.structural_site_matrix):
            required = ("gene_loci.tsv", "gene_loci.fasta", "transcript_paths.tsv") if str(args.model).startswith("exon-") else ("segment_occurrences.tsv", "segment_homology.tsv", "segment_sequences.fasta")
            for name in required:
                if not (directory / name).is_file():
                    raise FileNotFoundError(f"Prepared input lacks {name}")


def preflight(args):
    validate_arguments(args)
    validate_input_paths(args)
    missing = [name for name in required_tools(args) if not shutil.which(name)]
    if missing:
        raise RuntimeError("Required external tools are unavailable on PATH: " + ", ".join(missing)
                           + ". See docs/installation.md and run intraphy inspect-aligners.")
