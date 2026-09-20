"""Fail before expensive analysis when arguments or declared inputs are invalid."""
from collections import Counter
import math
from pathlib import Path
import shutil

from ..preparation.manifests import load_manifest
from ..storage.tabular import read_tsv
from ..topology import SpeciesTree


def required_tools(args):
    names = set()
    if args.command in {"build-case", "derive-tables", "check"}:
        names.add("mafft")  # Family protein alignment and exon-pair alignment.
        names.add(getattr(args, "context_aligner", "minimap2"))
        names.add(getattr(args, "aligner", "mafft"))
    if args.command == "run":
        names.update({"mafft", getattr(args, "evidence_aligner", "minimap2")})
        if getattr(args, "evidence_aligner", "minimap2") == "miniprot":
            names.add("minimap2")  # Nucleotide evidence is a separate channel.
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
    if args.command == "build-case" and not getattr(args, "species_tree", None):
        raise ValueError("--species-tree is required for a formal build-case; use extract-gene for extraction alone")
    if getattr(args, "manifest", None):
        rows = load_manifest(args.manifest)
        if args.command in {"build-case", "check"}:
            counts = Counter((r["family_id"], r["species"]) for r in rows)
            duplicates = ["/".join(key) for key, n in counts.items() if n > 1]
            if duplicates:
                raise ValueError("Formal input requires one gene per species and family: " + ", ".join(duplicates))
    if args.command in {"build-case", "check"} and getattr(args, "species_tree", None):
        import tempfile
        from ..orthofinder import _write_species_tree
        with tempfile.TemporaryDirectory(prefix="intraphy-tree-check-") as tmp:
            converted = Path(tmp) / "tree.tsv"
            _write_species_tree(args.species_tree, converted)
            tree = SpeciesTree(read_tsv(converted))
        panel = set(tree.leaf_by_label)
        if getattr(args, "manifest", None):
            families = {}
            for row in rows:
                families.setdefault(row["family_id"], set()).add(row["species"])
            for family, species in families.items():
                if species != panel:
                    raise ValueError(f"Tree/manifest species mismatch for {family}: "
                                     f"missing={sorted(panel-species)}, extra={sorted(species-panel)}")
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
        raise RuntimeError("Required alignment tools are unavailable on PATH: " + ", ".join(missing)
                           + ". See docs/installation.md and run intraphy inspect-aligners.")
