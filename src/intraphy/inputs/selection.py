"""A single validated input specification for direct files or an optional table."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections import Counter, defaultdict
import tempfile

from ..orthofinder import _write_species_tree
from ..preparation.annotation_index import load_annotation_index, read_annotation_for_gene
from ..preparation.manifests import load_manifest
from ..preparation.features import FeatureHierarchy
from ..storage.fasta import fasta_record_length, read_fasta_interval
from ..storage.tabular import read_tsv, write_tsv
from ..topology import SpeciesTree
from .orthologs import Target, targets_from_fastas
from .resources import FASTA_SUFFIXES, expand_files, pair_resources, resource_name

DNA_SYMBOLS = frozenset("ACGTRYSWKMBDHVN")


@dataclass(frozen=True)
class InputSelection:
    rows: tuple[dict, ...]
    source_mode: str
    tree_species: tuple[str, ...]

    def write(self, directory: Path) -> None:
        write_tsv(directory / "input_targets.tsv", self.rows, [
            "species", "family_id", "gene_id", "gene_copy_id", "genome_fasta",
            "annotation_file", "ortholog_fasta", "source_member_ids", "input_mode"])


def tree_panel(path: str | Path) -> tuple[str, ...]:
    if not path or not Path(path).is_file():
        raise ValueError("A rooted --species-tree is required (Newick or TSV)")
    with tempfile.TemporaryDirectory(prefix="intraphy-input-tree-") as directory:
        converted = Path(directory) / "tree.tsv"
        _write_species_tree(path, converted)
        tree = SpeciesTree(read_tsv(converted))
    return tuple(sorted(tree.leaf_by_label))


def _single_locus_targets(resources, family: str) -> tuple[Target, ...]:
    targets = []
    for resource in resources:
        features = load_annotation_index(resource.gff).rows
        genes = [r for r in features if r.get("type", "").lower() == "gene"]
        if len(genes) != 1 or not genes[0].get("id"):
            raise ValueError(
                f"{resource.gff} has {len(genes)} gene records. "
                "For whole-genome annotations, select the family with --orthologs FAMILY.fa. "
                "Without --orthologs, supply exactly one annotated gene locus per species. For GTF/GFF without explicit gene records, normalize the hierarchy first with AGAT."
            )
        targets.append(Target(family, resource.species, genes[0]["id"], (), "single_gene_gff"))
    return tuple(targets)


def _validate_rows(rows, panel):
    for row in rows:
        for field in ("family_id", "species"):
            name = row[field]
            if name in {".", ".."} or any(c in name for c in ("/", "\\", "\t", "\n", "\x00")):
                raise ValueError(f"Invalid {field}: {name!r}")
    reused = Counter((r["species"], r["gene_id"]) for r in rows)
    if any(n > 1 for n in reused.values()):
        raise ValueError("The same gene locus cannot belong to more than one input family")
    counts = Counter((r["family_id"], r["species"]) for r in rows)
    if any(count != 1 for count in counts.values()):
        raise ValueError("Formal input requires one gene per family and species")
    families = defaultdict(set)
    for row in rows:
        families[row["family_id"]].add(row["species"])
    for family, species in families.items():
        if species != set(panel):
            raise ValueError(f"Tree/input species mismatch for {family}: "
                             f"missing={sorted(set(panel)-species)}, extra={sorted(species-set(panel))}")
    record_owners = {}
    # Group by resource so the existing bounded annotation cache can be reused.
    for row in sorted(rows, key=lambda r: (r["annotation_file"], r["gene_id"])):
        features, gene, _, bounds = read_annotation_for_gene(row["annotation_file"], row["gene_id"])
        if gene.get("strand") not in {"+", "-"}:
            raise ValueError(f"Target gene {row['gene_id']} lacks a resolved genomic strand")
        target_features = FeatureHierarchy(features).transcript_features(gene["id"])
        if not any(f.get("type", "").lower() in {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr"}
                   for f in target_features):
            raise ValueError(f"Target gene {row['gene_id']} has no extractable exon, CDS or UTR annotation. "
                             "Genomic sequence alone does not define a complete gene model.")
        record = (row["genome_fasta"], gene["seqid"])
        previous_species = record_owners.setdefault(record, row["species"])
        if previous_species != row["species"]:
            raise ValueError("A combined genomic FASTA requires distinct sequence IDs for different species. "
                             f"Record {gene['seqid']} is assigned to {previous_species} and {row['species']}.")
        length = fasta_record_length(row["genome_fasta"], gene["seqid"])
        start, end = int(bounds["linked_start"]), int(bounds["linked_end"])
        if not 1 <= start <= end <= length:
            raise ValueError(f"GFF coordinates {gene['seqid']}:{start}-{end} exceed FASTA length {length}. "
                             "Use genomic DNA and matching coordinates; spliced CDS/protein cannot supply introns.")
        probe = read_fasta_interval(row["genome_fasta"], gene["seqid"], start, min(end, start+999))
        if set(probe) - DNA_SYMBOLS:
            raise ValueError(f"{row['genome_fasta']} contains non-DNA symbols at the target locus. "
                             "Use genomic DNA with --fasta; CDS/protein FASTA belongs in --orthologs.")


def resolve_inputs(args) -> InputSelection:
    panel = tree_panel(args.species_tree)
    if getattr(args, "manifest", None):
        if any(getattr(args, name, None) for name in ("fasta", "gff", "orthologs", "family_id")):
            raise ValueError("Choose direct FASTA/GFF inputs or --manifest, not both")
        rows = load_manifest(args.manifest)
        mode = "manifest"
    else:
        if not getattr(args, "fasta", None) or not getattr(args, "gff", None):
            raise ValueError("Provide --fasta and --gff files/directories, plus --species-tree")
        resources = pair_resources(args.fasta, args.gff)
        if {r.species for r in resources} != set(panel):
            raise ValueError("FASTA/GFF species filename stems must match all species-tree tips")
        orthologs = getattr(args, "orthologs", None)
        if orthologs:
            paths = expand_files(orthologs, FASTA_SUFFIXES + (".faa",))
            if getattr(args, "family_id", None) and len(paths) != 1:
                raise ValueError("--family-id applies only to a single ortholog FASTA")
            families = {}
            for path in paths:
                family = getattr(args, "family_id", None) or resource_name(path, FASTA_SUFFIXES + (".faa",))
                if family in families:
                    raise ValueError(f"Duplicate family filename: {family}")
                families[family] = path
            targets = targets_from_fastas(resources, families)
            mode = "genomic_fasta_gff_with_ortholog_fasta"
        else:
            targets = _single_locus_targets(resources, getattr(args, "family_id", None) or "target_gene")
            mode = "single_locus_genomic_fasta_gff"
        by_species = {r.species: r for r in resources}
        rows = [{"species": t.species, "family_id": t.family, "gene_id": t.gene,
                 "gene_copy_id": t.gene, "case_id": t.family,
                 "genome_fasta": str(by_species[t.species].fasta),
                 "annotation_file": str(by_species[t.species].gff),
                 "ortholog_fasta": t.source, "source_member_ids": ";".join(t.members),
                 "input_mode": mode} for t in targets]
    _validate_rows(rows, panel)
    return InputSelection(tuple(rows), mode, panel)
