"""Import one single-copy orthogroup without re-inferring gene homology."""

from __future__ import annotations

from pathlib import Path

from Bio import Phylo

from .io import read_tsv, write_tsv


MANIFEST_FIELDS = [
    "case_id",
    "species",
    "family_id",
    "gene_id",
    "gene_copy_id",
    "genome_fasta",
    "annotation_file",
    "assembly",
    "annotation",
    "source_url",
    "release",
    "notes",
]


def _orthogroups_table(orthofinder_dir):
    root = Path(orthofinder_dir)
    candidates = [
        root / "Orthogroups" / "Orthogroups.tsv",
        root / "Orthogroups.tsv",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise SystemExit(f"Orthogroups.tsv not found below {root}")


def _single_gene(value):
    genes = [part.strip() for part in str(value or "").split(",") if part.strip()]
    return genes


def _name_key(value):
    return str(value or "").strip().replace(" ", "_")


def _write_species_tree(tree_path, output_path):
    path = Path(tree_path)
    if path.suffix.lower() == ".tsv":
        rows = read_tsv(path, ["node_id", "parent_id", "label"])
        fields = ["node_id", "parent_id", "label", "branch_length"]
        write_tsv(output_path, rows, fields)
        return
    tree = Phylo.read(str(path), "newick")
    rows = []
    node_ids = {}
    counter = 0
    for clade in tree.find_clades(order="preorder"):
        counter += 1
        node_ids[clade] = f"node_{counter}"
    for clade in tree.find_clades(order="preorder"):
        parent = ""
        if clade is not tree.root:
            path_to_clade = tree.get_path(clade)
            parent_clade = tree.root if len(path_to_clade) == 1 else path_to_clade[-2]
            parent = node_ids[parent_clade]
        label = clade.name or node_ids[clade]
        rows.append(
            {
                "node_id": node_ids[clade],
                "parent_id": parent,
                "label": label,
                "branch_length": 0.0 if not parent else clade.branch_length if clade.branch_length is not None else "NA",
            }
        )
    write_tsv(output_path, rows, ["node_id", "parent_id", "label", "branch_length"])


def import_orthofinder(orthofinder_dir, orthogroup, genome_manifest, output_dir, species_tree=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    resource_rows = read_tsv(genome_manifest, ["species", "genome_fasta", "annotation_file"])
    resources = {_name_key(row["species"]): row for row in resource_rows}
    table = _orthogroups_table(orthofinder_dir)
    orthogroups = read_tsv(table)
    selected = [row for row in orthogroups if row.get("Orthogroup") == orthogroup]
    if len(selected) != 1:
        raise SystemExit(f"orthogroup {orthogroup} was not found exactly once in {table}")
    row = selected[0]
    manifest_rows = []
    excluded = []
    for species_column, value in row.items():
        if species_column == "Orthogroup":
            continue
        genes = _single_gene(value)
        if len(genes) != 1:
            excluded.append(
                {
                    "family_id": orthogroup,
                    "species": species_column,
                    "copy_count": len(genes),
                    "gene_copy_ids": ";".join(genes) or "NA",
                    "reason": "orthofinder_family_is_not_single_copy",
                }
            )
            continue
        resource = resources.get(_name_key(species_column))
        if resource is None:
            raise SystemExit(f"species {species_column} from Orthofinder is missing from genome manifest")
        manifest_rows.append(
            {
                "case_id": orthogroup,
                "species": resource["species"],
                "family_id": orthogroup,
                "gene_id": genes[0],
                "gene_copy_id": genes[0],
                "genome_fasta": resource["genome_fasta"],
                "annotation_file": resource["annotation_file"],
                "assembly": resource.get("assembly", "NA"),
                "annotation": resource.get("annotation", "NA"),
                "source_url": resource.get("source_url", "NA"),
                "release": resource.get("release", "NA"),
                "notes": "imported_from_orthofinder_single_copy_orthogroup",
            }
        )
    if excluded:
        write_tsv(
            output_dir / "excluded_families.tsv",
            excluded,
            ["family_id", "species", "copy_count", "gene_copy_ids", "reason"],
        )
        raise SystemExit(f"orthogroup {orthogroup} is not single-copy in every Orthofinder species")
    if not manifest_rows:
        raise SystemExit(f"orthogroup {orthogroup} contains no usable species")
    write_tsv(output_dir / "manifest.tsv", manifest_rows, MANIFEST_FIELDS)
    if species_tree:
        _write_species_tree(species_tree, output_dir / "species_tree.tsv")
    return manifest_rows
