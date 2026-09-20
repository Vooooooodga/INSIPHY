"""Gene manifests with paths resolved relative to the manifest file."""
from pathlib import Path
from ..storage.tabular import read_tsv

REQUIRED = ["species", "family_id", "gene_id", "genome_fasta", "annotation_file"]


def load_manifest(path):
    source = Path(path).resolve()
    rows = read_tsv(source, REQUIRED)
    if not rows:
        raise ValueError("Gene manifest contains no targets")
    for index, row in enumerate(rows, start=2):
        for key in REQUIRED:
            if row.get(key) in {None, "", "NA", "TBD"}:
                raise ValueError(f"{source}:{index}: missing {key}")
        for field, fallback in (("case_id", "family_id"), ("gene_copy_id", "gene_id")):
            if row.get(field) in {None, "", "NA"}:
                row[field] = row[fallback]
        for key in ("genome_fasta", "annotation_file"):
            resource = Path(row[key]).expanduser()
            resource = resource if resource.is_absolute() else source.parent / resource
            if not resource.is_file():
                raise FileNotFoundError(f"{source}:{index}: {key} does not exist: {resource}")
            row[key] = str(resource.resolve())
    return rows
