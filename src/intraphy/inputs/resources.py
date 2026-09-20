"""Match per-species genomic FASTA and GFF/GTF resources by filename."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

FASTA_SUFFIXES = (".fasta", ".fna", ".fa")
GFF_SUFFIXES = (".gff3", ".gff", ".gtf")


@dataclass(frozen=True)
class SpeciesFiles:
    species: str
    fasta: Path
    gff: Path


def resource_name(path: Path, suffixes: tuple[str, ...]) -> str:
    """Remove only the recognized extension; species labels stay unchanged."""
    name = path.name[:-3] if path.name.lower().endswith(".gz") else path.name
    for suffix in suffixes:
        if name.lower().endswith(suffix):
            label = name[:-len(suffix)]
            if label and not any(c.isspace() for c in label):
                return label
    raise ValueError(f"Unsupported filename or whitespace in species label: {path}")


def expand_files(values: Iterable[str], suffixes: tuple[str, ...]) -> tuple[Path, ...]:
    """Accept explicit paths or flat directories; never recurse unexpectedly."""
    files: list[Path] = []
    for value in values:
        path = Path(value).expanduser().resolve()
        if path.is_dir():
            for child in sorted(path.iterdir()):
                if not child.is_file():
                    continue
                try:
                    resource_name(child, suffixes)
                except ValueError:
                    continue
                files.append(child.resolve())
        elif path.is_file():
            resource_name(path, suffixes)
            files.append(path)
        else:
            raise FileNotFoundError(f"Input file or directory does not exist: {value}")
    if not files:
        raise ValueError("No input files with the required extensions were found")
    if len(files) != len(set(files)):
        raise ValueError("An input resource was supplied more than once")
    return tuple(files)


def _by_species(values: Iterable[str], suffixes: tuple[str, ...]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for path in expand_files(values, suffixes):
        label = resource_name(path, suffixes)
        if label in result:
            raise ValueError(f"Multiple files for species {label}: {result[label]} and {path}")
        result[label] = path
    return result


def pair_resources(fasta: Iterable[str], gff: Iterable[str]) -> tuple[SpeciesFiles, ...]:
    sequences = _by_species(fasta, FASTA_SUFFIXES)
    annotations = _by_species(gff, GFF_SUFFIXES)
    if len(sequences) == 1 and len(annotations) > 1:
        # A combined genomic FASTA is also valid when GFF sequence IDs identify
        # distinct species records. Selection validates IDs and coordinates.
        source = next(iter(sequences.values()))
        _check_combined_fasta(source)
        return tuple(SpeciesFiles(label, source, annotations[label])
                     for label in sorted(annotations))
    if set(sequences) != set(annotations):
        raise ValueError(
            "FASTA/GFF filenames must have the same species stem. "
            f"Without GFF: {sorted(set(sequences)-set(annotations))}; "
            f"without FASTA: {sorted(set(annotations)-set(sequences))}. "
            "Use, for example, Species_A.fa and Species_A.gff3."
        )
    return tuple(SpeciesFiles(label, sequences[label], annotations[label])
                 for label in sorted(sequences))


def _check_combined_fasta(path: Path) -> None:
    from ..storage.tabular import open_text
    identifiers = set()
    with open_text(path) as stream:
        for line in stream:
            if not line.startswith(">"):
                continue
            header = line[1:].split()
            if not header or header[0] in identifiers:
                raise ValueError(f"Combined genomic FASTA has an empty or duplicate record ID: {line.strip()}")
            identifiers.add(header[0])
    if not identifiers:
        raise ValueError("Combined genomic FASTA contains no sequence records")
