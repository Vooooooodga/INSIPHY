"""Prepared genomic sources to complete exon instances, without tree inference."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from ..storage.fasta import parse_fasta
from ..storage.tabular import read_tsv
from .types import ExonInstance


@dataclass(frozen=True)
class NativeLocus:
    family: str
    species: str
    locus: str
    contig: str
    strand: str
    search_start: int
    search_end: int
    sequence: str
    exons: tuple[ExonInstance, ...]
    paths: dict[str, tuple[str, ...]]
    partial: bool
    coding_by_transcript: dict[str, tuple[tuple[int, int, int | None], ...]] = field(default_factory=dict)
    translation_exceptions: tuple[str, ...] = ()

    def oriented_interval(self, exon: ExonInstance) -> tuple[int, int]:
        if self.strand == "+":
            return exon.start-(self.search_start-1), exon.end-(self.search_start-1)
        return self.search_end-exon.end, self.search_end-exon.start

    def source_position(self, oriented_base: int) -> int:
        """Return 0-based source-genomic base position."""
        return self.search_start-1+oriented_base if self.strand == "+" else self.search_end-1-oriented_base


def _cds_intervals(row):
    value = row.get("cds_intervals", "")
    phase = row.get("cds_phase", ".")
    phase = int(phase) if str(phase) in {"0", "1", "2"} else None
    result = []
    for part in value.split(";"):
        if "-" in part:
            left, right = part.split("-")
            result.append((int(left)-1, int(right), phase))
    return tuple(result)


def load_native(input_dir: str | Path) -> tuple[NativeLocus, ...]:
    root = Path(input_dir)
    loci = read_tsv(root / "gene_loci.tsv")
    paths = read_tsv(root / "transcript_paths.tsv")
    sequences = parse_fasta(root / "gene_loci.fasta")
    raw = read_tsv(root / "raw_gene_features.tsv", optional=True)
    metadata = defaultdict(list)
    for r in raw:
        metadata[(r["species"], r["gene_copy_id"])].append(r)
    exon_rows = [r for r in paths if r.get("role") == "exon"]
    by_locus = defaultdict(list)
    for row in exon_rows:
        by_locus[(row["species"], row["gene_copy_id"])].append(row)
    result = []
    for locus in loci:
        rows = by_locus.get((locus["species"], locus["gene_copy_id"]), [])
        raw_locus = metadata[(locus["species"], locus["gene_copy_id"])]
        families = {r["family_id"] for r in (rows or raw_locus)}
        if len(families) != 1:
            raise ValueError("Prepared locus belongs to multiple families")
        prefix = f"{locus['species']}|{locus['gene_copy_id']}|"
        keys = [key for key in sequences if key.startswith(prefix)]
        if len(keys) != 1:
            raise ValueError(f"Expected one prepared genomic sequence for {prefix}")
        by_span = defaultdict(list)
        for row in rows:
            by_span[(int(row["start"])-1, int(row["end"]))].append(row)
        instances = []
        by_tx = defaultdict(list)
        for (start, end), evidence in sorted(by_span.items()):
            occurrence = sorted({r["occurrence_id"] for r in evidence})[0]
            transcripts = tuple(sorted({r["transcript_id"] for r in evidence}))
            cds = tuple(sorted({p for r in evidence for p in _cds_intervals(r)}, key=lambda p: (p[0], p[1], str(p[2]))))
            instance = ExonInstance(occurrence, next(iter(families)), locus["species"], locus["gene_copy_id"],
                locus["contig"], start, end, locus["strand"], transcripts, cds,
                ";".join(sorted({r.get("annotation_source", "provided") for r in evidence})))
            instances.append(instance)
            for tx in transcripts:
                by_tx[tx].append(instance)
        reverse = locus["strand"] == "-"
        ordered_paths = {tx: tuple(e.id for e in sorted(es, key=lambda e: e.start, reverse=reverse)) for tx, es in by_tx.items()}
        coding = {}
        for tx in ordered_paths:
            coding[tx] = tuple(sorted({p for r in rows if r["transcript_id"] == tx
                                      for p in _cds_intervals(r)}, key=lambda p: (p[0], p[1], str(p[2])), reverse=reverse))
        exceptions = set()
        for r in raw_locus:
            attrs = r.get("attrs", "")
            if any(token in attrs for token in ("transl_except=", "translation_exception=", "exception=")):
                exceptions.add(attrs)
            for token in attrs.split(";"):
                if token.startswith(("transl_table=", "translation_table=")) and token.split("=", 1)[1] != "1":
                    exceptions.add(token)
        sequence = sequences[keys[0]]
        expected = int(locus["search_end"])-int(locus["search_start"])+1
        if len(sequence) != expected:
            raise ValueError("Prepared genomic sequence length does not match its coordinate record")
        partial = any(r.get("partial_start") == "1" or r.get("partial_end") == "1" for r in rows)
        result.append(NativeLocus(next(iter(families)), locus["species"], locus["gene_copy_id"],
            locus["contig"], locus["strand"], int(locus["search_start"]), int(locus["search_end"]),
            sequence, tuple(instances), ordered_paths, partial, coding, tuple(sorted(exceptions))))
    keys = [(r.family, r.species) for r in result]
    if len(keys) != len(set(keys)):
        raise ValueError("V19 exon configurations require one orthologous locus per family and species")
    return tuple(result)
