"""aligners / formats: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from insiphy.aligners.types import AlignmentBackendError
from pathlib import Path
import re


def _write_temp_fasta(path: Path, name: str, seq: str):
    path.write_text(f">{name}\n{(seq or '').upper()}\n")


def _write_temp_fasta_records(path: Path, records: dict[str, str]):
    lines = []
    for name, seq in records.items():
        safe_name = str(name).replace("\t", "_").replace("\n", "_").strip() or "sequence"
        lines.append(f">{safe_name}")
        lines.append((seq or "").upper())
    path.write_text("\n".join(lines) + "\n")


def _parse_fasta_records(text: str, backend: str) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if not name or name in records:
                raise AlignmentBackendError(f"{backend} returned invalid or duplicate FASTA identifiers")
            records[name] = ""
        elif name:
            records[name] += line.strip().upper()
    return records


def _parse_paf_tags(fields):
    tags = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def _paired_bases_from_cigar(cigar):
    if not cigar or cigar == "NA":
        return 0
    return sum(int(length) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in {"M", "=", "X"})
