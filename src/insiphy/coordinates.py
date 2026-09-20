"""Coordinate primitives shared by alignment and structural correspondence.

User-facing genome coordinates are 1-based and closed. Alignment coordinates
are 0-based and half-open. Conversion between the two conventions belongs at
the input/output boundary so that interval arithmetic has one internal form.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True, order=True)
class Interval0:
    """A 0-based, half-open interval."""

    start0: int
    end0: int

    def __post_init__(self):
        if self.start0 < 0 or self.end0 < self.start0:
            raise ValueError(f"invalid 0-based half-open interval: {self.start0}-{self.end0}")

    @property
    def length(self):
        return self.end0 - self.start0

    def overlaps(self, other: "Interval0"):
        return max(self.start0, other.start0) < min(self.end0, other.end0)

    def intersection(self, other: "Interval0") -> Optional["Interval0"]:
        start0 = max(self.start0, other.start0)
        end0 = min(self.end0, other.end0)
        return Interval0(start0, end0) if start0 < end0 else None


@dataclass(frozen=True, order=True)
class ClosedInterval1:
    """A 1-based, closed interval used by GFF and public TSV files."""

    start: int
    end: int

    def __post_init__(self):
        if self.start < 1 or self.end < self.start:
            raise ValueError(f"invalid 1-based closed interval: {self.start}-{self.end}")

    @property
    def length(self):
        return self.end - self.start + 1

    def to_interval0(self):
        return Interval0(self.start - 1, self.end)

    @classmethod
    def from_interval0(cls, interval):
        if interval.length == 0:
            raise ValueError("an empty half-open interval has no 1-based closed representation")
        return cls(interval.start0 + 1, interval.end0)


@dataclass(frozen=True)
class CoordinateBlock:
    """One ungapped query/target correspondence block in internal coordinates."""

    query: Interval0
    target: Interval0

    def __post_init__(self):
        if self.query.length != self.target.length:
            raise ValueError("aligned coordinate blocks must span equal query and target lengths")


def local_interval_to_genome(local, locus, strand):
    """Map a transcript-oriented local interval to forward genome coordinates."""

    if strand not in {"+", "-"}:
        raise ValueError(f"invalid strand for coordinate projection: {strand}")
    if local.end0 > locus.length:
        raise ValueError(
            f"local interval {local.start0}-{local.end0} exceeds locus length {locus.length}"
        )
    if strand == "+":
        return Interval0(locus.start0 + local.start0, locus.start0 + local.end0)
    return Interval0(locus.end0 - local.end0, locus.end0 - local.start0)


def genome_interval_to_local(genome, locus, strand):
    """Map a forward genome interval to transcript-oriented local coordinates."""

    if strand not in {"+", "-"}:
        raise ValueError(f"invalid strand for coordinate projection: {strand}")
    if genome.start0 < locus.start0 or genome.end0 > locus.end0:
        raise ValueError("genome interval is outside the containing locus")
    if strand == "+":
        return Interval0(genome.start0 - locus.start0, genome.end0 - locus.start0)
    return Interval0(locus.end0 - genome.end0, locus.end0 - genome.start0)


def parse_legacy_block(text):
    """Parse one legacy 1-based closed ``qstart-qend:tstart-tend`` block."""

    query_text, target_text = text.split(":", 1)
    query_start, query_end = (int(value) for value in query_text.split("-", 1))
    target_start, target_end = (int(value) for value in target_text.split("-", 1))
    return CoordinateBlock(
        ClosedInterval1(query_start, query_end).to_interval0(),
        ClosedInterval1(target_start, target_end).to_interval0(),
    )


def parse_legacy_blocks(text):
    if text in {None, "", "NA", "."}:
        return tuple()
    return tuple(parse_legacy_block(block) for block in str(text).split(";") if block)


def format_legacy_blocks(blocks):
    if not blocks:
        return "NA"
    formatted = []
    for block in blocks:
        query = ClosedInterval1.from_interval0(block.query)
        target = ClosedInterval1.from_interval0(block.target)
        formatted.append(f"{query.start}-{query.end}:{target.start}-{target.end}")
    return ";".join(formatted)
