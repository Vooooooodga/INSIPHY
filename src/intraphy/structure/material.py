"""Material-aware geometry; a deleted source is never silently resurrected."""
from __future__ import annotations

from .types import Catalogue, ExonConfiguration, ExonSpan


def retained_intervals(catalogue: Catalogue, states: tuple[int, ...],
                       start: int, end: int) -> tuple[ExonSpan, ...]:
    if len(states) != len(catalogue.material):
        raise ValueError("Material dimension mismatch")
    if start >= end:
        return ()
    missing = sorted((m.start, m.end) for m, s in zip(catalogue.material, states) if s != 1)
    cursor, parts = start, []
    for left, right in missing:
        if right <= cursor or left >= end:
            continue
        if cursor < left:
            parts.append(ExonSpan(cursor, min(left, end)))
        cursor = max(cursor, right)
    if cursor < end:
        parts.append(ExonSpan(cursor, end))
    return tuple(parts)


def retained_length(catalogue: Catalogue, states: tuple[int, ...], start: int, end: int) -> int:
    return sum(p.end - p.start for p in retained_intervals(catalogue, states, start, end))


def normalize_exons(catalogue: Catalogue, states: tuple[int, ...],
                    exons: tuple[ExonSpan, ...]) -> tuple[ExonSpan, ...]:
    """Trim deleted ends and merge only when *genomic* separating material is gone.

    A hole inside a retained exon is an indel, not an extra exon. Removing a
    whole intervening intron joins the retained exon pieces in one operation.
    """
    kept: list[ExonSpan] = []
    for exon in exons:
        parts = retained_intervals(catalogue, states, exon.start, exon.end)
        if not parts:
            continue
        e = ExonSpan(parts[0].start, parts[-1].end)
        if kept and retained_length(catalogue, states, kept[-1].end, e.start) == 0:
            kept[-1] = ExonSpan(kept[-1].start, e.end)
        else:
            kept.append(e)
    return tuple(kept)


def valid_configuration(catalogue: Catalogue, config: ExonConfiguration) -> bool:
    if len(config.material) != len(catalogue.material):
        return False
    if any(e.end > catalogue.length for e in config.exons):
        return False
    return normalize_exons(catalogue, config.material, config.exons) == config.exons
