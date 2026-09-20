"""aligners / columns: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from insiphy.aligners.runner import _backend_version
from insiphy.aligners.types import AlignmentBackendError
from insiphy.aligners.types import AlignmentStats
from insiphy.aligners.types import DNA_COMPLEMENT
from insiphy.aligners.types import KNOWN_NT
from typing import Optional
import re


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1].upper()


def ungapped_identity(seq_a: str, seq_b: str) -> float:
    if not seq_a or not seq_b:
        return 0.0
    n = min(len(seq_a), len(seq_b))
    matches = mismatches = 0
    for a, b in zip(seq_a[:n].upper(), seq_b[:n].upper()):
        a = "T" if a == "U" else a
        b = "T" if b == "U" else b
        if a not in KNOWN_NT or b not in KNOWN_NT:
            continue
        matches += int(a == b)
        mismatches += int(a != b)
    return matches / max(1, matches + mismatches)


def _compress_ops(ops):
    if not ops:
        return "NA"
    out = []
    last = ops[0]
    count = 1
    for op in ops[1:]:
        if op == last:
            count += 1
        else:
            out.append(f"{count}{last}")
            last = op
            count = 1
    out.append(f"{count}{last}")
    return "".join(out)


def _parse_cigar(cigar: str):
    if not cigar or cigar == "NA":
        return []
    return [(int(length or "1"), op) for length, op in re.findall(r"(\d*)([MIDNSHP=X])", cigar)]


def _nt(base: str) -> str:
    base = (base or "").upper()
    return "T" if base == "U" else base


def _known_base(base: str) -> bool:
    return _nt(base) in KNOWN_NT


def _empty_stats(
    query_len: int,
    target_len: int,
    backend: str,
    strand: str = "+",
    alignment_mode: str = "none",
    alignment_meaning: str = "no alignment reported",
) -> AlignmentStats:
    stats = AlignmentStats(
        0.0, 0.0, 0.0,
        query_end=query_len,
        target_end=target_len,
        backend=backend,
        strand=strand,
        alignment_mode=alignment_mode,
        alignment_meaning=alignment_meaning,
    )
    stats.backend_version = _backend_version(backend)
    stats.raw_score = 0.0
    stats.nt_identity = 0.0
    stats.query_length = query_len
    stats.target_length = target_len
    stats.relative_strand = strand
    return stats


def _stats_from_alignment_columns(
    query_aligned: str,
    target_aligned: str,
    query_len: int,
    target_len: int,
    *,
    query_start0: int = 0,
    query_end0: Optional[int] = None,
    target_start0: int = 0,
    strand: str = "+",
    score: float = 0.0,
    backend: str = "internal",
    mode: str = "global",
    alignment_mode: Optional[str] = None,
    alignment_meaning: str = "nucleotide alignment",
) -> AlignmentStats:
    if len(query_aligned) != len(target_aligned):
        raise AlignmentBackendError(f"{backend} returned unequal aligned sequence lengths")

    q_step = -1 if strand == "-" else 1
    q_pos = (query_end0 - 1) if strand == "-" and query_end0 is not None else query_start0
    t_pos = target_start0
    q_coords = []
    t_coords = []
    ops = []
    blocks: list[tuple[int, int, int, int]] = []
    active_block: Optional[list[int]] = None
    last_q = last_t = None
    matches = mismatches = gap_bases = unknown_bases = unknown_aligned_pairs = 0
    query_evidence = target_evidence = 0

    def flush_block():
        nonlocal active_block
        if active_block is not None:
            q1, q2, t1, t2 = active_block
            blocks.append((min(q1, q2), max(q1, q2), min(t1, t2), max(t1, t2)))
            active_block = None

    for q_char, t_char in zip(query_aligned, target_aligned):
        has_q = q_char != "-"
        has_t = t_char != "-"
        q_coord = t_coord = None
        if has_q:
            q_coord = q_pos + 1
            q_coords.append(q_coord)
            q_pos += q_step
        if has_t:
            t_coord = t_pos + 1
            t_coords.append(t_coord)
            t_pos += 1

        if has_q and has_t:
            q_known = _known_base(q_char)
            t_known = _known_base(t_char)
            if q_known and t_known:
                if active_block is None or last_q is None or last_t is None or q_coord != last_q + q_step or t_coord != last_t + 1:
                    flush_block()
                    active_block = [q_coord, q_coord, t_coord, t_coord]
                else:
                    active_block[1] = q_coord
                    active_block[3] = t_coord
                last_q, last_t = q_coord, t_coord
                query_evidence += 1
                target_evidence += 1
                if _nt(q_char) == _nt(t_char):
                    matches += 1
                    ops.append("=")
                else:
                    mismatches += 1
                    ops.append("X")
            else:
                flush_block()
                last_q = last_t = None
                unknown_bases += int(not q_known) + int(not t_known)
                unknown_aligned_pairs += 1
                ops.append("M")
        elif has_q:
            flush_block()
            last_q = last_t = None
            if _known_base(q_char):
                gap_bases += 1
            else:
                unknown_bases += 1
            ops.append("I")
        elif has_t:
            flush_block()
            last_q = last_t = None
            if _known_base(t_char):
                gap_bases += 1
            else:
                unknown_bases += 1
            ops.append("D")
    flush_block()

    denominator = matches + mismatches + gap_bases
    identity = matches / denominator if denominator else 0.0
    query_coverage = query_evidence / max(1, query_len)
    target_coverage = target_evidence / max(1, target_len)
    query_span_coverage = len(q_coords) / max(1, query_len)
    target_span_coverage = len(t_coords) / max(1, target_len)
    coverage = min(query_coverage, target_coverage) if mode == "global" else query_coverage
    stats = AlignmentStats(
        identity=identity,
        coverage=coverage,
        score=float(score),
        query_start=min(q_coords) if q_coords else 1,
        query_end=max(q_coords) if q_coords else 0,
        target_start=min(t_coords) if t_coords else 1,
        target_end=max(t_coords) if t_coords else 0,
        cigar=_compress_ops(ops),
        backend=backend,
        query_coverage=query_coverage,
        target_coverage=target_coverage,
        aligned_pairs=matches + mismatches,
        strand=strand,
        aligned_blocks=blocks,
        matches=matches,
        mismatches=mismatches,
        gap_bases=gap_bases,
        unknown_bases=unknown_bases,
        alignment_mode=alignment_mode or mode,
        alignment_meaning=alignment_meaning,
        query_span_coverage=query_span_coverage,
        target_span_coverage=target_span_coverage,
    )
    stats.backend_version = _backend_version(backend)
    stats.raw_score = float(score)
    stats.nt_identity = identity
    stats.known_aligned_pairs = matches + mismatches
    stats.unknown_aligned_pairs = unknown_aligned_pairs
    stats.query_covered_bases = query_evidence
    stats.target_covered_bases = target_evidence
    stats.query_length = query_len
    stats.target_length = target_len
    stats.relative_strand = strand
    return stats


def _columns_from_cigar(query_segment: str, target_segment: str, cigar: str):
    q_index = t_index = 0
    q_cols = []
    t_cols = []
    for length, op in _parse_cigar(cigar):
        if op in {"M", "=", "X"}:
            for _ in range(length):
                if q_index >= len(query_segment) or t_index >= len(target_segment):
                    raise AlignmentBackendError("CIGAR consumes beyond aligned sequence span")
                q_cols.append(query_segment[q_index])
                t_cols.append(target_segment[t_index])
                q_index += 1
                t_index += 1
        elif op == "I":
            for _ in range(length):
                if q_index >= len(query_segment):
                    raise AlignmentBackendError("CIGAR insertion consumes beyond query span")
                q_cols.append(query_segment[q_index])
                t_cols.append("-")
                q_index += 1
        elif op in {"D", "N"}:
            for _ in range(length):
                if t_index >= len(target_segment):
                    raise AlignmentBackendError("CIGAR deletion consumes beyond target span")
                q_cols.append("-")
                t_cols.append(target_segment[t_index])
                t_index += 1
        elif op == "S":
            q_index += length
        elif op in {"H", "P"}:
            continue
    return "".join(q_cols), "".join(t_cols)
