"""Alignment and splice-boundary helpers for INSIPHY."""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from urllib.parse import unquote


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
KNOWN_NT = {"A", "C", "G", "T"}
MAX_INTERNAL_DP_CELLS = 250_000
MAX_MAFFT_PAIR_CELLS = 25_000_000


@dataclass
class AlignmentStats:
    identity: float
    coverage: float
    score: float
    query_start: int = 1
    query_end: int = 0
    target_start: int = 1
    target_end: int = 0
    cigar: str = "NA"
    backend: str = "internal"
    query_coverage: float = 0.0
    target_coverage: float = 0.0
    aligned_pairs: int = 0
    strand: str = "+"
    aligned_blocks: list[tuple[int, int, int, int]] = field(default_factory=list)
    matches: int = 0
    mismatches: int = 0
    gap_bases: int = 0
    unknown_bases: int = 0
    alignment_mode: str = "global"
    alignment_meaning: str = "nucleotide alignment"
    query_span_coverage: float = 0.0
    target_span_coverage: float = 0.0


class AlignmentBackendError(RuntimeError):
    """Raised when an explicitly requested external alignment backend fails."""


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
    return AlignmentStats(
        0.0, 0.0, 0.0,
        query_end=query_len,
        target_end=target_len,
        backend=backend,
        strand=strand,
        alignment_mode=alignment_mode,
        alignment_meaning=alignment_meaning,
    )


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
    matches = mismatches = gap_bases = unknown_bases = 0
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
    return AlignmentStats(
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


def _write_temp_fasta(path: Path, name: str, seq: str):
    path.write_text(f">{name}\n{(seq or '').upper()}\n")


def _write_temp_fasta_records(path: Path, records: dict[str, str]):
    lines = []
    for name, seq in records.items():
        safe_name = str(name).replace("\t", "_").replace("\n", "_").strip() or "sequence"
        lines.append(f">{safe_name}")
        lines.append((seq or "").upper())
    path.write_text("\n".join(lines) + "\n")


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


def _external_minimap2_stats(query: str, target: str, mode: str, threads: int = 1) -> AlignmentStats:
    exe = shutil.which("minimap2")
    if not exe:
        raise AlignmentBackendError("minimap2 was requested but is not available on PATH")
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return _empty_stats(len(query), len(target), "minimap2")
    with tempfile.TemporaryDirectory(prefix="insiphy_minimap2_") as tmp:
        tmp = Path(tmp)
        query_path = tmp / "query.fa"
        target_path = tmp / "target.fa"
        _write_temp_fasta(query_path, "query", query)
        _write_temp_fasta(target_path, "target", target)
        cmd = [
            exe,
            "-c",
            "-x",
            "asm20",
            "-t",
            str(max(1, int(threads or 1))),
            str(target_path),
            str(query_path),
        ]
        proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise AlignmentBackendError(proc.stderr.strip() or "minimap2 failed")
    best = None
    for line in proc.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        tstart = int(fields[7])
        tend = int(fields[8])
        strand = fields[4]
        tags = _parse_paf_tags(fields)
        cigar = tags.get("cg", "NA")
        if cigar == "NA":
            continue
        query_segment = query[qstart:qend]
        if strand == "-":
            query_segment = revcomp(query_segment)
        target_segment = target[tstart:tend]
        query_cols, target_cols = _columns_from_cigar(query_segment, target_segment, cigar)
        score = float(tags.get("AS", fields[9]))
        stat = _stats_from_alignment_columns(
            query_cols,
            target_cols,
            qlen,
            len(target),
            query_start0=qstart,
            query_end0=qend,
            target_start0=tstart,
            strand=strand,
            score=score,
            backend="minimap2",
            mode=mode,
        )
        rank = (stat.coverage * stat.identity, stat.score)
        if best is None or rank > best[0]:
            best = (rank, stat)
    if best is None:
        return _empty_stats(len(query), len(target), "minimap2")
    return best[1]


def _external_miniprot_stats(query: str, target: str, threads: int = 1) -> AlignmentStats:
    query = (query or "").upper()
    target = (target or "").upper()
    if set(query) <= set("ACGTUN-"):
        raise AlignmentBackendError("miniprot requires an amino-acid query; a nucleotide query was supplied")
    exe = shutil.which("miniprot")
    if not exe:
        raise AlignmentBackendError("miniprot was requested but is not available on PATH")
    if not query or not target:
        return _empty_stats(len(query), len(target), "miniprot")
    with tempfile.TemporaryDirectory(prefix="insiphy_miniprot_") as tmp:
        tmp = Path(tmp)
        query_path = tmp / "query.fa"
        target_path = tmp / "target.fa"
        _write_temp_fasta(query_path, "query", query)
        _write_temp_fasta(target_path, "target", target)
        cmd = [exe, "-t", str(max(1, int(threads or 1))), str(target_path), str(query_path)]
        proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise AlignmentBackendError(proc.stderr.strip() or "miniprot failed")
    best = None
    for line in proc.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        strand = fields[4]
        tstart = int(fields[7])
        tend = int(fields[8])
        matching_nt = float(fields[9])
        aligned_nt = max(1.0, float(fields[10]))
        tags = _parse_paf_tags(fields)
        identity = matching_nt / aligned_nt
        query_coverage = (qend - qstart) / max(1, qlen)
        target_coverage = aligned_nt / max(1, len(target))
        score = float(tags.get("AS", matching_nt))
        stat = AlignmentStats(
            identity, query_coverage, score, qstart + 1, qend, min(tstart, tend) + 1,
            max(tstart, tend), tags.get("cg", "NA"), "miniprot",
            query_coverage, target_coverage, 0, strand,
        )
        stat.alignment_mode = "protein_to_genome"
        stat.alignment_meaning = "miniprot protein-to-genome alignment; query coordinates are amino-acid and target coordinates are nucleotide"
        rank = (stat.coverage * stat.identity, stat.score)
        if best is None or rank > best[0]:
            best = (rank, stat)
    if best is None:
        return _empty_stats(len(query), len(target), "miniprot")
    return best[1]


def _mafft_pair_alignment(query: str, target: str, threads: int = 1) -> tuple[str, str]:
    exe = shutil.which("mafft")
    if not exe:
        raise AlignmentBackendError("MAFFT was requested but is not available on PATH")
    with tempfile.TemporaryDirectory(prefix="insiphy_mafft_") as tmp:
        input_path = Path(tmp) / "pair.fa"
        input_path.write_text(f">query\n{query.upper()}\n>target\n{target.upper()}\n")
        proc = subprocess.run(
            [exe, "--quiet", "--thread", str(max(1, int(threads or 1))), "--nuc", "--auto", str(input_path)],
            text=True,
            capture_output=True,
            check=False,
        )
    if proc.returncode != 0:
        raise AlignmentBackendError(proc.stderr.strip() or "MAFFT failed")
    aligned = {}
    name = None
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            aligned[name] = ""
        elif name:
            aligned[name] += line.strip().upper()
    left, right = aligned.get("query", ""), aligned.get("target", "")
    if not left or len(left) != len(right):
        raise AlignmentBackendError("MAFFT did not return a valid two-sequence alignment")
    return left, right


def _external_mafft_stats(query: str, target: str, threads: int = 1) -> AlignmentStats:
    left, right = _mafft_pair_alignment(query, target, threads)
    return _stats_from_alignment_columns(
        left,
        right,
        len(query),
        len(target),
        score=0.0,
        backend="mafft",
        mode="global",
        alignment_mode="global",
        alignment_meaning="MAFFT global nucleotide alignment; terminal overhangs are counted as gaps",
    )


def _trim_terminal_overhangs(query_aligned: str, target_aligned: str) -> tuple[str, str, int, int]:
    paired = [idx for idx, (q_char, t_char) in enumerate(zip(query_aligned, target_aligned)) if q_char != "-" and t_char != "-"]
    if not paired:
        return "", "", 0, 0
    first = paired[0]
    last = paired[-1]
    query_offset = sum(1 for char in query_aligned[:first] if char != "-")
    target_offset = sum(1 for char in target_aligned[:first] if char != "-")
    return query_aligned[first : last + 1], target_aligned[first : last + 1], query_offset, target_offset


def overlap_alignment_stats(query: str, target: str, backend: str = "mafft", threads: int = 1) -> AlignmentStats:
    """Project a pairwise overlap by trimming terminal overhangs from an external global alignment."""
    backend = (backend or "mafft").lower()
    if backend != "mafft":
        raise AlignmentBackendError(f"unsupported overlap alignment backend: {backend}")
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return _empty_stats(
            len(query),
            len(target),
            "mafft_overlap",
            alignment_mode="overlap_projection",
            alignment_meaning="terminal-overhang-trimmed MAFFT overlap projection",
        )
    left, right = _mafft_pair_alignment(query, target, threads)
    trimmed_query, trimmed_target, query_offset, target_offset = _trim_terminal_overhangs(left, right)
    return _stats_from_alignment_columns(
        trimmed_query,
        trimmed_target,
        len(query),
        len(target),
        query_start0=query_offset,
        target_start0=target_offset,
        score=0.0,
        backend="mafft_overlap",
        mode="overlap",
        alignment_mode="overlap_projection",
        alignment_meaning="terminal-overhang-trimmed MAFFT overlap projection; not a true local algorithm",
    )


def _external_lastz_stats(query: str, target: str, mode: str, threads: int = 1) -> AlignmentStats:
    exe = shutil.which("lastz")
    if not exe:
        raise AlignmentBackendError("LASTZ was requested but is not available on PATH")
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return _empty_stats(len(query), len(target), "lastz")
    fields = "score,strand1,size1,zstart1,end1,strand2,size2,zstart2+,end2+,text1,text2,cigarx"
    with tempfile.TemporaryDirectory(prefix="insiphy_lastz_") as tmp:
        tmp = Path(tmp)
        query_path = tmp / "query.fa"
        target_path = tmp / "target.fa"
        _write_temp_fasta(query_path, "query", query)
        _write_temp_fasta(target_path, "target", target)
        cmd = [exe, str(target_path), str(query_path), "--strand=both", f"--format=general-:{fields}"]
        proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise AlignmentBackendError(proc.stderr.strip() or "LASTZ failed")
    best = None
    for line in proc.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 12:
            continue
        score, _strand1, _size1, tstart0, _tend, strand2, _size2, qstart0, qend, target_cols, query_cols, _cigarx = parts
        stat = _stats_from_alignment_columns(
            query_cols.upper(),
            target_cols.upper(),
            len(query),
            len(target),
            query_start0=int(qstart0),
            query_end0=int(qend),
            target_start0=int(tstart0),
            strand=strand2,
            score=float(score),
            backend="lastz",
            mode=mode,
        )
        rank = (stat.coverage * stat.identity, stat.score)
        if best is None or rank > best[0]:
            best = (rank, stat)
    if best is None:
        return _empty_stats(len(query), len(target), "lastz")
    return best[1]


def available_alignment_backends():
    rows = [{"aligner": "auto", "available": 1, "notes": "global uses MAFFT; local selects LASTZ/minimap2 when available"}]
    rows.append({"aligner": "internal", "available": 1, "notes": "Biopython PairwiseAligner, explicit small-pair backend"})
    rows.append({"aligner": "minimap2", "available": int(shutil.which("minimap2") is not None), "notes": "external nucleotide aligner"})
    rows.append({"aligner": "miniprot", "available": int(shutil.which("miniprot") is not None), "notes": "external protein-to-genome aligner"})
    rows.append({"aligner": "mafft", "available": int(shutil.which("mafft") is not None), "notes": "external global nucleotide aligner; local requests use terminal-overhang-trimmed overlap projection"})
    rows.append({"aligner": "lastz", "available": int(shutil.which("lastz") is not None), "notes": "external local nucleotide aligner"})
    return rows


def _pairwise_aligner_stats(seq_a: str, seq_b: str, mode: str) -> AlignmentStats:
    from Bio.Align import PairwiseAligner

    seq_a = (seq_a or "").upper()
    seq_b = (seq_b or "").upper()
    if not seq_a or not seq_b:
        return _empty_stats(len(seq_a), len(seq_b), "internal")
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(seq_a, seq_b)
    alignment = next(iter(alignments), None)
    if alignment is None:
        return _empty_stats(len(seq_a), len(seq_b), "internal")
    coordinates = alignment.coordinates
    query_cols = []
    target_cols = []
    for index in range(coordinates.shape[1] - 1):
        a_start, a_end = int(coordinates[0, index]), int(coordinates[0, index + 1])
        b_start, b_end = int(coordinates[1, index]), int(coordinates[1, index + 1])
        a_span = a_end - a_start
        b_span = b_end - b_start
        if a_span and b_span:
            span = min(a_span, b_span)
            query_cols.append(seq_a[a_start : a_start + span])
            target_cols.append(seq_b[b_start : b_start + span])
        elif a_span:
            query_cols.append(seq_a[a_start:a_end])
            target_cols.append("-" * a_span)
        elif b_span:
            query_cols.append("-" * b_span)
            target_cols.append(seq_b[b_start:b_end])
    return _stats_from_alignment_columns(
        "".join(query_cols),
        "".join(target_cols),
        len(seq_a),
        len(seq_b),
        query_start0=int(coordinates[0, 0]),
        target_start0=int(coordinates[1, 0]),
        score=float(alignment.score),
        backend="internal",
        mode=mode,
    )


def global_alignment_stats(seq_a: str, seq_b: str, backend: str = "internal", threads: int = 1) -> AlignmentStats:
    backend = (backend or "internal").lower()
    if backend == "auto":
        return global_alignment_stats(seq_a, seq_b, "mafft", threads)
    if backend == "minimap2":
        return _external_minimap2_stats(seq_a, seq_b, "global", threads)
    if backend == "miniprot":
        return _external_miniprot_stats(seq_a, seq_b, threads)
    if backend == "mafft":
        return _external_mafft_stats(seq_a, seq_b, threads)
    if backend == "lastz":
        raise AlignmentBackendError("LASTZ is a local nucleotide backend; use local_alignment_stats")
    if backend != "internal":
        raise AlignmentBackendError(f"unsupported alignment backend: {backend}")
    seq_a = (seq_a or "").upper()
    seq_b = (seq_b or "").upper()
    if not seq_a or not seq_b:
        return _empty_stats(len(seq_a), len(seq_b), "internal")
    if len(seq_a) * len(seq_b) > MAX_INTERNAL_DP_CELLS:
        raise AlignmentBackendError(
            f"internal global alignment requires {len(seq_a) * len(seq_b)} DP cells; "
            "select MAFFT or minimap2"
        )
    return _pairwise_aligner_stats(seq_a, seq_b, "global")


def local_alignment_stats(query: str, target: str, backend: str = "internal", threads: int = 1) -> AlignmentStats:
    backend = (backend or "internal").lower()
    if backend == "auto":
        if shutil.which("lastz"):
            return _external_lastz_stats(query, target, "local", threads)
        if shutil.which("minimap2"):
            return _external_minimap2_stats(query, target, "local", threads)
        raise AlignmentBackendError("local auto requires LASTZ or minimap2 on PATH; request internal explicitly for small Biopython alignment")
    if backend == "minimap2":
        return _external_minimap2_stats(query, target, "local", threads)
    if backend == "miniprot":
        return _external_miniprot_stats(query, target, threads)
    if backend == "mafft":
        return overlap_alignment_stats(query, target, backend="mafft", threads=threads)
    if backend == "lastz":
        return _external_lastz_stats(query, target, "local", threads)
    if backend != "internal":
        raise AlignmentBackendError(f"unsupported alignment backend: {backend}")
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return _empty_stats(len(query), len(target), "internal")
    if len(query) * len(target) > MAX_INTERNAL_DP_CELLS:
        raise AlignmentBackendError(
            f"internal local alignment requires {len(query) * len(target)} DP cells; select LASTZ or minimap2"
        )
    return _pairwise_aligner_stats(query, target, "local")


def best_ungapped_hit(query: str, target: str) -> AlignmentStats:
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return AlignmentStats(0.0, 0.0, 0.0)
    if len(target) < len(query):
        span = len(target)
        matches = sum(1 for a, b in zip(query[:span], target) if a == b and a != "N" and b != "N")
        return AlignmentStats(matches / max(1, span), span / max(1, len(query)), float(matches), 1, span, 1, span, f"{span}M")
    best = AlignmentStats(-1.0, 1.0, 0.0, 1, len(query), 1, len(query), f"{len(query)}M")
    qlen = len(query)
    for offset in range(0, len(target) - qlen + 1):
        window = target[offset : offset + qlen]
        matches = sum(1 for a, b in zip(query, window) if a == b and a != "N" and b != "N")
        identity = matches / qlen
        if identity > best.identity:
            best = AlignmentStats(identity, 1.0, float(matches), 1, qlen, offset + 1, offset + qlen, f"{qlen}M")
    return best


def _parse_gff_attributes(raw: str) -> dict[str, str]:
    attrs = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        attrs[unquote(key)] = unquote(value.strip().strip('"'))
    return attrs


def protein_locus_exons(proteins: dict[str, str], locus: str, threads=1) -> list[dict]:
    """Return miniprot CDS exon evidence rows for a batch of proteins against one locus."""
    exe = shutil.which("miniprot")
    if not exe:
        raise AlignmentBackendError("miniprot was requested but is not available on PATH")
    clean_proteins = {str(name): (seq or "").upper() for name, seq in proteins.items() if seq}
    locus = (locus or "").upper()
    if not clean_proteins or not locus:
        return []

    with tempfile.TemporaryDirectory(prefix="insiphy_miniprot_gff_") as tmp:
        tmp = Path(tmp)
        protein_path = tmp / "proteins.faa"
        locus_path = tmp / "locus.fa"
        _write_temp_fasta_records(protein_path, clean_proteins)
        _write_temp_fasta(locus_path, "locus", locus)
        cmd = [exe, "-t", str(max(1, int(threads or 1))), "--gff", str(locus_path), str(protein_path)]
        proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise AlignmentBackendError(proc.stderr.strip() or "miniprot --gff failed")

    parent_info: dict[str, dict] = {}
    cds_rows = []
    paf_cigars: dict[str, list[str]] = {}
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        if line.startswith("##PAF"):
            paf_fields = line.split("\t")[1:]
            if len(paf_fields) >= 12:
                tags = _parse_paf_tags(paf_fields)
                cigar = tags.get("cg")
                if cigar:
                    paf_cigars.setdefault(paf_fields[0], []).append(cigar)
            continue
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        seqid, _source, feature, start, end, _score, strand, phase, raw_attrs = fields
        attrs = _parse_gff_attributes(raw_attrs)
        if feature == "mRNA":
            parent_id = attrs.get("ID", "")
            target = attrs.get("Target", "")
            target_parts = target.split()
            protein_id = target_parts[0] if target_parts else attrs.get("ID", "")
            q_start = int(target_parts[1]) if len(target_parts) >= 3 and target_parts[1].isdigit() else 1
            q_end = int(target_parts[2]) if len(target_parts) >= 3 and target_parts[2].isdigit() else 0
            query_coverage = (q_end - q_start + 1) / max(1, len(clean_proteins.get(protein_id, ""))) if q_end else 0.0
            parent_info[parent_id] = {
                "protein_id": protein_id,
                "alignment_protein_id": protein_id,
                "query_start": q_start,
                "query_end": q_end,
                "identity": float(attrs.get("Identity", 0.0) or 0.0),
                "query_coverage": query_coverage,
                "cigar": "NA",
            }
        elif feature == "CDS":
            cds_rows.append((seqid, int(start), int(end), strand, phase, attrs))

    for info in parent_info.values():
        cigars = paf_cigars.get(info["protein_id"], [])
        if len(cigars) == 1:
            info["cigar"] = cigars[0]

    rows = []
    for _seqid, start, end, strand, phase, attrs in cds_rows:
        parent_id = attrs.get("Parent", "")
        info = parent_info.get(parent_id, {})
        cds_target = attrs.get("Target", "").split()
        has_cds_target = len(cds_target) >= 3 and cds_target[1].isdigit() and cds_target[2].isdigit()
        protein_id = cds_target[0] if cds_target else info.get("protein_id", parent_id)
        query_start = int(cds_target[1]) if has_cds_target else "NA"
        query_end = int(cds_target[2]) if has_cds_target else "NA"
        if attrs.get("Identity", "") not in {"", None}:
            identity = float(attrs.get("Identity", 0.0) or 0.0)
            identity_scope = "cds"
        else:
            identity = info.get("identity", 0.0)
            identity_scope = "parent_alignment"
        rows.append({
            "protein_id": protein_id,
            "alignment_protein_id": info.get("alignment_protein_id", protein_id),
            "query_start": query_start,
            "query_end": query_end,
            "parent_query_start": info.get("query_start", "NA"),
            "parent_query_end": info.get("query_end", "NA"),
            "target_start": start,
            "target_end": end,
            "strand": strand,
            "identity": identity,
            "identity_scope": identity_scope,
            "query_coverage": info.get("query_coverage", 0.0),
            "parent_query_coverage": info.get("query_coverage", 0.0),
            "parent_id": parent_id,
            "phase": phase,
            "cigar": info.get("cigar", "NA"),
            "projection_status": "cds_target" if has_cds_target else "missing_cds_target",
        })
    return rows


def splice_motif_score(oriented_intron_seq: str) -> tuple[float, str, str]:
    seq = (oriented_intron_seq or "").upper()
    if len(seq) < 4:
        return 0.0, "NA", "NA"
    donor = seq[:2]
    acceptor = seq[-2:]
    pair = f"{donor}-{acceptor}"
    if pair == "GT-AG":
        return 1.0, donor, acceptor
    if pair in {"GC-AG", "AT-AC"}:
        return 0.8, donor, acceptor
    if donor in {"GT", "GC", "AT"} or acceptor in {"AG", "AC"}:
        return 0.45, donor, acceptor
    return 0.1, donor, acceptor


def phase_compatibility(left_phase: str, right_phase: str, left_cds_length=None) -> str:
    if left_phase in {"", ".", "NA", None} or right_phase in {"", ".", "NA", None} or left_cds_length is None:
        return "unknown"
    try:
        left_phase = int(left_phase)
        right_phase = int(right_phase)
        coding_bases = int(left_cds_length) - left_phase
    except (TypeError, ValueError):
        return "unknown"
    if coding_bases < 0:
        return "incompatible"
    expected_right_phase = (3 - (coding_bases % 3)) % 3
    return "compatible" if right_phase == expected_right_phase else "incompatible"
