"""Alignment and splice-boundary helpers for INSIPHY.

The package declares Biopython as a runtime dependency, but these routines keep
small pure-Python fallbacks so source-tree tests can run before installation.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
import re
from dataclasses import dataclass
from pathlib import Path


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
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


class AlignmentBackendError(RuntimeError):
    """Raised when an explicitly requested external alignment backend fails."""


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1].upper()


def ungapped_identity(seq_a: str, seq_b: str) -> float:
    if not seq_a or not seq_b:
        return 0.0
    n = min(len(seq_a), len(seq_b))
    matches = sum(1 for a, b in zip(seq_a[:n].upper(), seq_b[:n].upper()) if a == b and a not in "-N" and b not in "-N")
    return matches / max(1, n)


def _fallback_global(seq_a: str, seq_b: str, match=2, mismatch=-1, gap=-2) -> AlignmentStats:
    seq_a = (seq_a or "").upper()
    seq_b = (seq_b or "").upper()
    m, n = len(seq_a), len(seq_b)
    if not m or not n:
        return AlignmentStats(0.0, 0.0, 0.0, query_end=m, target_end=n)
    if m * n > MAX_INTERNAL_DP_CELLS:
        raise AlignmentBackendError(
            f"internal global alignment requires {m * n} DP cells; "
            "select minimap2 or MAFFT for this sequence pair"
        )

    score = [[0] * (n + 1) for _ in range(m + 1)]
    trace = [[""] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        score[i][0] = i * gap
        trace[i][0] = "U"
    for j in range(1, n + 1):
        score[0][j] = j * gap
        trace[0][j] = "L"
    for i in range(1, m + 1):
        ai = seq_a[i - 1]
        for j in range(1, n + 1):
            bj = seq_b[j - 1]
            diag = score[i - 1][j - 1] + (match if ai == bj and ai != "N" and bj != "N" else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diag, up, left)
            score[i][j] = best
            trace[i][j] = "D" if best == diag else "U" if best == up else "L"

    i, j = m, n
    matches = aligned = q_used = t_used = 0
    ops = []
    while i > 0 or j > 0:
        step = trace[i][j] if i >= 0 and j >= 0 else ""
        if i > 0 and j > 0 and step == "D":
            a, b = seq_a[i - 1], seq_b[j - 1]
            matches += int(a == b and a != "N" and b != "N")
            aligned += 1
            q_used += 1
            t_used += 1
            ops.append("M")
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or step == "U"):
            q_used += 1
            ops.append("D")
            i -= 1
        else:
            t_used += 1
            ops.append("I")
            j -= 1
    ops.reverse()
    identity = matches / max(1, aligned)
    query_coverage = aligned / max(1, m)
    target_coverage = aligned / max(1, n)
    coverage = min(query_coverage, target_coverage)
    return AlignmentStats(
        identity, coverage, float(score[m][n]), query_end=m, target_end=n,
        cigar=_compress_ops(ops), query_coverage=query_coverage,
        target_coverage=target_coverage, aligned_pairs=aligned,
    )


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


def _write_temp_fasta(path: Path, name: str, seq: str):
    path.write_text(f">{name}\n{(seq or '').upper()}\n")


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
        return AlignmentStats(0.0, 0.0, 0.0, query_end=len(query), target_end=len(target), backend="minimap2")
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
        matches = float(fields[9])
        block = max(1.0, float(fields[10]))
        tags = _parse_paf_tags(fields)
        identity = matches / block
        paired = _paired_bases_from_cigar(tags.get("cg", "NA")) or int(min(qend - qstart, abs(tend - tstart)))
        query_coverage = paired / max(1, qlen)
        target_coverage = paired / max(1, len(target))
        coverage = min(query_coverage, target_coverage) if mode == "global" else query_coverage
        score = float(tags.get("AS", matches))
        stat = AlignmentStats(
            identity, coverage, score, qstart + 1, qend, min(tstart, tend) + 1,
            max(tstart, tend), tags.get("cg", "NA"), "minimap2",
            query_coverage, target_coverage, paired,
        )
        rank = (stat.coverage * stat.identity, stat.score)
        if best is None or rank > best[0]:
            best = (rank, stat)
    if best is None:
        return AlignmentStats(0.0, 0.0, 0.0, query_end=len(query), target_end=len(target), backend="minimap2")
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
        return AlignmentStats(0.0, 0.0, 0.0, query_end=len(query), target_end=len(target), backend="miniprot")
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
        tstart = int(fields[7])
        tend = int(fields[8])
        matches = float(fields[9])
        block = max(1.0, float(fields[10]))
        tags = _parse_paf_tags(fields)
        identity = matches / block
        paired = int(qend - qstart)
        coverage = paired / max(1, qlen)
        score = float(tags.get("AS", matches))
        stat = AlignmentStats(
            identity, coverage, score, qstart + 1, qend, min(tstart, tend) + 1,
            max(tstart, tend), tags.get("cg", "NA"), "miniprot",
            coverage, paired / max(1, len(target)), paired,
        )
        rank = (stat.coverage * stat.identity, stat.score)
        if best is None or rank > best[0]:
            best = (rank, stat)
    if best is None:
        return AlignmentStats(0.0, 0.0, 0.0, query_end=len(query), target_end=len(target), backend="miniprot")
    return best[1]


def _external_mafft_stats(query: str, target: str, threads: int = 1) -> AlignmentStats:
    exe = shutil.which("mafft")
    if not exe:
        raise AlignmentBackendError("MAFFT was requested but is not available on PATH")
    with tempfile.TemporaryDirectory(prefix="insiphy_mafft_") as tmp:
        input_path = Path(tmp) / "pair.fa"
        input_path.write_text(f">query\n{query.upper()}\n>target\n{target.upper()}\n")
        proc = subprocess.run(
            [exe, "--quiet", "--thread", str(max(1, int(threads or 1))), "--auto", str(input_path)],
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
    operations = []
    matches = paired = 0
    for a, b in zip(left, right):
        if a != "-" and b != "-":
            paired += 1
            matches += int(a == b and a != "N")
            operations.append("M")
        elif a != "-":
            operations.append("D")
        elif b != "-":
            operations.append("I")
    query_coverage = paired / max(1, len(query))
    target_coverage = paired / max(1, len(target))
    return AlignmentStats(
        matches / max(1, paired), min(query_coverage, target_coverage), float(matches),
        1, len(query), 1, len(target), _compress_ops(operations), "mafft",
        query_coverage, target_coverage, paired,
    )


def available_alignment_backends():
    rows = [{"aligner": "auto", "available": 1, "notes": "sequence-length-aware global alignment selection"}]
    rows.append({"aligner": "internal", "available": 1, "notes": "Biopython PairwiseAligner with pure-Python fallback"})
    rows.append({"aligner": "minimap2", "available": int(shutil.which("minimap2") is not None), "notes": "external nucleotide aligner"})
    rows.append({"aligner": "miniprot", "available": int(shutil.which("miniprot") is not None), "notes": "external protein-to-genome aligner"})
    rows.append({"aligner": "mafft", "available": int(shutil.which("mafft") is not None), "notes": "external progressive multiple-sequence aligner"})
    return rows


def _pairwise_aligner_stats(seq_a: str, seq_b: str, mode: str) -> AlignmentStats:
    from Bio.Align import PairwiseAligner

    seq_a = (seq_a or "").upper()
    seq_b = (seq_b or "").upper()
    if not seq_a or not seq_b:
        return AlignmentStats(0.0, 0.0, 0.0, query_end=len(seq_a), target_end=len(seq_b))
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(seq_a, seq_b)
    if len(alignments) == 0:
        return AlignmentStats(0.0, 0.0, 0.0)
    alignment = alignments[0]
    coordinates = alignment.coordinates
    operations = []
    matches = 0
    aligned = 0
    for index in range(coordinates.shape[1] - 1):
        a_start, a_end = int(coordinates[0, index]), int(coordinates[0, index + 1])
        b_start, b_end = int(coordinates[1, index]), int(coordinates[1, index + 1])
        a_span = a_end - a_start
        b_span = b_end - b_start
        if a_span and b_span:
            span = min(a_span, b_span)
            aligned += span
            matches += sum(
                1
                for left, right in zip(seq_a[a_start:a_end], seq_b[b_start:b_end])
                if left == right and left != "N"
            )
            operations.extend("M" for _ in range(span))
        elif a_span:
            operations.extend("D" for _ in range(a_span))
        elif b_span:
            operations.extend("I" for _ in range(b_span))
    query_start = int(coordinates[0, 0]) + 1
    query_end = int(coordinates[0, -1])
    target_start = int(coordinates[1, 0]) + 1
    target_end = int(coordinates[1, -1])
    query_used = aligned
    target_used = aligned
    if mode == "global":
        coverage = min(query_used / len(seq_a), target_used / len(seq_b))
    else:
        coverage = query_used / len(seq_a)
    return AlignmentStats(
        identity=matches / max(1, aligned),
        coverage=coverage,
        score=float(alignment.score),
        query_start=query_start,
        query_end=query_end,
        target_start=target_start,
        target_end=target_end,
        cigar=_compress_ops(operations),
        query_coverage=query_used / len(seq_a),
        target_coverage=target_used / len(seq_b),
        aligned_pairs=aligned,
    )


def global_alignment_stats(seq_a: str, seq_b: str, backend: str = "internal", threads: int = 1) -> AlignmentStats:
    backend = (backend or "internal").lower()
    if backend == "auto":
        cells = len(seq_a or "") * len(seq_b or "")
        if cells <= MAX_INTERNAL_DP_CELLS:
            return global_alignment_stats(seq_a, seq_b, "internal", threads)
        if cells <= MAX_MAFFT_PAIR_CELLS:
            return global_alignment_stats(seq_a, seq_b, "mafft", threads)
        return global_alignment_stats(seq_a, seq_b, "minimap2", threads)
    if backend == "minimap2":
        return _external_minimap2_stats(seq_a, seq_b, "global", threads)
    if backend == "miniprot":
        return _external_miniprot_stats(seq_a, seq_b, threads)
    if backend == "mafft":
        return _external_mafft_stats(seq_a, seq_b, threads)
    if backend != "internal":
        raise AlignmentBackendError(f"unsupported alignment backend: {backend}")
    try:
        seq_a = (seq_a or "").upper()
        seq_b = (seq_b or "").upper()
        if not seq_a or not seq_b:
            return AlignmentStats(0.0, 0.0, 0.0, query_end=len(seq_a), target_end=len(seq_b))
        if len(seq_a) * len(seq_b) > MAX_INTERNAL_DP_CELLS:
            raise AlignmentBackendError(
                f"internal global alignment requires {len(seq_a) * len(seq_b)} DP cells; "
                "select minimap2 or MAFFT"
            )
        return _pairwise_aligner_stats(seq_a, seq_b, "global")
    except ImportError:
        return _fallback_global(seq_a, seq_b)


def _fallback_local(query: str, target: str, match=2, mismatch=-1, gap=-2) -> AlignmentStats:
    query = (query or "").upper()
    target = (target or "").upper()
    m, n = len(query), len(target)
    if not m or not n:
        return AlignmentStats(0.0, 0.0, 0.0)
    if m * n > MAX_INTERNAL_DP_CELLS:
        raise AlignmentBackendError(
            f"internal local alignment requires {m * n} DP cells; select minimap2"
        )

    score = [[0] * (n + 1) for _ in range(m + 1)]
    trace = [[""] * (n + 1) for _ in range(m + 1)]
    best = (0, 0, 0)
    for i in range(1, m + 1):
        qi = query[i - 1]
        for j in range(1, n + 1):
            tj = target[j - 1]
            diag = score[i - 1][j - 1] + (match if qi == tj and qi != "N" and tj != "N" else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            val = max(0, diag, up, left)
            score[i][j] = val
            trace[i][j] = "" if val == 0 else "D" if val == diag else "U" if val == up else "L"
            if val > best[0]:
                best = (val, i, j)

    val, i, j = best
    end_i, end_j = i, j
    matches = aligned = q_used = 0
    ops = []
    while i > 0 and j > 0 and score[i][j] > 0:
        step = trace[i][j]
        if step == "D":
            a, b = query[i - 1], target[j - 1]
            matches += int(a == b and a != "N" and b != "N")
            aligned += 1
            q_used += 1
            ops.append("M")
            i -= 1
            j -= 1
        elif step == "U":
            q_used += 1
            ops.append("D")
            i -= 1
        else:
            ops.append("I")
            j -= 1
    ops.reverse()
    target_used = aligned
    return AlignmentStats(
        identity=matches / max(1, aligned),
        coverage=q_used / max(1, m),
        score=float(val),
        query_start=i + 1,
        query_end=end_i,
        target_start=j + 1,
        target_end=end_j,
        cigar=_compress_ops(ops),
        query_coverage=q_used / max(1, m),
        target_coverage=target_used / max(1, n),
        aligned_pairs=aligned,
    )


def local_alignment_stats(query: str, target: str, backend: str = "internal", threads: int = 1) -> AlignmentStats:
    backend = (backend or "internal").lower()
    if backend == "minimap2":
        return _external_minimap2_stats(query, target, "local", threads)
    if backend == "miniprot":
        return _external_miniprot_stats(query, target, threads)
    if backend == "mafft":
        raise AlignmentBackendError("MAFFT is a global correspondence backend and cannot perform locus-local searches")
    if backend != "internal":
        raise AlignmentBackendError(f"unsupported alignment backend: {backend}")
    try:
        query = (query or "").upper()
        target = (target or "").upper()
        if not query or not target:
            return AlignmentStats(0.0, 0.0, 0.0)
        if len(query) * len(target) > MAX_INTERNAL_DP_CELLS:
            raise AlignmentBackendError(
                f"internal local alignment requires {len(query) * len(target)} DP cells; select minimap2"
            )
        return _pairwise_aligner_stats(query, target, "local")
    except ImportError:
        return _fallback_local(query, target)


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
