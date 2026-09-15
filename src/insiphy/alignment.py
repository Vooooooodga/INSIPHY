"""Alignment and splice-boundary helpers for INSIPHY.

The package declares Biopython as a runtime dependency, but these routines keep
small pure-Python fallbacks so source-tree tests can run before installation.
"""

from __future__ import annotations

from dataclasses import dataclass


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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
    if m * n > 2_000_000:
        identity = ungapped_identity(seq_a, seq_b)
        coverage = min(m, n) / max(1, max(m, n))
        return AlignmentStats(identity, coverage, identity * coverage, query_end=m, target_end=n, cigar=f"{min(m, n)}M")

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
    coverage = min(q_used / max(1, m), t_used / max(1, n))
    return AlignmentStats(identity, coverage, float(score[m][n]), query_end=m, target_end=n, cigar=_compress_ops(ops))


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


def global_alignment_stats(seq_a: str, seq_b: str) -> AlignmentStats:
    try:
        from Bio import pairwise2  # type: ignore

        seq_a = (seq_a or "").upper()
        seq_b = (seq_b or "").upper()
        if not seq_a or not seq_b:
            return AlignmentStats(0.0, 0.0, 0.0, query_end=len(seq_a), target_end=len(seq_b))
        if len(seq_a) * len(seq_b) > 2_000_000:
            return _fallback_global(seq_a, seq_b)
        aln = pairwise2.align.globalms(seq_a, seq_b, 2, -1, -2, -0.5, one_alignment_only=True)[0]
        a, b, score, _, _ = aln
        aligned = matches = q_used = t_used = 0
        ops = []
        for ca, cb in zip(a, b):
            if ca != "-":
                q_used += 1
            if cb != "-":
                t_used += 1
            if ca != "-" and cb != "-":
                aligned += 1
                matches += int(ca == cb and ca != "N" and cb != "N")
                ops.append("M")
            elif ca != "-":
                ops.append("D")
            else:
                ops.append("I")
        return AlignmentStats(matches / max(1, aligned), min(q_used / max(1, len(seq_a)), t_used / max(1, len(seq_b))), float(score), query_end=len(seq_a), target_end=len(seq_b), cigar=_compress_ops(ops))
    except Exception:
        return _fallback_global(seq_a, seq_b)


def _fallback_local(query: str, target: str, match=2, mismatch=-1, gap=-2) -> AlignmentStats:
    query = (query or "").upper()
    target = (target or "").upper()
    m, n = len(query), len(target)
    if not m or not n:
        return AlignmentStats(0.0, 0.0, 0.0)
    if m * n > 2_000_000:
        return best_ungapped_hit(query, target)

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
    return AlignmentStats(
        identity=matches / max(1, aligned),
        coverage=q_used / max(1, m),
        score=float(val),
        query_start=i + 1,
        query_end=end_i,
        target_start=j + 1,
        target_end=end_j,
        cigar=_compress_ops(ops),
    )


def local_alignment_stats(query: str, target: str) -> AlignmentStats:
    try:
        from Bio import pairwise2  # type: ignore

        query = (query or "").upper()
        target = (target or "").upper()
        if not query or not target:
            return AlignmentStats(0.0, 0.0, 0.0)
        if len(query) * len(target) > 2_000_000:
            return best_ungapped_hit(query, target)
        aln = pairwise2.align.localms(query, target, 2, -1, -2, -0.5, one_alignment_only=True)[0]
        a, b, score, start, end = aln
        matches = aligned = q_used = 0
        t_used = 0
        target_before = 0
        target_after = 0
        in_aln = False
        ops = []
        for idx, (ca, cb) in enumerate(zip(a, b)):
            if idx == start:
                in_aln = True
            if idx == end:
                in_aln = False
            if cb != "-" and not in_aln and idx < start:
                target_before += 1
            if cb != "-" and not in_aln and idx >= end:
                target_after += 1
            if not in_aln:
                continue
            if ca != "-":
                q_used += 1
            if cb != "-":
                t_used += 1
            if ca != "-" and cb != "-":
                aligned += 1
                matches += int(ca == cb and ca != "N" and cb != "N")
                ops.append("M")
            elif ca != "-":
                ops.append("D")
            else:
                ops.append("I")
        return AlignmentStats(matches / max(1, aligned), q_used / max(1, len(query)), float(score), 1, q_used, target_before + 1, len(target) - target_after, _compress_ops(ops))
    except Exception:
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


def phase_compatibility(left_phase: str, right_phase: str) -> str:
    if left_phase in {"", ".", "NA", None} or right_phase in {"", ".", "NA", None}:
        return "unknown"
    return "compatible" if str(left_phase) == str(right_phase) else "incompatible"
