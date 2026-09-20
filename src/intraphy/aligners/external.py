"""aligners / external: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections.abc import Mapping
from intraphy.aligners.columns import _columns_from_cigar
from intraphy.aligners.columns import _empty_stats
from intraphy.aligners.columns import _stats_from_alignment_columns
from intraphy.aligners.columns import revcomp
from intraphy.aligners.formats import _parse_fasta_records
from intraphy.aligners.formats import _parse_paf_tags
from intraphy.aligners.formats import _write_temp_fasta
from intraphy.aligners.formats import _write_temp_fasta_records
from intraphy.aligners.types import AlignmentBackendError
from intraphy.aligners.types import AlignmentStats
from pathlib import Path
import shutil
import subprocess
import tempfile


def _external_minimap2_stats(query: str, target: str, mode: str, threads: int = 1) -> AlignmentStats:
    exe = shutil.which("minimap2")
    if not exe:
        raise AlignmentBackendError("minimap2 was requested but is not available on PATH")
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return _empty_stats(len(query), len(target), "minimap2")
    with tempfile.TemporaryDirectory(prefix="intraphy_minimap2_") as tmp:
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
    hits = []
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
        try:
            mapq = int(fields[11])
        except (TypeError, ValueError):
            mapq = 0
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
        stat.mapping_quality = mapq
        stat.is_secondary = tags.get("tp") == "S"
        rank = (stat.coverage * stat.identity, stat.score)
        hits.append((rank, stat, fields))
    if not hits:
        stats = _empty_stats(len(query), len(target), "minimap2")
        stats.score_scheme = "minimap2_AS_tag"
        stats.enumeration_complete = False
        stats.incomplete_reason = "backend_did_not_guarantee_complete_hit_enumeration"
        return stats
    hits.sort(key=lambda item: item[0], reverse=True)
    candidate_ids = [f"query->target.candidate_{index:03d}" for index in range(1, len(hits) + 1)]
    for index, (_rank, stat, _fields) in enumerate(hits):
        stat.score_scheme = "minimap2_AS_tag"
        stat.candidate_id = candidate_ids[index]
        stat.hit_count = len(hits)
        stat.alternative_candidate_ids = [
            candidate_id for candidate_id in candidate_ids if candidate_id != stat.candidate_id
        ]
        stat.is_secondary = stat.is_secondary or index > 0
        stat.search_interval = {
            "coordinate_system": "0-based-half-open",
            "start0": 0,
            "end0": len(target),
        }
        stat.enumeration_complete = False
        stat.incomplete_reason = "backend_did_not_guarantee_complete_hit_enumeration"
    best = hits[0][1]
    if len(hits) > 1:
        best_rank = hits[0][0][0]
        best.ambiguous_hit_count = sum(1 for rank, _stat, _fields in hits[1:] if rank[0] >= best_rank * 0.95)
        best.alternative_hits = [
            {
                "rank": index,
                "identity": f"{stat.identity:.6g}",
                "coverage": f"{stat.coverage:.6g}",
                "query_start": stat.query_start,
                "query_end": stat.query_end,
                "target_start": stat.target_start,
                "target_end": stat.target_end,
                "strand": stat.strand,
                "mapping_quality": stat.mapping_quality,
                "is_secondary": int(stat.is_secondary),
                "score": f"{stat.score:.6g}",
                "raw_score": f"{stat.raw_score:.6g}" if stat.raw_score is not None else "NA",
                "cigar": stat.cigar,
                "aligned_blocks": list(stat.aligned_blocks),
                "gap_blocks": list(stat.gap_blocks),
                "candidate_id": stat.candidate_id,
                "query_occurrence_id": stat.query_occurrence_id or "NA",
                "target_occurrence_id": stat.target_occurrence_id or "NA",
                "query_transcript_id": stat.query_transcript_id or "NA",
                "target_transcript_id": stat.target_transcript_id or "NA",
                "sequence_kind": stat.sequence_kind,
                "backend": stat.backend,
                "backend_version": stat.backend_version or "NA",
                "score_scheme": stat.score_scheme,
                "nt_identity": f"{stat.nt_identity:.6g}" if stat.nt_identity is not None else "NA",
                "aa_identity": "NA",
                "known_aligned_pairs": stat.known_aligned_pairs,
                "unknown_aligned_pairs": stat.unknown_aligned_pairs,
                "query_covered_bases": stat.query_covered_bases,
                "target_covered_bases": stat.target_covered_bases,
                "query_length": stat.query_length,
                "target_length": stat.target_length,
                "relative_strand": stat.relative_strand,
                "hit_count": stat.hit_count,
                "alternative_candidate_ids": list(stat.alternative_candidate_ids),
                "left_anchor_id": "NA",
                "right_anchor_id": "NA",
                "search_interval": stat.search_interval,
                "enumeration_complete": False,
                "incomplete_reason": "backend_did_not_guarantee_complete_hit_enumeration",
            }
            for index, (_rank, stat, _fields) in enumerate(hits[1:], start=2)
        ]
    return best


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
    with tempfile.TemporaryDirectory(prefix="intraphy_miniprot_") as tmp:
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


def _mafft_pair_alignment(query: str, target: str, threads: int = 1, *, amino: bool = False) -> tuple[str, str]:
    exe = shutil.which("mafft")
    if not exe:
        raise AlignmentBackendError("MAFFT was requested but is not available on PATH")
    with tempfile.TemporaryDirectory(prefix="intraphy_mafft_") as tmp:
        input_path = Path(tmp) / "pair.fa"
        input_path.write_text(f">query\n{query.upper()}\n>target\n{target.upper()}\n")
        proc = subprocess.run(
            [exe, "--quiet", "--thread", str(max(1, int(threads or 1))), "--amino" if amino else "--nuc", "--auto", str(input_path)],
            text=True,
            capture_output=True,
            check=False,
        )
    if proc.returncode != 0:
        raise AlignmentBackendError(proc.stderr.strip() or "MAFFT failed")
    aligned = _parse_fasta_records(proc.stdout, "MAFFT")
    left, right = aligned.get("query", ""), aligned.get("target", "")
    if not left or len(left) != len(right):
        raise AlignmentBackendError("MAFFT did not return a valid two-sequence alignment")
    return left, right


def protein_pair_alignment(query: str, target: str, threads: int = 1) -> tuple[str, str]:
    """Return MAFFT-aligned proteins in query/target order, preserving gaps and X.

    Inputs must be nonempty and ungapped. The caller normalizes internal stops
    to X and removes terminal stops together with their CDS coordinate entries.
    Returned strings contain amino acids; nucleotide scoring is not applied.
    """
    query, target = query.upper(), target.upper()
    if not query or not target or any(symbol in sequence for sequence in (query, target) for symbol in "-*"):
        raise AlignmentBackendError("protein_pair_alignment requires nonempty ungapped proteins with stops normalized by the caller")
    left, right = _mafft_pair_alignment(query, target, threads, amino=True)
    if left.replace("-", "") != query or right.replace("-", "") != target:
        raise AlignmentBackendError("MAFFT protein alignment did not preserve the input residues")
    return left, right


def protein_multiple_alignment(records, mode: str = "linsi", threads: int = 1) -> dict[str, str]:
    """Align a protein family with a deterministic MAFFT L-INS-i/E-INS-i command.

    ``records`` may be a mapping or an iterable of ``(identifier, sequence)``
    pairs. Identifiers are sorted before writing FASTA and the returned mapping
    has the same stable order. Inputs must already have stops normalized by the
    caller, matching :func:`protein_pair_alignment`.
    """
    mode = (mode or "linsi").lower()
    if mode not in {"linsi", "einsi"}:
        raise AlignmentBackendError(f"unsupported protein MSA mode: {mode}")

    items = records.items() if isinstance(records, Mapping) else records
    normalized: dict[str, str] = {}
    for raw_name, raw_sequence in items:
        name = str(raw_name)
        sequence = (raw_sequence or "").upper()
        if not name or any(char.isspace() for char in name) or ">" in name:
            raise AlignmentBackendError("protein MSA identifiers must be nonempty FASTA-safe tokens")
        if name in normalized:
            raise AlignmentBackendError(f"duplicate protein MSA identifier: {name}")
        if not sequence or any(symbol in sequence for symbol in "-*"):
            raise AlignmentBackendError(
                "protein_multiple_alignment requires nonempty ungapped proteins with stops normalized by the caller"
            )
        normalized[name] = sequence
    if not normalized:
        raise AlignmentBackendError("protein_multiple_alignment requires at least one protein")

    ordered = {name: normalized[name] for name in sorted(normalized)}
    if len(ordered) == 1:
        return dict(ordered)

    exe = shutil.which("mafft")
    if not exe:
        raise AlignmentBackendError("MAFFT was requested but is not available on PATH")
    with tempfile.TemporaryDirectory(prefix="intraphy_mafft_msa_") as tmp:
        input_path = Path(tmp) / "proteins.fa"
        _write_temp_fasta_records(input_path, ordered)
        command = [
            exe,
            "--quiet",
            "--thread",
            str(max(1, int(threads or 1))),
            "--threadit",
            "0",
            "--amino",
        ]
        if mode == "linsi":
            command.extend(["--localpair", "--maxiterate", "1000"])
        else:
            command.extend(["--genafpair", "--ep", "0", "--maxiterate", "1000"])
        command.append(str(input_path))
        proc = subprocess.run(command, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise AlignmentBackendError(proc.stderr.strip() or f"MAFFT {mode} failed")

    aligned = _parse_fasta_records(proc.stdout, "MAFFT")
    if set(aligned) != set(ordered):
        raise AlignmentBackendError("MAFFT protein MSA did not return exactly the input identifiers")
    lengths = {len(sequence) for sequence in aligned.values()}
    if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
        raise AlignmentBackendError("MAFFT protein MSA returned unequal or empty aligned sequences")
    for name, sequence in ordered.items():
        if aligned[name].replace("-", "") != sequence:
            raise AlignmentBackendError(f"MAFFT protein MSA did not preserve input residues for {name}")
    return {name: aligned[name] for name in ordered}


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


def _external_lastz_stats(query: str, target: str, mode: str, threads: int = 1) -> AlignmentStats:
    exe = shutil.which("lastz")
    if not exe:
        raise AlignmentBackendError("LASTZ was requested but is not available on PATH")
    query = (query or "").upper()
    target = (target or "").upper()
    if not query or not target:
        return _empty_stats(len(query), len(target), "lastz")
    fields = "score,strand1,size1,zstart1,end1,strand2,size2,zstart2+,end2+,text1,text2,cigarx"
    with tempfile.TemporaryDirectory(prefix="intraphy_lastz_") as tmp:
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
