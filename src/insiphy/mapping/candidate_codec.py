"""mapping / candidate_codec: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from insiphy.aligners.types import NT_BLASTN_V1_GAP_EXTEND
from insiphy.aligners.types import NT_BLASTN_V1_GAP_OPEN
from insiphy.aligners.types import NT_BLASTN_V1_MATCH
from insiphy.aligners.types import NT_BLASTN_V1_MISMATCH
from insiphy.coordinates import ClosedInterval1
from insiphy.coordinates import CoordinateBlock
from insiphy.coordinates import Interval0
from insiphy.coordinates import format_legacy_blocks
from insiphy.coordinates import local_interval_to_genome
from insiphy.mapping.fields import STRUCTURAL_ROLES
import re


def _coordinate_block0(block):
    if isinstance(block, CoordinateBlock):
        return block
    if isinstance(block, dict):
        return CoordinateBlock(
            Interval0(int(block["query_start0"]), int(block["query_end0"])),
            Interval0(int(block["target_start0"]), int(block["target_end0"])),
        )
    query_start, query_end, target_start, target_end = block
    return CoordinateBlock(
        ClosedInterval1(int(query_start), int(query_end)).to_interval0(),
        ClosedInterval1(int(target_start), int(target_end)).to_interval0(),
    )


def _alignment_blocks0(aln):
    blocks = getattr(aln, "aligned_blocks", None)
    if blocks:
        return tuple(_coordinate_block0(block) for block in blocks)
    return tuple()


def _alignment_blocks(aln):
    """Compatibility view for callers that still consume public 1-based blocks."""

    public = []
    for block in _alignment_blocks0(aln):
        query = ClosedInterval1.from_interval0(block.query)
        target = ClosedInterval1.from_interval0(block.target)
        public.append((query.start, query.end, target.start, target.end))
    return public


def _format_alignment_blocks(blocks):
    if not blocks:
        return "NA"
    return format_legacy_blocks(tuple(_coordinate_block0(block) for block in blocks))


def _block_signature(block):
    block = _coordinate_block0(block)
    return (
        block.query.start0, block.query.end0,
        block.target.start0, block.target.end0,
    )


def _valid_local_boundary_range(row, sequence):
    try:
        start = int(row["start"])
        end = int(row["end"])
    except (KeyError, TypeError, ValueError):
        return False
    return (
        row.get("contig") not in {None, "", "NA"}
        and row.get("strand") in {"+", "-"}
        and start <= end
        and len(sequence) == end - start + 1
    )


def _explicit_bounded_target(row, sequence):
    explicit = row.get("target_interval_bounded", row.get("search_interval_bounded"))
    if explicit not in {None, "", "NA"}:
        return explicit in {True, 1, "1", "true", "True", "yes"} and _valid_local_boundary_range(row, sequence)
    if str(row.get("boundary_class", "")).lower() in {
        "whole_locus",
        "whole_locus_search_interval",
        "unbounded_locus_search",
    }:
        return False
    return _valid_local_boundary_range(row, sequence) and (
        row.get("role") in STRUCTURAL_ROLES
        or row.get("source_feature_id") not in {None, "", "NA"}
        or row.get("boundary_class") not in {None, "", "NA"}
    )


def _genomic_matched_blocks(row, blocks, side):
    mapped = []
    for block in blocks:
        block = _coordinate_block0(block)
        local = block.query if side == "query" else block.target
        try:
            locus = ClosedInterval1(int(row["start"]), int(row["end"])).to_interval0()
            genome = local_interval_to_genome(local, locus, row.get("strand"))
            public = ClosedInterval1.from_interval0(genome)
        except (KeyError, TypeError, ValueError):
            continue
        mapped.append(
            f"{row.get('contig', 'NA')}:{public.start}-{public.end}:{row.get('strand', 'NA')}"
        )
    return ";".join(mapped) if mapped else "NA"


def _genomic_blocks0(row, blocks, side):
    mapped = []
    try:
        locus = ClosedInterval1(int(row["start"]), int(row["end"])).to_interval0()
    except (KeyError, TypeError, ValueError):
        return tuple()
    for block in blocks:
        block = _coordinate_block0(block)
        local = block.query if side == "query" else block.target
        try:
            mapped.append(local_interval_to_genome(local, locus, row.get("strand")))
        except ValueError:
            continue
    return tuple(mapped)


def _format_genomic_blocks(contig, strand, intervals):
    blocks = []
    for interval in intervals or ():
        try:
            public = ClosedInterval1.from_interval0(interval)
        except ValueError:
            continue
        blocks.append(f"{contig}:{public.start}-{public.end}:{strand}")
    return ";".join(blocks) if blocks else "NA"


def _candidate_value(candidate, key, default="NA"):
    if isinstance(candidate, dict):
        return candidate.get(key, default)
    return getattr(candidate, key, default)


def _alignment_gap_blocks(candidate):
    gaps = []
    for gap in _candidate_value(candidate, "gap_blocks", ()) or ():
        if isinstance(gap, dict):
            if gap.get("gap_in") == "query":
                gaps.append({
                    "gap_in": "query",
                    "query_cut0": int(gap.get("query_cut0", gap.get("query_start0", 0))),
                    "target_start0": int(gap["target_start0"]),
                    "target_end0": int(gap["target_end0"]),
                })
            elif gap.get("gap_in") == "target":
                gaps.append({
                    "gap_in": "target",
                    "query_start0": int(gap["query_start0"]),
                    "query_end0": int(gap["query_end0"]),
                    "target_cut0": int(gap.get("target_cut0", gap.get("target_start0", 0))),
                })
            continue
        query = getattr(gap, "query", None)
        target = getattr(gap, "target", None)
        if query is None or target is None:
            continue
        gaps.append(
            ({
                "gap_in": "query",
                "query_cut0": int(query.start0),
                "target_start0": int(target.start0),
                "target_end0": int(target.end0),
            } if query.length == 0 else {
                "gap_in": "target",
                "query_start0": int(query.start0),
                "query_end0": int(query.end0),
                "target_cut0": int(target.start0),
            })
        )
    return gaps


def _covered_bases(blocks, side):
    intervals = sorted(
        (
            (_coordinate_block0(block).query if side in {"query", 0} else _coordinate_block0(block).target)
        )
        for block in blocks
    )
    if not intervals:
        return 0
    covered = 0
    start, end = intervals[0].start0, intervals[0].end0
    for interval in intervals[1:]:
        if interval.start0 <= end:
            end = max(end, interval.end0)
        else:
            covered += end - start
            start, end = interval.start0, interval.end0
    return covered + end - start


def _unknown_pair_count(candidate):
    explicit = _candidate_value(candidate, "unknown_aligned_pairs", None)
    if explicit not in {None, "", "NA"}:
        return int(explicit)
    return sum(
        int(length)
        for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", str(_candidate_value(candidate, "cigar", "")))
        if operation == "M"
    )


def _candidate_record(candidate, rank, fallback_backend, fallback_scheme):
    if isinstance(candidate, dict):
        record = dict(candidate)
        blocks = tuple(_coordinate_block0(block) for block in record.get("aligned_blocks", ()))
        record["aligned_blocks"] = blocks
        if blocks:
            record["query_start0"] = min(block.query.start0 for block in blocks)
            record["query_end0"] = max(block.query.end0 for block in blocks)
            record["target_start0"] = min(block.target.start0 for block in blocks)
            record["target_end0"] = max(block.target.end0 for block in blocks)
        elif record.get("query_start") not in {None, "", "NA"}:
            query = ClosedInterval1(
                int(record["query_start"]), int(record["query_end"]),
            ).to_interval0()
            target = ClosedInterval1(
                int(record["target_start"]), int(record["target_end"]),
            ).to_interval0()
            record["query_start0"], record["query_end0"] = query.start0, query.end0
            record["target_start0"], record["target_end0"] = target.start0, target.end0
        for field in ("query_start", "query_end", "target_start", "target_end"):
            record.pop(field, None)
        record.setdefault("rank", rank)
        record.setdefault("backend", fallback_backend)
        record.setdefault("query_coverage", record.get("coverage", "NA"))
        record.setdefault("target_coverage", record.get("coverage", "NA"))
        record.setdefault("aligned_pairs", sum(block.query.length for block in blocks))
        record["gap_blocks"] = _alignment_gap_blocks(record)
        if record.get("score_scheme") in {None, "", "unspecified"}:
            record["score_scheme"] = fallback_scheme
        sequence_kind = record.setdefault("sequence_kind", "nucleotide")
        record.setdefault("backend_version", "NA")
        record.setdefault("raw_score", record.get("score", "NA"))
        record.setdefault("nt_identity", record.get("identity", "NA") if sequence_kind == "nucleotide" else "NA")
        record.setdefault("aa_identity", record.get("identity", "NA") if sequence_kind == "amino_acid" else "NA")
        unknown = _unknown_pair_count(record)
        record.setdefault("unknown_aligned_pairs", unknown)
        record.setdefault("known_aligned_pairs", int(record.get("aligned_pairs", 0) or 0))
        record.setdefault("query_covered_bases", _covered_bases(blocks, "query"))
        record.setdefault("target_covered_bases", _covered_bases(blocks, "target"))
        record.setdefault("query_length", "NA")
        record.setdefault("target_length", "NA")
        record.setdefault("relative_strand", record.get("strand", "+"))
        record.setdefault("mapq", record.get("mapping_quality", "NA"))
        if record.get("mapq") is None:
            record["mapq"] = "NA"
        record.setdefault("search_interval_side", "target")
        return record
    adapter_fields = dict(vars(candidate)) if hasattr(candidate, "__dict__") else {}
    for field in (
        "aligned_blocks", "gap_blocks", "alternative_hits",
        "query_interval", "target_interval",
    ):
        adapter_fields.pop(field, None)
    blocks = _alignment_blocks0(candidate)
    aligned_pairs = int(_candidate_value(candidate, "aligned_pairs", 0) or 0)
    unknown_pairs = _unknown_pair_count(candidate)
    sequence_kind = _candidate_value(candidate, "sequence_kind", "nucleotide")
    mapq = _candidate_value(
        candidate, "mapq", _candidate_value(candidate, "mapping_quality", "NA"),
    )
    if mapq is None:
        mapq = "NA"
    return {
        **adapter_fields,
        "candidate_id": _candidate_value(candidate, "candidate_id", "NA"),
        "rank": rank,
        "identity": _candidate_value(candidate, "identity", "NA"),
        "coverage": _candidate_value(candidate, "coverage", "NA"),
        "query_coverage": _candidate_value(candidate, "query_coverage", "NA"),
        "target_coverage": _candidate_value(candidate, "target_coverage", "NA"),
        "aligned_pairs": aligned_pairs,
        "query_start0": min((block.query.start0 for block in blocks), default="NA"),
        "query_end0": max((block.query.end0 for block in blocks), default="NA"),
        "target_start0": min((block.target.start0 for block in blocks), default="NA"),
        "target_end0": max((block.target.end0 for block in blocks), default="NA"),
        "strand": _candidate_value(candidate, "strand", "+"),
        "mapping_quality": _candidate_value(candidate, "mapping_quality", "NA"),
        "is_secondary": int(bool(_candidate_value(candidate, "is_secondary", rank > 1))),
        "score": _candidate_value(candidate, "score", "NA"),
        "cigar": _candidate_value(candidate, "cigar", "NA"),
        "aligned_blocks": blocks,
        "gap_blocks": _alignment_gap_blocks(candidate),
        "sequence_kind": sequence_kind,
        "backend": _candidate_value(candidate, "backend", fallback_backend),
        "backend_version": _candidate_value(candidate, "backend_version", "NA"),
        "score_scheme": _candidate_value(candidate, "score_scheme", fallback_scheme),
        "raw_score": _candidate_value(candidate, "raw_score", _candidate_value(candidate, "score", "NA")),
        "nt_identity": _candidate_value(
            candidate,
            "nt_identity",
            _candidate_value(candidate, "identity", "NA") if sequence_kind == "nucleotide" else "NA",
        ),
        "aa_identity": _candidate_value(
            candidate,
            "aa_identity",
            _candidate_value(candidate, "identity", "NA") if sequence_kind == "amino_acid" else "NA",
        ),
        "known_aligned_pairs": _candidate_value(
            candidate, "known_aligned_pairs", aligned_pairs,
        ),
        "unknown_aligned_pairs": unknown_pairs,
        "query_covered_bases": _candidate_value(
            candidate, "query_covered_bases", _covered_bases(blocks, "query"),
        ),
        "target_covered_bases": _candidate_value(
            candidate, "target_covered_bases", _covered_bases(blocks, "target"),
        ),
        "query_length": _candidate_value(candidate, "query_length", "NA"),
        "target_length": _candidate_value(candidate, "target_length", "NA"),
        "relative_strand": _candidate_value(
            candidate, "relative_strand", _candidate_value(candidate, "strand", "+"),
        ),
        "mapq": mapq,
        "left_anchor_id": _candidate_value(candidate, "left_anchor_id", "NA"),
        "right_anchor_id": _candidate_value(candidate, "right_anchor_id", "NA"),
        "search_interval": _candidate_value(candidate, "search_interval", "NA"),
        "search_interval_side": "target",
        "enumeration_complete": _candidate_value(
            candidate, "enumeration_complete", "NA",
        ),
        "incomplete_reason": _candidate_value(candidate, "incomplete_reason", "NA"),
    }


def _public_interval(interval):
    """The sole adapter for nonempty half-open intervals written to TSV JSON."""

    if not isinstance(interval, dict) or interval.get("start0") in {None, "", "NA"}:
        return interval
    internal = Interval0(int(interval["start0"]), int(interval["end0"]))
    public = ClosedInterval1.from_interval0(internal)
    return {
        key: value
        for key, value in interval.items()
        if key not in {"coordinate_system", "start0", "end0"}
    } | {
        "coordinate_system": "1-based-closed",
        "start": public.start,
        "end": public.end,
    }


def _public_gap_blocks(gaps):
    public = []
    for gap in gaps or ():
        if gap.get("gap_in") == "query":
            interval = ClosedInterval1.from_interval0(Interval0(
                int(gap["target_start0"]), int(gap["target_end0"]),
            ))
            public.append({
                "gap_in": "query",
                "query_cut0": int(gap["query_cut0"]),
                "target_start": interval.start,
                "target_end": interval.end,
            })
        elif gap.get("gap_in") == "target":
            interval = ClosedInterval1.from_interval0(Interval0(
                int(gap["query_start0"]), int(gap["query_end0"]),
            ))
            public.append({
                "gap_in": "target",
                "query_start": interval.start,
                "query_end": interval.end,
                "target_cut0": int(gap["target_cut0"]),
            })
    return public


def _public_candidate_record(record):
    public = dict(record)
    blocks = tuple(_coordinate_block0(block) for block in record.get("aligned_blocks", ()))
    public["aligned_blocks"] = []
    for block in blocks:
        query = ClosedInterval1.from_interval0(block.query)
        target = ClosedInterval1.from_interval0(block.target)
        public["aligned_blocks"].append(
            (query.start, query.end, target.start, target.end)
        )
    for side in ("query", "target"):
        start0 = record.get(f"{side}_start0")
        end0 = record.get(f"{side}_end0")
        if start0 not in {None, "", "NA"} and int(end0) > int(start0):
            interval = ClosedInterval1.from_interval0(Interval0(int(start0), int(end0)))
            public[f"{side}_start"] = interval.start
            public[f"{side}_end"] = interval.end
        public.pop(f"{side}_start0", None)
        public.pop(f"{side}_end0", None)
    public["gap_blocks"] = _public_gap_blocks(record.get("gap_blocks", ()))
    public["search_interval"] = _public_interval(record.get("search_interval"))
    for side in ("query", "target"):
        genomic = record.get(f"{side}_genomic_blocks0")
        if genomic:
            public[f"{side}_genomic_blocks"] = [
                _public_interval({"start0": interval.start0, "end0": interval.end0})
                for interval in genomic
            ]
        public.pop(f"{side}_genomic_blocks0", None)
    return public


def _nt_column_score(aln):
    score = (
        float(getattr(aln, "matches", 0) or 0) * NT_BLASTN_V1_MATCH
        + float(getattr(aln, "mismatches", 0) or 0) * NT_BLASTN_V1_MISMATCH
    )
    for length_text, operation in re.findall(r"(\d+)([ID])", str(getattr(aln, "cigar", ""))):
        length = int(length_text)
        score += NT_BLASTN_V1_GAP_OPEN
        score += max(0, length - 1) * NT_BLASTN_V1_GAP_EXTEND
    return score


def _set_explicit_alignment_score(aln):
    if getattr(aln, "backend", "") == "mafft_overlap":
        aln.score = _nt_column_score(aln)
        aln.score_scheme = "nt_blastn_v1/sum_of_column_scores"


def _alignment_candidate_records(aln, candidate_set=None):
    backend = getattr(aln, "backend", "internal")
    score_scheme = getattr(aln, "score_scheme", "unspecified")
    if candidate_set is not None:
        return [
            _candidate_record(candidate, rank, backend, candidate_set.score_scheme)
            for rank, candidate in enumerate(candidate_set.candidates, start=1)
        ]
    records = [_candidate_record(aln, 1, backend, score_scheme)] if _alignment_blocks(aln) else []
    records.extend(
        _candidate_record(candidate, rank, backend, score_scheme)
        for rank, candidate in enumerate(getattr(aln, "alternative_hits", ()) or (), start=2)
    )
    return records


def _transpose_cigar(cigar):
    return str(cigar).translate(str.maketrans({"I": "D", "D": "I"}))


def _transpose_gap_blocks(gaps):
    transposed = []
    for gap in gaps or ():
        if gap.get("gap_in") == "query":
            transposed.append({
                "gap_in": "target",
                "query_start0": int(gap["target_start0"]),
                "query_end0": int(gap["target_end0"]),
                "target_cut0": int(gap["query_cut0"]),
            })
        elif gap.get("gap_in") == "target":
            transposed.append({
                "gap_in": "query",
                "query_cut0": int(gap["target_cut0"]),
                "target_start0": int(gap["query_start0"]),
                "target_end0": int(gap["query_end0"]),
            })
    return transposed


def _transpose_candidate_record(record):
    transposed = dict(record)
    transposed["query_start0"], transposed["target_start0"] = (
        record.get("target_start0", "NA"), record.get("query_start0", "NA"),
    )
    transposed["query_end0"], transposed["target_end0"] = (
        record.get("target_end0", "NA"), record.get("query_end0", "NA"),
    )
    transposed["query_coverage"], transposed["target_coverage"] = (
        record.get("target_coverage", "NA"), record.get("query_coverage", "NA"),
    )
    transposed["query_covered_bases"], transposed["target_covered_bases"] = (
        record.get("target_covered_bases", "NA"), record.get("query_covered_bases", "NA"),
    )
    transposed["query_length"], transposed["target_length"] = (
        record.get("target_length", "NA"), record.get("query_length", "NA"),
    )
    transposed["query_occurrence_id"], transposed["target_occurrence_id"] = (
        record.get("target_occurrence_id", "NA"),
        record.get("query_occurrence_id", "NA"),
    )
    transposed["query_transcript_id"], transposed["target_transcript_id"] = (
        record.get("target_transcript_id", "NA"),
        record.get("query_transcript_id", "NA"),
    )
    transposed["aligned_blocks"] = tuple(
        CoordinateBlock(
            _coordinate_block0(block).target,
            _coordinate_block0(block).query,
        )
        for block in record.get("aligned_blocks", ())
    )
    transposed["gap_blocks"] = _transpose_gap_blocks(record.get("gap_blocks", ()))
    transposed["cigar"] = _transpose_cigar(record.get("cigar", "NA"))
    transposed["short_sequence_coverage"] = record.get("query_coverage", "NA")
    transposed["alignment_input_transposed"] = 1
    transposed["search_interval_side"] = "query"
    return transposed


def _mapped_genomic_interval(row, rel_start, rel_end):
    if rel_start in {"NA", None, ""} or rel_end in {"NA", None, ""}:
        return "NA", "NA", "NA"
    rel_start = int(rel_start)
    rel_end = int(rel_end)
    if rel_end < rel_start:
        return "NA", "NA", "NA"
    start = int(row["start"])
    end = int(row["end"])
    if row.get("strand") == "-":
        genomic_start = end - rel_end + 1
        genomic_end = end - rel_start + 1
    else:
        genomic_start = start + rel_start - 1
        genomic_end = start + rel_end - 1
    genomic_start, genomic_end = min(genomic_start, genomic_end), max(genomic_start, genomic_end)
    return genomic_start, genomic_end, genomic_end - genomic_start + 1


def _mapped_genomic_interval0(row, local):
    if local is None:
        return "NA", "NA", "NA"
    try:
        locus = ClosedInterval1(int(row["start"]), int(row["end"])).to_interval0()
        genome = local_interval_to_genome(local, locus, row.get("strand"))
        public = ClosedInterval1.from_interval0(genome)
    except (KeyError, TypeError, ValueError):
        return "NA", "NA", "NA"
    return public.start, public.end, public.length
