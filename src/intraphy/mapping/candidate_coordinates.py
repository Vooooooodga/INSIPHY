"""mapping / candidate coordinates: explicit implementation ownership."""
from __future__ import annotations

from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import CoordinateBlock
from intraphy.coordinates import Interval0
from intraphy.coordinates import format_legacy_blocks
from intraphy.coordinates import local_interval_to_genome
from intraphy.mapping.fields import STRUCTURAL_ROLES
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
