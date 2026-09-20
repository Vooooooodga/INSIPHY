"""mapping / candidate serialization: explicit implementation ownership."""
from __future__ import annotations

from intraphy.aligners.types import NT_BLASTN_V1_GAP_EXTEND
from intraphy.aligners.types import NT_BLASTN_V1_GAP_OPEN
from intraphy.aligners.types import NT_BLASTN_V1_MATCH
from intraphy.aligners.types import NT_BLASTN_V1_MISMATCH
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import CoordinateBlock
from intraphy.coordinates import Interval0
from intraphy.coordinates import local_interval_to_genome
from intraphy.mapping.candidate_coordinates import _alignment_blocks
from intraphy.mapping.candidate_coordinates import _coordinate_block0
from intraphy.mapping.candidate_records import _candidate_record
import re


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
