"""coding / projection: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.coding.types import _ANCHOR_WINDOW
from intraphy.coding.types import _BLOSUM62
from intraphy.coding.types import _MIN_ANCHOR_PAIRS
from intraphy.coordinates import CoordinateBlock
from intraphy.coordinates import Interval0


def _blosum62_score(query_aa, target_aa):
    return float(_BLOSUM62[query_aa, target_aa])


def _aligned_occurrence_pairs(query, target, aligned_query, aligned_target):
    """Index aligned codons in 0-based MSA, occurrence and genome coordinates."""
    pairs = {}
    known_columns = []
    query_index = target_index = 0
    for column0, (query_aa, target_aa) in enumerate(zip(aligned_query, aligned_target)):
        query_residue = query_index if query_aa != "-" else None
        target_residue = target_index if target_aa != "-" else None
        query_codon = query.codon_sources[query_residue] if query_residue is not None else None
        target_codon = target.codon_sources[target_residue] if target_residue is not None else None
        query_known = query_residue is not None and query.codons[query_residue] is not None
        target_known = target_residue is not None and target.codons[target_residue] is not None
        query_index += query_aa != "-"
        target_index += target_aa != "-"
        if not query_known or not target_known:
            continue
        known_columns.append((column0, query_aa, target_aa, query_codon, target_codon))
        touched = set()
        for query_base, target_base in zip(query_codon, target_codon):
            key = (query_base.occurrence_id, target_base.occurrence_id)
            record = pairs.setdefault(key, {
                "positions0": [], "aa_pairs": 0, "aa_matches": 0,
                "blosum62_score": 0.0, "msa_columns0": [],
            })
            record["positions0"].append((
                query_base.occurrence_interval.start0,
                target_base.occurrence_interval.start0,
            ))
            touched.add(key)
        for key in touched:
            pairs[key]["aa_pairs"] += 1
            pairs[key]["aa_matches"] += query_aa == target_aa
            pairs[key]["blosum62_score"] += _blosum62_score(query_aa, target_aa)
            pairs[key]["msa_columns0"].append(column0)
    return pairs, tuple(known_columns)


def _coordinate_blocks0(positions0):
    blocks = []
    for query0, target0 in sorted(positions0):
        if (
            blocks
            and query0 == blocks[-1].query.end0
            and target0 == blocks[-1].target.end0
        ):
            blocks[-1] = CoordinateBlock(
                Interval0(blocks[-1].query.start0, query0 + 1),
                Interval0(blocks[-1].target.start0, target0 + 1),
            )
        else:
            blocks.append(CoordinateBlock(
                Interval0(query0, query0 + 1), Interval0(target0, target0 + 1),
            ))
    return tuple(blocks)


def _anchor_metrics(
    record, known_columns, aligned_query, aligned_target,
    query_occurrence, target_occurrence,
):
    columns = sorted(set(record["msa_columns0"]))
    first, last = columns[0], columns[-1]

    def independent_of_focal_occurrences(item):
        query_codon, target_codon = item[3], item[4]
        return (
            all(base.occurrence_id != query_occurrence for base in query_codon)
            and all(base.occurrence_id != target_occurrence for base in target_codon)
        )

    independent_columns = [
        item for item in known_columns if independent_of_focal_occurrences(item)
    ]
    left = [item for item in independent_columns if item[0] < first][-_ANCHOR_WINDOW:]
    right = [item for item in independent_columns if item[0] > last][:_ANCHOR_WINDOW]
    left_score = sum(_blosum62_score(item[1], item[2]) for item in left)
    right_score = sum(_blosum62_score(item[1], item[2]) for item in right)
    left_supported = len(left) >= _MIN_ANCHOR_PAIRS and left_score > 0
    right_supported = len(right) >= _MIN_ANCHOR_PAIRS and right_score > 0
    left_terminal = first == 0
    right_terminal = last + 1 == len(aligned_query) == len(aligned_target)
    span = max(1, last - first + 1)
    gap_columns = sum(
        aligned_query[column0] == "-" or aligned_target[column0] == "-"
        for column0 in range(first, last + 1)
    )
    resolved = (
        (left_supported and right_supported)
        or (left_terminal and right_supported)
        or (right_terminal and left_supported)
    )
    terminal_side = (
        "both" if left_terminal and right_terminal
        else "left" if left_terminal
        else "right" if right_terminal
        else "none"
    )
    return {
        "left_pairs": len(left), "right_pairs": len(right),
        "left_score": left_score, "right_score": right_score,
        "left_supported": left_supported, "right_supported": right_supported,
        "terminal_side": terminal_side, "gap_fraction": gap_columns / span,
        "resolved": resolved,
        "column_interval": Interval0(first, last + 1),
    }


def _candidate_json(candidate):
    """Return one stable, scalar-only record for table serialization."""
    interval = candidate.msa_column_interval
    return {
        "aa_identity": candidate.aa_identity,
        "anchor_resolved": candidate.anchor_resolved,
        "blosum62_score": candidate.blosum62_score,
        "candidate_id": candidate.candidate_id,
        "competing_occurrences": list(candidate.competing_occurrences),
        "known_aa_pairs": candidate.known_aa_pairs,
        "left_anchor_pairs": candidate.left_anchor_pairs,
        "left_anchor_score": candidate.left_anchor_score,
        "left_anchor_supported": candidate.left_anchor_supported,
        "msa_column_end0": interval.end0,
        "msa_column_start0": interval.start0,
        "position_eligible": candidate.position_eligible,
        "position_monotonic": candidate.position_monotonic,
        "projected_blocks": [
            {
                "query_end0": block.query.end0,
                "query_start0": block.query.start0,
                "target_end0": block.target.end0,
                "target_start0": block.target.start0,
            }
            for block in candidate.coordinate_blocks
        ],
        "query_cds_coverage": candidate.query_cds_coverage,
        "query_transcript_id": candidate.query_transcript_id,
        "right_anchor_pairs": candidate.right_anchor_pairs,
        "right_anchor_score": candidate.right_anchor_score,
        "right_anchor_supported": candidate.right_anchor_supported,
        "score_scheme": "BLOSUM62",
        "target_cds_coverage": candidate.target_cds_coverage,
        "target_transcript_id": candidate.target_transcript_id,
        "terminal_side": candidate.terminal_side,
    }
