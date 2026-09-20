"""mapping / coverage: explicit implementation ownership."""
from __future__ import annotations

from intraphy.coordinates import Interval0
from intraphy.mapping.chain_types import CoverageClassification
from typing import Iterable
from typing import Optional


def classify_reference_coverage(
    intervals: Iterable[Interval0],
    parent_interval: Optional[Interval0] = None,
):
    """Describe whether several projections are complementary or repetitive."""

    intervals = tuple(sorted(intervals))
    if not intervals:
        uncovered = parent_interval.length if parent_interval is not None else None
        return CoverageClassification("uncovered", 0, 0, uncovered)
    if parent_interval is not None and any(
        interval.start0 < parent_interval.start0 or interval.end0 > parent_interval.end0
        for interval in intervals
    ):
        raise ValueError("reference projection lies outside the parent interval")

    events = []
    for interval in intervals:
        events.append((interval.start0, 1))
        events.append((interval.end0, -1))
    covered_bases = 0
    overlap_bases = 0
    depth = 0
    previous = events[0][0]
    for coordinate, delta in sorted(events, key=lambda item: (item[0], item[1])):
        span = coordinate - previous
        if depth > 0:
            covered_bases += span
        if depth > 1:
            overlap_bases += span
        depth += delta
        previous = coordinate
    uncovered = None if parent_interval is None else parent_interval.length - covered_bases

    if len(intervals) == 1:
        if parent_interval is None:
            relation = "single_projection"
        else:
            relation = "single_complete" if uncovered == 0 else "single_partial"
    elif overlap_bases:
        relation = "repeated_overlap"
    else:
        if parent_interval is None:
            relation = "complementary_disjoint"
        else:
            relation = "complementary_complete" if uncovered == 0 else "complementary_partial"
    return CoverageClassification(relation, covered_bases, overlap_bases, uncovered)
