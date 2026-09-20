"""inference / formatting: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
import math


def _fmt(value):
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):.8g}"


def _bh_adjust(rows):
    by_test = defaultdict(list)
    for index, row in enumerate(rows):
        try:
            by_test[row["test_id"]].append((index, float(row["p_value"])))
        except (TypeError, ValueError):
            row["q_value"] = "NA"
            row["q_value_method"] = "not_available"
    for entries in by_test.values():
        if len(entries) == 1:
            rows[entries[0][0]]["q_value"] = _fmt(entries[0][1])
            rows[entries[0][0]]["q_value_method"] = "Benjamini-Hochberg_m_equals_1"
            continue
        ordered = sorted(entries, key=lambda item: item[1])
        adjusted = [0.0] * len(ordered)
        running = 1.0
        for rank in range(len(ordered), 0, -1):
            value = min(running, ordered[rank - 1][1] * len(ordered) / rank)
            adjusted[rank - 1] = value
            running = value
        for (index, _p), value in zip(ordered, adjusted):
            rows[index]["q_value"] = _fmt(value)
            rows[index]["q_value_method"] = "Benjamini-Hochberg_within_test_type"
