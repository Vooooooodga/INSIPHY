"""storage / values: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations




UNKNOWN = {
    "",
    ".",
    "NA",
    "?",
    "unknown",
    "ambiguous",
    "truncated",
    "unavailable",
    "unresolved",
}


def norm_state(value):
    value = (value or "unknown").strip()
    return "unknown" if value in UNKNOWN else value


def to_float(value, default=0.0):
    try:
        if value in ("", "NA", None):
            return default
        return float(value)
    except ValueError:
        return default


def uniq(values):
    return len({v for v in values if v not in ("", "NA", None)})
