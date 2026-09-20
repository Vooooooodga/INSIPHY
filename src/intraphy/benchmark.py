"""Compatibility exports. Implementations live in the documented submodules."""


# Public entry points; implementations have a single owner.
from intraphy.verification.benchmark import (
    NON_BIOLOGICAL_BENCHMARK_CLASSES,
    COPY_CONTEXT_CLASSES,
    AMBIGUOUS_CLASSES,
    event_key,
    call_scope,
    benchmark_events,
)
