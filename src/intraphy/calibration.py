"""Compatibility exports. Implementations live in the documented submodules."""


# Public entry points; implementations have a single owner.
from intraphy.verification.calibration import (
    DEFAULT_SCENARIOS,
    _run_case,
    _mean,
    _fmt,
    calibrate_simulations,
)
