"""inference / parameters: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from scipy.stats import chi2
import math
import numpy as np


RATE_MIN = 1e-8


RATE_MAX = 100.0


MULTIPLIER_MIN = 1e-3


MULTIPLIER_MAX = 1e3


PROFILE_DROP_95 = 0.5 * float(chi2.ppf(0.95, 1))


COMPLETE_UNIVERSE_RULES = {"independent_catalogue", "curated_complete_universe"}


GENERATED_ALL_ZERO_MASKS = {"generated_all_zero", "explicit_all_zero"}


def _decode_parameters(theta, model, root_frequency="estimated", root_presence=0.5):
    values = np.exp(np.asarray(theta, dtype=float))
    if model == "ER":
        gain, loss, multiplier, offset = float(values[0]), float(values[0]), 1.0, 1
    elif model == "ARD":
        gain, loss, multiplier, offset = float(values[0]), float(values[1]), 1.0, 2
    elif model == "ARD_FOREGROUND":
        gain, loss, multiplier, offset = float(values[0]), float(values[1]), float(values[2]), 3
    else:
        raise ValueError(f"unknown model: {model}")
    if root_frequency == "stationary":
        rho = gain / max(gain + loss, 1e-300)
    elif root_frequency == "fixed":
        rho = float(root_presence)
    else:
        logit = float(theta[offset])
        rho = 1.0 / (1.0 + math.exp(-logit))
    return gain, loss, multiplier, rho


def _model_bounds(model, root_frequency="estimated"):
    rate_bounds = [(math.log(RATE_MIN), math.log(RATE_MAX))]
    if model in {"ARD", "ARD_FOREGROUND"}:
        rate_bounds.append((math.log(RATE_MIN), math.log(RATE_MAX)))
    if model == "ARD_FOREGROUND":
        rate_bounds.append((math.log(MULTIPLIER_MIN), math.log(MULTIPLIER_MAX)))
    if root_frequency == "estimated":
        rate_bounds.append((-13.8155095579, 13.8155095579))
    return rate_bounds


def _comparison_models(model):
    if model == "foreground":
        return "ARD", "ARD_FOREGROUND", "homogeneous_vs_foreground"
    return "ER", "ARD", "equal_rates_vs_gain_loss"


def _model_starts(model, root_frequency="estimated", root_presence=0.5):
    if model == "ER":
        starts = [[math.log(value)] for value in (0.01, 0.1, 1.0)]
    elif model == "ARD":
        starts = [[math.log(a), math.log(b)] for a, b in ((0.01, 0.1), (0.1, 0.1), (0.1, 1.0), (1.0, 0.1))]
    else:
        starts = [
        [math.log(a), math.log(b), math.log(m)]
        for a, b, m in ((0.01, 0.1, 0.5), (0.1, 0.1, 1.0), (0.1, 1.0, 2.0), (1.0, 0.1, 5.0))
        ]
    if root_frequency == "estimated":
        rho = min(1.0 - 1e-6, max(1e-6, float(root_presence)))
        root_logit = math.log(rho / (1.0 - rho))
        starts = [start + [root_logit] for start in starts]
    return starts
