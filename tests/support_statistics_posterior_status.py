import math

import tempfile

import unittest

from pathlib import Path

from types import SimpleNamespace

from unittest.mock import patch

import numpy as np

from scipy.optimize import brentq as scipy_brentq

from intraphy.io import (
    read_structural_site_matrix,
    read_tsv,
    structural_site_observed_state,
    write_structural_site_matrix,
)

from intraphy.parsimony import (
    _event_type,
    _parsimony_tables,
    _placement_status,
    _structural_relation,
    infer_single_copy_parsimony,
)

from intraphy.structural_phylogeny import (
    PROFILE_DROP_95,
    _ascertainment_log_probability,
    _dataset_log_likelihood,
    _expected_transition_count,
    _fit_valid_for_posterior,
    _pattern_log_likelihood,
    _posterior_messages,
    _profile_contrast_interval,
    _profile_interval,
    _selected_for_ascertainment,
    _transition_matrix,
    _transition_event_type,
    _transition_structural_relation,
    _validated_tree_rows,
    fit_model,
    infer_single_copy_phylogeny,
)

from intraphy.tree import SpeciesTree

from intraphy.structural_sites import structural_matrix_annotation_view


class PosteriorStatusTestsSupport:
    def _fit(self, model, status, identifiable):
        intervals = [{"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None}]
        if model in {"ARD", "ARD_FOREGROUND"}:
            intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        if model == "ARD_FOREGROUND":
            intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        intervals.append({"low": None, "high": None, "status": status, "raw_low": None, "raw_high": None})
        theta_length = 2 if model == "ER" else 3 if model == "ARD" else 4
        if model == "ARD":
            contrast = {
                "name": "gain_loss_rate_ratio", "low": 0.5, "high": 2.0,
                "status": "two_sided", "raw_low": math.log(0.5), "raw_high": math.log(2.0),
                "theta_low": np.zeros(theta_length), "theta_high": np.zeros(theta_length),
            }
        elif model == "ARD_FOREGROUND":
            contrast = {
                "name": "foreground_multiplier", "low": 0.5, "high": 2.0,
                "status": "two_sided", "raw_low": math.log(0.5), "raw_high": math.log(2.0),
                "theta_low": np.zeros(theta_length), "theta_high": np.zeros(theta_length),
            }
        else:
            contrast = {"name": "NA", "status": "not_applicable"}
        return {
            "model": model,
            "theta": np.zeros(theta_length),
            "gain_rate": 0.1,
            "loss_rate": 0.1,
            "foreground_multiplier": 1.0,
            "root_presence": 0.5,
            "log_likelihood": -2.0 if model == "ER" else -1.9,
            "parameter_count": theta_length,
            "aic": 8.0 if model == "ER" else 9.8,
            "converged": True,
            "fit_status": status,
            "optimizer_message": status,
            "boundary": False,
            "intervals": intervals,
            "start_count": 1,
            "identifiable": identifiable,
            "information_condition": math.inf,
            "root_frequency": "estimated",
            "profile_support_thetas": [],
            "contrast_profile": contrast,
        }
