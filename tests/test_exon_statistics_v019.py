"""Configuration likelihood fitting, observation masking and resampling contracts."""
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from intraphy.topology import SpeciesTree
from intraphy.structure.types import Catalogue, ExonSpan, ExonConfiguration, ObservationEvidence
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.space import enumerate_space
from intraphy.structure.observations import observation_scenarios
from intraphy.structure.serialization import write_catalogues
from intraphy.inference.configuration_model import RateModel
from intraphy.inference.exon_rates import InferenceUnit, load_units, fit_scale, gene_bootstrap
from intraphy.inference.exon_resampling import foreground_bootstrap, simulate_unit


def fixture(count=12, discovery="independent_catalogue"):
    rows = [{"node_id": "root", "parent_id": "", "label": "root", "branch_length": 0}]
    rows += [{"node_id": s, "parent_id": "root", "label": s, "branch_length": .5} for s in "ABCD"]
    tree = SpeciesTree(rows)
    exon = ExonSpan(0, 30)
    full, empty = ExonConfiguration((exon,)), ExonConfiguration(())
    units = []
    for i in range(count):
        obs = tuple(ObservationEvidence(s, (empty if s == "D" and i < 2 else full,)) for s in "ABCD")
        c = Catalogue(f"g{i}", "u", 30, (exon,), (), observations=obs,
                      discovery=discovery, boundary_candidates=(exon,))
        space = enumerate_space(c)
        tips = observation_scenarios(space, tuple("ABCD"), "annotation")[0][1]
        units.append(InferenceUnit(c.family, c.unit, space, tree, tips))
    model = RateModel({kind: .2 if kind in {"exonization", "exon_inactivation"} else 0. for kind in EDIT_KINDS})
    return tuple(units), model


class ExonRateTests(unittest.TestCase):
    def test_pooled_scale_estimates_only_one_parameter(self):
        units, base = fixture()
        fit = fit_scale(units, base)
        self.assertEqual(fit["status"], "estimated_conditional_scale")
        self.assertEqual(fit["estimated_parameter_count"], 1)
        self.assertIsNone(fit["formal_asymptotic_P"])
        self.assertAlmostEqual(fit["scale"], .4368194, places=4)

    def test_one_gene_does_not_fit_all_edit_rates(self):
        units, model = fixture(1)
        self.assertEqual(fit_scale(units, model)["status"], "fewer_than_two_genes")

    def test_gene_multiplicity_is_applied_to_all_local_units(self):
        from intraphy.inference.exon_rates import collection_log_likelihood
        units, model = fixture()
        a = collection_log_likelihood(units, model, {"g0": 2})
        b = collection_log_likelihood(units[:1], model)
        self.assertAlmostEqual(a, 2*b, places=10)
        with self.assertRaises(ValueError):
            collection_log_likelihood(units, model, {"g0": -1})

    def test_bootstrap_reproducible_and_resamples_genes(self):
        units, base = fixture()
        a = gene_bootstrap(units, base, 2, 17)
        b = gene_bootstrap(units, base, 2, 17)
        self.assertEqual(a, b)
        for draw in a["replicates"]:
            self.assertEqual(sum(draw["gene_multiplicities"].values()), 12)

    def test_failed_draws_are_not_discarded(self):
        units, model = fixture()
        with patch("intraphy.inference.exon_rates.fit_scale", return_value={"status": "boundary", "valid_for_resampling": False}):
            result = gene_bootstrap(units, model, 3, 11)
        self.assertEqual(len(result["replicates"]), 3)
        self.assertIsNone(result["percentile_interval"])

    def test_no_asymptotic_shortcut_for_discovered_catalogues(self):
        units, base = fixture(discovery="annotation_discovered")
        fitted = {"log_likelihood": -10., "valid_for_resampling": True, "scale": 1., "status": "valid"}
        with patch("intraphy.inference.exon_resampling.fit_scale", return_value=fitted):
            result = foreground_bootstrap(units, base, {"D"}, 10, 1)
        self.assertEqual(result["status"], "annotation_discovered_catalogue_not_calibrated")
        self.assertIsNone(result["conditional_monte_carlo_P"])
        self.assertEqual(result["parametric_replicates"], [])

    def test_parametric_test_retains_all_simulation_draws(self):
        units, base = fixture()
        fitted = {"log_likelihood": -10., "valid_for_resampling": True, "scale": 1., "status": "valid"}
        with patch("intraphy.inference.exon_resampling.fit_scale", return_value=fitted):
            result = foreground_bootstrap(units, base, {"D"}, 3, 1)
        self.assertEqual(len(result["parametric_replicates"]), 3)
        self.assertEqual(result["conditional_monte_carlo_P"], 1.)
        self.assertFalse(result["discovery_pipeline_calibrated"])

    def test_simulation_does_not_reveal_hidden_origin_states(self):
        units, base = fixture()
        a = simulate_unit(units[0], base, np.random.default_rng(2))
        b = simulate_unit(units[0], base, np.random.default_rng(2))
        for s in a.tips:
            np.testing.assert_array_equal(a.tips[s], b.tips[s])
            self.assertEqual(a.tips[s].sum(), 1)

    def test_foreground_must_not_be_all_edges(self):
        units, base = fixture()
        with self.assertRaises(ValueError):
            fit_scale(units, base, foreground=frozenset("ABCD"))

    def test_duplicate_or_overlapping_units_rejected(self):
        units, _ = fixture()
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"c.jsonl"
            with self.assertRaises(ValueError):
                write_catalogues(p, (units[0].space.catalogue, units[0].space.catalogue))
            first = units[0].space.catalogue
            with self.assertRaises(ValueError):
                write_catalogues(p, (first, replace(first, unit="other")))

    def test_partial_observation_is_not_exact_parametric_validation(self):
        units, base = fixture()
        # A three-state partially resolved observation is checked separately by
        # the calibration wrapper; ordinary unknown 2-state tips are allowed.
        fitted = {"log_likelihood": -10., "valid_for_resampling": True, "scale": 1., "status": "valid"}
        with patch("intraphy.inference.exon_resampling.fit_scale", return_value=fitted):
            result = foreground_bootstrap(units, base, {"D"}, 0, 1)
        self.assertIsNone(result["conditional_monte_carlo_P"])
        self.assertEqual(result["status"], "bootstrap_not_run_or_contains_failed_draws")


if __name__ == "__main__": unittest.main()
