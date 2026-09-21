"""Regression tests for the independently identified 0.19.0 audit failures."""
from dataclasses import replace
from itertools import combinations
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from intraphy.structure.types import (Catalogue, ExonSpan as E, ExonConfiguration as C,
    Material, ObservationEvidence as O, ConfigurationAlternative as A, InsertionPayload)
from intraphy.structure.space import enumerate_space, estimate_geometry_count, close_span_catalogue
from intraphy.structure.paths import make_graph
from intraphy.structure.observations import compatibility
from intraphy.structure.serialization import write_catalogues, read_catalogues, encode_catalogue
from intraphy.structure.validation import validate_collection
from intraphy.structure.tree_context import canonical_tree
from intraphy.inference.configuration_model import RateModel, evaluate_model, generator
from intraphy.inference.exon_rates import fit_scale, load_units, gene_bootstrap
from intraphy.inference.exon_resampling import nested_statistic
from intraphy.structure.edits import EDIT_KINDS
from intraphy.topology import SpeciesTree
from test_exon_statistics_v019 import fixture


def trees():
    base = [{"node_id": "r", "parent_id": "", "label": "r"},
        {"node_id": "A", "parent_id": "r", "label": "A", "branch_length": 1.},
        {"node_id": "B", "parent_id": "r", "label": "B", "branch_length": 1.}]
    subdivided = [base[0], {"node_id": "u", "parent_id": "r", "label": "u", "branch_length": .5},
        {"node_id": "A", "parent_id": "u", "label": "A", "branch_length": .5}, base[2]]
    return SpeciesTree(base), SpeciesTree(subdivided)


def source_example():
    c = Catalogue("g", "u", 10, (E(0, 10),), (), (Material("m", 0, 10),),
                  boundary_candidates=(E(0, 10),))
    space = enumerate_space(c)
    observations = [O("A", (C((E(0, 10),), (1,)),), material_presence=(1,)),
                    O("B", (C((), (0,)),), material_presence=(0,))]
    tips = {o.species: compatibility(space, o, "annotation") for o in observations}
    model = RateModel({k: .3 if k == "dna_insertion" else .2 if k == "dna_deletion" else 0. for k in EDIT_KINDS})
    return space, tips, model


class AuditModelTests(unittest.TestCase):
    def test_relabelled_same_region_rejected_by_writer(self):
        c = Catalogue("g", "u", 10, (E(0, 10),), ())
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "Overlapping physical"):
                write_catalogues(Path(tmp)/"c.jsonl", (c, replace(c, unit="renamed")))

    def test_relabelled_same_region_rejected_by_reader(self):
        c = Catalogue("g", "u", 10, (E(0, 10),), ())
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"c.jsonl"
            p.write_text("\n".join(json.dumps(encode_catalogue(x)) for x in (c, replace(c, unit="other"))))
            with self.assertRaisesRegex(ValueError, "Overlapping physical"):
                read_catalogues(p)

    def test_adjacent_regions_are_not_overlap(self):
        c = Catalogue("g", "u", 10, (E(0, 10),), ())
        self.assertEqual(len(validate_collection((c, replace(c, unit="v", alignment_offset=10)))), 2)

    def test_overlap_across_files_checked_together(self):
        c = Catalogue("g", "u", 10, (E(0, 10),), ())
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp)/"a.jsonl", Path(tmp)/"b.jsonl"]
            write_catalogues(paths[0], (c,)); write_catalogues(paths[1], (replace(c, unit="v"),))
            with self.assertRaisesRegex(ValueError, "Overlapping physical"):
                load_units(paths, trees()[0])

    def test_unary_subdivision_does_not_change_origin_likelihood_or_root(self):
        space, tips, model = source_example()
        a, b = [evaluate_model(space, t, tips, model, counts=True) for t in trees()]
        self.assertAlmostEqual(a["log_likelihood"], b["log_likelihood"], places=12)
        np.testing.assert_allclose(a["nodes"]["r"], b["nodes"]["r"], atol=1e-12)
        self.assertNotIn("u", b["nodes"])
        self.assertEqual(b["tree_normalization"]["normalized_edge_to_input_children"]["A"], ("u", "A"))

    def test_real_foreground_boundary_not_erased(self):
        context = canonical_tree(trees()[1], frozenset({"A"}))
        self.assertIn("u", context.tree.parent)
        all_foreground = canonical_tree(trees()[1], frozenset({"u", "A"}))
        self.assertNotIn("u", all_foreground.tree.parent)
        self.assertEqual(all_foreground.foreground, frozenset({"A"}))

    def test_unknown_lengths_are_not_invented(self):
        _, tree = trees(); tree.length["u"] = None
        self.assertIsNone(canonical_tree(tree).tree.length["A"])

    def test_origin_weight_prior_is_normalized_before_data(self):
        space, _, model = source_example()
        for w in (.1, 1., 10.):
            result = evaluate_model(space, trees()[0], {s: np.ones(len(space.states)) for s in "AB"},
                                    replace(model, origin_root_weight=w), counts=False)
            self.assertAlmostEqual(result["log_likelihood"], 0., places=11)

    def test_invalid_origin_weight_rejected(self):
        _, _, model = source_example()
        for w in (0., -1., float("nan")):
            with self.assertRaises(ValueError): replace(model, origin_root_weight=w)

    def test_probability_any_edit_available_without_marked_counts(self):
        space, tips, model = source_example()
        result = evaluate_model(space, trees()[0], tips, model, counts=False)
        self.assertIn("probability_at_least_one_edit", result["branches"][0])
        self.assertNotIn("expected_edits", result["branches"][0])
        _, marks = generator(space, model, {"m": "A"}, "A", include_marks=False)
        self.assertEqual(marks, {})

    def test_geometry_preflight_matches_bruteforce(self):
        spans = (E(0, 2), E(0, 5), E(3, 5), E(6, 9))
        count = sum(all(a.end < b.start for a, b in zip(xs, xs[1:]))
                    for n in range(len(spans)+1) for xs in combinations(sorted(spans), n))
        self.assertEqual(estimate_geometry_count(spans), count)

    def test_four_exons_are_usable_without_discarding_states(self):
        exons = tuple(E(i*30, i*30+20) for i in range(4)); full = E(0, exons[-1].end)
        junctions = tuple((a.end, b.start) for a, b in zip(exons, exons[1:]))
        material = tuple(Material(str(i), a, b) for i, (a, b) in enumerate(junctions))
        c = Catalogue("g", "u", full.end, (full,)+exons, junctions, material,
                      boundary_candidates=(full,)+exons)
        space = enumerate_space(c)
        self.assertTrue(space.complete); self.assertEqual(len(space.states), 556)
        self.assertFalse(enumerate_space(c, 256).complete)
        self.assertEqual(space.diagnostics["geometry_count_exact"], 34)

    def test_sparse_shortest_paths_match_floyd_warshall(self):
        c = Catalogue("g", "u", 100, (E(0, 100), E(0, 40), E(60, 100)), ((40, 60),))
        space = enumerate_space(c); graph = make_graph(space, {}, "A")
        d = np.full_like(graph.distance, np.inf); np.fill_diagonal(d, 0)
        for i, j, _, w in graph.edges: d[i, j] = min(w, d[i, j])
        for k in range(len(d)): d = np.minimum(d, d[:, k, None]+d[None, k, :])
        np.testing.assert_allclose(d, graph.distance)

    def payload_catalogue(self):
        a, b = E(10, 30), E(60, 90)
        payload = InsertionPayload("m", (a, b), "explicit_source_locus", "source_sequence_and_exon_annotation")
        c = Catalogue("g", "u", 100, (a, b), (), (Material("m", 0, 100),),
                      boundary_candidates=(a, b), insertion_payloads=(payload,))
        return c, C((), (0,)), C((a, b), (1,))

    def test_source_supported_module_arrives_in_one_edit(self):
        c, source, target = self.payload_catalogue(); space = enumerate_space(c)
        graph = make_graph(space, {"m": "A"}, "A")
        self.assertEqual(graph.distance[space.index[source], space.index[target]], 1.)
        edits = graph.witness_path(space.index[source], space.index[target])
        self.assertEqual([e.kind for e in edits], ["dna_insertion"])
        self.assertIn("introduced_exons:2", edits[0].consequences)

    def test_module_evidence_is_required_and_no_deleted_source_revival(self):
        c, _, target = self.payload_catalogue(); space = enumerate_space(c)
        self.assertFalse(any(e.source.material == (2,) and e.target.material == (1,) for e in space.edits))
        with self.assertRaises(ValueError): InsertionPayload("m", target.exons, "", "")
        with self.assertRaises(ValueError): replace(c, insertion_payloads=(replace(c.insertion_payloads[0], material_id="absent"),))

    def test_duplicate_payload_does_not_increase_generator(self):
        c, _, _ = self.payload_catalogue()
        first = enumerate_space(c); second = enumerate_space(replace(c, insertion_payloads=c.insertion_payloads*2))
        model = RateModel({k: .2 for k in EDIT_KINDS})
        np.testing.assert_allclose(generator(first, model, {"m": "A"}, "A")[0],
                                   generator(second, model, {"m": "A"}, "A")[0])

    def test_atomic_alternatives_do_not_create_a_chimeric_combination(self):
        native = C((E(10, 90),)); left = C((E(0, 90),)); right = C((E(10, 100),))
        alternatives = tuple(A(v, native.key, "B", f"tx{i}", "test", 1.) for i, v in enumerate((left, right)))
        c = Catalogue("g", "u", 110, (E(10, 90), E(0, 90), E(10, 100)), ())
        space = enumerate_space(c); o = O("A", (native,), "partial", alternatives=alternatives)
        weights = compatibility(space, o)
        self.assertTrue(weights[space.index[left]]); self.assertTrue(weights[space.index[right]])
        self.assertFalse(weights[space.index[C((E(0, 100),))]])
        self.assertEqual(compatibility(space, o, "annotation").sum(), 1)

    def test_payload_and_alternatives_schema_roundtrip(self):
        c, _, target = self.payload_catalogue()
        o = O("A", (target,), material_presence=(1,), alternatives=(A(C((), (1,)), target.key, "B", "tx", "test", 1.),))
        c = replace(c, observations=(o,))
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"c.jsonl"; write_catalogues(p, (c,))
            self.assertEqual(read_catalogues(p), (c,))

    def test_bootstrap_duplicate_draw_reaches_numeric_fit(self):
        units, model = fixture(2)
        result = fit_scale(units, model, weights={"g0": 2, "g1": 0})
        self.assertNotEqual(result["status"], "fewer_than_two_genes")
        self.assertIn("scale", result)
        self.assertEqual(result["positive_weight_gene_count"], 1)
        self.assertEqual(result["weighted_gene_count"], 2)

    def test_negative_and_unknown_gene_weights_rejected(self):
        units, model = fixture(2)
        for weights in ({"g0": -1}, {"unknown": 2}):
            with self.assertRaises(ValueError): fit_scale(units, model, weights=weights)

    def test_nested_fit_failure_is_not_zero_statistic(self):
        self.assertIsNone(nested_statistic({"log_likelihood": -10}, {"log_likelihood": -11}))
        self.assertEqual(nested_statistic({"log_likelihood": -10}, {"log_likelihood": -10-1e-9}), 0)
        self.assertEqual(nested_statistic({"log_likelihood": -10}, {"log_likelihood": -9}), 2)

    def test_mixed_unary_rate_regime_cannot_change_the_origin_model(self):
        space, tips, model = source_example()
        with self.assertRaisesRegex(ValueError, "Mixed rate regimes"):
            evaluate_model(space, trees()[1], tips, replace(model, foreground=frozenset({"A"})))

    def test_bounded_cache_evicts_without_changing_values(self):
        from intraphy.inference.kernel_cache import KernelCache
        cache = KernelCache(32)
        a = cache.get_or_compute("a", lambda: np.arange(4, dtype=float))
        b = cache.get_or_compute("b", lambda: np.arange(4, dtype=float)+1)
        again = cache.get_or_compute("a", lambda: np.arange(4, dtype=float))
        self.assertLessEqual(cache.bytes, 32)
        np.testing.assert_array_equal(a, again)
        self.assertEqual(cache.misses, 3)
