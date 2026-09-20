import shutil

import tempfile

import unittest

from pathlib import Path

from unittest.mock import patch

from intraphy.alignment import AlignmentBackendError, AlignmentStats, global_alignment_stats, phase_compatibility

from intraphy.annotation import _overlapping_annotation_role

from intraphy.cli import run_all

from intraphy.correspondence import simple_identity

from intraphy.io import read_tsv

from intraphy.benchmark import benchmark_events

from intraphy.calibration import calibrate_simulations

from intraphy.case import build_case, inspect_annotation, scan_hidden_segments

from intraphy.orthofinder import import_orthofinder

from intraphy.preprocess import derive_tables, extract_gene, graph_components

from intraphy.simulate import simulate_dataset

from intraphy.structural_phylogeny import build_structural_site_matrix, fit_model, infer_single_copy_phylogeny

from intraphy.tree import SpeciesTree

from intraphy.visualize import visualize_results

ROOT = Path(__file__).resolve().parents[1]

from support_intraphy_fixture import FixtureTestsSupport

class FixtureTests(FixtureTestsSupport, unittest.TestCase):
    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_jingwei_case_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_all(ROOT / "demos" / "jingwei", tmp, analysis_scope="experimental-multicopy")
            rows = read_tsv(Path(tmp) / "case_summary.tsv")
            scores = read_tsv(Path(tmp) / "character_model_scores.tsv")
            tests = read_tsv(Path(tmp) / "hypothesis_tests.tsv")
            baselines = read_tsv(Path(tmp) / "baseline_comparison.tsv")
            progressive = read_tsv(Path(tmp) / "progressive_correspondence.tsv")
            progressive_elements = read_tsv(Path(tmp) / "progressive_element_correspondence.tsv")
            ancestral_graph = read_tsv(Path(tmp) / "observed_element_tree_coverage.tsv")
            paths = read_tsv(Path(tmp) / "observed_intragenic_paths.tsv")
            events = read_tsv(Path(tmp) / "candidate_structural_events.tsv")
            hints = read_tsv(Path(tmp) / "interpretation_hints.tsv")
            elements = read_tsv(Path(tmp) / "element_correspondence.tsv")
            element_coverage = read_tsv(Path(tmp) / "element_phylogenetic_coverage.tsv")
            coverage = read_tsv(Path(tmp) / "internal_homology_phylogenetic_coverage.tsv")
            self.assertEqual(rows[0]["family_id"], "jingwei")
            self.assertEqual(rows[0]["hidden_segment_candidates"], "0")
            self.assertEqual(rows[0]["best_compound_model"], "compound_chimeric_or_copy_event")
            self.assertTrue(scores)
            self.assertTrue(elements)
            self.assertTrue(all(row["element_id"].startswith("EG_") for row in elements))
            self.assertTrue(element_coverage)
            self.assertIn("p_value", tests[0])
            self.assertIn("q_value", tests[0])
            self.assertTrue(progressive)
            self.assertTrue(progressive_elements)
            self.assertTrue(ancestral_graph)
            self.assertTrue(paths)
            self.assertTrue(events)
            self.assertFalse(("mechanism_" + "hypothesis") in events[0])
            self.assertIn("structural_change_type", events[0])
            self.assertTrue(hints)
            self.assertIn("possible_interpretation", hints[0])
            self.assertTrue(coverage)
            self.assertIn("coverage_class", coverage[0])
            self.assertEqual({row["baseline_model"] for row in baselines}, {"not_evaluated"})

    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_sdic_case_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_all(ROOT / "demos" / "sdic", tmp, analysis_scope="experimental-multicopy")
            rows = read_tsv(Path(tmp) / "case_summary.tsv")
            self.assertEqual(rows[0]["family_id"], "sdic")
            self.assertEqual(rows[0]["hidden_segment_candidates"], "0")
            self.assertEqual(rows[0]["best_annotation_model"], "strict_annotation")

    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_visualize_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            out = tmp / "out"
            fig = tmp / "fig"
            run_all(
                ROOT / "demos" / "jingwei",
                out,
                bootstrap_replicates=2,
                stochastic_maps=2,
                analysis_scope="experimental-multicopy",
            )
            visualize_results(ROOT / "demos" / "jingwei", out, fig, layout="legacy-overview")
            manifest = read_tsv(fig / "visualization_manifest.tsv")
            self.assertEqual({row["description"] for row in manifest}, {"Exon-like gene-internal synteny by species and copy", "Phylogenetic structural event map", "Integrated phylogenetic and exon-like synteny map"})
            synteny_svg = (fig / "intragenic_synteny.svg").read_text()
            self.assertIn("EG_", synteny_svg)
            self.assertIn("#0072B2", synteny_svg)
            self.assertTrue((fig / "integrated_phylo_synteny.svg").exists())
            fig_pattern = tmp / "fig_pattern"
            visualize_results(ROOT / "demos" / "jingwei", out, fig_pattern, correspondence_encoding="pattern", layout="legacy-overview")
            synteny_pattern_svg = (fig_pattern / "intragenic_synteny.svg").read_text()
            self.assertIn("<pattern", synteny_pattern_svg)
            self.assertIn('stroke-dasharray="2 2"', synteny_pattern_svg)
