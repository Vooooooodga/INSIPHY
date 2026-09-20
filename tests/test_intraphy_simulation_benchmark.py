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

from support_intraphy_simulation_benchmark import SimulationBenchmarkTestsSupport

class SimulationBenchmarkTests(SimulationBenchmarkTestsSupport, unittest.TestCase):
    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_simulation_benchmark(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sim = tmp / "sim"
            out = tmp / "out"
            simulate_dataset(sim, seed=7)
            run_all(sim, out, analysis_scope="experimental-multicopy")
            benchmark_events(sim, out)
            truth = read_tsv(sim / "truth_events.tsv")
            bench = read_tsv(out / "benchmark_summary.tsv")
            fit = read_tsv(out / "model_fit.tsv")
            element_coverage = read_tsv(out / "element_phylogenetic_coverage.tsv")
            self.assertTrue(truth)
            self.assertGreaterEqual(float(bench[0]["recall"]), 0.0)
            self.assertIn("precision", bench[0])
            self.assertTrue(fit)
            self.assertTrue(element_coverage)
            self.assertEqual(fit[0]["model"], "ctmc_mk_branch_length")

    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_named_simulation_scenarios(self):
        scenarios = [
            "exonization",
            "te_exonization",
            "splice_boundary_shift",
            "segment_split",
            "segment_fusion",
            "source_join",
            "tandem_duplication",
            "processed_copy_or_intron_loss",
            "annotation_dropout",
            "negative_control",
            "gene_conversion",
        ]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for scenario in scenarios:
                sim = tmp / f"sim_{scenario}"
                out = tmp / f"out_{scenario}"
                simulate_dataset(sim, seed=19, scenario=scenario)
                run_all(sim, out, analysis_scope="experimental-multicopy")
                benchmark_events(sim, out)
                self.assertTrue(read_tsv(out / "element_correspondence.tsv"))
                self.assertTrue(read_tsv(out / "candidate_structural_events.tsv", optional=True) or scenario in {"negative_control", "annotation_dropout"})
                self.assertTrue(read_tsv(out / "benchmark_summary.tsv"))

    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_multicopy_structural_characters_use_copy_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sim = tmp / "sim_compound"
            out = tmp / "out_compound"
            simulate_dataset(sim, seed=23, scenario="compound")
            run_all(sim, out, analysis_scope="experimental-multicopy")
            benchmark_events(sim, out)
            fig = tmp / "fig_compound"
            visualize_results(sim, out, fig, layout="legacy-overview")
            scope = read_tsv(out / "phylogeny_scope.tsv")
            branches = read_tsv(out / "branch_event_probabilities.tsv")
            bench = read_tsv(out / "benchmark_summary.tsv")
            integrated_svg = (fig / "integrated_phylo_synteny.svg").read_text()
            self.assertEqual(bench[0]["precision"], "1")
            self.assertEqual(bench[0]["recall"], "1")
            self.assertEqual(bench[0]["branch_accuracy"], "1")
            self.assertEqual(
                [row for row in scope if row["scope"] == "structural_characters"][0]["tree_scope"],
                "copy_tree",
            )
            self.assertIn("copy_tree.tsv", integrated_svg)
            self.assertIn("clade34_copy1", integrated_svg)
            self.assertTrue(
                any(
                    row["layer"] == "element_presence"
                    and row["object_id"] == "EG_sim_C"
                    and row["parent_label"] == "clade34_copy1"
                    and row["child_label"] in {"Sp3:copy1", "Sp4:copy1"}
                    for row in branches
                )
            )
            self.assertTrue(
                any(
                    row["layer"] == "element_presence"
                    and row["object_id"] == "EG_sim_C"
                    and row["parent_label"] == "clade34_dup"
                    and row["child_label"] == "clade34_copy2"
                    for row in branches
                )
            )

    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_bootstrap_and_stochastic_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sim = tmp / "sim_stats"
            out = tmp / "out_stats"
            simulate_dataset(sim, seed=13, scenario="exonization")
            run_all(
                sim,
                out,
                bootstrap_replicates=5,
                stochastic_maps=5,
                seed=13,
                analysis_scope="experimental-multicopy",
            )
            benchmark_events(sim, out)
            boot = read_tsv(out / "hypothesis_bootstrap.tsv")
            histories = read_tsv(out / "branch_history_posteriors.tsv")
            calibration = read_tsv(out / "benchmark_calibration.tsv")
            self.assertTrue(boot)
            self.assertIn("empirical_p_value", boot[0])
            self.assertTrue(histories)
            self.assertIn("posterior_pr_any_change", histories[0])
            self.assertEqual(calibration[0]["bootstrap_tests"], str(len(boot)))

    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_simulation_calibration_operating_characteristics(self):
        with tempfile.TemporaryDirectory() as tmp:
            calibrate_simulations(tmp, scenarios=["negative_control"], replicates=1, seed=29)
            summary = read_tsv(Path(tmp) / "calibration_operating_characteristics.tsv")
            reps = read_tsv(Path(tmp) / "calibration_replicates.tsv")
            self.assertEqual(summary[0]["scenario"], "negative_control")
            self.assertEqual(summary[0]["replicates"], "1")
            self.assertTrue(reps)

    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_annotation_dropout_negative_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sim = tmp / "sim_dropout"
            out = tmp / "out_dropout"
            simulate_dataset(sim, seed=11, scenario="annotation_dropout")
            run_all(sim, out, analysis_scope="experimental-multicopy")
            benchmark_events(sim, out)
            annot = read_tsv(out / "annotation_completion_candidates.tsv")
            bench = read_tsv(out / "benchmark_summary.tsv")
            self.assertIn("ambiguous_evidence", {row["completion_call"] for row in annot})
            self.assertEqual(bench[0]["truth_events"], "0")
            self.assertEqual(bench[0]["called_events"], "0")

    @unittest.skipUnless((shutil.which("mafft") and shutil.which("minimap2")), "MAFFT and minimap2 are required")
    def test_paralogous_similarity_is_ambiguous_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sim = tmp / "sim_gene_conversion"
            out = tmp / "out_gene_conversion"
            simulate_dataset(sim, seed=17, scenario="gene_conversion")
            run_all(sim, out, analysis_scope="experimental-multicopy")
            benchmark_events(sim, out)
            events = read_tsv(out / "candidate_structural_events.tsv")
            bench = read_tsv(out / "benchmark_summary.tsv")
            ambiguous = [row for row in events if row.get("structural_pattern") == "ambiguous_paralogous_similarity"]
            self.assertTrue(ambiguous)
            self.assertEqual({row["call_scope"] for row in ambiguous}, {"ambiguous_evidence"})
            self.assertGreaterEqual(int(bench[0]["ambiguous_evidence_called"]), 1)
