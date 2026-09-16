import tempfile
import unittest
from pathlib import Path

from insiphy.cli import run_all
from insiphy.io import read_tsv
from insiphy.benchmark import benchmark_events
from insiphy.case import build_case, inspect_annotation, scan_hidden_segments
from insiphy.preprocess import derive_tables, extract_gene
from insiphy.simulate import simulate_dataset
from insiphy.visualize import visualize_results


ROOT = Path(__file__).resolve().parents[1]


class FixtureTests(unittest.TestCase):
    def test_jingwei_case_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_all(ROOT / "demos" / "jingwei", tmp)
            rows = read_tsv(Path(tmp) / "case_summary.tsv")
            scores = read_tsv(Path(tmp) / "character_model_scores.tsv")
            tests = read_tsv(Path(tmp) / "hypothesis_tests.tsv")
            baselines = read_tsv(Path(tmp) / "baseline_comparison.tsv")
            progressive = read_tsv(Path(tmp) / "progressive_correspondence.tsv")
            self.assertEqual(rows[0]["family_id"], "jingwei")
            self.assertEqual(rows[0]["hidden_segment_candidates"], "1")
            self.assertEqual(rows[0]["best_compound_model"], "compound_chimeric_or_copy_event")
            self.assertTrue(scores)
            self.assertIn("p_value", tests[0])
            self.assertIn("q_value", tests[0])
            self.assertTrue(progressive)
            self.assertEqual({row["baseline_model"] for row in baselines}, {"annotation_only", "sequence_only", "synteny_aware_phylogenetic"})

    def test_sdic_case_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_all(ROOT / "demos" / "sdic", tmp)
            rows = read_tsv(Path(tmp) / "case_summary.tsv")
            self.assertEqual(rows[0]["family_id"], "sdic")
            self.assertEqual(rows[0]["hidden_segment_candidates"], "1")
            self.assertEqual(rows[0]["best_annotation_model"], "annotation_error")

    def test_visualize_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            out = tmp / "out"
            fig = tmp / "fig"
            run_all(ROOT / "demos" / "jingwei", out, bootstrap_replicates=2, stochastic_maps=2)
            visualize_results(ROOT / "demos" / "jingwei", out, fig)
            manifest = read_tsv(fig / "visualization_manifest.tsv")
            self.assertEqual({row["description"] for row in manifest}, {"Gene-internal synteny by species and copy", "Species-tree structural event map"})
            self.assertIn("<pattern", (fig / "intragenic_synteny.svg").read_text())


class PreprocessTests(unittest.TestCase):
    def test_extract_gene_and_derive_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genome = tmp / "genome.fa"
            annot = tmp / "annotation.gff3"
            out = tmp / "work"
            genome.write_text(">chr1\n" + "ACGT" * 80 + "\n")
            annot.write_text(
                "\n".join(
                    [
                        "chr1\tINSIPHY\tgene\t10\t120\t.\t+\t.\tID=geneA;Name=GeneA",
                        "chr1\tINSIPHY\tmRNA\t10\t120\t.\t+\t.\tID=txA;Parent=geneA",
                        "chr1\tINSIPHY\texon\t10\t40\t.\t+\t.\tID=ex1;Parent=txA",
                        "chr1\tINSIPHY\tCDS\t15\t40\t.\t+\t0\tID=cds1;Parent=txA",
                        "chr1\tINSIPHY\texon\t80\t120\t.\t+\t.\tID=ex2;Parent=txA",
                        "chr1\tINSIPHY\tCDS\t80\t110\t.\t+\t2\tID=cds2;Parent=txA",
                    ]
                )
                + "\n"
            )
            extract_gene(genome, annot, "GeneA", "family_a", "SpeciesA", "SpeciesA_geneA", out)
            derive_tables(out, identity_threshold=0.5)
            occ = read_tsv(out / "segment_occurrences.tsv")
            adj = read_tsv(out / "physical_adjacencies.tsv")
            hom = read_tsv(out / "segment_homology.tsv")
            matches = read_tsv(out / "segment_matches.tsv")
            paths = read_tsv(out / "transcript_paths.tsv")
            introns = read_tsv(out / "intron_sites.tsv")
            copy_context = read_tsv(out / "copy_context.tsv")
            self.assertEqual(len(occ), 3)
            self.assertEqual(len(adj), 2)
            self.assertEqual({row["role"] for row in occ}, {"CDS", "intron"})
            self.assertTrue(hom)
            self.assertIn("phase_score", matches[0])
            self.assertTrue(paths)
            self.assertEqual(introns[0]["phase_compatibility"], "incompatible")
            self.assertIn("copy_subclass", copy_context[0])


class SimulationBenchmarkTests(unittest.TestCase):
    def test_simulation_benchmark(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sim = tmp / "sim"
            out = tmp / "out"
            simulate_dataset(sim, seed=7)
            run_all(sim, out)
            benchmark_events(sim, out)
            truth = read_tsv(sim / "truth_events.tsv")
            bench = read_tsv(out / "benchmark_summary.tsv")
            fit = read_tsv(out / "model_fit.tsv")
            self.assertTrue(truth)
            self.assertGreaterEqual(float(bench[0]["recall"]), 0.0)
            self.assertIn("precision", bench[0])
            self.assertTrue(fit)
            self.assertEqual(fit[0]["model"], "ctmc_mk_branch_length")

    def test_bootstrap_and_stochastic_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sim = tmp / "sim_stats"
            out = tmp / "out_stats"
            simulate_dataset(sim, seed=13, scenario="exonization")
            run_all(sim, out, bootstrap_replicates=5, stochastic_maps=5, seed=13)
            benchmark_events(sim, out)
            boot = read_tsv(out / "hypothesis_bootstrap.tsv")
            histories = read_tsv(out / "branch_history_posteriors.tsv")
            calibration = read_tsv(out / "benchmark_calibration.tsv")
            self.assertTrue(boot)
            self.assertIn("empirical_p_value", boot[0])
            self.assertTrue(histories)
            self.assertIn("posterior_pr_any_change", histories[0])
            self.assertEqual(calibration[0]["bootstrap_tests"], str(len(boot)))

    def test_annotation_dropout_negative_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sim = tmp / "sim_dropout"
            out = tmp / "out_dropout"
            simulate_dataset(sim, seed=11, scenario="annotation_dropout")
            run_all(sim, out)
            benchmark_events(sim, out)
            annot = read_tsv(out / "annotation_completion_candidates.tsv")
            bench = read_tsv(out / "benchmark_summary.tsv")
            self.assertIn("hidden_segment_candidate", {row["completion_call"] for row in annot})
            self.assertEqual(bench[0]["truth_events"], "0")
            self.assertEqual(bench[0]["called_events"], "0")


class RealCasePreparationTests(unittest.TestCase):
    def write_case_files(self, tmp, species):
        genome = tmp / f"{species}.fa"
        annot = tmp / f"{species}.gff3"
        genome.write_text(">chr1\n" + "ACGT" * 120 + "\n")
        annot.write_text(
            "\n".join(
                [
                    f"chr1\tINSIPHY\tgene\t10\t160\t.\t+\t.\tID={species}_geneA;Name=GeneA;Alias=jgw,jingwei",
                    f"chr1\tINSIPHY\tmRNA\t10\t160\t.\t+\t.\tID={species}_txA;Parent={species}_geneA",
                    f"chr1\tINSIPHY\tCDS\t20\t50\t.\t+\t0\tID={species}_cds1;Parent={species}_txA",
                    f"chr1\tINSIPHY\tCDS\t100\t140\t.\t+\t1\tID={species}_cds2;Parent={species}_txA",
                ]
            )
            + "\n"
        )
        return genome, annot

    def test_inspect_annotation_aliases(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genome, annot = self.write_case_files(tmp, "SpA")
            out = tmp / "inspect"
            inspect_annotation(annot, out, queries=["jingwei"], species="SpA", case_id="case1")
            rows = read_tsv(out / "gene_candidate_report.tsv")
            self.assertTrue(rows)
            self.assertEqual(rows[0]["feature_id"], "SpA_geneA")
            self.assertEqual(rows[0]["match_rank"], "exact_token")

    def test_build_case_from_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            g1, a1 = self.write_case_files(tmp, "SpA")
            g2, a2 = self.write_case_files(tmp, "SpB")
            manifest = tmp / "manifest.tsv"
            tree = tmp / "species_tree.tsv"
            out = tmp / "case"
            manifest.write_text(
                "\n".join(
                    [
                        "case_id\tspecies\tfamily_id\tgene_id\tgene_copy_id\tgenome_fasta\tannotation_file\tassembly\tannotation\tsource_url\trelease\tnotes",
                        f"case1\tSpA\tfam1\tSpA_geneA\tSpA_geneA\t{g1}\t{a1}\tasmA\tannA\tlocal\tv1\tok",
                        f"case1\tSpB\tfam1\tSpB_geneA\tSpB_geneA\t{g2}\t{a2}\tasmB\tannB\tlocal\tv1\tok",
                    ]
                )
                + "\n"
            )
            tree.write_text("node_id\tparent_id\tlabel\nroot\t\troot\nspa\troot\tSpA\nspb\troot\tSpB\n")
            build_case(manifest, out, species_tree=tree)
            report = read_tsv(out / "case_build_report.tsv")
            prov = read_tsv(out / "case_provenance.tsv")
            occ = read_tsv(out / "segment_occurrences.tsv")
            self.assertEqual({row["status"] for row in report}, {"extracted"})
            self.assertEqual(len(prov), 2)
            self.assertTrue(occ)
            self.assertTrue((out / "species_tree.tsv").exists())

    def test_scan_hidden_segments(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "source.fa"
            target = tmp / "target.fa"
            out = tmp / "scan"
            source.write_text(">hidden_seg\nAACCGGTTAACC\n")
            target.write_text(">gene_interval\nTTTTTAACCGGTTAACCGGGGG\n")
            scan_hidden_segments(source, target, out, family_id="fam", species="SpA", gene_copy_id="copy1", min_identity=0.9)
            rows = read_tsv(out / "hidden_segment_scan.tsv")
            self.assertEqual(rows[0]["support_call"], "hidden_segment_candidate")
            self.assertGreaterEqual(float(rows[0]["identity"]), 0.9)
            self.assertIn("alignment_cigar", rows[0])
            self.assertIn("frame_status", rows[0])


if __name__ == "__main__":
    unittest.main()
