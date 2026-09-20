import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from insiphy.alignment import AlignmentBackendError, AlignmentStats, global_alignment_stats, phase_compatibility
from insiphy.annotation import _overlapping_annotation_role
from insiphy.cli import run_all
from insiphy.correspondence import simple_identity
from insiphy.io import read_tsv
from insiphy.benchmark import benchmark_events
from insiphy.calibration import calibrate_simulations
from insiphy.case import build_case, inspect_annotation, scan_hidden_segments
from insiphy.orthofinder import import_orthofinder
from insiphy.preprocess import derive_tables, extract_gene, graph_components
from insiphy.simulate import simulate_dataset
from insiphy.structural_phylogeny import build_structural_site_matrix, fit_model, infer_single_copy_phylogeny
from insiphy.tree import SpeciesTree
from insiphy.visualize import visualize_results


ROOT = Path(__file__).resolve().parents[1]


class SingleCopyPhylogenyTests(unittest.TestCase):
    def write_structural_case(self, root):
        input_dir = root / "input"
        result_dir = root / "result"
        input_dir.mkdir()
        result_dir.mkdir()
        (input_dir / "species_tree.tsv").write_text(
            "node_id\tparent_id\tlabel\tbranch_length\n"
            "root\t\troot\t0\n"
            "ab\troot\tab\t0.5\n"
            "a\tab\tA\t0.2\n"
            "b\tab\tB\t0.2\n"
            "cd\troot\tcd\t0.5\n"
            "c\tcd\tC\t0.2\n"
            "d\tcd\tD\t0.2\n"
        )
        occurrence_rows = [
            "occurrence_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tsegment_id\tcontig\tstart\tend\tstrand\trole\tpresence_status\tboundary_state\tevidence"
        ]
        element_rows = [
            "element_id\tfamily_id\thomology_id\toccurrence_id\tspecies\tgene_copy_id\telement_class\tdisplay_role\tsource_label\tsupport_type\tconfidence\tmembership_score\tmembership_call"
        ]
        path_rows = [
            "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\tcontig\tstart\tend\tstrand\tphase\tpath_status"
        ]
        for species in "ABCD":
            roles = [
                "CDS",
                "CDS" if species in "AB" else "intron",
                "CDS" if species in "AC" else "intron",
                "CDS" if species in "AD" else "intron",
            ]
            for index, role in enumerate(roles, start=1):
                occurrence = f"{species}_e{index}"
                occurrence_rows.append(
                    f"{occurrence}\tfam\t{species}\t{species}_gene\t{species}_tx\te{index}\tchr1\t{index * 100}\t{index * 100 + 50}\t+\t{role}\tpresent\tconserved\ttest_fixture"
                )
                element_class = "exon_like" if role == "CDS" else "candidate_source"
                element_rows.append(
                    f"EG_{index}\tfam\tH_{index}\t{occurrence}\t{species}\t{species}_gene\t{element_class}\t{role}\tfam\ttest_fixture\thigh\t0.95\tcore_member"
                )
                path_rows.append(
                    f"{species}_p{index}\tfam\t{species}\t{species}_gene\t{species}_tx\t{index}\t{occurrence}\t{role}\tchr1\t{index * 100}\t{index * 100 + 50}\t+\t0\tcanonical_transcript_path"
                )
        (input_dir / "segment_occurrences.tsv").write_text("\n".join(occurrence_rows) + "\n")
        (input_dir / "transcript_paths.tsv").write_text("\n".join(path_rows) + "\n")
        (result_dir / "element_correspondence.tsv").write_text("\n".join(element_rows) + "\n")
        return input_dir, result_dir

    def test_single_copy_models_and_posteriors(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_structural_case(Path(tmp))
            infer_single_copy_phylogeny(input_dir, result_dir, threads=2)
            matrix = read_tsv(result_dir / "structural_site_matrix.tsv")
            tests = read_tsv(result_dir / "model_tests.tsv")
            fits = read_tsv(result_dir / "model_fits.tsv")
            nodes = read_tsv(result_dir / "node_state_posteriors.tsv")
            branches = read_tsv(result_dir / "branch_transition_posteriors.tsv")
            self.assertEqual({row["layer"] for row in matrix}, {"exon_presence", "exon_role", "splice_junction"})
            self.assertEqual({row["model"] for row in fits}, {"ER", "ARD"})
            self.assertTrue(
                all(
                    row["test_status"] in {
                        "tested",
                        "parameters_not_estimable",
                        "optimization_failure_alternative_below_null",
                    }
                    for row in tests
                )
            )
            invariant = [row for row in tests if row["layer"] == "exon_presence"][0]
            self.assertEqual(invariant["p_value"], "NA")
            self.assertTrue(all(0 <= float(row["posterior_probability"]) <= 1 for row in nodes))
            self.assertTrue(all("fit_status=success" in row["conditioning"] for row in nodes))
            self.assertTrue(all("fit_status=success" in row["conditioning"] for row in branches))

    def test_foreground_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_structural_case(Path(tmp))
            foreground = Path(tmp) / "foreground.tsv"
            foreground.write_text("parent_id\tchild_id\nab\ta\n")
            infer_single_copy_phylogeny(
                input_dir,
                result_dir,
                model="foreground",
                foreground_branches=foreground,
                threads=2,
            )
            fits = read_tsv(result_dir / "model_fits.tsv")
            self.assertIn("ARD_FOREGROUND", {row["model"] for row in fits})
            self.assertEqual(
                {row["test_id"] for row in read_tsv(result_dir / "model_tests.tsv")},
                {"homogeneous_vs_foreground"},
            )

    def test_foreground_rejects_unmatched_and_all_tree_branches(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_structural_case(Path(tmp))
            unmatched = Path(tmp) / "unmatched.tsv"
            unmatched.write_text("parent_id\tchild_id\nmissing\tbranch\n")
            with self.assertRaises(SystemExit):
                infer_single_copy_phylogeny(
                    input_dir, result_dir, model="foreground", foreground_branches=unmatched
                )
            all_branches = Path(tmp) / "all.tsv"
            all_branches.write_text(
                "parent_id\tchild_id\nroot\tab\nab\ta\nab\tb\nroot\tcd\ncd\tc\ncd\td\n"
            )
            with self.assertRaises(SystemExit):
                infer_single_copy_phylogeny(
                    input_dir, result_dir, model="foreground", foreground_branches=all_branches
                )

    def test_predicted_completion_preserves_explicit_nonexonic_observation(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir = self.write_structural_case(Path(tmp))
            (result_dir / "annotation_completion_candidates.tsv").write_text(
                "evidence_id\tfamily_id\tspecies\tgene_copy_id\thomology_id\tinterval\tannotation_status\tinferred_role\tsupport_score\tcompletion_call\tevidence_status\tinferred_event\tframe_status\n"
                "ev1\tfam\tC\tC_gene\tH_2\tchr1:200-250:+\tmissing_annotation\tCDS\t0.9\thidden_segment_candidate\tsupports_hidden_segment\thidden_exon\tcoding_frame_preserved\n"
            )
            infer_single_copy_phylogeny(input_dir, result_dir)
            states = [
                row
                for row in read_tsv(result_dir / "structural_site_matrix.tsv")
                if row["layer"] == "exon_role" and row["site_id"] == "EG_2" and row["species"] == "C"
            ]
            self.assertEqual(states[0]["state"], "not_exonic")
            self.assertIn("homologous_non_exonic_sequence", states[0]["evidence"])

    def test_orthofinder_single_copy_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            orthofinder = root / "orthofinder" / "Orthogroups"
            orthofinder.mkdir(parents=True)
            (orthofinder / "Orthogroups.tsv").write_text(
                "Orthogroup\tA\tB\nOG0001\tgene_a\tgene_b\n"
            )
            resources = root / "genomes.tsv"
            gff_a = root / "a.gff3"
            gff_b = root / "b.gff3"
            gff_a.write_text("chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=gene_a\n")
            gff_b.write_text("chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=gene_b\n")
            resources.write_text(
                "species\tgenome_fasta\tannotation_file\n"
                f"A\t/a.fa\t{gff_a}\n"
                f"B\t/b.fa\t{gff_b}\n"
            )
            tree = root / "tree.nwk"
            tree.write_text("(A:0.1,B:0.1)root;\n")
            out = root / "imported"
            import_orthofinder(root / "orthofinder", "OG0001", resources, out, tree)
            manifest = read_tsv(out / "manifest.tsv")
            tree_rows = read_tsv(out / "species_tree.tsv")
            self.assertEqual({row["gene_id"] for row in manifest}, {"gene_a", "gene_b"})
            self.assertEqual({row["label"] for row in tree_rows if row["parent_id"]}, {"A", "B"})

    def test_one_to_many_correspondence_and_within_element_junction(self):
        components = graph_components(
            ["A_left", "A_right", "B_whole"],
            [("A_left", "B_whole", 0.9), ("A_right", "B_whole", 0.9)],
        )
        self.assertEqual(components, [["A_left", "A_right", "B_whole"]])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "input"
            result_dir = root / "result"
            input_dir.mkdir()
            result_dir.mkdir()
            (input_dir / "segment_occurrences.tsv").write_text(
                "occurrence_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\trole\tpresence_status\tcontig\tstart\tend\tstrand\n"
                "A_left\tfam\tA\tA_gene\tA_tx\tCDS\tpresent\tchrA\t1\t50\t+\n"
                "A_intron\tfam\tA\tA_gene\tA_tx\tintron\tpresent\tchrA\t51\t100\t+\n"
                "A_right\tfam\tA\tA_gene\tA_tx\tCDS\tpresent\tchrA\t101\t150\t+\n"
                "B_whole\tfam\tB\tB_gene\tB_tx\tCDS\tpresent\tchrB\t1\t100\t+\n"
            )
            (input_dir / "segment_matches.tsv").write_text(
                "match_id\tquery_occurrence_id\tsubject_occurrence_id\tmatch_status\tprojected_reference_blocks\tprojected_reference_strand\n"
                "m1\tA_left\tB_whole\tmapped\t1-50:1-50\t+\n"
                "m2\tA_right\tB_whole\tmapped\t1-50:51-100\t+\n"
            )
            (input_dir / "transcript_paths.tsv").write_text(
                "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\n"
                "p1\tfam\tA\tA_gene\tA_tx\t1\tA_left\tCDS\n"
                "p2\tfam\tA\tA_gene\tA_tx\t2\tA_intron\tintron\n"
                "p3\tfam\tA\tA_gene\tA_tx\t3\tA_right\tCDS\n"
                "p4\tfam\tB\tB_gene\tB_tx\t1\tB_whole\tCDS\n"
            )
            (result_dir / "element_correspondence.tsv").write_text(
                "element_id\tfamily_id\thomology_id\toccurrence_id\tspecies\tgene_copy_id\telement_class\tmembership_call\n"
                "EG_1\tfam\tH_1\tA_left\tA\tA_gene\texon_like\tcore_member\n"
                "EG_1\tfam\tH_1\tA_right\tA\tA_gene\texon_like\tcore_member\n"
                "EG_1\tfam\tH_1\tB_whole\tB\tB_gene\texon_like\tcore_member\n"
            )
            matrix, excluded = build_structural_site_matrix(input_dir, result_dir)
            junctions = {
                row["species"]: row["state"]
                for row in matrix
                if row["site_id"] == "JG_EG_1_REF_B_whole_D50_A51"
            }
            self.assertFalse(excluded)
            self.assertEqual(junctions, {"A": "present", "B": "absent"})

    def test_variable_site_ascertainment(self):
        tree = SpeciesTree(
            [
                {"node_id": "root", "parent_id": "", "label": "root", "branch_length": 0},
                {"node_id": "a", "parent_id": "root", "label": "A", "branch_length": 1},
                {"node_id": "b", "parent_id": "root", "label": "B", "branch_length": 1},
                {"node_id": "c", "parent_id": "root", "label": "C", "branch_length": 1},
                {"node_id": "d", "parent_id": "root", "label": "D", "branch_length": 1},
            ]
        )
        patterns = [
            {"A": 0, "B": 0, "C": 1, "D": 1},
            {"A": 0, "B": 1, "C": 0, "D": 1},
            {"A": 0, "B": 1, "C": 1, "D": 0},
        ]
        fit = fit_model(tree, patterns, "ER", ascertainment="variable-only")
        self.assertTrue(fit["converged"])
        with self.assertRaises(SystemExit):
            fit_model(tree, [{"A": 1, "B": 1, "C": 1, "D": 1}], "ER", ascertainment="variable-only")


class FixtureTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
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
            self.assertEqual(rows[0]["hidden_segment_candidates"], "1")
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

    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
    def test_sdic_case_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_all(ROOT / "demos" / "sdic", tmp, analysis_scope="experimental-multicopy")
            rows = read_tsv(Path(tmp) / "case_summary.tsv")
            self.assertEqual(rows[0]["family_id"], "sdic")
            self.assertEqual(rows[0]["hidden_segment_candidates"], "1")
            self.assertEqual(rows[0]["best_annotation_model"], "annotation_error")

    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
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


class PreprocessTests(unittest.TestCase):
    def test_alignment_coverage_counts_paired_bases(self):
        stats = global_alignment_stats("ACGT", "ACGTTTTT")
        self.assertAlmostEqual(stats.query_coverage, 1.0)
        self.assertAlmostEqual(stats.target_coverage, 0.5)
        self.assertAlmostEqual(stats.coverage, 0.5)
        self.assertEqual(stats.aligned_pairs, 4)

    def test_auto_global_alignment_requires_mafft_for_short_exons(self):
        sequence = "ACGT" * 25
        with patch("insiphy.alignment.shutil.which", return_value="/mock/mafft") as available, patch(
            "insiphy.alignment.subprocess.run"
        ) as run:
            run.return_value.returncode = 0
            run.return_value.stdout = f">query\n{sequence}\n>target\n{sequence}\n"
            stats = global_alignment_stats(sequence, sequence, backend="auto")
        available.assert_called_once_with("mafft")
        run.assert_called_once()
        self.assertEqual(run.call_args.args[0][0], "/mock/mafft")
        self.assertEqual(stats.backend, "mafft")
        self.assertEqual(stats.identity, 1.0)
        self.assertEqual(stats.coverage, 1.0)
        with patch("insiphy.alignment.shutil.which", return_value=None), patch(
            "insiphy.alignment.subprocess.run"
        ) as run:
            with self.assertRaisesRegex(AlignmentBackendError, "MAFFT was requested but is not available"):
                global_alignment_stats(sequence, sequence, backend="auto")
        run.assert_not_called()

    def test_internal_alignment_does_not_count_all_optimal_paths(self):
        import numpy as np

        class Alignment:
            coordinates = np.array([[0, 4], [0, 4]])
            score = 8.0

        class Alignments:
            def __len__(self):
                raise OverflowError("too many optimal alignments")

            def __iter__(self):
                yield Alignment()

        class Aligner:
            def align(self, _left, _right):
                return Alignments()

        with patch("Bio.Align.PairwiseAligner", return_value=Aligner()):
            stats = global_alignment_stats("ACGT", "ACGT")
        self.assertEqual(stats.identity, 1.0)
        self.assertEqual(stats.aligned_pairs, 4)

    def test_correspondence_identity_uses_auto_alignment(self):
        expected = AlignmentStats(0.8, 1.0, 8.0, backend="mafft")
        with patch("insiphy.correspondence.global_alignment_stats", return_value=expected) as align:
            identity = simple_identity("ACGT", "ACGA")
        self.assertEqual(identity, 0.8)
        align.assert_called_once_with("ACGT", "ACGA", backend="auto")

    def test_locus_hit_on_annotated_exon_is_not_hidden_exon_evidence(self):
        plus = [{"role": "exon", "contig": "chr1", "start": "120", "end": "150", "strand": "+"}]
        minus = [{"role": "exon", "contig": "chr1", "start": "250", "end": "280", "strand": "-"}]
        self.assertEqual(_overlapping_annotation_role(20, 50, "Sp|gene|chr1:101-300:+", plus, "+"), ("exon", True))
        self.assertEqual(_overlapping_annotation_role(60, 80, "Sp|gene|chr1:101-300:+", plus, "+"), ("unknown", True))
        self.assertEqual(_overlapping_annotation_role(20, 50, "Sp|gene|chr1:101-300:-", minus, "-"), ("exon", True))
        self.assertEqual(_overlapping_annotation_role(60, 80, "Sp|gene|chr1:101-300:-", minus, "-"), ("unknown", True))
        self.assertEqual(_overlapping_annotation_role(20, 50, "Sp|gene|chr1:101-300:+", plus, "-"), ("exon", False))
        self.assertEqual(_overlapping_annotation_role(20, 50, "Sp|gene|chr1:101-300:-", minus, "+"), ("exon", False))

    def test_internal_alignment_refuses_oversized_dynamic_program(self):
        with self.assertRaises(AlignmentBackendError):
            global_alignment_stats("A" * 501, "A" * 501)

    def test_phase_compatibility_uses_cds_length(self):
        self.assertEqual(phase_compatibility("0", "1", 26), "compatible")
        self.assertEqual(phase_compatibility("0", "2", 26), "incompatible")

    def test_gene_without_exon_cds_or_utr_is_unresolved(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genome = tmp / "genome.fa"
            annotation = tmp / "annotation.gff3"
            genome.write_text(">chr1\n" + "ACGT" * 50 + "\n")
            annotation.write_text(
                "chr1\ttest\tgene\t10\t100\t.\t+\t.\tID=g\n"
                "chr1\ttest\tmRNA\t10\t100\t.\t+\t.\tID=t;Parent=g\n"
            )
            with self.assertRaises(SystemExit):
                extract_gene(genome, annotation, "g", "fam", "Sp", "g", tmp / "out")

    def test_negative_strand_is_written_in_transcript_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genome = tmp / "genome.fa"
            annotation = tmp / "annotation.gff3"
            output = tmp / "out"
            genome.write_text(">chr1\n" + "ACGT" * 100 + "\n")
            annotation.write_text(
                "chr1\ttest\tgene\t20\t180\t.\t-\t.\tID=g\n"
                "chr1\ttest\tmRNA\t20\t180\t.\t-\t.\tID=t;Parent=g\n"
                "chr1\ttest\texon\t20\t60\t.\t-\t.\tID=e1;Parent=t\n"
                "chr1\ttest\tCDS\t30\t60\t.\t-\t1\tID=c1;Parent=t\n"
                "chr1\ttest\texon\t140\t180\t.\t-\t.\tID=e2;Parent=t\n"
                "chr1\ttest\tCDS\t140\t170\t.\t-\t0\tID=c2;Parent=t\n"
            )
            extract_gene(genome, annotation, "g", "fam", "Sp", "g", output)
            exons = [row for row in read_tsv(output / "segment_occurrences.tsv") if row["role"] == "exon"]
            ordered = sorted(exons, key=lambda row: int(row["transcript_order"]))
            self.assertEqual([int(row["start"]) for row in ordered], [140, 20])

    def test_gff_identifier_with_pipe_preserves_parent_relationship(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            genome = tmp / "genome.fa"
            annotation = tmp / "annotation.gff3"
            output = tmp / "out"
            genome.write_text(">chr1\n" + "ACGT" * 100 + "\n")
            transcript_id = "rna-gnl|WGS:TEST|tx1"
            annotation.write_text(
                "chr1\ttest\tgene\t20\t180\t.\t+\t.\tID=gene-g1;Name=g1\n"
                f"chr1\ttest\tmRNA\t20\t180\t.\t+\t.\tID={transcript_id};Parent=gene-g1\n"
                f"chr1\ttest\texon\t20\t60\t.\t+\t.\tID=e1;Parent={transcript_id}\n"
                f"chr1\ttest\tCDS\t30\t60\t.\t+\t0\tID=c1;Parent={transcript_id}\n"
                f"chr1\ttest\texon\t140\t180\t.\t+\t.\tID=e2;Parent={transcript_id}\n"
                f"chr1\ttest\tCDS\t140\t170\t.\t+\t0\tID=c2;Parent={transcript_id}\n"
            )
            extract_gene(genome, annotation, "g1", "fam", "Sp", "g1", output)
            exons = [row for row in read_tsv(output / "segment_occurrences.tsv") if row["role"] == "exon"]
            self.assertEqual(len(exons), 2)
            self.assertEqual({row["transcript_id"] for row in exons}, {transcript_id})

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
            self.assertEqual({row["role"] for row in occ}, {"exon", "intron"})
            exons = [row for row in occ if row["role"] == "exon"]
            self.assertEqual([(int(row["start"]), int(row["end"])) for row in exons], [(10, 40), (80, 120)])
            self.assertEqual({row["coding_status"] for row in exons}, {"coding"})
            self.assertTrue(hom)
            self.assertEqual(matches, [])
            self.assertTrue(paths)
            self.assertEqual(introns[0]["phase_compatibility"], "incompatible")
            self.assertIn("copy_subclass", copy_context[0])


class SimulationBenchmarkTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
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

    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
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

    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
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

    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
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

    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
    def test_simulation_calibration_operating_characteristics(self):
        with tempfile.TemporaryDirectory() as tmp:
            calibrate_simulations(tmp, scenarios=["negative_control"], replicates=1, seed=29)
            summary = read_tsv(Path(tmp) / "calibration_operating_characteristics.tsv")
            reps = read_tsv(Path(tmp) / "calibration_replicates.tsv")
            self.assertEqual(summary[0]["scenario"], "negative_control")
            self.assertEqual(summary[0]["replicates"], "1")
            self.assertTrue(reps)

    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
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
            self.assertIn("hidden_segment_candidate", {row["completion_call"] for row in annot})
            self.assertEqual(bench[0]["truth_events"], "0")
            self.assertEqual(bench[0]["called_events"], "0")

    @unittest.skipUnless(shutil.which("mafft"), "MAFFT external integration dependency is unavailable")
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
            copy_tree = tmp / "copy_tree.tsv"
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
            copy_tree.write_text("node_id\tparent_id\tlabel\nroot\t\troot\nspa_copy\troot\tSpA:SpA_geneA\nspb_copy\troot\tSpB:SpB_geneA\n")
            build_case(manifest, out, species_tree=tree, copy_tree=copy_tree)
            report = read_tsv(out / "case_build_report.tsv")
            prov = read_tsv(out / "case_provenance.tsv")
            occ = read_tsv(out / "segment_occurrences.tsv")
            self.assertEqual({row["status"] for row in report}, {"extracted"})
            self.assertEqual(len(prov), 2)
            self.assertTrue(occ)
            self.assertTrue((out / "species_tree.tsv").exists())
            self.assertTrue((out / "copy_tree.tsv").exists())

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
