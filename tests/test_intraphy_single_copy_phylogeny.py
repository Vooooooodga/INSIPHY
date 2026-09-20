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

from support_intraphy_single_copy_phylogeny import SingleCopyPhylogenyTestsSupport

class SingleCopyPhylogenyTests(SingleCopyPhylogenyTestsSupport, unittest.TestCase):
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
