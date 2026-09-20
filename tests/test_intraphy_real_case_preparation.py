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

from support_intraphy_real_case_preparation import RealCasePreparationTestsSupport

class RealCasePreparationTests(RealCasePreparationTestsSupport, unittest.TestCase):
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
