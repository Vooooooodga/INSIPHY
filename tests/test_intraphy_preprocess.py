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

from support_intraphy_preprocess import PreprocessTestsSupport

class PreprocessTests(PreprocessTestsSupport, unittest.TestCase):
    def test_alignment_coverage_counts_paired_bases(self):
        stats = global_alignment_stats("ACGT", "ACGTTTTT")
        self.assertAlmostEqual(stats.query_coverage, 1.0)
        self.assertAlmostEqual(stats.target_coverage, 0.5)
        self.assertAlmostEqual(stats.coverage, 0.5)
        self.assertEqual(stats.aligned_pairs, 4)

    def test_auto_global_alignment_requires_mafft_for_short_exons(self):
        sequence = "ACGT" * 25
        with patch("intraphy.alignment.shutil.which", return_value="/mock/mafft") as available, patch(
            "intraphy.alignment.subprocess.run"
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
        with patch("intraphy.alignment.shutil.which", return_value=None), patch(
            "intraphy.alignment.subprocess.run"
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
        with patch("intraphy.correspondence.global_alignment_stats", return_value=expected) as align:
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
                        "chr1\tIntraPhy\tgene\t10\t120\t.\t+\t.\tID=geneA;Name=GeneA",
                        "chr1\tIntraPhy\tmRNA\t10\t120\t.\t+\t.\tID=txA;Parent=geneA",
                        "chr1\tIntraPhy\texon\t10\t40\t.\t+\t.\tID=ex1;Parent=txA",
                        "chr1\tIntraPhy\tCDS\t15\t40\t.\t+\t0\tID=cds1;Parent=txA",
                        "chr1\tIntraPhy\texon\t80\t120\t.\t+\t.\tID=ex2;Parent=txA",
                        "chr1\tIntraPhy\tCDS\t80\t110\t.\t+\t2\tID=cds2;Parent=txA",
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
