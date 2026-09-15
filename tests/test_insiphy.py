import tempfile
import unittest
from pathlib import Path

from insiphy.cli import run_all
from insiphy.io import read_tsv
from insiphy.preprocess import derive_tables, extract_gene


ROOT = Path(__file__).resolve().parents[1]


class DemoTests(unittest.TestCase):
    def test_jingwei_demo_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_all(ROOT / "demos" / "jingwei", tmp)
            rows = read_tsv(Path(tmp) / "demo_summary.tsv")
            self.assertEqual(rows[0]["family_id"], "jingwei")
            self.assertEqual(rows[0]["hidden_segment_candidates"], "1")
            self.assertEqual(rows[0]["best_compound_model"], "compound_chimeric_or_copy_event")

    def test_sdic_demo_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_all(ROOT / "demos" / "sdic", tmp)
            rows = read_tsv(Path(tmp) / "demo_summary.tsv")
            self.assertEqual(rows[0]["family_id"], "sdic")
            self.assertEqual(rows[0]["hidden_segment_candidates"], "1")
            self.assertEqual(rows[0]["best_annotation_model"], "annotation_error")


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
            self.assertEqual(len(occ), 3)
            self.assertEqual(len(adj), 2)
            self.assertEqual({row["role"] for row in occ}, {"CDS", "intron"})
            self.assertTrue(hom)


if __name__ == "__main__":
    unittest.main()
