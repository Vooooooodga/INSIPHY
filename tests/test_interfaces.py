import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from insiphy.cli import main
from insiphy.io import read_tsv
from insiphy.orthofinder import _resolve_members_to_locus, import_orthofinder


class InterfaceTests(unittest.TestCase):
    def _orthofinder_root(self, root, header, row):
        og = root / "orthofinder" / "Orthogroups"
        og.mkdir(parents=True)
        (og / "Orthogroups.tsv").write_text(header + "\n" + row + "\n")
        return root / "orthofinder"

    def _orthofinder_txt_root(self, root, line):
        og = root / "orthofinder" / "Orthogroups"
        og.mkdir(parents=True)
        (og / "Orthogroups.txt").write_text(line + "\n")
        return root / "orthofinder"

    def _manifest(self, root, rows):
        path = root / "genomes.tsv"
        path.write_text(
            "species\tgenome_fasta\tannotation_file\n"
            + "\n".join(f"{species}\t{root / (species + '.fa')}\t{gff}" for species, gff in rows)
            + "\n"
        )
        return path

    def test_cli_defaults_keep_all_annotated_transcripts(self):
        with patch("insiphy.cli.extract_gene") as extract:
            main(
                [
                    "extract-gene",
                    "--genome",
                    "genome.fa",
                    "--annotation",
                    "annotation.gff3",
                    "--gene-id",
                    "gene1",
                    "--family-id",
                    "fam",
                    "--species",
                    "sp",
                    "--gene-copy-id",
                    "copy1",
                    "--output-dir",
                    "out",
                ]
            )
        self.assertEqual(extract.call_args.args[8], "all")

        with patch("insiphy.cli.build_case") as build:
            main(["build-case", "--manifest", "manifest.tsv", "--output-dir", "case"])
        self.assertEqual(build.call_args.args[4], "all")

    def test_simulate_and_calibrate_help_mark_legacy_scope(self):
        for command in ["simulate", "calibrate"]:
            out = StringIO()
            with self.assertRaises(SystemExit) as ctx, redirect_stdout(out):
                main([command, "--help"])
            self.assertEqual(ctx.exception.code, 0)
            help_text = " ".join(out.getvalue().split()).replace("- ", "-")
            self.assertIn("Legacy experimental-multicopy", help_text)
            self.assertIn("single-copy", help_text)

    def test_orthofinder_headers_map_through_gff_parent_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ncbi_gff = root / "A.gff3"
            ncbi_gff.write_text(
                "chr1\tRefSeq\tgene\t1\t1000\t.\t+\t.\tID=gene-LOC410353;Name=LOC410353\n"
                "chr1\tRefSeq\tmRNA\t1\t1000\t.\t+\t.\tID=rna-XM_016915455.2;Parent=gene-LOC410353\n"
                "chr1\tRefSeq\tCDS\t10\t900\t.\t+\t0\tID=cds-XP_016770944.1;Parent=rna-XM_016915455.2;protein_id=XP_016770944.1\n"
            )
            ensembl_gff = root / "B.gff3"
            ensembl_gff.write_text(
                "LG1\tEnsembl\tgene\t1\t1000\t.\t+\t.\tID=gene:ENSBIGG00000000001;Name=ENSBIGG00000000001;gene_id=ENSBIGG00000000001\n"
                "LG1\tEnsembl\tmRNA\t1\t1000\t.\t+\t.\tID=transcript:ENSBIGT00000000002;Parent=gene:ENSBIGG00000000001;gene_id=ENSBIGG00000000001\n"
                "LG1\tEnsembl\tCDS\t20\t920\t.\t+\t0\tID=CDS:ENSBIGP00000000012;Parent=transcript:ENSBIGT00000000002;protein_id=ENSBIGP00000000012;version=1;gene_id=ENSBIGG00000000001\n"
            )
            manifest = self._manifest(root, [("A", ncbi_gff), ("B", ensembl_gff)])
            orthofinder = self._orthofinder_root(
                root,
                "Orthogroup\tA\tB\tunused_extra_species",
                "OG0001\trna-XM_016915455.2 gene=gene-LOC410353\tENSBIGP00000000012.1\textra_copy",
            )
            out = root / "out"

            import_orthofinder(orthofinder, "OG0001", manifest, out)

            rows = {row["species"]: row for row in read_tsv(out / "manifest.tsv")}
            self.assertEqual(rows["A"]["gene_id"], "gene-LOC410353")
            self.assertEqual(rows["B"]["gene_id"], "gene:ENSBIGG00000000001")
            self.assertEqual(rows["A"]["orthofinder_member_ids"], "rna-XM_016915455.2 gene=gene-LOC410353")
            mapping = read_tsv(out / "orthofinder_member_mapping.tsv")
            self.assertTrue(any(row["source_member_id"] == "ENSBIGP00000000012.1" for row in mapping))

    def test_orthogroups_txt_supports_official_pipe_colon_normalization(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bee_gff = root / "bee.gff3"
            bee_gff.write_text(
                "LG1\tEnsembl\tgene\t1\t1000\t.\t+\t.\tID=gene:ENSBIGG00000018832;gene_id=ENSBIGG00000018832\n"
                "LG1\tEnsembl\tmRNA\t1\t1000\t.\t+\t.\tID=transcript:ENSBIGT00000018832;Parent=gene:ENSBIGG00000018832\n"
                "LG1\tEnsembl\tCDS\t10\t900\t.\t+\t0\tID=CDS:ENSBIGP00000018832;Parent=transcript:ENSBIGT00000018832;protein_id=ENSBIGP00000018832;version=1\n"
            )
            tetra_gff = root / "tetra.gff3"
            tetra_gff.write_text(
                "scaf1\tmaker\tgene\t1\t800\t.\t+\t.\tID=g10177\n"
                "scaf1\tmaker\tmRNA\t1\t800\t.\t+\t.\tID=rna-gnl|WGS:JAWNGG|g10177.t1;Parent=g10177\n"
            )
            manifest = self._manifest(root, [("Bombus_ignitus", bee_gff), ("Tetragonisca_angustula", tetra_gff)])
            orthofinder = self._orthofinder_txt_root(
                root,
                "dsxOG0006023: ENSBIGP00000018832.1 rna-gnl|WGS_JAWNGG|g10177.t1",
            )
            out = root / "out"

            import_orthofinder(orthofinder, "dsxOG0006023", manifest, out)

            rows = {row["species"]: row for row in read_tsv(out / "manifest.tsv")}
            self.assertEqual(rows["Bombus_ignitus"]["gene_id"], "gene:ENSBIGG00000018832")
            self.assertEqual(rows["Tetragonisca_angustula"]["gene_id"], "g10177")
            self.assertEqual(rows["Tetragonisca_angustula"]["orthofinder_member_ids"], "rna-gnl|WGS_JAWNGG|g10177.t1")

    def test_orthogroups_txt_prefers_sequenceids_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "A.gff3"
            gff.write_text(
                "chr1\tRefSeq\tgene\t1\t500\t.\t+\t.\tID=geneA\n"
                "chr1\tRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=raw_tx_A;Parent=geneA\n"
            )
            manifest = self._manifest(root, [("A", gff)])
            orthofinder = self._orthofinder_txt_root(root, "OG0001: 0_0")
            work = orthofinder / "WorkingDirectory"
            work.mkdir()
            (work / "SequenceIDs.txt").write_text("0_0: raw_tx_A\n")
            out = root / "out"

            import_orthofinder(orthofinder, "OG0001", manifest, out)

            rows = read_tsv(out / "manifest.tsv")
            self.assertEqual(rows[0]["gene_id"], "geneA")

    def test_exact_member_match_precedes_derived_namespace_tokens(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "A.gff3"
            gff.write_text(
                "chr1\tEnsembl\tgene\t1\t500\t.\t+\t.\tID=gene:X;gene_id=X\n"
                "chr1\tEnsembl\tmRNA\t1\t500\t.\t+\t.\tID=txX;Parent=gene:X\n"
                "chr1\tEnsembl\tgene\t700\t1200\t.\t+\t.\tID=gene:Y;Name=X\n"
                "chr1\tEnsembl\tmRNA\t700\t1200\t.\t+\t.\tID=txY;Parent=gene:Y\n"
            )
            manifest = self._manifest(root, [("A", gff)])
            orthofinder = self._orthofinder_root(root, "Orthogroup\tA", "OG0001\tgene:X")
            out = root / "out"

            import_orthofinder(orthofinder, "OG0001", manifest, out)

            rows = read_tsv(out / "manifest.tsv")
            self.assertEqual(rows[0]["gene_id"], "gene:X")

    def test_dangling_parent_does_not_create_locus(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "A.gff3"
            gff.write_text(
                "chr1\tRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=tx_dangling;Parent=missing_gene\n"
            )
            manifest = self._manifest(root, [("A", gff)])
            orthofinder = self._orthofinder_root(root, "Orthogroup\tA", "OG0001\ttx_dangling")
            out = root / "out"

            with self.assertRaises(SystemExit):
                import_orthofinder(orthofinder, "OG0001", manifest, out)

            excluded = read_tsv(out / "excluded_families.tsv")
            self.assertEqual(excluded[0]["reason"], "unresolved_member_id")

    def test_gene_id_alias_and_parent_take_priority_over_gene_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "A.gff3"
            gff.write_text(
                "chr1\tEnsembl\tgene\t1\t500\t.\t+\t.\tID=gene:X;gene_id=X\n"
                "chr1\tEnsembl\tmRNA\t1\t500\t.\t+\t.\tID=txX;Parent=gene:X;gene_id=X\n"
                "chr1\tEnsembl\tmRNA\t1\t500\t.\t+\t.\tID=txAlias;gene_id=X\n"
                "chr1\tEnsembl\tgene\t700\t1200\t.\t+\t.\tID=gene:Y;Name=X;Alias=gene:X\n"
                "chr1\tEnsembl\tmRNA\t700\t1200\t.\t+\t.\tID=txDangling;Parent=missing;gene_id=X\n"
            )
            for member in ("gene:X", "X", "txX", "txAlias"):
                with self.subTest(member=member):
                    locus, status, _ = _resolve_members_to_locus([member], gff)
                    self.assertEqual((locus, status), ("gene:X", "only_gene_unique"))
            locus, status, _ = _resolve_members_to_locus(["txDangling"], gff)
            self.assertEqual((locus, status), (None, "unresolved_member_id"))

    def test_exact_transcript_id_precedes_derived_prefix_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "A.gff3"
            gff.write_text(
                "chr1\tEnsembl\tgene\t1\t500\t.\t+\t.\tID=geneA\n"
                "chr1\tEnsembl\tmRNA\t1\t500\t.\t+\t.\tID=tx1;Parent=geneA\n"
                "chr1\tEnsembl\tgene\t700\t1200\t.\t+\t.\tID=geneB\n"
                "chr1\tEnsembl\tmRNA\t700\t1200\t.\t+\t.\tID=transcript:tx1;Parent=geneB\n"
            )
            locus, status, _ = _resolve_members_to_locus(["tx1 gene=geneB"], gff)
            self.assertEqual((locus, status), ("geneA", "only_gene_unique"))
            locus, status, _ = _resolve_members_to_locus(["unknown:tx1"], gff)
            self.assertEqual((locus, status), (None, "unresolved_member_id"))

    def test_protein_version_requires_explicit_matching_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "A.gff3"
            gff.write_text(
                "chr1\tEnsembl\tgene\t1\t500\t.\t+\t.\tID=gene:X;gene_id=X;version=1\n"
                "chr1\tEnsembl\tmRNA\t1\t500\t.\t+\t.\tID=transcript:T;Parent=gene:X;transcript_id=T;version=1\n"
                "chr1\tEnsembl\tCDS\t10\t400\t.\t+\t0\tID=CDS:P;Parent=transcript:T;protein_id=P;version=1\n"
                "chr1\tEnsembl\tCDS\t10\t400\t.\t+\t0\tID=CDS:Q;Parent=transcript:T;protein_id=Q\n"
            )
            for member in ("P.1", "T.1", "unknown_protein gene:X.1"):
                with self.subTest(member=member):
                    locus, status, _ = _resolve_members_to_locus([member], gff)
                    self.assertEqual((locus, status), ("gene:X", "only_gene_unique"))
            for member in ("P.2", "Q.1", "unknown_protein seq_id=gene:X"):
                with self.subTest(member=member):
                    locus, status, _ = _resolve_members_to_locus([member], gff)
                    self.assertEqual((locus, status), (None, "unresolved_member_id"))

    def test_orthofinder_normalization_collision_requires_sequenceids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "A.gff3"
            raw = "rna-gnl|WGS:JAWNGG|g10177.t1"
            normalized = raw.replace(":", "_")
            gff.write_text(
                "chr1\tGenBank\tgene\t1\t500\t.\t+\t.\tID=geneA\n"
                f"chr1\tGenBank\tmRNA\t1\t500\t.\t+\t.\tID={raw};Parent=geneA\n"
                "chr1\tGenBank\tgene\t700\t1200\t.\t+\t.\tID=geneB\n"
                f"chr1\tGenBank\tmRNA\t700\t1200\t.\t+\t.\tID={normalized};Parent=geneB\n"
            )
            locus, status, _ = _resolve_members_to_locus([normalized], gff)
            self.assertEqual((locus, status), (None, "ambiguous_member_id"))

            manifest = self._manifest(root, [("A", gff)])
            orthofinder = self._orthofinder_root(root, "Orthogroup\tA", f"OG0001\t{normalized}")
            work = orthofinder / "WorkingDirectory"
            work.mkdir()
            sequence_ids = work / "SequenceIDs.txt"
            sequence_ids.write_text(f"0_0: {raw} gene=geneA\n")
            rows = import_orthofinder(orthofinder, "OG0001", manifest, root / "mapped")
            self.assertEqual(rows[0]["gene_id"], "geneA")

            sequence_ids.write_text(f"0_0: {raw} gene=geneA\n0_1: {normalized} gene=geneB\n")
            with self.assertRaises(SystemExit):
                import_orthofinder(orthofinder, "OG0001", manifest, root / "ambiguous")
            excluded = read_tsv(root / "ambiguous" / "excluded_families.tsv")
            self.assertEqual(excluded[0]["locus_count"], "2")

    def test_orthofinder_multiple_isoforms_same_locus_is_singlecopy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "A.gff3"
            gff.write_text(
                "chr1\tRefSeq\tgene\t1\t1000\t.\t+\t.\tID=geneA\n"
                "chr1\tRefSeq\tmRNA\t1\t1000\t.\t+\t.\tID=tx1;Parent=geneA\n"
                "chr1\tRefSeq\tmRNA\t1\t1000\t.\t+\t.\tID=tx2;Parent=geneA\n"
            )
            manifest = self._manifest(root, [("A", gff)])
            orthofinder = self._orthofinder_root(root, "Orthogroup\tA", "OG0001\trna-tx1, rna-tx2")
            out = root / "out"

            import_orthofinder(orthofinder, "OG0001", manifest, out)

            rows = read_tsv(out / "manifest.tsv")
            self.assertEqual(rows[0]["gene_id"], "geneA")
            self.assertEqual(rows[0]["orthofinder_member_count"], "2")
            self.assertEqual(rows[0]["orthofinder_mapping_status"], "only_gene_unique")

    def test_orthofinder_members_from_multiple_loci_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "A.gff3"
            gff.write_text(
                "chr1\tRefSeq\tgene\t1\t500\t.\t+\t.\tID=geneA\n"
                "chr1\tRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=tx1;Parent=geneA\n"
                "chr1\tRefSeq\tgene\t700\t1200\t.\t+\t.\tID=geneB\n"
                "chr1\tRefSeq\tmRNA\t700\t1200\t.\t+\t.\tID=tx2;Parent=geneB\n"
            )
            manifest = self._manifest(root, [("A", gff)])
            orthofinder = self._orthofinder_root(root, "Orthogroup\tA", "OG0001\trna-tx1, rna-tx2")
            out = root / "out"

            with self.assertRaises(SystemExit):
                import_orthofinder(orthofinder, "OG0001", manifest, out)

            excluded = read_tsv(out / "excluded_families.tsv")
            self.assertEqual(excluded[0]["reason"], "orthofinder_members_map_to_multiple_loci")

    def test_orthofinder_shared_protein_id_ambiguity_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "A.gff3"
            gff.write_text(
                "chr1\tRefSeq\tgene\t1\t500\t.\t+\t.\tID=geneA\n"
                "chr1\tRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=tx1;Parent=geneA\n"
                "chr1\tRefSeq\tCDS\t1\t300\t.\t+\t0\tID=cds1;Parent=tx1;protein_id=sharedP\n"
                "chr1\tRefSeq\tgene\t700\t1200\t.\t+\t.\tID=geneB\n"
                "chr1\tRefSeq\tmRNA\t700\t1200\t.\t+\t.\tID=tx2;Parent=geneB\n"
                "chr1\tRefSeq\tCDS\t700\t1000\t.\t+\t0\tID=cds2;Parent=tx2;protein_id=sharedP\n"
            )
            manifest = self._manifest(root, [("A", gff)])
            orthofinder = self._orthofinder_root(root, "Orthogroup\tA", "OG0001\tsharedP")
            out = root / "out"

            with self.assertRaises(SystemExit):
                import_orthofinder(orthofinder, "OG0001", manifest, out)

            excluded = read_tsv(out / "excluded_families.tsv")
            self.assertEqual(excluded[0]["reason"], "ambiguous_member_id")
            self.assertEqual(excluded[0]["locus_count"], "2")


if __name__ == "__main__":
    unittest.main()
