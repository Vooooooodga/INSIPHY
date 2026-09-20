import tempfile

import unittest

from contextlib import redirect_stdout, nullcontext

from io import StringIO

from pathlib import Path

from intraphy.candidate_chain import (
    DEFAULT_CHAIN_CONFIGURATION,
    ChainCandidate,
    ChainPathMembership,
    classify_reference_coverage,
    ordered_candidate_chain,
)

from intraphy.coordinates import (
    ClosedInterval1,
    CoordinateBlock,
    Interval0,
    format_legacy_blocks,
    genome_interval_to_local,
    local_interval_to_genome,
    parse_legacy_blocks,
)

from intraphy.io import (
    LEGACY_STRUCTURAL_SITE_SCHEMA_VERSION,
    normalize_structural_site_row,
    read_structural_site_matrix,
    validate_structural_site_tip_rows,
    write_structural_site_matrix,
)

from unittest.mock import patch

from intraphy.cli import main

from intraphy.io import read_tsv

from intraphy.orthofinder import _resolve_members_to_locus, import_orthofinder

from support_interfaces_interface import InterfaceTestsSupport

class InterfaceTests(InterfaceTestsSupport, unittest.TestCase):
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
