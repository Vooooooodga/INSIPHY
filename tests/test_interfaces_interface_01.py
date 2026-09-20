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
    def test_coordinate_contract_round_trips_on_both_strands(self):
        public = ClosedInterval1(101, 150)
        locus = public.to_interval0()
        self.assertEqual(locus, Interval0(100, 150))
        self.assertEqual(ClosedInterval1.from_interval0(locus), public)

        local = Interval0(5, 12)
        for strand, expected in (("+", Interval0(105, 112)), ("-", Interval0(138, 145))):
            with self.subTest(strand=strand):
                genome = local_interval_to_genome(local, locus, strand)
                self.assertEqual(genome, expected)
                self.assertEqual(genome_interval_to_local(genome, locus, strand), local)

    def test_legacy_alignment_blocks_have_one_boundary_conversion(self):
        text = "1-5:11-15;8-10:20-22"
        blocks = parse_legacy_blocks(text)
        self.assertEqual(
            blocks,
            (
                CoordinateBlock(Interval0(0, 5), Interval0(10, 15)),
                CoordinateBlock(Interval0(7, 10), Interval0(19, 22)),
            ),
        )
        self.assertEqual(format_legacy_blocks(blocks), text)

    def test_structural_site_compatibility_keeps_unknown_unobserved(self):
        row = normalize_structural_site_row(
            {
                "family_id": "fam",
                "layer": "exon_role",
                "site_id": "site",
                "species": "A",
                "state": "unknown",
                "state_0": "not_exonic",
                "state_1": "exonic",
                "evidence": "not_covered",
            }
        )
        self.assertEqual(row["schema_version"], "3")
        self.assertEqual(row["applicability"], "undetermined")
        self.assertEqual(row["observation_mask"], "missing")
        self.assertEqual(row["transcript_scope"], "annotated_transcript_repertoire")

    def test_legacy_matrix_is_not_relabelled_as_schema_v2(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "legacy.tsv"
            path.write_text(
                "family_id\tlayer\tsite_id\tspecies\tstate\tstate_0\tstate_1\tevidence\n"
                "fam\texon_presence\tE1\tA\tpresent\tabsent\tpresent\tlegacy\n"
            )
            rows = read_structural_site_matrix(path)
            self.assertEqual(
                rows[0]["schema_version"], LEGACY_STRUCTURAL_SITE_SCHEMA_VERSION
            )

    def test_matrix_rejects_duplicate_observation_keys(self):
        row = {
            "family_id": "fam",
            "layer": "exon_presence",
            "site_id": "E1",
            "species": "A",
            "state": "present",
            "state_0": "absent",
            "state_1": "present",
            "evidence": "fixture",
        }
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "duplicate structural observation"):
                write_structural_site_matrix(
                    Path(tmp) / "matrix.tsv", [row, dict(row)]
                )

    def test_tip_contract_rejects_duplicate_rows_even_without_writer(self):
        row = {
            "family_id": "fam", "layer": "exon_presence", "site_id": "E1",
            "species": "A", "state": "present", "state_0": "absent",
            "state_1": "present",
        }
        with self.assertRaisesRegex(SystemExit, "duplicate=Ax2"):
            validate_structural_site_tip_rows([row, dict(row)], {"A"})

    def test_ordered_candidate_chain_retains_near_optimal_coordinate_paths(self):
        candidates = [
            ChainCandidate("left", Interval0(0, 5), Interval0(10, 15), 10.0, "nt"),
            ChainCandidate("middle_best", Interval0(5, 10), Interval0(15, 20), 8.0, "nt"),
            ChainCandidate("middle_near", Interval0(5, 10), Interval0(16, 21), 7.5, "nt"),
            ChainCandidate("reversed", Interval0(5, 10), Interval0(2, 7), 20.0, "nt"),
            ChainCandidate("right", Interval0(10, 15), Interval0(22, 27), 10.0, "nt"),
        ]
        result = ordered_candidate_chain(
            candidates,
            score_delta=0.5,
            start_ids={"left"},
            end_ids={"right"},
        )
        self.assertEqual(result.best_score, 28.0)
        self.assertEqual(
            result.best_path_member_ids,
            frozenset({"left", "middle_best", "right"}),
        )
        self.assertEqual(
            result.retained_ids,
            frozenset({"left", "middle_best", "middle_near", "right"}),
        )

    def test_reference_coverage_separates_split_and_repeat_patterns(self):
        parent = Interval0(0, 20)
        split = classify_reference_coverage([Interval0(0, 8), Interval0(8, 20)], parent)
        repeat = classify_reference_coverage([Interval0(0, 12), Interval0(6, 18)], parent)
        partial = classify_reference_coverage([Interval0(0, 8), Interval0(12, 20)], parent)
        self.assertEqual(split.relation, "complementary_complete")
        self.assertEqual(repeat.relation, "repeated_overlap")
        self.assertGreater(repeat.overlap_bases, 0)
        self.assertEqual(partial.relation, "complementary_partial")
        self.assertEqual(partial.uncovered_bases, 4)

    def test_fixed_anchor_chain_requires_one_compatible_transcript_path(self):
        def path(query_path, target_path, order):
            return ChainPathMembership(
                query_path,
                target_path,
                order,
                order,
                "chrQ",
                "chrT",
                "+",
                "+",
            )

        candidates = [
            ChainCandidate(
                "start", Interval0(0, 5), Interval0(0, 5), 2.0, "nt_blastn_v1",
                path_memberships=(path("q_tx1", "t_tx1", 1),),
            ),
            ChainCandidate(
                "compatible", Interval0(5, 10), Interval0(5, 10), 4.0, "nt_blastn_v1",
                path_memberships=(path("q_tx1", "t_tx1", 2),),
            ),
            ChainCandidate(
                "other_path", Interval0(5, 10), Interval0(5, 10), 40.0, "nt_blastn_v1",
                path_memberships=(path("q_tx2", "t_tx2", 2),),
            ),
            ChainCandidate(
                "end", Interval0(10, 15), Interval0(10, 15), 2.0, "nt_blastn_v1",
                path_memberships=(path("q_tx1", "t_tx1", 3),),
            ),
        ]
        delta = DEFAULT_CHAIN_CONFIGURATION.score_delta(8.0, "nt_blastn_v1")
        result = ordered_candidate_chain(
            candidates,
            delta,
            start_ids={"start"},
            end_ids={"end"},
            configuration_name=DEFAULT_CHAIN_CONFIGURATION.name,
        )

        self.assertEqual(result.best_score, 8.0)
        self.assertEqual(result.retained_ids, frozenset({"start", "compatible", "end"}))
        self.assertNotIn("other_path", result.retained_ids)
        self.assertFalse(result.local_mode)
        self.assertEqual(result.start_anchor_ids, frozenset({"start"}))
        self.assertEqual(result.end_anchor_ids, frozenset({"end"}))
        self.assertEqual(result.configuration_name, DEFAULT_CHAIN_CONFIGURATION.name)

    def test_arbitrary_complementary_blocks_and_exact_duplicates_use_actual_intervals(self):
        parent = Interval0(0, 30)
        one_to_three = classify_reference_coverage(
            [Interval0(0, 7), Interval0(7, 19), Interval0(19, 30)],
            parent,
        )
        exact_repeat = classify_reference_coverage(
            [Interval0(0, 12), Interval0(0, 12)],
            parent,
        )

        self.assertEqual(one_to_three.relation, "complementary_complete")
        self.assertEqual(one_to_three.covered_bases, 30)
        self.assertEqual(exact_repeat.relation, "repeated_overlap")
        self.assertEqual(exact_repeat.covered_bases, 12)
        self.assertEqual(exact_repeat.overlap_bases, 12)

    @patch("intraphy.cli.command_session", lambda args: nullcontext())
    def test_cli_defaults_keep_all_annotated_transcripts(self):
        with patch("intraphy.cli.extract_gene") as extract:
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

        with patch("intraphy.cli.build_case") as build:
            main(["build-case", "--manifest", "manifest.tsv", "--output-dir", "case"])
        self.assertEqual(build.call_args.args[4], "all")

    @patch("intraphy.cli.command_session", lambda args: nullcontext())
    def test_cli_passes_frozen_matrix_and_annotation_view(self):
        with tempfile.TemporaryDirectory() as tmp, patch(
            "intraphy.cli.infer_phylogeny"
        ) as infer:
            main([
                "infer-phylogeny",
                "--input-dir", tmp,
                "--output-dir", str(Path(tmp) / "out"),
                "--structural-site-matrix", "matrix.tsv",
                "--annotation-view", "canonical",
            ])
        self.assertEqual(infer.call_args.kwargs["structural_site_matrix_path"], "matrix.tsv")
        self.assertEqual(infer.call_args.kwargs["annotation_view"], "canonical")

    @patch("intraphy.cli.command_session", lambda args: nullcontext())
    def test_run_defaults_to_repertoire_annotation_view(self):
        with tempfile.TemporaryDirectory() as tmp, patch("intraphy.cli.run_all") as run:
            main([
                "run",
                "--input-dir", tmp,
                "--output-dir", str(Path(tmp) / "out"),
                "--structural-site-matrix", "matrix.tsv",
            ])
        self.assertEqual(run.call_args.kwargs["annotation_view"], "repertoire")
        self.assertEqual(
            run.call_args.kwargs["structural_site_matrix_path"], "matrix.tsv"
        )

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
