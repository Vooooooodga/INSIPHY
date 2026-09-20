import json

import tempfile

import unittest

from pathlib import Path

from unittest.mock import patch

from intraphy.correspondence import _membership_match_details, infer_correspondence, match_total, tree_distances

from intraphy.coordinates import parse_legacy_blocks

from intraphy.elements import collect_element_profiles

from intraphy.io import read_tsv, write_tsv

from intraphy.cli import main as cli_main

from intraphy.preprocess import MATCH_FIELDS, _apply_ordered_candidate_chains, _species_tree_distances, assess_short_candidate_thresholds, cheap_match_evidence, cluster_segments, copy_order_context, extract_gene, graph_components, introns_from_path, match_evidence, read_annotation_for_gene, transcript_cds_length

from intraphy.alignment import AlignmentStats, phase_compatibility

from intraphy.structural_sites import (
    _completion_presence,
    _element_site_rows,
    _genomically_contiguous,
    _junction_site_rows,
)

from support_observations_observation_semantics import ObservationSemanticsTestsSupport

class ObservationSemanticsTests(ObservationSemanticsTestsSupport, unittest.TestCase):
    def test_exon_prediction_preserves_confirmed_and_alternative_exonic_usage(self):
        for alternative_usage in (False, True):
            with self.subTest(alternative_usage=alternative_usage):
                occurrences, elements, candidate = self._intronic_cds_prediction_case()
                exonic = {
                    **occurrences[1], "occurrence_id": "Bter_exon", "role": "CDS",
                    "start": "2476398", "end": "2476441",
                }
                exonic_element = {
                    **elements[1], "occurrence_id": "Bter_exon", "element_class": "exon_like",
                }
                if alternative_usage:
                    occurrences.append(exonic)
                    elements.append(exonic_element)
                else:
                    occurrences[1] = exonic
                    elements[1] = exonic_element
                rows = _element_site_rows(
                    occurrences, elements, [candidate], {"OG0006454"},
                    {"OG0006454": {"Bombus_ignitus", "Bombus_terrestris"}},
                )
                by_layer = {row["layer"]: row for row in rows if row["species"] == "Bombus_terrestris"}
                self.assertEqual(by_layer["exon_presence"]["state"], "present")
                self.assertEqual(by_layer["exon_role"]["state"], "exonic")
                self.assertNotIn("role_conflict_unresolved", by_layer["exon_role"]["conclusion_flag"])
                if alternative_usage:
                    self.assertIn("alternative_usage", by_layer["exon_role"]["evidence"])

    def test_exon_prediction_preserves_contradictory_sequence_absence_as_unknown(self):
        occurrences, elements, candidate = self._intronic_cds_prediction_case()
        absence = {
            **candidate,
            "completion_call": "supports_true_absence",
            "homologous_dna_presence": "absent",
            "homologous_dna_evidence": "ordered_flank_deletion",
            "absence_evidence": "ordered_flank_deletion",
        }
        rows = _element_site_rows(
            occurrences, elements, [candidate, absence], {"OG0006454"},
            {"OG0006454": {"Bombus_ignitus", "Bombus_terrestris"}},
        )
        by_layer = {row["layer"]: row for row in rows if row["species"] == "Bombus_terrestris"}
        self.assertEqual(by_layer["exon_presence"]["state"], "unknown")
        self.assertIn("conflicting_presence_states", by_layer["exon_presence"]["evidence"])
        self.assertEqual(by_layer["exon_role"]["state"], "unknown")
        self.assertNotIn("role_conflict_unresolved", by_layer["exon_role"]["conclusion_flag"])

    def test_transcript_cds_length_uses_union_of_cds_metadata(self):
        tx_features = [
            {"type": "exon", "start": 1, "end": 500, "cds_intervals": [(10, 90), (80, 120)]},
            {"type": "exon", "start": 600, "end": 700, "cds_intervals": [(620, 650)]},
        ]
        self.assertEqual(transcript_cds_length(tx_features), 111 + 31)

    def test_completion_with_explicit_unresolved_correspondence_does_not_assign_eg_presence(self):
        occurrences = [
            {"occurrence_id": "A_ex", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present"},
        ]
        elements = [
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "A_ex", "element_class": "exon_like", "membership_call": "core_member"},
        ]
        candidate = {
            "homology_id": "H1", "family_id": "fam", "species": "B", "gene_copy_id": "gB",
            "completion_call": "homologous_sequence_candidate",
            "evidence_status": "homologous_sequence_candidate", "inferred_role": "unknown",
        }
        for status in ("unknown", "ambiguous", "resolved", "", "NA", None):
            with self.subTest(correspondence_status=status):
                completion = dict(candidate)
                if status is not None:
                    completion["correspondence_status"] = status
                unresolved = status != "resolved"
                if unresolved:
                    completion["primary_mapping_status"] = "ambiguous_repeated_mapping"
                    completion["interval_candidates"] = '[{"target_start":10,"target_end":30},{"target_start":50,"target_end":70}]'
                rows = _element_site_rows(occurrences, elements, [completion], {"fam"}, {"fam": {"A", "B"}})
                by_key = {(row["layer"], row["species"]): row for row in rows}
                presence = by_key[("exon_presence", "B")]
                role = by_key[("exon_role", "B")]
                self.assertEqual(presence["state"], "unknown")
                self.assertEqual(role["state"], "unknown")
                self.assertEqual(by_key[("exon_presence", "A")]["state"], "present")
                if unresolved:
                    self.assertIn("completion_correspondence_unresolved", presence["evidence"])
                    self.assertIn("completion_correspondence_unresolved", role["evidence"])
                    self.assertEqual(presence["confidence_flag"], "low")
                else:
                    self.assertNotIn("completion_correspondence_unresolved", presence["evidence"])

    def test_candidate_label_without_resolved_presence_evidence_stays_unknown(self):
        occurrences = [{
            "occurrence_id": "A_ex", "family_id": "fam", "species": "A",
            "gene_copy_id": "gA", "role": "CDS", "presence_status": "present",
        }]
        elements = [{
            "element_id": "EG_1", "homology_id": "H1", "occurrence_id": "A_ex",
            "element_class": "exon_like", "membership_call": "core_member",
        }]
        candidate = {
            "homology_id": "H1", "family_id": "fam", "species": "B",
            "gene_copy_id": "gB", "completion_call": "hidden_segment_candidate",
            "evidence_status": "supports_hidden_segment",
            "candidate_resolution_status": "resolved",
            "correspondence_status": "resolved",
        }
        rows = _element_site_rows(
            occurrences, elements, [candidate], {"fam"}, {"fam": {"A", "B"}},
        )
        by_layer = {
            row["layer"]: row for row in rows if row["species"] == "B"
        }
        self.assertEqual(by_layer["exon_presence"]["state"], "unknown")
        self.assertEqual(by_layer["exon_role"]["state"], "unknown")

    def test_legacy_true_absence_requires_resolved_ordered_flank_evidence(self):
        accepted = {
            "completion_call": "supports_true_absence",
            "candidate_resolution_status": "resolved",
            "absence_evidence": "ordered_flank_deletion",
        }
        self.assertEqual(_completion_presence(accepted)[0], "absent")
        self.assertIsNone(
            _completion_presence({
                **accepted,
                "candidate_resolution_status": "ambiguous",
            })[0]
        )
        self.assertIsNone(
            _completion_presence({
                **accepted,
                "correspondence_status": "unknown",
            })[0]
        )
        self.assertIsNone(
            _completion_presence({**accepted, "absence_evidence": "no_sequence_hit"})[0]
        )

    def test_transcript_paths_and_raw_features_preserve_each_cds_phase(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff, out = root / "genome.fa", root / "gene.gff3", root / "out"
            fasta.write_text(">chr1\n" + "A" * 100 + "\n")
            gff.write_text(
                "chr1\ttest\tgene\t10\t90\t.\t+\t.\tID=g1\n"
                "chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=tx1;Parent=g1\n"
                "chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=tx2;Parent=g1\n"
                "chr1\ttest\tncRNA\t10\t90\t.\t+\t.\tID=tx3;Parent=g1\n"
                "chr1\ttest\texon\t10\t39\t.\t+\t.\tID=ex;Parent=tx1,tx2\n"
                "chr1\ttest\texon\t10\t39\t.\t+\t.\tID=ex_nc;Parent=tx3\n"
                "chr1\ttest\tCDS\t10\t39\t.\t+\t0\tID=c1;Parent=tx1\n"
                "chr1\ttest\tCDS\t13\t39\t.\t+\t1\tID=c2;Parent=tx2\n"
            )
            extract_gene(fasta, gff, "g1", "fam", "A", "gA", out, flank=0, max_extension=0)
            paths = {row["transcript_id"]: row for row in read_tsv(out / "transcript_paths.tsv")}
            self.assertEqual((paths["tx1"]["cds_intervals"], paths["tx1"]["cds_phase"]), ("10-39", "0"))
            self.assertEqual((paths["tx2"]["cds_intervals"], paths["tx2"]["cds_phase"]), ("13-39", "1"))
            self.assertEqual(paths["tx1"]["cds_length"], "30")
            self.assertEqual(paths["tx2"]["cds_length"], "27")
            self.assertEqual(paths["tx3"]["cds_length"], "0")
            self.assertEqual(paths["tx1"]["path_role"], "exonic")
            self.assertEqual(paths["tx1"]["coding_role"], "CDS")
            self.assertEqual(paths["tx1"]["position_role"], "single")
            self.assertEqual(paths["tx1"]["annotation_source"], "test")
            self.assertIn("ID=ex", paths["tx1"]["original_attributes"])
            # v0.17: the gene extent does not establish transcript incompleteness.
            self.assertEqual((paths["tx1"]["partial_start"], paths["tx1"]["partial_end"]), ("0", "0"))
            self.assertEqual(paths["tx3"]["coding_role"], "noncoding_exon")
            raw = {row["id"]: row for row in read_tsv(out / "raw_gene_features.tsv")}
            self.assertEqual((raw["c1"]["parent"], raw["c1"]["phase"]), ("tx1", "0"))
            self.assertEqual((raw["c2"]["parent"], raw["c2"]["phase"]), ("tx2", "1"))

    def test_negative_strand_context_uses_each_transcript_without_first_path_rank(self):
        common = {"family_id": "fam", "species": "A", "gene_copy_id": "gA", "contig": "chr1", "strand": "-"}
        occurrences = [
            {**common, "occurrence_id": "shared", "start": "300", "end": "330", "role": "CDS", "transcript_id": "tx1;tx2", "transcript_order": "2"},
            {**common, "occurrence_id": "low1", "start": "100", "end": "130", "role": "exon", "transcript_id": "tx1", "transcript_order": "1"},
            {**common, "occurrence_id": "low2", "start": "200", "end": "230", "role": "UTR", "transcript_id": "tx2", "transcript_order": "1"},
        ]
        context = copy_order_context(occurrences)
        self.assertEqual(context["shared"]["left_role"], "terminal")
        self.assertEqual(context["shared"]["right_role"], "alternative")
        self.assertEqual(context["shared"]["scaled_index"], 0.0)
        self.assertEqual(context["shared"]["copy_size"], 2)
        self.assertEqual(context["low1"]["left_role"], "CDS")
        self.assertEqual(context["low2"]["left_role"], "CDS")

    def test_direct_junction_absence_requires_known_contig_strand_and_order(self):
        left = {"contig": "chr1", "strand": "+", "start": "1", "end": "10"}
        right = {"contig": "chr1", "strand": "+", "start": "11", "end": "20"}
        self.assertTrue(_genomically_contiguous(left, right))
        self.assertFalse(_genomically_contiguous(right, left))
        self.assertTrue(_genomically_contiguous({**right, "strand": "-"}, {**left, "strand": "-"}))
        self.assertFalse(_genomically_contiguous({**left, "contig": "NA"}, {**right, "contig": "NA"}))
        self.assertFalse(_genomically_contiguous({**left, "strand": "."}, {**right, "strand": "."}))

    def test_match_score_separates_sequence_and_structural_context(self):
        left = {"occurrence_id": "left", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "phase": "0", "start": "1", "end": "9", "splice_motif_score": "0.5"}
        right = {"occurrence_id": "right", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "role": "CDS", "phase": "1", "start": "1", "end": "9", "splice_motif_score": "0.5"}
        evidence = cheap_match_evidence(left, right, {"left": {}, "right": {}})

        self.assertEqual(evidence["sequence_score"], 0.0)
        self.assertEqual(evidence["boundary_score"], 1.0)
        self.assertEqual(evidence["phase_score"], 0.1)
        self.assertAlmostEqual(match_total(evidence), evidence["total_score"])

    def test_between_element_gap_without_intron_record_stays_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "transcript_paths.tsv").write_text(
                "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\tcontig\tstart\tend\tstrand\tphase\ttranscript_order\tcoding_status\tcds_phase\tpath_status\n"
                "p1\tfam\tA\tgA\tt1\t1\to1\tCDS\tchr\t1\t10\t+\t0\t1\tcoding\t0\tannotated_transcript_path\n"
                "p2\tfam\tA\tgA\tt1\t2\to2\tCDS\tchr\t21\t30\t+\t0\t2\tcoding\t0\tannotated_transcript_path\n"
            )
            occurrences = [
                {"occurrence_id": "o1", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present", "contig": "chr", "strand": "+", "start": "1", "end": "10"},
                {"occurrence_id": "o2", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present", "contig": "chr", "strand": "+", "start": "21", "end": "30"},
            ]
            element_rows = [
                {"element_id": "EG_1", "occurrence_id": "o1", "element_class": "exon_like", "membership_call": "core_member"},
                {"element_id": "EG_2", "occurrence_id": "o2", "element_class": "exon_like", "membership_call": "core_member"},
            ]

            rows = _junction_site_rows(input_dir, output_dir, occurrences, element_rows, {"fam"}, {"fam": {"A"}})
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["state"], "unknown")
            self.assertIn("missing_intron_record_with_genomic_gap", rows[0]["evidence"])

    def test_within_element_repertoire_presence_overrides_continuous_unsplit_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "transcript_paths.tsv").write_text(
                "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\tcontig\tstart\tend\tstrand\tphase\ttranscript_order\tcoding_status\tcds_phase\tpath_status\n"
                "p1\tfam\tA\tgA\tsplit\t1\tleft\tCDS\tchr\t1\t40\t+\t0\t1\tcoding\t0\tannotated_transcript_path\n"
                "p2\tfam\tA\tgA\tsplit\t2\tint1\tintron\tchr\t41\t50\t+\t.\t2\tnot_applicable\t.\tannotated_transcript_path\n"
                "p3\tfam\tA\tgA\tsplit\t3\tright\tCDS\tchr\t51\t110\t+\t0\t3\tcoding\t0\tannotated_transcript_path\n"
                "p4\tfam\tA\tgA\tcontinuous\t1\tref\tCDS\tchr\t1\t110\t+\t0\t1\tcoding\t0\tannotated_transcript_path\n"
                "p5\tfam\tB\tgB\tpartial\t1\tpartial\tCDS\tchr\t1\t40\t+\t0\t1\tcoding\t0\tannotated_transcript_path\n"
                "p6\tfam\tC\tgC\tcontinuous\t1\tcontinuous\tCDS\tchr\t1\t105\t+\t0\t1\tcoding\t0\tannotated_transcript_path\n"
            )
            (output_dir / "segment_matches.tsv").write_text(
                "match_id\tquery_occurrence_id\tsubject_occurrence_id\tmatch_status\tprojected_reference_blocks\talignment_strand\n"
                "m1\tleft\tref\tmapped\t1-40:1-40\t+\n"
                "m2\tright\tref\tmapped\t1-60:51-110\t+\n"
                "m3\tpartial\tref\tmapped\t1-40:1-40\t+\n"
                "m4\tcontinuous\tref\tmapped\t1-105:1-105\t+\n"
            )
            occurrences = [
                {"occurrence_id": "left", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present", "contig": "chr", "strand": "+", "start": "1", "end": "40"},
                {"occurrence_id": "int1", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "intron", "presence_status": "present", "contig": "chr", "strand": "+", "start": "41", "end": "50"},
                {"occurrence_id": "right", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present", "contig": "chr", "strand": "+", "start": "51", "end": "110"},
                {"occurrence_id": "ref", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present", "contig": "chr", "strand": "+", "start": "1", "end": "110"},
                {"occurrence_id": "partial", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "role": "CDS", "presence_status": "present", "contig": "chr", "strand": "+", "start": "1", "end": "40"},
                {"occurrence_id": "continuous", "family_id": "fam", "species": "C", "gene_copy_id": "gC", "role": "CDS", "presence_status": "present", "contig": "chr", "strand": "+", "start": "1", "end": "105"},
            ]
            element_rows = [
                {"element_id": "EG_1", "occurrence_id": "left", "element_class": "exon_like", "membership_call": "core_member"},
                {"element_id": "EG_1", "occurrence_id": "right", "element_class": "exon_like", "membership_call": "core_member"},
                {"element_id": "EG_1", "occurrence_id": "ref", "element_class": "exon_like", "membership_call": "core_member"},
                {"element_id": "EG_1", "occurrence_id": "partial", "element_class": "exon_like", "membership_call": "core_member"},
                {"element_id": "EG_1", "occurrence_id": "continuous", "element_class": "exon_like", "membership_call": "core_member"},
            ]

            rows = _junction_site_rows(input_dir, output_dir, occurrences, element_rows, {"fam"}, {"fam": {"A", "B", "C"}})
            by_species = {row["species"]: row for row in rows}
            self.assertEqual(by_species["A"]["state"], "present")
            self.assertEqual(by_species["A"]["conclusion_flag"], "repertoire_present")
            self.assertIn("alternative_transcript_without_this_junction", by_species["A"]["evidence"])
            self.assertIn("reference_unsplit_exon_spans_boundary", by_species["A"]["evidence"])
            self.assertEqual(by_species["B"]["state"], "unknown")
            self.assertEqual(by_species["C"]["state"], "absent")
            self.assertIn("single_transcript_continuous_alignment_block_spans_both_boundary_anchors", by_species["C"]["evidence"])

    def test_correspondence_writes_observed_not_ancestral_summaries(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            input_dir = root / "input"
            output_dir = root / "output"
            input_dir.mkdir()
            output_dir.mkdir()
            (input_dir / "segment_homology.tsv").write_text(
                "homology_id\toccurrence_id\tsupport_type\tconfidence\tsource_label\n"
                "H1\to1\tfixture\t0.95\tknown\n"
            )
            (input_dir / "segment_occurrences.tsv").write_text(
                "occurrence_id\tfamily_id\tspecies\tgene_copy_id\trole\tpresence_status\n"
                "o1\tfam\tA\tgA\tCDS\tpresent\n"
            )
            (input_dir / "segment_matches.tsv").write_text(
                "match_id\tquery_occurrence_id\tsubject_occurrence_id\tmatch_status\n"
            )
            (input_dir / "segment_sequences.fasta").write_text(">o1\nATG\n")

            infer_correspondence(input_dir, output_dir)
            self.assertTrue((output_dir / "observed_element_tree_coverage.tsv").exists())
            self.assertTrue((output_dir / "observed_intragenic_paths.tsv").exists())
            self.assertFalse((output_dir / "ancestral_element_graph.tsv").exists())
            self.assertFalse((output_dir / "ancestral_intragenic_paths.tsv").exists())
            self.assertTrue(read_tsv(output_dir / "element_correspondence.tsv"))

    def test_conservation_summary_describes_observed_species_not_fragment_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output_dir = root / "output"
            output_dir.mkdir()
            (root / "segment_homology.tsv").write_text(
                "homology_id\toccurrence_id\tsupport_type\tconfidence\n"
                "H1\ta\tfixture\t0.95\nH1\tb\tfixture\t0.95\n"
                "H2\tc\tfixture\t0.95\nH2\td\tfixture\t0.95\n"
                "H3\te\tfixture\t0.95\nH3\tf\tfixture\t0.95\n"
                "H4\tg\tfixture\t0.95\n"
            )
            (root / "segment_occurrences.tsv").write_text(
                "occurrence_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\trole\tpresence_status\n"
                "a\tfam\tA\tgA\ttx1\tCDS\tpresent\n"
                "b\tfam\tA\tgA\ttx2\tCDS\tpresent\n"
                "c\tfam\tA\tgA\ttx1\tCDS\tpresent\n"
                "d\tfam\tB\tgB\ttx1\tCDS\tpresent\n"
                "e\tfam\tA\tgA\ttx1\tCDS\tpresent\n"
                "f\tfam\tB\tgB\ttx1\tCDS\tabsent\n"
                "g\tfam\tNA\tunknown_copy\ttx1\tCDS\tpresent\n"
            )
            (root / "segment_sequences.fasta").write_text("".join(f">{name}\nATG\n" for name in "abcdefg"))
            with patch("intraphy.correspondence.simple_identity", side_effect=[0.95, 0.25, 0.95]) as identity:
                infer_correspondence(root, output_dir)
            self.assertEqual(identity.call_count, 3)
            rows = {row["homology_id"]: row for row in read_tsv(output_dir / "segment_conservation.tsv")}
            self.assertEqual(rows["H1"]["occurrence_count"], "2")
            self.assertEqual(rows["H1"]["species_count"], "1")
            self.assertEqual(rows["H1"]["conservation_call"], "observed_single_species")
            self.assertEqual(rows["H1"]["mean_pairwise_identity"], "0.95")
            self.assertEqual(rows["H2"]["species_count"], "2")
            self.assertEqual(rows["H2"]["conservation_call"], "observed_in_multiple_species")
            self.assertEqual(rows["H2"]["mean_pairwise_identity"], "0.25")
            self.assertEqual(rows["H3"]["species_count"], "1")
            self.assertEqual(rows["H3"]["conservation_call"], "observed_single_species")
            self.assertEqual(rows["H4"]["species_count"], "0")
            self.assertEqual(rows["H4"]["conservation_call"], "no_species_presence_observed")

    def test_graph_components_allows_ordered_many_fragment_split_but_keeps_repeats_ambiguous(self):
        occurrence_by_id = {
            "ref": {"family_id": "fam", "species": "A", "gene_copy_id": "gA"},
            "frag1": {"family_id": "fam", "species": "B", "gene_copy_id": "gB"},
            "frag2": {"family_id": "fam", "species": "B", "gene_copy_id": "gB"},
            "frag3": {"family_id": "fam", "species": "B", "gene_copy_id": "gB"},
            "repeat": {"family_id": "fam", "species": "B", "gene_copy_id": "gB"},
        }
        edges = [("ref", "frag1", 0.9), ("ref", "frag2", 0.9), ("ref", "frag3", 0.9), ("ref", "repeat", 0.9)]
        complementary = {
            frozenset(("frag1", "frag2")),
            frozenset(("frag1", "frag3")),
            frozenset(("frag2", "frag3")),
        }
        components = graph_components(
            list(occurrence_by_id),
            edges,
            occurrence_by_id=occurrence_by_id,
            same_copy_compatible=complementary,
            cross_copy_compatible=set(),
        )
        component_sets = {frozenset(component) for component in components}

        self.assertIn(frozenset(("ref", "frag1", "frag2", "frag3")), component_sets)
        self.assertIn(frozenset(("repeat",)), component_sets)
