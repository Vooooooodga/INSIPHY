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
    def test_protein_rescued_short_exon_does_not_enter_dna_sensitivity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            segment_matches = root / "segment_matches.tsv"
            output = root / "candidate_sensitivity.tsv"
            protein_candidate = {
                "candidate_id": "match_00001.candidate_001",
                "source": "protein_msa_projection",
                "score_scheme": "blosum62_cds_projection",
                "identity": 0.99,
                "query_coverage": 1.0,
                "aligned_pairs": 20,
                "acceptance_threshold": 0.7,
            }
            dna_candidate = {
                "candidate_id": "match_00001.dna_candidate_001",
                "source": "nucleotide_alignment",
                "score_scheme": "nt_blastn_v1",
                "identity": 0.72,
                "coverage": 0.90,
                "query_coverage": 0.85,
                "aligned_pairs": 20,
                "accepted": 1,
                "acceptance_threshold": 0.70,
            }
            segment_matches.write_text(
                "match_id\tquery_occurrence_id\tsubject_occurrence_id\tshort_context_route\t"
                "candidate_assessments\tdna_candidate_assessments\n"
                "match_00001\tq\ts\tbounded_local\t"
                f"{json.dumps([protein_candidate], separators=(',', ':'))}\t"
                f"{json.dumps([dna_candidate], separators=(',', ':'))}\n"
            )

            rows = assess_short_candidate_thresholds(
                segment_matches,
                output,
                short_identity_thresholds=(0.70,),
                short_query_coverage_thresholds=(0.80,),
            )

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["candidate_id"], "match_00001.dna_candidate_001")
            self.assertEqual(rows[0]["acceptance"], 1)
            self.assertNotIn("blosum62", output.read_text())

    def test_short_dna_sensitivity_retains_saved_pair_threshold(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            segment_matches = root / "segment_matches.tsv"
            output = root / "candidate_sensitivity.tsv"

            def candidate(candidate_id, threshold):
                return {
                    "candidate_id": candidate_id,
                    "source": "nucleotide_alignment",
                    "score_scheme": "nt_blastn_v1",
                    "identity": 0.75,
                    "coverage": 0.90,
                    "query_coverage": 0.90,
                    "aligned_pairs": 25,
                    "accepted": int(threshold <= 0.75),
                    "acceptance_threshold": threshold,
                }

            segment_matches.write_text(
                "match_id\tquery_occurrence_id\tsubject_occurrence_id\tshort_context_route\t"
                "dna_candidate_assessments\n"
                f"m1\tq1\ts1\tbounded_local\t{json.dumps([candidate('m1.dna_candidate_001', 0.70)], separators=(',', ':'))}\n"
                f"m2\tq2\ts2\tbounded_local\t{json.dumps([candidate('m2.dna_candidate_001', 0.80)], separators=(',', ':'))}\n"
            )

            rows = assess_short_candidate_thresholds(
                segment_matches,
                output,
                short_identity_thresholds=(0.60,),
                short_query_coverage_thresholds=(0.80,),
            )
            by_match = {row["match_id"]: row for row in rows}

            self.assertEqual(by_match["m1"]["saved_pair_threshold"], "0.7")
            self.assertEqual(by_match["m1"]["effective_identity_cutoff"], "0.7")
            self.assertEqual(by_match["m1"]["acceptance"], 1)
            self.assertEqual(by_match["m2"]["saved_pair_threshold"], "0.8")
            self.assertEqual(by_match["m2"]["effective_identity_cutoff"], "0.8")
            self.assertEqual(by_match["m2"]["acceptance"], 0)

    def test_candidate_sensitivity_cli_keeps_public_flags(self):
        with patch("intraphy.cli.assess_short_candidate_thresholds") as assess:
            cli_main([
                "candidate-sensitivity",
                "--segment-matches", "segment_matches.tsv",
                "--output", "candidate_sensitivity.tsv",
                "--identity", "0.65",
                "--identity", "0.75",
                "--coverage", "0.60",
                "--coverage", "0.80",
            ])

        assess.assert_called_once_with(
            "segment_matches.tsv",
            "candidate_sensitivity.tsv",
            [0.65, 0.75],
            [0.60, 0.80],
        )

    def test_negative_strand_intron_phase_uses_transcriptional_upstream_cds(self):
        gene = {"seqid": "chr", "strand": "-"}
        exons = [
            {"id": "first", "seqid": "chr", "type": "exon", "start": 3000, "end": 3157, "strand": "-", "phase": "0", "cds_length": 158},
            {"id": "second", "seqid": "chr", "type": "exon", "start": 2000, "end": 2088, "strand": "-", "phase": "1", "cds_length": 89},
            {"id": "third", "seqid": "chr", "type": "exon", "start": 1, "end": 1070, "strand": "-", "phase": "2", "cds_length": 1070},
        ]
        introns = introns_from_path(exons, gene, "tx")
        by_upstream = {row["left_feature_id"]: row for row in introns}
        self.assertEqual(by_upstream["first"]["right_feature_id"], "second")
        self.assertEqual(by_upstream["second"]["right_feature_id"], "third")
        self.assertEqual((by_upstream["first"]["start"], by_upstream["first"]["end"]), (2089, 2999))
        for upstream, expected in (("first", 1), ("second", 2)):
            row = by_upstream[upstream]
            status = phase_compatibility(row["left_phase"], row["right_phase"], row["left_cds_length"])
            self.assertEqual(status, "compatible")
            self.assertEqual(int(row["right_phase"]), expected)

    def test_correspondence_distances_use_topology_for_any_missing_nonroot_length(self):
        topology_distances = {("A", "B"): 2.0, ("A", "C"): 3.0, ("B", "C"): 3.0}
        cases = [
            ("topology_only", ("NA", "NA", "NA", "NA"), topology_distances),
            ("missing_internal", ("NA", "2", "0", "7"), topology_distances),
            ("missing_leaf", ("4", "NA", "0", "7"), topology_distances),
            ("complete_nonroot", ("4", "2", "0", "7"), {("A", "B"): 2.0, ("A", "C"): 13.0, ("B", "C"): 11.0}),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "species_tree.tsv"
            for name, (internal, a, b, c), expected in cases:
                with self.subTest(case=name):
                    path.write_text(
                        "node_id\tparent_id\tlabel\tbranch_length\n"
                        "root\t\troot\tNA\n"
                        f"internal\troot\tinternal\t{internal}\n"
                        f"a\tinternal\tA\t{a}\n"
                        f"b\tinternal\tB\t{b}\n"
                        f"c\troot\tC\t{c}\n"
                    )
                    preprocessing = _species_tree_distances(path)
                    correspondence = tree_distances(root)
                    self.assertEqual(preprocessing, correspondence)
                    for (left, right), distance in expected.items():
                        self.assertEqual(preprocessing[(left, right)], distance)
                        self.assertEqual(preprocessing[(right, left)], distance)
                    for label in ("A", "B", "C"):
                        self.assertEqual(preprocessing[(label, label)], 0.0)

    def test_element_matrix_retains_absent_nonexonic_and_unknown_roles(self):
        occurrences = [
            {"occurrence_id": "A_ex", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present"},
            {"occurrence_id": "B_in", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "role": "intron", "presence_status": "present"},
            {"occurrence_id": "C_abs", "family_id": "fam", "species": "C", "gene_copy_id": "gC", "role": "CDS", "presence_status": "absent"},
            {"occurrence_id": "D_unk", "family_id": "fam", "species": "D", "gene_copy_id": "gD", "role": "unknown", "presence_status": "present"},
            {"occurrence_id": "E_noncoding", "family_id": "fam", "species": "E", "gene_copy_id": "gE", "role": "noncoding", "presence_status": "present"},
            {"occurrence_id": "F_ncexon", "family_id": "fam", "species": "F", "gene_copy_id": "gF", "role": "noncoding_exon", "presence_status": "present"},
            {"occurrence_id": "G_reg", "family_id": "fam", "species": "G", "gene_copy_id": "gG", "role": "regulatory", "presence_status": "present"},
            {"occurrence_id": "H_uncertain", "family_id": "fam", "species": "H", "gene_copy_id": "gH", "role": "intron_or_noncoding", "presence_status": "present"},
            {"occurrence_id": "I_intergenic", "family_id": "fam", "species": "I", "gene_copy_id": "gI", "role": "intergenic", "presence_status": "present"},
        ]
        element_rows = [
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "A_ex", "element_class": "exon_like", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "B_in", "element_class": "candidate_source", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "C_abs", "element_class": "absent", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "D_unk", "element_class": "candidate_source", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "E_noncoding", "element_class": "candidate_source", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "F_ncexon", "element_class": "exon_like", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "G_reg", "element_class": "candidate_source", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "H_uncertain", "element_class": "candidate_source", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "I_intergenic", "element_class": "candidate_source", "membership_call": "core_member"},
        ]
        rows = _element_site_rows(occurrences, element_rows, [], {"fam"}, {"fam": {"A", "B", "C", "D", "E", "F", "G", "H", "I"}})
        by_key = {(row["layer"], row["species"]): row for row in rows}

        self.assertEqual(by_key[("exon_presence", "B")]["state"], "present")
        self.assertEqual(by_key[("exon_role", "B")]["state"], "not_exonic")
        self.assertEqual(by_key[("exon_presence", "C")]["state"], "absent")
        self.assertEqual(by_key[("exon_role", "C")]["state"], "unknown")
        self.assertIn("sequence_absent_role_not_applicable", by_key[("exon_role", "C")]["evidence"])
        self.assertEqual(by_key[("exon_role", "D")]["state"], "unknown")
        self.assertIn("homologous_sequence_role_unknown", by_key[("exon_role", "D")]["evidence"])
        self.assertEqual(by_key[("exon_role", "E")]["state"], "unknown")
        self.assertEqual(by_key[("exon_role", "F")]["state"], "exonic")
        self.assertEqual(by_key[("exon_role", "G")]["state"], "unknown")
        self.assertEqual(by_key[("exon_role", "H")]["state"], "unknown")
        self.assertEqual(by_key[("exon_role", "I")]["state"], "not_exonic")

    def test_raw_regulatory_feature_is_retained_without_hard_role_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "genome.fa"
            gff = root / "annotation.gff3"
            out = root / "out"
            fasta.write_text(">chr1\n" + "A" * 120 + "\n")
            gff.write_text(
                "chr1\ttest\tgene\t10\t90\t.\t+\t.\tID=g1;Name=g1\n"
                "chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=tx1;Parent=g1;partial=true\n"
                "chr1\ttest\texon\t10\t30\t.\t+\t.\tID=ex1;Parent=tx1\n"
                "chr1\ttest\tCDS\t10\t30\t.\t+\t0\tID=cds1;Parent=tx1\n"
                "chr1\ttest\tregulatory\t40\t45\t.\t+\t.\tID=reg1;Parent=g1\n"
            )

            extract_gene(fasta, gff, "g1", "fam", "A", "gA", out, flank=0, max_extension=0)
            raw = read_tsv(out / "raw_gene_features.tsv")
            segments = read_tsv(out / "segment_occurrences.tsv")

            regulatory = [row for row in raw if row["type"] == "regulatory"]
            self.assertEqual(len(regulatory), 1)
            self.assertEqual(regulatory[0]["ownership"], "target_gene_descendant")
            self.assertNotIn("regulatory", {row["role"] for row in segments})

    def test_alternative_exonic_usage_requires_same_locus_without_absence_conflict(self):
        exonic = {
            "occurrence_id": "ex", "family_id": "fam", "species": "A",
            "gene_copy_id": "gA", "transcript_id": "exonic_tx", "role": "exon",
            "presence_status": "present", "contig": "chr1", "strand": "+",
            "start": "40", "end": "60",
        }
        intronic = {
            **exonic, "occurrence_id": "in", "role": "intron",
            "transcript_id": "intronic_tx", "start": "30", "end": "70",
        }
        memberships = [
            {"element_id": "EG_1", "occurrence_id": "ex", "element_class": "exon_like", "membership_call": "core_member"},
            {"element_id": "EG_1", "occurrence_id": "in", "element_class": "candidate_source", "membership_call": "core_member"},
        ]
        cases = [
            ("same_locus", intronic, "exonic"),
            ("different_interval", {**intronic, "start": "130", "end": "170"}, "unknown"),
            ("different_gene", {**intronic, "gene_copy_id": "gA_other"}, "unknown"),
            ("absence_conflict", {**intronic, "presence_status": "absent"}, "unknown"),
        ]
        for name, other, expected in cases:
            with self.subTest(case=name):
                rows = _element_site_rows([exonic, other], memberships, [], {"fam"}, {"fam": {"A"}})
                by_layer = {row["layer"]: row for row in rows}
                self.assertEqual(by_layer["exon_role"]["state"], expected)
                if name == "same_locus":
                    self.assertIn("alternative_usage", by_layer["exon_role"]["evidence"])
                    self.assertNotIn("conflicting_transcript_states", by_layer["exon_role"]["evidence"])
                if name == "absence_conflict":
                    self.assertEqual(by_layer["exon_presence"]["state"], "unknown")

    def test_gene_boundary_prefers_unique_exact_id_and_rejects_ambiguous_alias_loci(self):
        with tempfile.TemporaryDirectory() as tmp:
            gff = Path(tmp) / "annotation.gff3"
            gff.write_text(
                "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1;Name=shared_name;Alias=shared,unique_g1\n"
                "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=g1\n"
                "chr1\ttest\texon\t1\t20\t.\t+\t.\tID=ex1;Parent=tx1\n"
                "chr1\ttest\texon\t120\t140\t.\t+\t.\tID=ex1_extended;Parent=tx1\n"
                "chr1\ttest\tgene\t500\t600\t.\t+\t.\tID=g2;Name=shared_name;Alias=g1,shared\n"
                "chr1\ttest\tmRNA\t500\t600\t.\t+\t.\tID=tx2;Parent=g2;gene_name=g1\n"
                "chr1\ttest\texon\t500\t520\t.\t+\t.\tID=ex2;Parent=tx2\n"
            )

            rows, gene, gene_ids, bounds = read_annotation_for_gene(gff, "g1")
            self.assertEqual(gene["id"], "g1")
            self.assertNotIn("g2", {row.get("id") for row in rows})
            self.assertNotIn("tx2", {row.get("id") for row in rows})
            self.assertIn("ex1_extended", {row.get("id") for row in rows})
            self.assertEqual(bounds["annotation_end"], 100)
            self.assertEqual(bounds["linked_end"], 140)
            self.assertTrue({"g1", "tx1", "ex1", "ex1_extended"} <= gene_ids)
            self.assertFalse({"g2", "tx2", "shared", "shared_name", "unique_g1"} & gene_ids)
            alias_rows, alias_gene, alias_ids, alias_bounds = read_annotation_for_gene(gff, "unique_g1")
            self.assertEqual(alias_gene["id"], "g1")
            self.assertEqual({row.get("id") for row in alias_rows}, {row.get("id") for row in rows})
            self.assertEqual(alias_ids, gene_ids)
            self.assertEqual(alias_bounds, bounds)

            with self.assertRaises(SystemExit):
                read_annotation_for_gene(gff, "shared")

    def test_ambiguous_singleton_membership_keeps_own_annotation(self):
        occurrences = [
            {"occurrence_id": "A_ex", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present"},
        ]
        element_rows = [
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "A_ex", "element_class": "exon_like", "membership_call": "ambiguous_member"},
        ]
        rows = _element_site_rows(occurrences, element_rows, [], {"fam"}, {"fam": {"A"}})
        by_layer = {row["layer"]: row for row in rows}

        self.assertEqual(by_layer["exon_presence"]["state"], "present")
        self.assertEqual(by_layer["exon_role"]["state"], "exonic")
        self.assertIn("ambiguous_correspondence", by_layer["exon_role"]["evidence"])

    def test_ambiguous_crossspecies_membership_does_not_harden_eg_state(self):
        occurrences = [
            {"occurrence_id": "A_ex", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present"},
            {"occurrence_id": "B_ex", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "role": "CDS", "presence_status": "present"},
        ]
        element_rows = [
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "A_ex", "element_class": "exon_like", "membership_call": "core_member"},
            {"element_id": "EG_1", "homology_id": "H1", "occurrence_id": "B_ex", "element_class": "exon_like", "membership_call": "ambiguous_member"},
        ]
        rows = _element_site_rows(occurrences, element_rows, [], {"fam"}, {"fam": {"A", "B"}})
        by_key = {(row["layer"], row["species"]): row for row in rows}

        self.assertEqual(by_key[("exon_presence", "B")]["state"], "unknown")
        self.assertEqual(by_key[("exon_role", "B")]["state"], "unknown")
        self.assertIn("own_annotation_not_used_as_crossspecies_EG_state", by_key[("exon_role", "B")]["evidence"])

    def test_predicted_exon_candidate_keeps_role_unknown(self):
        occurrences = [
            {"occurrence_id": "A_ex", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "role": "CDS", "presence_status": "present"},
        ]
        element_rows = [
            {"element_id": "EG_1", "homology_id": "H_1", "occurrence_id": "A_ex", "element_class": "exon_like", "membership_call": "core_member"},
        ]
        completion_rows = [
            {
                "homology_id": "H_1",
                "family_id": "fam",
                "species": "B",
                "gene_copy_id": "gB",
                "homologous_dna_presence": "present",
                "homologous_dna_evidence": "resolved_nucleotide_alignment",
                "candidate_resolution_status": "resolved",
                "completion_call": "predicted_exon_candidate",
                "inferred_role": "predicted_CDS",
                "predicted_role": "CDS",
                "confidence_flag": "medium",
            }
        ]
        _profiles, element_by_homology = collect_element_profiles(
            [{"homology_id": "H_1", "occurrence_id": "A_ex", "support_type": "fixture", "confidence": "0.95"}],
            occurrences,
            completion_rows,
        )
        self.assertEqual(element_by_homology["H_1"], "EG_1")
        rows = _element_site_rows(occurrences, element_rows, completion_rows, {"fam"}, {"fam": {"A", "B"}})
        by_key = {(row["layer"], row["species"]): row for row in rows}

        self.assertEqual(by_key[("exon_presence", "B")]["state"], "present")
        self.assertEqual(by_key[("exon_role", "B")]["state"], "unknown")
        self.assertIn("sequence_supported_exon_completion_role_unresolved", by_key[("exon_role", "B")]["evidence"])

    def test_supported_exon_prediction_downgrades_overlapping_nonexonic_role(self):
        for call in ("predicted_exon_candidate", "boundary_conflict_candidate"):
            with self.subTest(completion_call=call):
                occurrences, elements, candidate = self._intronic_cds_prediction_case()
                candidate["completion_call"] = call
                if call == "boundary_conflict_candidate":
                    candidate.update(
                        evidence_status="conflicts_annotation",
                        annotation_status="protein_projection_boundary_conflict",
                    )
                rows = _element_site_rows(
                    occurrences, elements, [candidate], {"OG0006454"},
                    {"OG0006454": {"Bombus_ignitus", "Bombus_terrestris"}},
                )
                by_key = {(row["layer"], row["species"]): row for row in rows}
                presence = by_key[("exon_presence", "Bombus_terrestris")]
                role = by_key[("exon_role", "Bombus_terrestris")]
                self.assertEqual(presence["state"], "present")
                self.assertNotIn("role_conflict_unresolved", presence["conclusion_flag"])
                self.assertNotIn("low", presence["confidence_flag"])
                self.assertEqual(role["state"], "unknown")
                self.assertIn("homologous_non_exonic_sequence", role["evidence"])
                self.assertIn("nonexonic_annotation_conflicts_with_exon_prediction", role["evidence"])
                self.assertIn("role_conflict_unresolved", role["conclusion_flag"])
                self.assertIn("low", role["confidence_flag"])
                self.assertEqual(by_key[("exon_role", "Bombus_ignitus")]["state"], "exonic")

    def test_unrelated_or_unresolved_exon_prediction_keeps_nonexonic_role(self):
        cases = [
            ("different_family", {"family_id": "other_family"}),
            ("different_species", {"species": "Bombus_ignitus"}),
            ("different_copy", {"gene_copy_id": "other_gene"}),
            ("different_element", {"homology_id": "HC_other"}),
            ("different_contig", {"interval": "other_contig:2476398-2476441:-"}),
            ("different_strand", {"interval": "NC_063273.1:2476398-2476441:+"}),
            ("nonoverlapping", {"interval": "NC_063273.1:2507266-2507309:-"}),
            ("missing_interval", {"interval": "NA"}),
            ("ambiguous_correspondence", {"correspondence_status": "ambiguous"}),
            ("unknown_correspondence", {"correspondence_status": "unknown"}),
            ("legacy_unspecified_correspondence", {"correspondence_status": "unspecified"}),
            ("low_support", {"confidence_flag": "low"}),
            ("unknown_support", {"confidence_flag": "unknown"}),
            ("unknown_predicted_role", {"predicted_role": "unknown"}),
            ("ordinary_dna_candidate", {
                "completion_call": "homologous_sequence_candidate",
                "primary_mapping_status": "unique_nucleotide_mapping",
                "interval_scope": "aligned_sequence", "predicted_role": "unknown",
            }),
        ]
        for name, overrides in cases:
            with self.subTest(case=name):
                occurrences, elements, candidate = self._intronic_cds_prediction_case()
                occurrences.append({**occurrences[0], "occurrence_id": "other_element_exon"})
                elements.append({
                    "element_id": "EG_other", "homology_id": "HC_other",
                    "occurrence_id": "other_element_exon", "element_class": "exon_like",
                    "membership_call": "core_member",
                })
                candidate.update(overrides)
                rows = _element_site_rows(
                    occurrences, elements, [candidate], {"OG0006454", "other_family"},
                    {"OG0006454": {"Bombus_ignitus", "Bombus_terrestris"}},
                )
                by_layer = {
                    row["layer"]: row for row in rows
                    if row["site_id"] == "EG_0029" and row["species"] == "Bombus_terrestris"
                }
                self.assertEqual(by_layer["exon_presence"]["state"], "present")
                self.assertEqual(by_layer["exon_role"]["state"], "not_exonic")
                self.assertNotIn("nonexonic_annotation_conflicts_with_exon_prediction", by_layer["exon_role"]["evidence"])
