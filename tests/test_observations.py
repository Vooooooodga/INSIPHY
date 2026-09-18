import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from insiphy.correspondence import infer_correspondence, match_total, tree_distances
from insiphy.elements import collect_element_profiles
from insiphy.io import read_tsv
from insiphy.preprocess import _species_tree_distances, cheap_match_evidence, cluster_segments, copy_order_context, extract_gene, graph_components, introns_from_path, read_annotation_for_gene, transcript_cds_length
from insiphy.alignment import phase_compatibility
from insiphy.structural_sites import _element_site_rows, _genomically_contiguous, _junction_site_rows


class ObservationSemanticsTests(unittest.TestCase):
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
                "chr1\ttest\tmRNA\t10\t90\t.\t+\t.\tID=tx1;Parent=g1\n"
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

    def _intronic_cds_prediction_case(self):
        occurrences = [
            {
                "occurrence_id": "Bign_exon", "family_id": "OG0006454",
                "species": "Bombus_ignitus", "gene_copy_id": "gene:ENSBIGG00000010441",
                "role": "exon", "presence_status": "present", "contig": "LG5",
                "start": "1947147", "end": "1947446", "strand": "-",
            },
            {
                "occurrence_id": "Bter_intron", "family_id": "OG0006454",
                "species": "Bombus_terrestris", "gene_copy_id": "gene-LOC100645322",
                "role": "intron", "presence_status": "present", "contig": "NC_063273.1",
                "start": "2472297", "end": "2507265", "strand": "-",
            },
        ]
        elements = [
            {"element_id": "EG_0029", "homology_id": "HC_0029", "occurrence_id": "Bign_exon", "element_class": "exon_like", "membership_call": "core_member"},
            {"element_id": "EG_0029", "homology_id": "HC_0029", "occurrence_id": "Bter_intron", "element_class": "candidate_source", "membership_call": "core_member"},
        ]
        candidate = {
            "family_id": "OG0006454", "homology_id": "HC_0029",
            "species": "Bombus_terrestris", "gene_copy_id": "gene-LOC100645322",
            "interval": "NC_063273.1:2476398-2476441:-",
            "completion_call": "predicted_exon_candidate",
            "evidence_status": "supports_hidden_segment",
            "annotation_status": "protein_projection_supports_missing_cds",
            "inferred_event": "protein_cds_projection", "inferred_role": "predicted_CDS",
            "predicted_role": "CDS", "support_score": "0.8667",
            "primary_mapping_status": "protein_projection", "correspondence_status": "resolved",
            "interval_scope": "projected_cds_container", "confidence_flag": "high",
        }
        return occurrences, elements, candidate

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
        absence = {**candidate, "completion_call": "supports_true_absence"}
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
                unresolved = status in {"unknown", "ambiguous"}
                if unresolved:
                    completion["primary_mapping_status"] = "ambiguous_repeated_mapping"
                    completion["interval_candidates"] = '[{"target_start":10,"target_end":30},{"target_start":50,"target_end":70}]'
                rows = _element_site_rows(occurrences, elements, [completion], {"fam"}, {"fam": {"A", "B"}})
                by_key = {(row["layer"], row["species"]): row for row in rows}
                presence = by_key[("exon_presence", "B")]
                role = by_key[("exon_role", "B")]
                self.assertEqual(presence["state"], "unknown" if unresolved else "present")
                self.assertEqual(role["state"], "unknown")
                self.assertEqual(by_key[("exon_presence", "A")]["state"], "present")
                if unresolved:
                    self.assertIn("completion_correspondence_unresolved", presence["evidence"])
                    self.assertIn("completion_correspondence_unresolved", role["evidence"])
                    self.assertEqual(presence["confidence_flag"], "low")
                else:
                    self.assertNotIn("completion_correspondence_unresolved", presence["evidence"])

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
            with patch("insiphy.correspondence.simple_identity", side_effect=[0.95, 0.25, 0.95]) as identity:
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

    def test_same_locus_alternative_exons_share_element_without_merging_tandem_repeat(self):
        occurrences = [
            {
                "occurrence_id": "alt_short",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx1",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "3535611",
                "end": "3536164",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
            {
                "occurrence_id": "alt_long",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx2",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "3535611",
                "end": "3536182",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
            {
                "occurrence_id": "nearby_repeat",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx3",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "3537000",
                "end": "3537553",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
        ]

        homology, matches = cluster_segments(occurrences, {}, threads=1)
        by_occ = {row["occurrence_id"]: row["homology_id"] for row in homology}
        overlap_match = [row for row in matches if {row["query_occurrence_id"], row["subject_occurrence_id"]} == {"alt_short", "alt_long"}]

        self.assertEqual(by_occ["alt_short"], by_occ["alt_long"])
        self.assertNotEqual(by_occ["alt_short"], by_occ["nearby_repeat"])
        self.assertEqual(overlap_match[0]["alignment_backend"], "genomic_overlap")
        self.assertNotEqual(overlap_match[0]["projected_reference_blocks"], "NA")

    def test_same_locus_overlap_with_shared_transcript_membership_is_not_alternative_evidence(self):
        occurrences = [
            {
                "occurrence_id": "left",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx1;tx_shared",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "100",
                "end": "180",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
            {
                "occurrence_id": "right",
                "family_id": "fam",
                "species": "A",
                "gene_copy_id": "gA",
                "transcript_id": "tx2;tx_shared",
                "role": "CDS",
                "presence_status": "present",
                "contig": "chr1",
                "start": "120",
                "end": "200",
                "strand": "+",
                "phase": "0",
                "splice_motif_score": "0.5",
            },
        ]

        homology, matches = cluster_segments(occurrences, {}, threads=1)
        by_occ = {row["occurrence_id"]: row["homology_id"] for row in homology}

        self.assertNotEqual(by_occ["left"], by_occ["right"])
        self.assertFalse(matches)


if __name__ == "__main__":
    unittest.main()
