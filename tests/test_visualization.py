import re
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree

from insiphy.io import read_tsv, write_tsv
from insiphy.visualize import draw_integrated_phylo_synteny, draw_phylogeny, draw_synteny, visualize_results


class VisualizationSemanticTests(unittest.TestCase):
    def write_matches(self, directory, rows):
        write_tsv(
            directory / "segment_matches.tsv",
            rows,
            ["match_id", "query_occurrence_id", "subject_occurrence_id", "match_status",
             "alignment_strand", "projected_reference_blocks", "alignment_backend",
             "correspondence_basis", "protein_projected_blocks"],
        )

    def write_case(self, root):
        input_dir = root / "input"
        result_dir = root / "result"
        output_dir = root / "fig"
        input_dir.mkdir()
        result_dir.mkdir()
        output_dir.mkdir()
        (input_dir / "species_tree.tsv").write_text(
            "node_id\tparent_id\tlabel\tbranch_length\n"
            "root\t\troot\t0\n"
            "a\troot\tA\t0.1\n"
            "b\troot\tB\t0.1\n"
        )
        (input_dir / "segment_occurrences.tsv").write_text(
            "occurrence_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tsegment_id\tcontig\tstart\tend\tstrand\trole\tpresence_status\tboundary_state\tevidence\n"
            "A_e1\tfam\tA\tA_gene\tA_tx1\te1\tchr1\t100\t150\t+\tCDS\tpresent\tconserved\tannotated\n"
            "A_e2\tfam\tA\tA_gene\tA_tx2\te2\tchr1\t220\t260\t+\tCDS\tpresent\tconserved\tannotated\n"
            "A_e3\tfam\tA\tA_gene\tA_tx3\te3\tchr1\t320\t360\t+\tCDS\tpresent\tconserved\tannotated\n"
            "A_e4\tfam\tA\tA_gene\tA_tx3\te4\tchr1\t400\t440\t+\tCDS\tpresent\tconserved\tannotated\n"
            "B_e1\tfam\tB\tB_gene\tB_tx1\te1\tchr1\t100\t150\t+\tCDS\tpresent\tconserved\tannotated\n"
            "B_cand\tfam\tB\tB_gene\tB_tx1\tcand\tchr1\t220\t260\t+\tintron\tpresent\tunknown\tsequence_candidate\n"
            "B_pred\tfam\tB\tB_gene\tB_tx2\tpred\tchr1\t320\t360\t+\tCDS\tpresent\tpredicted\tsupports_hidden_segment\n"
            "B_unknown\tfam\tB\tB_gene\tB_tx2\te4\tchr1\t400\t440\t+\tunknown\tpresent\tconserved\tsequence_supported\n"
        )
        (input_dir / "transcript_paths.tsv").write_text(
            "path_id\tfamily_id\tspecies\tgene_copy_id\ttranscript_id\tpath_rank\toccurrence_id\trole\tcontig\tstart\tend\tstrand\tphase\tpath_status\n"
            "A_tx1_p1\tfam\tA\tA_gene\tA_tx1\t1\tA_e1\tCDS\tchr1\t100\t150\t+\t0\tobserved_transcript_path\n"
            "A_tx2_p0\tfam\tA\tA_gene\tA_tx2\t1\tA_e1\tCDS\tchr1\t100\t150\t+\t0\tobserved_transcript_path\n"
            "A_tx2_p1\tfam\tA\tA_gene\tA_tx2\t1\tA_e2\tCDS\tchr1\t220\t260\t+\t0\tobserved_transcript_path\n"
            "A_tx3_p1\tfam\tA\tA_gene\tA_tx3\t1\tA_e3\tCDS\tchr1\t320\t360\t+\t0\tobserved_transcript_path\n"
            "A_tx3_p2\tfam\tA\tA_gene\tA_tx3\t2\tA_e4\tCDS\tchr1\t400\t440\t+\t0\tobserved_transcript_path\n"
            "B_tx1_p1\tfam\tB\tB_gene\tB_tx1\t1\tB_e1\tCDS\tchr1\t100\t150\t+\t0\tobserved_transcript_path\n"
            "B_tx1_p2\tfam\tB\tB_gene\tB_tx1\t2\tB_cand\tintron\tchr1\t220\t260\t+\t.\tobserved_transcript_path\n"
            "B_tx2_p1\tfam\tB\tB_gene\tB_tx2\t1\tB_pred\tCDS\tchr1\t320\t360\t+\t0\tpredicted_transcript_path\n"
            "B_tx2_p2\tfam\tB\tB_gene\tB_tx2\t2\tB_unknown\tunknown\tchr1\t400\t440\t+\t.\tobserved_transcript_path\n"
        )
        (result_dir / "element_correspondence.tsv").write_text(
            "element_id\tfamily_id\thomology_id\toccurrence_id\tspecies\tgene_copy_id\telement_class\tdisplay_role\tsource_label\tsupport_type\tconfidence\tmembership_score\tmembership_call\tinferred_role\tpredicted_role\n"
            "EG_1\tfam\tH_1\tA_e1\tA\tA_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_1\tfam\tH_1\tB_e1\tB\tB_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_2\tfam\tH_2\tA_e2\tA\tA_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_2\tfam\tH_2\tB_cand\tB\tB_gene\tcandidate_source\tintron\tfam\tsequence_candidate\t0.8\t0.8\tambiguous_member\tnon_exonic_source\tNA\n"
            "EG_3\tfam\tH_3\tA_e3\tA\tA_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_3\tfam\tH_3\tB_pred\tB\tB_gene\texon_like\tCDS\tfam\tpredicted_exon_candidate\t0.8\t0.8\tambiguous_member\tpredicted_CDS\tCDS\n"
            "EG_4\tfam\tH_4\tA_e4\tA\tA_gene\texon_like\tCDS\tfam\tannotation\t1\t1\tcore_member\tCDS\tCDS\n"
            "EG_4\tfam\tH_4\tB_unknown\tB\tB_gene\texon_like\tunknown\tfam\tsequence_supported\t1\t1\tcore_member\tNA\tNA\n"
        )
        self.write_matches(input_dir, [
            {
                "match_id": f"m{i}", "query_occurrence_id": f"A_e{i}",
                "subject_occurrence_id": target, "match_status": "mapped",
                "alignment_strand": "+", "projected_reference_blocks": f"1-{length}:1-{length}",
                "alignment_backend": "mafft_overlap",
            }
            for i, target, length in ((1, "B_e1", 51), (2, "B_cand", 41), (3, "B_pred", 41), (4, "B_unknown", 41))
        ])
        return input_dir, result_dir, output_dir

    def test_transcript_paths_are_drawn_as_separate_lanes(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            draw_synteny(input_dir, result_dir, output_dir)
            svg = (output_dir / "intragenic_synteny.svg").read_text()
            self.assertIn('data-lane-id="A_tx1"', svg)
            self.assertIn('data-lane-id="A_tx2"', svg)
            self.assertIn('data-lane-id="B_tx2"', svg)

    def test_candidate_source_is_visible_but_not_established_ribbon(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            draw_synteny(input_dir, result_dir, output_dir)
            svg = (output_dir / "intragenic_synteny.svg").read_text()
            self.assertIn('data-unit-class="candidate_source"', svg)
            self.assertIn('data-visual-status="sequence_candidate"', svg)
            self.assertIn('data-connector="exon_correspondence" data-element-id="EG_1"', svg)
            self.assertNotIn('data-connector="exon_correspondence" data-element-id="EG_2"', svg)
            self.assertNotIn('data-connector="exon_correspondence" data-element-id="EG_3"', svg)

    def test_shared_transcript_members_do_not_create_same_species_ribbons(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            for renderer in (draw_synteny, draw_integrated_phylo_synteny):
                with self.subTest(renderer=renderer.__name__):
                    path = renderer(input_dir, result_dir, output_dir)
                    svg = ElementTree.parse(path).getroot()
                    ribbons = [
                        node for node in svg.iter()
                        if node.get("data-connector") == "exon_correspondence"
                        and node.get("data-element-id") == "EG_1"
                    ]
                    self.assertEqual(len(ribbons), 1)
                    source = next(
                        node for node in svg.iter()
                        if node.get("data-occurrence-id") == "A_e1"
                        and node.get("data-lane-id") == "A_tx2"
                    )
                    source_y = float(source.get("y")) + float(source.get("height")) / 2
                    self.assertTrue(ribbons[0].get("d").startswith(f'M{float(source.get("x")):.2f},{source_y:.2f} '))

    def test_unknown_role_does_not_create_confirmed_ribbon(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            for renderer in (draw_synteny, draw_integrated_phylo_synteny):
                with self.subTest(renderer=renderer.__name__):
                    path = renderer(input_dir, result_dir, output_dir)
                    svg = ElementTree.parse(path).getroot()
                    unknown = [node for node in svg.iter() if node.get("data-occurrence-id") == "B_unknown"]
                    self.assertEqual(len(unknown), 1)
                    self.assertEqual(unknown[0].get("data-visual-status"), "unknown")
                    self.assertFalse(any(
                        node.get("data-connector") == "exon_correspondence"
                        and node.get("data-element-id") == "EG_4"
                        for node in svg.iter()
                    ))

    def test_predicted_exon_is_visually_distinct(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            draw_synteny(input_dir, result_dir, output_dir)
            svg = (output_dir / "intragenic_synteny.svg").read_text()
            self.assertIn('data-occurrence-id="B_pred"', svg)
            self.assertIn('data-visual-status="predicted"', svg)

    def test_partial_protein_alignment_clips_both_endpoints_without_covering_utr(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            self.write_matches(input_dir, [{
                "match_id": "protein_partial", "query_occurrence_id": "B_e1",
                "subject_occurrence_id": "A_e1", "match_status": "mapped",
                "alignment_strand": "-", "projected_reference_blocks": "1-51:1-51",
                "alignment_backend": "mafft_overlap", "correspondence_basis": "annotated_CDS_protein",
                "protein_projected_blocks": "10-39:7-36",
            }])
            for strand in ("+", "-"):
                for table in ("segment_occurrences.tsv", "transcript_paths.tsv"):
                    rows = read_tsv(input_dir / table)
                    for row in rows:
                        if row["species"] == "B":
                            row["strand"] = strand
                    write_tsv(input_dir / table, rows, list(rows[0]))
                for renderer in (draw_synteny, draw_integrated_phylo_synteny):
                    with self.subTest(renderer=renderer.__name__, strand=strand):
                        svg = ElementTree.parse(renderer(input_dir, result_dir, output_dir)).getroot()
                        ribbons = [node for node in svg.iter() if node.get("data-connector") == "exon_correspondence"]
                        self.assertEqual(len(ribbons), 1)
                        ribbon = ribbons[0]
                        self.assertEqual(ribbon.get("data-source-start"), "7")
                        self.assertEqual(ribbon.get("data-source-end"), "36")
                        self.assertEqual(ribbon.get("data-target-start"), "10")
                        self.assertEqual(ribbon.get("data-target-end"), "39")
                        self.assertEqual(ribbon.get("data-match-id"), "protein_partial")
                        self.assertEqual(ribbon.get("data-correspondence-basis"), "annotated_CDS_protein")
                        self.assertEqual(ribbon.get("data-projection-field"), "protein_projected_blocks")
                        coords = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", ribbon.get("d"))]
                        for occurrence, lane, start, left, right in (
                            ("A_e1", "A_tx2", 7, coords[0], coords[14]),
                            ("B_e1", "B_tx1", 10, coords[6], coords[8]),
                        ):
                            box = next(node for node in svg.iter()
                                       if node.get("data-occurrence-id") == occurrence and node.get("data-lane-id") == lane)
                            x, width = float(box.get("x")), float(box.get("width"))
                            self.assertAlmostEqual(left, x + (start - 1) * width / 51, delta=0.02)
                            self.assertAlmostEqual(right - left, 30 * width / 51, delta=0.02)
                            self.assertGreater(left, x)
                            self.assertLess(right, x + width)

    def test_split_members_attach_to_complementary_reference_intervals(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            occurrences = []
            for occurrence, species, start, end in (
                ("reference", "A", 100, 219), ("left", "B", 100, 159), ("right", "B", 220, 279),
            ):
                occurrences.append({
                    "occurrence_id": occurrence, "family_id": "fam", "species": species,
                    "gene_copy_id": f"{species}_gene", "transcript_id": f"{species}_tx",
                    "contig": "chr1", "start": start, "end": end, "strand": "+",
                    "role": "exon", "presence_status": "present", "boundary_state": "conserved",
                    "evidence": "annotated",
                })
            write_tsv(input_dir / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
            paths = [{**row, "path_rank": i, "path_status": "observed_transcript_path"}
                     for i, row in enumerate(occurrences, 1)]
            write_tsv(input_dir / "transcript_paths.tsv", paths, list(paths[0]))
            elements = [{"element_id": "EG_split", "occurrence_id": row["occurrence_id"],
                         "element_class": "exon_like", "membership_call": "core_member",
                         "support_type": "annotation"} for row in occurrences]
            write_tsv(result_dir / "element_correspondence.tsv", elements, list(elements[0]))
            self.write_matches(input_dir, [
                {"match_id": "split_left", "query_occurrence_id": "left", "subject_occurrence_id": "reference",
                 "match_status": "mapped", "alignment_strand": "+", "projected_reference_blocks": "1-60:1-60",
                 "alignment_backend": "mafft_overlap", "correspondence_basis": "annotated_CDS_protein",
                 "protein_projected_blocks": "1-40:21-60"},
                {"match_id": "split_right", "query_occurrence_id": "reference", "subject_occurrence_id": "right",
                 "match_status": "mapped", "alignment_strand": "+", "projected_reference_blocks": "1-60:1-60",
                 "alignment_backend": "mafft_overlap", "correspondence_basis": "annotated_CDS_protein",
                 "protein_projected_blocks": "61-100:1-40"},
            ])
            for renderer in (draw_synteny, draw_integrated_phylo_synteny):
                with self.subTest(renderer=renderer.__name__):
                    svg = ElementTree.parse(renderer(input_dir, result_dir, output_dir)).getroot()
                    ribbons = [node for node in svg.iter() if node.get("data-connector") == "exon_correspondence"]
                    self.assertEqual(len(ribbons), 2)
                    self.assertEqual(
                        {(node.get("data-source-start"), node.get("data-source-end"), node.get("data-target-occurrence"))
                         for node in ribbons},
                        {("21", "60", "left"), ("61", "100", "right")},
                    )
                    reference = next(node for node in svg.iter() if node.get("data-occurrence-id") == "reference")
                    endpoints = []
                    for node in ribbons:
                        coords = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", node.get("d"))]
                        endpoints.append((coords[0], coords[14]))
                        self.assertAlmostEqual(coords[14] - coords[0], float(reference.get("width")) / 3, delta=0.02)
                    endpoints.sort()
                    self.assertLessEqual(endpoints[0][1], endpoints[1][0])

    def test_missing_ambiguous_or_reverse_matches_do_not_create_full_span_ribbons(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            base = {"match_id": "m", "query_occurrence_id": "A_e1", "subject_occurrence_id": "B_e1",
                    "match_status": "mapped", "alignment_strand": "+", "projected_reference_blocks": "1-30:1-30"}
            cases = (
                [],
                [{**base, "match_status": "ambiguous"}],
                [{**base, "alignment_strand": "-"}],
                [{**base, "alignment_strand": "NA"}],
                [{**base, "projected_reference_blocks": "NA"}],
                [{**base, "correspondence_basis": "annotated_CDS_protein", "protein_projected_blocks": "NA"}],
                [{**base, "correspondence_basis": "annotated_CDS_protein", "protein_projected_blocks": "1-30:30-1"}],
                [base, {**base, "match_id": "conflict", "projected_reference_blocks": "1-30:2-31"}],
            )
            for rows in cases:
                self.write_matches(result_dir, rows)
                for renderer in (draw_synteny, draw_integrated_phylo_synteny):
                    with self.subTest(renderer=renderer.__name__, matches=rows):
                        svg = ElementTree.parse(renderer(input_dir, result_dir, output_dir)).getroot()
                        self.assertFalse(any(node.get("data-connector") for node in svg.iter()))
                        self.assertTrue(any(node.get("data-occurrence-id") == "A_e1" for node in svg.iter()))

    def test_legacy_group_membership_without_base_matches_has_no_ribbon(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            (input_dir / "segment_matches.tsv").unlink()
            for renderer in (draw_synteny, draw_integrated_phylo_synteny):
                with self.subTest(renderer=renderer.__name__):
                    svg = ElementTree.parse(renderer(input_dir, result_dir, output_dir)).getroot()
                    self.assertFalse(any(node.get("data-connector") for node in svg.iter()))
                    self.assertTrue(any(node.get("data-element-id") == "EG_1" for node in svg.iter()))

    def test_invalid_statistical_fit_does_not_draw_probability_event(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            (result_dir / "structural_changes.tsv").write_text(
                "family_id\tlayer\tsite_id\tparent_node\tchild_node\tbranch_scope\tstructural_change_type\tstructural_pattern\tendpoint_change_probability\trate_test_status\tconditioning\n"
                "fam\texon_presence\tEG_bad\troot\ta\troot->A\tposterior_not_reported\tabsent->present\t0.9\tparameters_not_estimable\tfit_status=success\n"
                "fam\texon_presence\tEG_fail\troot\ta\troot->A\tgain\tabsent->present\t0.8\tparameters_not_estimable\tfit_status=nonidentifiable\n"
            )
            draw_phylogeny(input_dir, result_dir, output_dir)
            svg = (output_dir / "phylogenetic_event_map.svg").read_text()
            self.assertIn("No valid branch event display", svg)
            self.assertNotIn("Pr=0.900", svg)
            self.assertNotIn("EG_bad", svg)
            self.assertNotIn("EG_fail", svg)

    def test_valid_selected_fit_posterior_can_draw_without_lrt_test(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            (result_dir / "structural_changes.tsv").write_text(
                "family_id\tlayer\tsite_id\tparent_node\tchild_node\tbranch_scope\tstructural_change_type\tstructural_pattern\tendpoint_change_probability\trate_test_status\tconditioning\n"
                "fam\texon_presence\tEG_er\troot\ta\troot->A\tbidirectional_transition_probabilities\tabsent<->present\t0.4\tparameters_not_estimable\tconditional_MLE;fit_status=success;uncertainty=sensitivity_not_estimated\n"
            )
            draw_phylogeny(input_dir, result_dir, output_dir)
            svg = (output_dir / "phylogenetic_event_map.svg").read_text()
            self.assertIn("Conditional branch probabilities", svg)
            self.assertNotIn("Tested conditional branch probabilities", svg)
            self.assertIn("EG_er", svg)
            self.assertIn("Pr=0.400", svg)

    def test_invalid_endpoint_probability_is_not_replaced_by_another_probability(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            for probability in ("NA", "nan", "inf", "-0.1", "1.1"):
                with self.subTest(probability=probability):
                    (result_dir / "structural_changes.tsv").write_text(
                        "site_id\tbranch_scope\tstructural_change_type\tendpoint_change_probability\tctmc_change_probability\trate_test_status\tconditioning\n"
                        f"EG_bad\troot->A\tgain\t{probability}\t0.9\ttested\tconditional_MLE;fit_status=success\n"
                    )
                    draw_phylogeny(input_dir, result_dir, output_dir)
                    svg = (output_dir / "phylogenetic_event_map.svg").read_text()
                    self.assertIn("No valid branch event display", svg)
                    self.assertNotIn("EG_bad", svg)

    def test_empty_current_posterior_table_does_not_use_legacy_events(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            (result_dir / "structural_changes.tsv").write_text(
                "site_id\tbranch_scope\tstructural_change_type\tendpoint_change_probability\tconditioning\n"
            )
            (result_dir / "event_support_summary.tsv").write_text(
                "site_id\tbranch_scope\tevent_type\tctmc_change_probability\n"
                "EG_legacy\troot->A\tgain\t0.9\n"
            )
            draw_phylogeny(input_dir, result_dir, output_dir)
            svg = (output_dir / "phylogenetic_event_map.svg").read_text()
            self.assertIn("No valid branch event display", svg)
            self.assertNotIn("EG_legacy", svg)

    def test_visualize_keeps_colorblind_palette_and_integrated_view(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            visualize_results(input_dir, result_dir, output_dir)
            synteny = (output_dir / "intragenic_synteny.svg").read_text()
            integrated = (output_dir / "integrated_phylo_synteny.svg").read_text()
            self.assertIn("#0072B2", synteny)
            self.assertIn('data-lane-id="B_tx2"', integrated)
            self.assertIn("No valid branch event display", integrated)

    def test_tree_drawings_use_units_for_missing_lengths_and_preserve_known_lengths(self):
        cases = (
            ("0", "", "0.75", "unit branches (missing lengths)", 1.0),
            ("0", "NA", "0", "unit branches (missing lengths)", 1.0),
            ("NA", "0", "0.75", "supplied branch lengths", 0.0),
            ("", "0.25", "0.75", "supplied branch lengths", 1.0 / 3.0),
            ("", "0", "0", "supplied branch lengths", None),
        )
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            for root_length, a_length, b_length, layout, ratio in cases:
                (input_dir / "species_tree.tsv").write_text(
                    "node_id\tparent_id\tlabel\tbranch_length\n"
                    f"root\t\troot\t{root_length}\n"
                    f"a\troot\tA\t{a_length}\n"
                    f"b\troot\tB\t{b_length}\n"
                )
                for renderer in (draw_phylogeny, draw_integrated_phylo_synteny):
                    with self.subTest(renderer=renderer.__name__, lengths=(root_length, a_length, b_length)):
                        path = renderer(input_dir, result_dir, output_dir)
                        svg = ElementTree.parse(path).getroot()
                        self.assertIn(layout, "".join(svg.itertext()))
                        branch_widths = [
                            float(node.get("x2")) - float(node.get("x1"))
                            for node in svg.iter("{http://www.w3.org/2000/svg}line")
                            if node.get("stroke-width") == "1.2" and node.get("y1") == node.get("y2")
                        ]
                        self.assertEqual(len(branch_widths), 2)
                        if ratio is None:
                            self.assertEqual(branch_widths, [0.0, 0.0])
                        else:
                            self.assertGreater(branch_widths[1], 0)
                            self.assertAlmostEqual(branch_widths[0] / branch_widths[1], ratio, places=4)


if __name__ == "__main__":
    unittest.main()
