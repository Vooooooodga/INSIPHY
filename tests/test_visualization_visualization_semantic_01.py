import re

import tempfile

import unittest

from pathlib import Path

from xml.etree import ElementTree

from intraphy.io import read_tsv, write_tsv

from intraphy.visualize import draw_integrated_phylo_synteny, draw_phylogeny, draw_synteny, visualize_results

from support_visualization_visualization_semantic import VisualizationSemanticTestsSupport

class VisualizationSemanticTests(VisualizationSemanticTestsSupport, unittest.TestCase):
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
                "protein_projected_blocks": "10-39:7-36", "matched_blocks": "10-39:7-36",
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
                        self.assertEqual(ribbon.get("data-block-source"), "matched_blocks")
                        self.assertEqual(ribbon.get("data-block-field"), "matched_blocks")
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
                 "protein_projected_blocks": "1-40:21-60", "matched_blocks": "1-40:21-60"},
                {"match_id": "split_right", "query_occurrence_id": "reference", "subject_occurrence_id": "right",
                 "match_status": "mapped", "alignment_strand": "+", "projected_reference_blocks": "1-60:1-60",
                 "alignment_backend": "mafft_overlap", "correspondence_basis": "annotated_CDS_protein",
                 "protein_projected_blocks": "61-100:1-40", "matched_blocks": "61-100:1-40"},
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

    def test_legacy_projection_is_labeled_and_never_used_as_ribbon_geometry(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            write_tsv(input_dir / "segment_matches.tsv", [{
                "match_id": "legacy", "query_occurrence_id": "A_e1", "subject_occurrence_id": "B_e1",
                "match_status": "mapped", "alignment_strand": "+", "projected_reference_blocks": "1-51:1-51",
                "alignment_backend": "old_backend",
            }], ["match_id", "query_occurrence_id", "subject_occurrence_id", "match_status",
                "alignment_strand", "projected_reference_blocks", "alignment_backend"])
            svg = ElementTree.parse(draw_synteny(input_dir, result_dir, output_dir)).getroot()
            self.assertFalse(any(node.get("data-connector") == "exon_correspondence" for node in svg.iter()))
            self.assertTrue(any(node.get("data-match-state") == "actual_blocks_missing" for node in svg.iter()))
            self.assertIn("legacy projection fields were not used as ribbon geometry", "".join(svg.itertext()))

    def test_matched_alignment_blocks_show_both_sides_of_an_indel(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            self.write_matches(result_dir, [{
                "match_id": "indel",
                "query_occurrence_id": "A_e1",
                "subject_occurrence_id": "B_e1",
                "match_status": "mapped",
                "alignment_strand": "+",
                "projected_reference_blocks": "7-20:10-23;21-36:25-40",
                "matched_blocks": "7-20:10-23;21-36:25-40",
                "alignment_backend": "mafft_overlap",
            }])
            svg = ElementTree.parse(draw_synteny(input_dir, result_dir, output_dir)).getroot()
            ribbons = [node for node in svg.iter() if node.get("data-connector") == "exon_correspondence"]
            self.assertEqual(len(ribbons), 2)
            self.assertEqual(
                {(node.get("data-source-start"), node.get("data-target-start")) for node in ribbons},
                {("7", "10"), ("21", "25")},
            )

    def test_membership_blocks_can_draw_when_pairwise_table_is_legacy(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            membership_rows = read_tsv(result_dir / "element_correspondence.tsv")
            for row in membership_rows:
                if row["occurrence_id"] == "A_e1":
                    row["matched_blocks"] = "legacy:query:7-20:10-23;30-36:30-36"
                elif row["occurrence_id"] == "B_e1":
                    row["matched_blocks"] = "legacy:subject:7-20:10-23;30-36:30-36"
                else:
                    row["matched_blocks"] = ""
            write_tsv(result_dir / "element_correspondence.tsv", membership_rows, list(membership_rows[0]))
            write_tsv(input_dir / "segment_matches.tsv", [{
                "match_id": "legacy", "query_occurrence_id": "A_e1", "subject_occurrence_id": "B_e1",
                "match_status": "mapped", "alignment_strand": "+", "projected_reference_blocks": "1-51:1-51",
                "alignment_backend": "old_backend",
            }], ["match_id", "query_occurrence_id", "subject_occurrence_id", "match_status",
                "alignment_strand", "projected_reference_blocks", "alignment_backend"])
            svg = ElementTree.parse(draw_synteny(input_dir, result_dir, output_dir)).getroot()
            ribbons = [node for node in svg.iter() if node.get("data-connector") == "membership_correspondence"]
            self.assertEqual(len(ribbons), 2)
            self.assertTrue(all(node.get("data-block-source") == "membership_blocks" for node in ribbons))
            self.assertEqual(
                {(node.get("data-source-start"), node.get("data-target-start")) for node in ribbons},
                {("7", "10"), ("30", "30")},
            )

    def test_tree_is_left_of_species_tracks_and_candidates_have_a_separate_lane(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            svg = ElementTree.parse(draw_synteny(input_dir, result_dir, output_dir)).getroot()
            tree_x = [float(node.get("x1")) for node in svg.iter()
                      if node.get("data-tree-edge") == "horizontal"]
            track_groups = [node for node in svg.iter()
                            if node.get("data-lane-id") == "A_tx1" and node.tag.endswith("g")]
            track_x = [float(line.get("x1")) for group in track_groups for line in group.iter()
                       if line.get("y1") == line.get("y2")]
            self.assertTrue(tree_x)
            self.assertTrue(track_x)
            self.assertLess(max(tree_x), min(track_x))
            candidate = next(node for node in svg.iter()
                             if node.get("data-occurrence-id") == "B_cand" and node.tag.endswith("rect"))
            self.assertEqual(candidate.get("data-lane-id"), "candidate evidence")
            self.assertEqual(candidate.get("data-visual-status"), "sequence_candidate")

    def test_cds_utr_and_noncoding_features_have_distinct_shapes(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            rows = read_tsv(input_dir / "segment_occurrences.tsv")
            utr = dict(rows[0], occurrence_id="A_utr", transcript_id="A_tx1", start="170", end="185", role="three_prime_UTR")
            noncoding = dict(rows[0], occurrence_id="A_nc", transcript_id="A_tx1", start="190", end="205", role="noncoding_exon")
            rows.extend((utr, noncoding))
            write_tsv(input_dir / "segment_occurrences.tsv", rows, list(rows[0]))
            svg = ElementTree.parse(draw_synteny(input_dir, result_dir, output_dir)).getroot()
            shapes = {node.get("data-occurrence-id"): node.get("data-feature-shape") for node in svg.iter()
                      if node.tag.endswith("rect") and node.get("data-occurrence-id")}
            self.assertEqual(shapes["A_e1"], "cds")
            self.assertEqual(shapes["A_utr"], "utr")
            self.assertEqual(shapes["A_nc"], "noncoding")

    def test_short_feature_gets_a_labeled_local_zoom_panel(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            rows = read_tsv(input_dir / "segment_occurrences.tsv")
            micro = dict(rows[0], occurrence_id="A_micro", transcript_id="A_tx_micro", start="50000", end="50029", role="CDS")
            rows.append(micro)
            write_tsv(input_dir / "segment_occurrences.tsv", rows, list(rows[0]))
            svg = ElementTree.parse(draw_synteny(input_dir, result_dir, output_dir)).getroot()
            panels = [node for node in svg.iter() if node.get("data-zoom-panel")]
            self.assertTrue(panels)
            self.assertIn("Local detail, magnified", "".join(svg.itertext()))
            self.assertTrue(any(node.get("data-zoom-start") and node.get("data-zoom-end") for node in panels))

    def test_each_event_branch_has_counts_in_both_phylogenetic_views(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            self.write_branch_events(result_dir, [
                {"family_id": "fam", "layer": "exon_role", "site_id": "EG_A1",
                 "branch_scope": "root->A", "event_type": "exon_role_gain",
                 "placement_status": "required", "call_scope": "core_structural_event"},
                {"family_id": "fam", "layer": "exon_role", "site_id": "EG_A2",
                 "branch_scope": "root->A", "event_type": "exon_role_loss",
                 "placement_status": "required", "call_scope": "core_structural_event"},
                {"family_id": "fam", "layer": "splice_junction", "site_id": "J_A3",
                 "branch_scope": "root->A", "event_type": "intron_gain",
                 "placement_status": "possible", "call_scope": "core_structural_event"},
                {"family_id": "fam", "layer": "splice_junction", "site_id": "J_B1",
                 "branch_scope": "root->B", "event_type": "intron_loss",
                 "placement_status": "possible", "call_scope": "core_structural_event"},
            ])
            expected_a = (
                'data-branch-scope="root-&gt;A" data-event-total="3" '
                'data-event-required="2" data-event-possible="1"'
            )
            expected_b = (
                'data-branch-scope="root-&gt;B" data-event-total="1" '
                'data-event-required="0" data-event-possible="1"'
            )
            for renderer in (draw_phylogeny, draw_integrated_phylo_synteny):
                with self.subTest(renderer=renderer.__name__):
                    svg = renderer(input_dir, result_dir, output_dir).read_text()
                    self.assertIn(expected_a, svg)
                    self.assertIn(expected_b, svg)
                    self.assertIn("3 total (R2/P1)", svg)
                    self.assertIn("1 total (R0/P1)", svg)
                    self.assertIn("EG_A1", svg)
                    self.assertIn("J_B1", svg)

    def test_event_sidebar_reports_exact_overflow_from_parsimony_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            self.write_branch_events(result_dir, [
                {"family_id": "fam", "layer": "exon_role", "site_id": f"EG_{index:02d}",
                 "branch_scope": "root->A", "event_type": "exon_role_gain",
                 "placement_status": "required", "call_scope": "core_structural_event"}
                for index in range(14)
            ])
            svg = draw_phylogeny(input_dir, result_dir, output_dir).read_text()
            self.assertIn(
                'data-event-sidebar-total="14" data-event-sidebar-shown="12" '
                'data-event-sidebar-additional="2" data-event-source="branch_structural_events.tsv"',
                svg,
            )
            self.assertIn(
                "Valid events: total=14; shown=12; additional=2.",
                svg,
            )
            self.assertIn(
                "Additional events remaining in branch_structural_events.tsv: 2.",
                svg,
            )
            self.assertIn("bounded preview", svg)

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
            visualize_results(input_dir, result_dir, output_dir, layout="legacy-overview")
            synteny = (output_dir / "intragenic_synteny.svg").read_text()
            integrated = (output_dir / "integrated_phylo_synteny.svg").read_text()
            self.assertIn("#0072B2", synteny)
            self.assertIn('data-lane-id="B_tx2"', integrated)
            self.assertIn("No valid branch event display", integrated)

    def test_pattern_encoding_is_an_alternative_to_color(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_dir, result_dir, output_dir = self.write_case(Path(tmp))
            draw_synteny(input_dir, result_dir, output_dir, correspondence_encoding="pattern")
            svg = (output_dir / "intragenic_synteny.svg").read_text()
            self.assertIn('data-encoding="pattern"', svg)
            self.assertIn('fill="url(#eg_', svg)
            self.assertNotIn("Color: homologous exon membership", svg)
