import re

import tempfile

import unittest

from pathlib import Path

from xml.etree import ElementTree

from intraphy.io import read_tsv, write_tsv

from intraphy.visualize import draw_integrated_phylo_synteny, draw_phylogeny, draw_synteny, visualize_results

from support_visualization_visualization_semantic import VisualizationSemanticTestsSupport

class VisualizationSemanticTests(VisualizationSemanticTestsSupport, unittest.TestCase):
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
