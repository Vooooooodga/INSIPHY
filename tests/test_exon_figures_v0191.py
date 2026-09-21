"""Figure semantics, not just XML parseability or unchanged screenshot filenames."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
import numpy as np
from intraphy.reporting.exon_guide import render_guide
from intraphy.reporting.exon_results import draw_unit, render_exon_results
from intraphy.topology import SpeciesTree


class ExonFigureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.root=Path(cls.tmp.name)
        cls.records=render_guide(cls.root)
        cls.data=json.loads((cls.root/"illustrative_model.json").read_text())
        cls.tree=SpeciesTree(cls.data["tree"])
        cls.unit=next(u for u in cls.data["units"] if u["unit_id"]=="split")

    @classmethod
    def tearDownClass(cls): cls.tmp.cleanup()

    def test_old_model_history_cannot_be_rendered_as_current(self):
        source = self.root/"example_results"
        history = source/"exon_history.json"
        original = history.read_text()
        changed = json.loads(original)
        changed["model"] = "exon_configuration_v1"
        try:
            history.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "current configuration model"):
                render_exon_results(source, self.root/"rejected_legacy_result")
        finally:
            history.write_text(original)

    def test_all_current_plates_are_present(self):
        self.assertEqual({r["svg"] for r in self.records},{"method_overview.svg","homology_inference.svg","phylogenetic_model.svg","architecture.svg"})
        for r in self.records:
            ET.parse(self.root/r["svg"])
            self.assertEqual(r["model"],"exon_configuration_v2")
            self.assertFalse(r["data_result"])

    def test_overview_shows_more_than_exon_splitting(self):
        svg=(self.root/"method_overview.svg").read_text()
        for label in ("Sequence deletion","Boundary difference","Two exons / one exon","rooted species tree"):
            self.assertIn(label,svg)

    def test_model_plate_values_are_actual_model_outputs(self):
        xml=ET.parse(self.root/"phylogenetic_model.svg")
        for element in xml.iter():
            node=element.get("data-node")
            if node:
                self.assertAlmostEqual(float(element.get("data-probability")),max(self.unit["ctmc"]["nodes"][node]),places=12)

    def test_annotation_example_keeps_both_conditions(self):
        row=next(r for r in self.data["summary"] if r["family_id"]=="Annotation_control")
        self.assertEqual(row["annotation_conditional_minima"],[1])
        self.assertEqual(row["evidence_compatible_minima"],[0])

    def test_ctmc_values_change_the_result_svg(self):
        a=deepcopy(self.unit);b=deepcopy(self.unit)
        n=len(a["states"])
        a["ctmc"]["nodes"]["root"]=[.8,.2]+[0.]*(n-2)
        b["ctmc"]["nodes"]["root"]=[.1,.9]+[0.]*(n-2)
        first=draw_unit(a,self.tree,[]);second=draw_unit(b,self.tree,[])
        self.assertNotEqual(first,second)
        self.assertIn('data-probability="0.90000000000000002"',second)

    def test_branch_probability_changes_corresponding_glyph(self):
        a=deepcopy(self.unit);b=deepcopy(self.unit)
        a["ctmc"]["branches"][0]["probability_at_least_one_edit"]=.2
        b["ctmc"]["branches"][0]["probability_at_least_one_edit"]=.7
        self.assertNotEqual(draw_unit(a,self.tree,[]),draw_unit(b,self.tree,[]))

    def test_renderer_does_not_refit_or_reconstruct(self):
        with patch("intraphy.inference.configuration_model.evaluate_model",side_effect=AssertionError("refit")), \
             patch("intraphy.inference.configuration_history.reconstruct",side_effect=AssertionError("reinfer")):
            svg=draw_unit(self.unit,self.tree,[],mode="ctmc")
        self.assertIn("conditional CTMC",svg)

    def test_ctmc_and_parsimony_are_distinct_files(self):
        records=json.loads((self.root/"example_figures/figure_manifest.json").read_text())
        self.assertEqual({r["mode"] for r in records},{"parsimony","ctmc"})
        self.assertEqual(len(records),10)
        a=(self.root/"example_figures/exon_structure_00001.svg").read_text()
        b=(self.root/"example_figures/exon_structure_00001_ctmc.svg").read_text()
        self.assertIn("data-event-id",a)
        self.assertNotIn("data-event-id",b)
        self.assertIn("data-any-edit-probability",b)

    def test_invalid_saved_probability_is_not_silently_normalized(self):
        a=deepcopy(self.unit);a["ctmc"]["nodes"]["root"]=[.5]*len(a["states"])
        with self.assertRaises(ValueError): draw_unit(a,self.tree,[])

    def test_missing_ctmc_not_replaced_with_parsimony(self):
        a=deepcopy(self.unit);a.pop("ctmc")
        svg=draw_unit(a,self.tree,[],mode="ctmc")
        self.assertIn("CTMC probabilities unavailable",svg)
        self.assertNotIn("data-probability",svg)

    def test_origin_sensitivity_is_saved_not_selected_as_best_model(self):
        unit=next(u for u in self.data["units"] if u["unit_id"]=="deletion")
        sensitivity=unit["origin_prior_sensitivity"]
        self.assertEqual([r["origin_root_weight"] for r in sensitivity],[.25,4.])
        self.assertTrue(all(r["conditional_sensitivity_not_model_selection"] for r in sensitivity))

    def test_text_is_xml_escaped(self):
        a=deepcopy(self.unit);a["catalogue"]["family"]='g<&"'
        ET.fromstring(draw_unit(a,self.tree,[]))

    def test_no_exon_usage_model_in_guide_inputs(self):
        self.assertTrue(self.data["synthetic"])
        self.assertEqual(self.data["engine"],"exon-ctmc")
        self.assertFalse(self.data["biological_accuracy_validation"])
        self.assertTrue((self.root/"example_inputs/catalogue.jsonl").exists())
