"""Scientific invariants behind the illustrative diagrams, not validation data."""
from contextlib import redirect_stderr
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

from intraphy.cli import main
from intraphy.inference.sankoff import _parsimony_tables
from intraphy.reporting.methods import render_guide, FIGURES
from intraphy.reporting.methods.overview import TOY_TREE, TOY_STATES
from intraphy.topology import SpeciesTree


class MethodGuideTests(unittest.TestCase):
    def test_overview_has_one_required_junction_change(self):
        tree=SpeciesTree(TOY_TREE)
        minimum,states,pairs=_parsimony_tables(tree,TOY_STATES)
        self.assertEqual(minimum,1)
        self.assertEqual(pairs[('AB','B')],{(1,0)})
        self.assertEqual(states['root'],{1})

    def test_two_tip_alternatives_are_not_two_events(self):
        tree=SpeciesTree([{'node_id':'root','parent_id':'','label':'root','branch_length':'0'},
            {'node_id':'A','parent_id':'root','label':'A','branch_length':'1'},
            {'node_id':'B','parent_id':'root','label':'B','branch_length':'1'}])
        minimum,states,pairs=_parsimony_tables(tree,{'A':1,'B':0})
        self.assertEqual(minimum,1)
        self.assertEqual(states['root'],{0,1})
        self.assertIn((0,1),pairs[('root','A')])
        self.assertIn((1,0),pairs[('root','B')])

    def test_same_character_can_change_twice(self):
        minimum,_,_=_parsimony_tables(SpeciesTree(TOY_TREE),{'A':1,'B':0,'C':1,'D':0})
        self.assertEqual(minimum,2)

    def test_all_figures_are_labelled_editable_schematics(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows=render_guide(tmp)
            self.assertEqual(len(rows),3)
            for name,_,_ in FIGURES:
                path=Path(tmp)/f'{name}.svg'
                root=ET.parse(path).getroot()
                self.assertEqual(root.get('data-purpose'),'synthetic-teaching-schematic')
                self.assertIn('SYNTHETIC EXAMPLE',path.read_text())
            self.assertTrue(all(not row['data_result'] for row in rows))
            self.assertEqual(json.loads((Path(tmp)/'figure_manifest.json').read_text()),rows)

    def test_cli_explain_without_genomic_input_and_guards_output(self):
        with tempfile.TemporaryDirectory() as tmp,redirect_stderr(StringIO()):
            output=Path(tmp)/'guide'
            self.assertEqual(main(['explain','--output-dir',str(output)]),0)
            self.assertTrue((output/'index.html').is_file())
            self.assertEqual(main(['explain','--output-dir',str(output)]),2)
            self.assertEqual(main(['explain','--output-dir',str(output),'--force']),0)
            self.assertEqual(len(list(Path(tmp).glob('guide.previous-*'))),1)


if __name__=='__main__': unittest.main()
