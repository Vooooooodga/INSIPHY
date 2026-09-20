"""The illustrated probabilities must match the actual phylogenetic algorithm."""
import math
import tempfile
import unittest
from pathlib import Path
from intraphy.inference.ctmc import _pattern_log_likelihood
from intraphy.inference.posterior import _posterior_messages
from intraphy.inference.sankoff import _parsimony_tables
from intraphy.reporting.methods import render_guide
from intraphy.reporting.methods.example import (
    CHARACTERS, GAIN, LOSS, ROOT_PRESENCE, TREE, probability_example,
)
from intraphy.topology import SpeciesTree


class ModelFigureTests(unittest.TestCase):
    def test_displayed_probabilities_equal_pruning_results(self):
        tree=SpeciesTree(TREE)
        obs=CHARACTERS['intron_presence']
        enumerated=probability_example(obs)
        log_likelihood=_pattern_log_likelihood(tree,obs,GAIN,LOSS,1.0,set(),ROOT_PRESENCE)
        nodes,_=_posterior_messages(tree,obs,GAIN,LOSS,1.0,set(),ROOT_PRESENCE)
        self.assertAlmostEqual(log_likelihood,enumerated['log_likelihood'],places=12)
        for node,value in enumerated['node_presence'].items():
            self.assertAlmostEqual(nodes[node][1],value,places=12)
        self.assertTrue(enumerated['rates_are_fixed_for_illustration'])

    def test_unknown_tip_sums_both_possible_states(self):
        unknown={**CHARACTERS['intron_presence'],'C':'unknown'}
        zero={**unknown,'C':0}
        one={**unknown,'C':1}
        likelihood=lambda obs: math.exp(probability_example(obs)['log_likelihood'])
        self.assertAlmostEqual(likelihood(unknown),likelihood(zero)+likelihood(one),places=12)
        self.assertGreater(likelihood(unknown),likelihood(zero))
        self.assertAlmostEqual(probability_example({})['log_likelihood'],0,places=12)

    def test_other_displayed_loss_is_a_different_character(self):
        minimum,states,pairs=_parsimony_tables(SpeciesTree(TREE),CHARACTERS['sequence_presence'])
        self.assertEqual(minimum,1)
        self.assertEqual(pairs[('CD','D')],{(1,0)})
        self.assertEqual(states['root'],{1})
        self.assertEqual(CHARACTERS['exonic_status']['C'],'unknown')
        self.assertEqual(CHARACTERS['exonic_status']['D'],'inapplicable')

    def test_figures_use_biological_labels_and_real_algorithm_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            render_guide(tmp)
            text='\n'.join(p.read_text() for p in Path(tmp).glob('*.svg'))
            for stale in ('Supplied tree','annotated paths','freeze observations',
                          'DNA b','exon role b','junction j'):
                self.assertNotIn(stale,text)
            for expected in ('MAFFT','minimap2','Ordered-chain dynamic programming',
                             'Homologous DNA','Q =','Unknown','ancestral'):
                self.assertIn(expected,text)
            self.assertTrue((Path(tmp)/'illustrative_model.json').is_file())

if __name__=='__main__': unittest.main()
