"""Regression tests for elementary state changes and conservative likelihood scope."""
import itertools
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from helpers_v018 import observation, small_tree, write_case
from intraphy.inference.history import history_rows, optimal_history
from intraphy.inference.sankoff import _parsimony_tables
from intraphy.inference.inputs import _validate_complete_universe
from intraphy.inference.likelihood_run import infer_single_copy_phylogeny
from intraphy.observations.schema import read_structural_site_matrix
from intraphy.parsimony import infer_single_copy_parsimony
from intraphy.storage.tabular import read_tsv


class ElementaryInferenceTests(unittest.TestCase):
    def test_witness_attains_optimum_for_every_three_tip_pattern(self):
        tree=small_tree()
        for values in itertools.product((0,1,"unknown"),repeat=3):
            with self.subTest(values=values):
                obs=dict(zip("ABC",values)); result=_parsimony_tables(tree,obs)
                states=optimal_history(tree,obs,*result)
                self.assertEqual(sum(states[p]!=states[c] for p,c in tree.edges()),result[0])

    def test_all_unknown_has_no_witness(self):
        tree=small_tree(); obs=dict.fromkeys("ABC","unknown")
        self.assertEqual(history_rows(("f","exon_presence","x"),tree,obs,*_parsimony_tables(tree,obs)),[])

    def test_alternative_placements_are_not_added(self):
        rows=[observation('j',sp,state) for sp,state in zip("ABC",(0,1,"unknown"))]
        with tempfile.TemporaryDirectory() as tmp:
            source,out,matrix=write_case(tmp,rows)
            infer_single_copy_parsimony(source,out,structural_site_matrix_path=matrix)
            summary=read_tsv(out/'gene_change_summary.tsv')[0]
            self.assertEqual(summary['minimum_character_changes'],'1')
            self.assertEqual(summary['mutation_event_count'],'not_estimated')
            changes=[r for r in read_tsv(out/'branch_structural_events.tsv') if r['event_type']!='no_change']
            self.assertGreater(len(changes),1)
            self.assertEqual(len(changes),len({r['event_id'] for r in changes}))

    def test_no_compound_event_file_even_after_low_level_rerun(self):
        rows=[observation(site,sp,state,group='linked') for site in ('j1','j2')
              for sp,state in zip('ABC',(0,0,1))]
        with tempfile.TemporaryDirectory() as tmp:
            source,out,matrix=write_case(tmp,rows); out.mkdir()
            (out/'compound_structural_events.tsv').write_text('obsolete\n')
            infer_single_copy_parsimony(source,out,structural_site_matrix_path=matrix)
            self.assertFalse((out/'compound_structural_events.tsv').exists())
            self.assertEqual(read_tsv(out/'gene_change_summary.tsv')[0]['minimum_character_changes'],'2')

    def test_linked_characters_never_reach_independent_fit(self):
        rows=[observation(site,sp,state,group='linked') for site in ('j1','j2')
              for sp,state in zip('ABC',(0,0,1))]
        with tempfile.TemporaryDirectory() as tmp:
            source,out,matrix=write_case(tmp,rows)
            with patch('intraphy.inference.layer_fitting.fit_model') as fit:
                infer_single_copy_phylogeny(source,out,structural_site_matrix_path=matrix)
            fit.assert_not_called()
            for record in read_tsv(out/'model_fits.tsv'):
                self.assertEqual(record['aic'],'NA')
                self.assertEqual(record['gain_rate'],'NA')
                self.assertEqual(record['posterior_available'],'false')
            test=read_tsv(out/'model_tests.tsv')[0]
            self.assertEqual(test['p_value'],'NA')
            self.assertEqual(test['posterior_available'],'false')
            self.assertEqual(test['lrt_unavailable_reason'],'correlated_linked_sites_not_modelled')
            self.assertEqual(read_tsv(out/'node_state_posteriors.tsv'),[])
            self.assertEqual(read_tsv(out/'branch_transition_posteriors.tsv'),[])

    def test_complete_catalogue_accepts_explicit_unknown(self):
        rows=[observation('j',sp,value) for sp,value in zip('ABC',(0,1,'unknown'))]
        with tempfile.TemporaryDirectory() as tmp:
            _,_,matrix=write_case(tmp,rows)
            _validate_complete_universe(read_structural_site_matrix(matrix),small_tree())

    def test_nonindependent_catalogue_remains_invalid(self):
        rows=[dict(observation('j',sp,1),discovery_rule='observed_at_least_one') for sp in 'ABC']
        with tempfile.TemporaryDirectory() as tmp:
            _,_,matrix=write_case(tmp,rows)
            with self.assertRaises(SystemExit):
                _validate_complete_universe(read_structural_site_matrix(matrix),small_tree())
