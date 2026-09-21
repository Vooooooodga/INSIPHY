"""Independent enumeration and analytic calculations, not precomputed answers."""
import itertools
import math
import unittest
import numpy as np
from scipy.linalg import expm
from intraphy.inference.configuration_dp import sankoff
from intraphy.inference.configuration_ctmc import likelihood,transition_matrix,expected_count,probability_any_change
from intraphy.topology import SpeciesTree


def tree():
    return SpeciesTree([{'node_id':'r','parent_id':'','label':'r'},
        {'node_id':'x','parent_id':'r','label':'x','branch_length':.4},
        {'node_id':'A','parent_id':'x','label':'A','branch_length':.5},
        {'node_id':'B','parent_id':'x','label':'B','branch_length':.6},
        {'node_id':'C','parent_id':'r','label':'C','branch_length':.7}])


class NumericConfigurationTests(unittest.TestCase):
    def test_sankoff_matches_full_enumeration_all_optima(self):
        t=tree();c=np.array([[0,1,3],[2,0,1],[1,2,0]],float)
        costs={child:c for _,child in t.edges()}
        for values in itertools.product(range(3),repeat=3):
            tips={s:np.eye(3)[v] for s,v in zip('ABC',values)}
            result=sankoff(t,tips,costs)
            all_states=[]
            for r,x in itertools.product(range(3),repeat=2):
                states=dict(zip('ABC',values));states.update(r=r,x=x)
                score=sum(c[states[p],states[ch]] for p,ch in t.edges())
                all_states.append((score,states))
            best=min(s for s,_ in all_states)
            optimal=[states for s,states in all_states if s==best]
            self.assertEqual(result.cost,best)
            for node in t.parent:self.assertEqual(set(result.nodes[node]),{s[node] for s in optimal})
            for p,ch in t.edges():self.assertEqual(set(result.pairs[p,ch]),{(s[p],s[ch]) for s in optimal})
            self.assertEqual(sum(c[result.witness[p],result.witness[ch]] for p,ch in t.edges()),best)

    def test_pruning_matches_enumeration_with_partial_observation(self):
        t=tree();q=np.array([[-.5,.3,.2],[.1,-.3,.2],[0,0,0]])
        ps={c:transition_matrix(q,t.branch_length(c)) for _,c in t.edges()}
        tips={'A':np.array([1.,0,0]),'B':np.array([0.,1.,1.]),'C':np.ones(3)};prior=np.array([.5,.3,.2])
        result=likelihood(t,tips,ps,prior)
        total=0.;counts={v:np.zeros(3) for v in t.parent}
        for values in itertools.product(range(3),repeat=5):
            states=dict(zip(t.parent,values));prob=prior[states['r']]
            for p,c in t.edges():prob*=ps[c][states[p],states[c]]
            for label,weights in tips.items():prob*=weights[states[label]]
            total+=prob
            for node in states:counts[node][states[node]]+=prob
        self.assertAlmostEqual(result.log_likelihood,math.log(total))
        for node in counts:np.testing.assert_allclose(result.nodes[node],counts[node]/total,atol=1e-12)
        for endpoint in result.endpoints.values():self.assertAlmostEqual(endpoint.sum(),1)

    def test_transition_semigroup_and_zero(self):
        q=np.array([[-.3,.2,.1],[.4,-.7,.3],[.5,.6,-1.1]])
        np.testing.assert_allclose(transition_matrix(q,0),np.eye(3))
        np.testing.assert_allclose(transition_matrix(q,.3)@transition_matrix(q,.7),transition_matrix(q,1),atol=1e-12)

    def test_expected_counts_differ_from_endpoint_changes(self):
        rate,length=2.,1.3;q=np.array([[-rate,rate],[rate,-rate]])
        p=transition_matrix(q,length);endpoints=.5*p
        self.assertAlmostEqual(expected_count(q,p,endpoints,length),rate*length,places=10)
        self.assertAlmostEqual(probability_any_change(q,p,endpoints,length),1-math.exp(-rate*length),places=10)
        self.assertLess(1-np.trace(endpoints),probability_any_change(q,p,endpoints,length))
        b=np.array([[0.,rate],[0.,0.]])
        self.assertAlmostEqual(expected_count(q,p,endpoints,length,b),rate*length/2,places=10)

    def test_absorbing_deletion_marked_count(self):
        q=np.array([[-.7,.7],[0.,0.]])
        p=transition_matrix(q,2);endpoint=np.zeros((2,2));endpoint[0,:]=p[0,:]
        self.assertAlmostEqual(expected_count(q,p,endpoint,2),1-math.exp(-1.4),places=10)

    def test_unknown_likelihood_is_one_not_half(self):
        t=tree();q=np.array([[-.3,.3],[.2,-.2]]);ps={c:transition_matrix(q,t.branch_length(c)) for _,c in t.edges()}
        result=likelihood(t,{s:np.ones(2) for s in 'ABC'},ps,np.array([.5,.5]))
        self.assertAlmostEqual(result.log_likelihood,0)

    def test_zero_length_impossible_observation(self):
        t=SpeciesTree([{'node_id':'r','parent_id':'','label':'r'},
                       {'node_id':'A','parent_id':'r','label':'A','branch_length':0},
                       {'node_id':'B','parent_id':'r','label':'B','branch_length':0}])
        r=likelihood(t,{'A':np.array([1.,0.]),'B':np.array([0.,1.])},{'A':np.eye(2),'B':np.eye(2)},np.array([.5,.5]))
        self.assertEqual(r.log_likelihood,-math.inf);self.assertFalse(r.nodes)

    def test_invalid_q_not_renormalized(self):
        for q in [np.array([[-.5,.4],[0,0]]),np.array([[.1,-.1],[0,0]]),np.array([[math.nan,0],[0,0]])]:
            with self.assertRaises(ValueError):transition_matrix(q,1)


if __name__=='__main__':unittest.main()
