"""Explicit, reproducible examples; no fitted or calibrated biological numbers."""
from __future__ import annotations
from itertools import product
import math
from intraphy.inference.ctmc import _transition_matrix

TREE = [
    {'node_id':'root','parent_id':'','label':'root','branch_length':'0'},
    {'node_id':'AB','parent_id':'root','label':'AB','branch_length':'0.4'},
    {'node_id':'CD','parent_id':'root','label':'CD','branch_length':'0.4'},
    *[{'node_id':s,'parent_id':'AB' if s in 'AB' else 'CD','label':s,
       'branch_length':'1.0'} for s in 'ABCD']]
CHARACTERS = {
    'sequence_presence':dict(zip('ABCD',[1,1,1,0])),
    'exonic_status':dict(zip('ABCD',[1,1,'unknown','inapplicable'])),
    'intron_presence':dict(zip('ABCD',[1,0,1,1]))}
GAIN, LOSS, ROOT_PRESENCE = 0.4, 0.2, 2/3

def probability_example(observations=None):
    """Enumerate node states to cross-check the production pruning algorithm."""
    observed=CHARACTERS['intron_presence'] if observations is None else observations
    free=['root','AB','CD']+[s for s in 'ABCD' if observed.get(s) not in (0,1)]
    weighted=[]
    for values in product((0,1),repeat=len(free)):
        states={s:v for s,v in observed.items() if v in (0,1)}
        states.update(zip(free,values))
        weight=ROOT_PRESENCE if states['root'] else 1-ROOT_PRESENCE
        for row in TREE[1:]:
            p=_transition_matrix(GAIN,LOSS,float(row['branch_length']))
            weight*=p[states[row['parent_id']],states[row['node_id']]]
        weighted.append((states,float(weight)))
    likelihood=sum(w for _,w in weighted)
    return {'gain':GAIN,'loss':LOSS,'root_presence':ROOT_PRESENCE,
            'log_likelihood':math.log(likelihood),
            'node_presence':{n:sum(w for st,w in weighted if st[n])/likelihood
                             for n in ['root','AB','CD']},
            'rates_are_fixed_for_illustration':True}
