"""Exon objects, finite closure, elementary edits, origins and observation contracts."""
import unittest
from dataclasses import replace
import numpy as np
from intraphy.structure.types import Catalogue, ExonSpan as E, ExonConfiguration as C, Material, ObservationEvidence as O
from intraphy.structure.space import enumerate_space
from intraphy.structure.edits import elementary_edits, EDIT_KINDS
from intraphy.structure.material import normalize_exons
from intraphy.structure.observations import compatibility, observation_scenarios
from intraphy.structure.paths import make_graph
from intraphy.structure.origins import origin_scenarios
from intraphy.inference.configuration_history import reconstruct
from intraphy.inference.configuration_model import RateModel, generator, evaluate_model
from intraphy.topology import SpeciesTree


def tree():
    return SpeciesTree([{"node_id":"r","parent_id":"","label":"r"},
        {"node_id":"x","parent_id":"r","label":"x","branch_length":.2},
        {"node_id":"A","parent_id":"x","label":"A","branch_length":.3},
        {"node_id":"B","parent_id":"x","label":"B","branch_length":.3},
        {"node_id":"C","parent_id":"r","label":"C","branch_length":.5}])


def split_catalogue():
    full, split = C((E(0,100),)), C((E(0,40),E(60,100)))
    return Catalogue("g","u",100,(E(0,100),E(0,40),E(60,100)),((40,60),),
        observations=(O("A",(full,)),O("B",(full,)),O("C",(split,))),
        boundary_candidates=(E(0,100),E(0,40),E(60,100)))


class ConfigurationObjects(unittest.TestCase):
    def test_invalid_coordinate_and_overlap(self):
        for args in [(1,1),(-1,4),(1.,5),(True,5)]:
            with self.assertRaises(ValueError): E(*args)
        with self.assertRaises(ValueError): C((E(0,10),E(5,20)))

    def test_split_and_fusion_are_single_edits(self):
        sp=enumerate_space(split_catalogue()); index=sp.index
        a,b=C((E(0,100),)),C((E(0,40),E(60,100)))
        graph=make_graph(sp,{},"A")
        self.assertEqual(graph.distance[index[a],index[b]],1)
        self.assertEqual(graph.distance[index[b],index[a]],1)
        self.assertEqual([e.kind for e in graph.witness_path(index[a],index[b])],["split"])
        self.assertEqual([e.kind for e in graph.witness_path(index[b],index[a])],["fusion"])

    def test_one_to_three_is_two_separator_changes(self):
        c=Catalogue("g","u",100,(E(0,100),E(0,20),E(30,60),E(70,100)),((20,30),(60,70)),
                    boundary_candidates=(E(0,100),E(0,20),E(30,60),E(70,100)))
        sp=enumerate_space(c); graph=make_graph(sp,{},"A")
        a,b=C((E(0,100),)),C((E(0,20),E(30,60),E(70,100)))
        self.assertEqual(graph.distance[sp.index[a],sp.index[b]],2)
        self.assertEqual(len(graph.witness_path(sp.index[a],sp.index[b])),2)

    def test_contiguous_deletion_affecting_two_exons_is_one(self):
        c=Catalogue("g","u",100,(E(10,30),E(60,90)),(),(Material("m",0,100),),
                    boundary_candidates=(E(10,30),E(60,90)))
        sp=enumerate_space(c)
        a,b=C((E(10,30),E(60,90)),(1,)),C((),(2,))
        path=make_graph(sp,{"m":"r"},"A").witness_path(sp.index[a],sp.index[b])
        self.assertEqual(len(path),1)
        self.assertEqual(path[0].kind,"dna_deletion")
        self.assertIn("deleted_exons:2",path[0].consequences)

    def test_exact_intron_deletion_fuses_without_duplicate_count(self):
        c=replace(split_catalogue(),material=(Material("m",40,60),),observations=())
        sp=enumerate_space(c)
        a,b=C((E(0,40),E(60,100)),(1,)),C((E(0,100),),(2,))
        path=make_graph(sp,{"m":"r"},"A").witness_path(sp.index[a],sp.index[b])
        self.assertEqual([e.kind for e in path],["dna_deletion"])
        self.assertIn("exon_fusion",path[0].consequences)

    def test_deleted_material_has_no_reverse_introduction(self):
        c=replace(split_catalogue(),material=(Material("m",40,60),),observations=())
        sp=enumerate_space(c)
        for edit in sp.edits:
            if edit.source.material==(2,): self.assertEqual(edit.target.material,(2,))

    def test_introduction_split_is_single_edit(self):
        c=replace(split_catalogue(),material=(Material("m",40,60),),observations=())
        sp=enumerate_space(c)
        a,b=C((E(0,100),),(0,)),C((E(0,40),E(60,100)),(1,))
        graph=make_graph(sp,{"m":"A"},"A")
        self.assertEqual(graph.distance[sp.index[a],sp.index[b]],1)
        forbidden=make_graph(sp,{"m":"B"},"A")
        self.assertTrue(np.isinf(forbidden.distance[sp.index[a],sp.index[b]]))

    def test_intermediates_not_restricted_to_observed_states(self):
        sp=enumerate_space(split_catalogue())
        self.assertIn(C((E(0,40),)),sp.states)
        self.assertIn(C(()),sp.states)

    def test_state_cap_gates_every_numeric_model(self):
        sp=enumerate_space(split_catalogue(),max_states=2)
        self.assertFalse(sp.complete)
        with self.assertRaises(ValueError): make_graph(sp,{},"A")
        with self.assertRaises(ValueError): generator(sp,RateModel({k:.1 for k in EDIT_KINDS}),{},"A")

    def test_repeated_candidate_does_not_change_q(self):
        c=split_catalogue(); other=replace(c,spans=c.spans+c.spans,junctions=c.junctions*3,
                                        boundary_candidates=c.boundary_candidates*2)
        a,b=enumerate_space(c),enumerate_space(other)
        self.assertEqual(a.states,b.states)
        model=RateModel({k:.1 for k in EDIT_KINDS})
        np.testing.assert_allclose(generator(a,model,{},"A")[0],generator(b,model,{},"A")[0])

    def test_opportunity_mass_is_not_number_of_destinations(self):
        sp=enumerate_space(split_catalogue())
        for state in sp.states:
            masses={}
            for edit in elementary_edits(sp.catalogue,state):
                k=(edit.kind,edit.opportunity);masses[k]=masses.get(k,0)+edit.weight
            for value in masses.values():self.assertAlmostEqual(value,1.)

    def test_coexistence_is_not_unknown(self):
        c=split_catalogue();o=O("A",(C((E(0,100),)),C((E(0,40),E(60,100)))),"coexisting")
        sp=enumerate_space(replace(c,observations=(o,*c.observations[1:])))
        with self.assertRaises(ValueError): compatibility(sp,o)
        self.assertEqual(len(observation_scenarios(sp,tuple("ABC"),"annotation")),2)

    def test_predicted_exon_does_not_harden_absence(self):
        c=Catalogue("g","u",100,(E(0,100),),(),boundary_candidates=(E(0,100),))
        sp=enumerate_space(c)
        o=O("A",(C(()),),"partial",alternative_exons=(E(0,100),))
        self.assertEqual(sum(compatibility(sp,o,"annotation")),1)
        self.assertEqual(sum(compatibility(sp,o,"evidence")),2)

    def test_unknown_not_assigned_one_half(self):
        sp=enumerate_space(split_catalogue());o=O("A",(),"unknown")
        np.testing.assert_equal(compatibility(sp,o),np.ones(len(sp.states)))

    def test_single_origin_mixture_is_normalized_for_unknown_data(self):
        c=replace(split_catalogue(),material=(Material("m",40,60),),observations=())
        sp=enumerate_space(c);t=tree();tips={s:np.ones(len(sp.states)) for s in "ABC"}
        result=evaluate_model(sp,t,tips,RateModel({k:.1 for k in EDIT_KINDS}),counts=False)
        self.assertAlmostEqual(result['log_likelihood'],0,places=10)
        self.assertAlmostEqual(sum(r['posterior_weight'] for r in result['origins']),1)

    def test_one_structural_character_can_change_on_several_branches(self):
        c=split_catalogue();sp=enumerate_space(c);t=tree()
        obs={s:compatibility(sp,o,"annotation") for s,o in zip("ABC",c.observations)}
        result=reconstruct(sp,t,obs)
        self.assertEqual(result['minimum_structural_edits'],1)
        self.assertEqual(len(result['witness']),1)
        self.assertTrue(result['nodes'])


if __name__=='__main__':unittest.main()
