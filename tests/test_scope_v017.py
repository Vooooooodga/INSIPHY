"""0.17 scientific-scope and target-figure contracts, using explicit synthetic data.

These test known/unknown semantics and computations. They are not empirical
calibration of alignment accuracy, CTMC intervals or evolutionary truth.
"""
from pathlib import Path
import sys
import json
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from build_v017_example import build_example, observation, write_rows, FAMILY, SPECIES
from intraphy.observations.eligibility import ScopePolicy,assess_scope
from intraphy.observations.matrix import prepare_observation_matrix
from intraphy.observations.schema import write_structural_site_matrix
from intraphy.storage.tabular import read_tsv
from intraphy.phylogeny import infer_phylogeny
from intraphy.visualize import visualize_results
from intraphy.candidate_chain import (ChainCandidate,ChainPathMembership,ordered_candidate_chain,genomic_candidate_chain)
from intraphy.coordinates import Interval0
from intraphy.preparation.transcripts import _path_role_record
from intraphy.mapping.policies import _candidate_sequence_accepted
from intraphy.observations.junctions import _explicit_intron_between
from intraphy.reporting.target_model import load_targets,FigureTarget,resolve_intervals,target_ribbons
from intraphy.reporting.target_scene import build_scene

SVG='{http://www.w3.org/2000/svg}'


def membership(qrank,trank,qp='qex',tp='tex',qid='qt',tid='tt'):
    return ChainPathMembership(qid,tid,qrank,trank,'qchr','tchr','+','+',qp,tp)


def cand(name,qa,qb,ta,tb,pm=(),strand='+'):
    return ChainCandidate(name,Interval0(qa,qb),Interval0(ta,tb),100.,'nt_test',strand,tuple(pm))


class CoordinateAndCandidateV017(unittest.TestCase):
    def test_shorter_alternative_first_is_not_partial(self):
        f={'id':'E','type':'exon','role':'exon','start':200,'end':400,'strand':'+','attrs':{}}
        tx={'id':'T','strand':'+','attrs':{},'_gene_start':100,'_gene_end':500}
        row=_path_role_record('T',1,[f],0,f,tx)
        self.assertEqual((row['partial_start'],row['partial_end']),('0','0'))

    def test_shorter_alternative_negative_is_not_partial(self):
        f={'id':'E','type':'exon','start':200,'end':400,'strand':'-','attrs':{}}
        row=_path_role_record('T',1,[f],0,f,{'strand':'-','attrs':{},'_gene_start':100,'_gene_end':500})
        self.assertEqual((row['partial_start'],row['partial_end']),('0','0'))

    def test_explicit_partial_is_retained(self):
        f={'id':'E','type':'exon','start':200,'end':400,'strand':'+','attrs':{'partial':'true'}}
        row=_path_role_record('T',1,[f],0,f,{'strand':'+','attrs':{}})
        self.assertEqual(row['partial_start'],'1')

    def test_same_parent_split_complements_chain(self):
        a=cand('a',0,50,0,50,[membership(1,1,'Q','T1')]);b=cand('b',50,100,65,115,[membership(1,2,'Q','T2')])
        result=ordered_candidate_chain([a,b],0)
        self.assertEqual(result.best_score,200)
        self.assertEqual(result.retained_ids,{'a','b'})

    def test_same_parent_merge_complements_chain(self):
        a=cand('a',0,50,0,50,[membership(1,1,'Q1','T')]);b=cand('b',65,115,50,100,[membership(2,1,'Q2','T')])
        self.assertEqual(ordered_candidate_chain([a,b],0).best_score,200)

    def test_parent_ids_required_for_equal_ranks(self):
        a=cand('a',0,50,0,50,[membership(1,1,'','')]);b=cand('b',50,100,50,100,[membership(1,1,'','')])
        self.assertEqual(ordered_candidate_chain([a,b],0).best_score,100)

    def test_repeated_coverage_is_not_complementary_split(self):
        a=cand('a',0,60,0,60,[membership(1,1)]);b=cand('b',20,80,80,140,[membership(1,1)])
        self.assertEqual(ordered_candidate_chain([a,b],0).best_score,100)

    def test_pathless_DNA_and_annotated_path_are_separate(self):
        candidates=[cand('a',0,20,0,20,[membership(1,1)]),cand('u',30,50,30,50),cand('c',60,80,60,80,[membership(3,3)])]
        physical=genomic_candidate_chain(candidates,0);rna=ordered_candidate_chain(candidates,0)
        self.assertEqual(physical.best_score,300)
        self.assertEqual(rna.best_score,200)
        self.assertIn('u',physical.retained_ids);self.assertNotIn('u',rna.retained_ids)

    def test_reversed_piece_is_not_an_ordered_chain(self):
        a=cand('a',0,20,0,20);b=cand('b',30,50,30,50,strand='-')
        self.assertNotIn(('a','b'),genomic_candidate_chain([a,b],0).retained_edges)

    def test_microexon_needs_context(self):
        r={'identity':1.,'coverage':1.,'query_coverage':1.,'known_aligned_pairs':9}
        self.assertFalse(_candidate_sequence_accepted(r,.7,True))
        self.assertTrue(_candidate_sequence_accepted(r,.7,True,anchored_microexon=True))

    def test_microexon_mismatch_not_rescued(self):
        r={'identity':.8,'coverage':1.,'query_coverage':1.,'known_aligned_pairs':9}
        self.assertFalse(_candidate_sequence_accepted(r,.7,True,anchored_microexon=True))

    def test_dinucleotide_not_a_microexon_match(self):
        r={'identity':1.,'coverage':1.,'query_coverage':1.,'known_aligned_pairs':2}
        self.assertFalse(_candidate_sequence_accepted(r,.7,True,anchored_microexon=True))

    def test_unmapped_middle_exon_cannot_create_junction(self):
        path=[{'role':'exon'},{'role':'intron'},{'role':'exon'},{'role':'intron'},{'role':'exon'}]
        self.assertFalse(_explicit_intron_between(path,0,4,{}, {},[], 'T'))
        self.assertTrue(_explicit_intron_between(path,0,2,{}, {},[], 'T'))


class CallabilityV017(unittest.TestCase):
    def rows(self,values,layer='exon_presence'):
        return [observation('X',str(i),x,layer) for i,x in enumerate(values)]

    def test_known_zeros_are_calls_not_alignment_failure(self):
        _,selected,_,summary=assess_scope(self.rows([1,1]+[0]*8),policy=ScopePolicy('high-coverage'))
        self.assertEqual(len(selected),10)
        self.assertEqual(summary[0]['callable_fraction_panel'],1)
        self.assertEqual(summary[0]['state1_fraction_called'],.2)

    def test_unknowns_not_counted_as_absence(self):
        full,selected,_,summary=assess_scope(self.rows([1,1]+[None]*8),policy=ScopePolicy('high-coverage'))
        self.assertEqual(len(full),10);self.assertEqual(selected,[])
        self.assertEqual(summary[0]['state0_n'],0)
        self.assertEqual(summary[0]['callable_fraction_panel'],.2)

    def test_role_panel_and_applicable_denominators(self):
        rows=self.rows([1,1]+[None]*8,'exon_role')
        for row in rows[2:]:row['applicability']='inapplicable'
        _,selected,_,summary=assess_scope(rows,policy=ScopePolicy('high-coverage'))
        self.assertEqual(summary[0]['callable_fraction_applicable'],1)
        self.assertEqual(summary[0]['callable_fraction_panel'],.2)
        self.assertFalse(selected)

    def test_all_range_retains_local_low_coverage_sites(self):
        rows=self.rows([1,None,None,None]);full,selected,_,summary=assess_scope(rows)
        self.assertEqual(full,selected);self.assertEqual(summary[0]['callable_n'],1)

    def test_masked_known_state_preserved_only_in_audit(self):
        rows=self.rows([1,0]);rows[0]['observation_mask']='masked'
        full,_,audit,_=assess_scope(rows)
        self.assertEqual(full[0]['state'],'unknown');self.assertEqual(audit[0]['input_state'],'present')
        self.assertEqual(rows[0]['state'],'present')

    def test_threshold_uses_greater_equal(self):
        rows=self.rows([1]+[0]*6+[None]*3)
        self.assertEqual(len(assess_scope(rows,policy=ScopePolicy('high-coverage',.7))[1]),10)

    def test_thresholds_validate(self):
        for x in (-.1,1.1,float('nan'),float('inf')):
            with self.subTest(value=x),self.assertRaises(ValueError):ScopePolicy(min_callable_fraction=x)
        with self.assertRaises(ValueError):ScopePolicy('anything')

    def test_exact_panel_required(self):
        with self.assertRaises(ValueError):assess_scope(self.rows([1,0]),panel_species=['0','1','2'])
        with self.assertRaises(ValueError):assess_scope(self.rows([1,0]),panel_species=['0','0'])

    def test_duplicate_row_rejected(self):
        r=self.rows([1,0]);r.append(dict(r[0]))
        with self.assertRaises(ValueError):assess_scope(r)

    def test_scope_does_not_filter_by_variability(self):
        full,selected,_,_=assess_scope(self.rows([0,0,0,0]),policy=ScopePolicy('high-coverage'))
        self.assertEqual(full,selected)

    def test_channels_not_masked_together(self):
        a=self.rows([None,None]);b=self.rows([1,0],layer='splice_junction')
        _,selected,_,_=assess_scope(a+b,policy=ScopePolicy('high-coverage'))
        self.assertEqual({r['layer'] for r in selected},{'splice_junction'})

    def test_reported_tree_scope_is_local(self):
        from intraphy.topology import SpeciesTree
        tree=SpeciesTree([{'node_id':n,'parent_id':p,'label':n,'branch_length':'1'} for n,p in [('root',''),('AB','root'),('CD','root'),('A','AB'),('B','AB'),('C','CD'),('D','CD')]])
        rows=[observation('X',sp,1 if sp in 'AB' else None) for sp in 'ABCD']
        summary=assess_scope(rows,tree=tree)[3][0]
        self.assertEqual(summary['called_mrca'],'AB');self.assertEqual(summary['represented_root_children'],1)

    def test_matrix_snapshot_and_full_catalogue(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp,prepared=build_example(tmp);out=Path(tmp)/'out'
            matrix=prepare_observation_matrix(inp,out,prepared/'structural_site_matrix.tsv',analysis_range='high-coverage')
            self.assertEqual(len(matrix.full_rows),24);self.assertEqual(len(matrix.rows),16)
            self.assertEqual(len(read_tsv(out/'structural_site_matrix.tsv')),24)
            self.assertEqual(len(read_tsv(out/'analysis_structural_site_matrix.tsv')),16)
            with self.assertRaises(TypeError):matrix.rows[0]['state']='absent'
            metadata=json.loads((out/'analysis_scope.json').read_text())
            self.assertFalse(metadata['coordinates_and_adjacencies_changed'])
            self.assertEqual(metadata['selected_characters'],4)

    def test_same_matrix_ownership_for_formal_engines(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp,prepared=build_example(tmp)
            captures=[]
            def capture(*args,**kwargs):captures.append(kwargs['observation_matrix'])
            with patch('intraphy.parsimony.infer_single_copy_parsimony',side_effect=capture):
                infer_phylogeny(inp,Path(tmp)/'p',structural_site_matrix_path=prepared/'structural_site_matrix.tsv',analysis_range='high-coverage')
            with patch('intraphy.structural_phylogeny.infer_single_copy_phylogeny',side_effect=capture):
                infer_phylogeny(inp,Path(tmp)/'c',model='er-ard',structural_site_matrix_path=prepared/'structural_site_matrix.tsv',analysis_range='high-coverage')
            self.assertEqual(captures[0].rows,captures[1].rows)


class TargetFiguresV017(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.inp,self.prepared=build_example(self.root)
        self.result=self.root/'out';self.result.mkdir()
        import shutil
        for name in ('element_correspondence.tsv','splice_boundary_correspondence.tsv'):
            shutil.copyfile(self.prepared/name,self.result/name)
        infer_phylogeny(self.inp,self.result,model='parsimony',structural_site_matrix_path=self.prepared/'structural_site_matrix.tsv')
        self.fig=self.root/'fig'
        self.manifest=visualize_results(self.inp,self.result,self.fig,target_manifest=self.root/'targets.tsv')

    def tearDown(self):self.tmp.cleanup()

    def svg(self,target='exon_3',view='phylogeny'):
        row=next(x for x in self.manifest if x['target_id']==target and x['view']==view)
        return ET.parse(self.fig/row['path']).getroot()

    def test_three_views_for_each_target(self):
        self.assertEqual(len(self.manifest),12)
        self.assertEqual({r['view'] for r in self.manifest},{'structure','ribbons','phylogeny'})

    def test_background_identical_across_all_twelve_figures(self):
        backgrounds=[]
        for r in self.manifest:
            root=ET.parse(self.fig/r['path']).getroot()
            backgrounds.append(ET.tostring(root.find(f"{SVG}g[@id='native-background']")))
        self.assertEqual(len(set(backgrounds)),1)

    def test_background_is_grey_not_target_colored(self):
        bg=self.svg().find(f"{SVG}g[@id='native-background']")
        for n in bg.iter():
            color=n.get('fill','')
            if color.startswith('#'):
                self.assertEqual(color[1:3],color[3:5]);self.assertEqual(color[1:3],color[5:7])

    def test_other_owner_and_RNA_annotations_retained(self):
        bg=self.svg().find(f"{SVG}g[@id='native-background']")
        fs=[n for n in bg.iter() if n.get('data-native')=='true']
        self.assertIn('ncRNA',{n.get('data-kind') for n in fs})
        self.assertIn('repeat_region',{n.get('data-kind') for n in fs})
        anti=next(n for n in fs if n.get('data-feature')=='A_antisense')
        self.assertEqual(anti.get('data-strand'),'-')

    def test_target_colors_unique_and_consistent(self):
        colors={}
        for r in self.manifest:
            self.assertEqual(colors.setdefault(r['target_id'],r['color']),r['color'])
            root=ET.parse(self.fig/r['path']).getroot()
            self.assertEqual(root.get('data-target-color'),r['color'])
            self.assertEqual(root.find(f"{SVG}g[@id='target-legend']").get('data-color'),r['color'])
        self.assertEqual(len(set(colors.values())),4)

    def test_ribbons_only_selected_target(self):
        for target in ('exon_3','exon_4','intron_segment'):
            root=self.svg(target);ribbons=root.find(f"{SVG}g[@id='target-ribbons']")
            self.assertTrue(len(ribbons)>0)
            for n in ribbons:
                self.assertEqual(n.get('data-ribbon-target'),target)
                oid=n.get('data-query-occurrence')
                self.assertIn(('_e3' if target=='exon_3' else '_e4' if target=='exon_4' else '_q'),oid)

    def test_split_ribbons_do_not_cover_intervening_DNA(self):
        nodes=self.svg().find(f"{SVG}g[@id='target-ribbons']")
        a_ranges={(int(n.get('data-query-start0')),int(n.get('data-query-end0'))) for n in nodes if n.get('data-query-occurrence')=='A_e3'}
        self.assertEqual(a_ranges,{(700,750),(765,820)})

    def test_unknown_species_does_not_have_confirmed_ribbon(self):
        for n in self.svg().find(f"{SVG}g[@id='target-ribbons']"):
            self.assertFalse(n.get('data-query-occurrence').startswith('D_'))
            self.assertFalse(n.get('data-target-occurrence').startswith('D_'))
        ds=[n for n in self.svg().find(f"{SVG}g[@id='target-highlights']") if n.get('data-species')=='D']
        self.assertTrue(ds);self.assertEqual(ds[0].get('fill'),'none')

    def test_missing_E4_is_not_drawn_present(self):
        root=self.svg('exon_4')
        self.assertFalse(any(n.get('data-species')=='D' for n in root.find(f"{SVG}g[@id='target-highlights']")))
        self.assertIn('absent',' '.join(root.itertext()))

    def test_negative_strand_coordinates_retained(self):
        root=self.svg();n=next(n for n in root.find(f"{SVG}g[@id='target-highlights']") if n.get('data-species')=='C')
        self.assertEqual((int(n.get('data-start0')),int(n.get('data-end0'))),(580,700))

    def test_junction_is_cuts_not_intron_body_ribbon(self):
        root=self.svg('splice_position')
        self.assertEqual(len(root.find(f"{SVG}g[@id='target-ribbons']")),0)
        highlights=root.find(f"{SVG}g[@id='target-highlights']")
        self.assertEqual(len(highlights),6)
        self.assertTrue(all(n.tag==SVG+'line' for n in highlights))

    def test_structure_view_has_no_ribbons(self):
        self.assertEqual(len(self.svg(view='structure').find(f"{SVG}g[@id='target-ribbons']")),0)

    def test_selected_targets_keep_catalogue_color(self):
        all_targets=load_targets(self.result,manifest=self.root/'targets.tsv')[0]
        one=load_targets(self.result,selectors=['exon_4'],manifest=self.root/'targets.tsv')[0]
        self.assertEqual(one[0].color,next(t.color for t in all_targets if t.target_id=='exon_4'))

    def test_excluded_scope_still_uses_identical_raw_background(self):
        import shutil
        high=self.root/'high';high.mkdir()
        for name in ('element_correspondence.tsv','splice_boundary_correspondence.tsv'):shutil.copyfile(self.prepared/name,high/name)
        infer_phylogeny(self.inp,high,structural_site_matrix_path=self.prepared/'structural_site_matrix.tsv',analysis_range='high-coverage')
        rows=visualize_results(self.inp,high,self.root/'high_fig',target_manifest=self.root/'targets.tsv')
        r=next(r for r in rows if r['target_id']=='intron_segment' and r['view']=='phylogeny')
        self.assertEqual(r['selected_for_analysis'],'0')
        root=ET.parse(self.root/'high_fig'/r['path']).getroot()
        self.assertIn('excluded from this fit',' '.join(root.itertext()))
        self.assertEqual(ET.tostring(root.find(f"{SVG}g[@id='native-background']")),ET.tostring(self.svg().find(f"{SVG}g[@id='native-background']")))

    def test_plot_does_not_reinfer(self):
        with patch('intraphy.phylogeny.infer_phylogeny',side_effect=AssertionError('must not infer')):
            visualize_results(self.inp,self.result,self.root/'again',targets=[FAMILY+'/exon_presence/E3'])

    def test_manifest_rejects_mixed_targets(self):
        p=self.root/'bad.tsv'
        write_rows(p,[{'target_id':'X','family_id':FAMILY,'layer':'exon_presence','site_id':'E3'}, {'target_id':'X','family_id':FAMILY,'layer':'exon_presence','site_id':'E4'}])
        with self.assertRaises(ValueError):load_targets(self.result,manifest=p)

    def test_manifest_rejects_duplicate_colors(self):
        p=self.root/'bad.tsv'
        write_rows(p,[{'target_id':s,'family_id':FAMILY,'layer':'exon_presence','site_id':s,'color':'#008877'} for s in ('E3','E4')])
        with self.assertRaises(ValueError):load_targets(self.result,manifest=p)

    def test_default_catalogue_requires_unambiguous_selection(self):
        with self.assertRaises(ValueError):load_targets(self.result,selectors=['E3'])
        self.assertEqual(len(load_targets(self.result,selectors=[FAMILY+'/exon_presence/E3'])[0]),1)

    def test_label_is_XML_escaped(self):
        rows=read_tsv(self.root/'targets.tsv');rows[0]['label']='E3 <special> & "quoted"'
        write_rows(self.root/'special.tsv',rows)
        out=self.root/'special';manifest=visualize_results(self.inp,self.result,out,target_manifest=self.root/'special.tsv')
        root=ET.parse(out/manifest[0]['path']).getroot();self.assertIn(rows[0]['label'],' '.join(root.itertext()))

    def test_manual_subinterval_not_whole_exon(self):
        p=self.root/'manual.tsv'
        write_rows(p,[{'target_id':'portion','family_id':FAMILY,'layer':'exon_presence','site_id':'E3','element_id':'E3',
            'species':'A','gene_copy_id':'A_gene','contig':'chrA','start':721,'end':740,'strand':'+','occurrence_id':'A_e3'}])
        targets,manual,_,members=load_targets(self.result,manifest=p)
        ints=resolve_intervals(targets[0],read_tsv(self.inp/'segment_occurrences.tsv'),members,[],manual['portion'])
        self.assertEqual([i.interval for i in ints],[Interval0(720,740)])

    def test_layer_does_not_claim_function(self):
        from intraphy.reporting.target_model import LAYER_LABELS
        self.assertEqual(LAYER_LABELS['exon_role'],'Annotated exonic use')
        self.assertNotIn('function',' '.join(LAYER_LABELS.values()).lower())



class IntegratedContractsV017(unittest.TestCase):
    def test_pathless_actual_mapping_stage_preserves_DNA(self):
        from intraphy.mapping.chains import _apply_ordered_candidate_chains
        occurrences=[];paths=[];rows=[]
        for sp in ('A','B'):
            for i,name in enumerate(('a','u','c')):
                occ={'occurrence_id':sp+name,'family_id':'f','species':sp,'gene_copy_id':'g'+sp,
                     'contig':'chr'+sp,'start':str(i*200+1),'end':str(i*200+100),'strand':'+','role':'unknown' if name=='u' else 'exon'}
                occurrences.append(occ)
                if name!='u':paths.append({**occ,'transcript_id':'tx'+sp,'path_rank':i+1})
        for name in ('a','u','c'):
            rows.append({'match_id':name,'query_occurrence_id':'A'+name,'subject_occurrence_id':'B'+name,
                  'match_status':'mapped','enumeration_complete':1,'short_context_route':'not_used',
                  '_candidate_records':[{'candidate_id':name,'accepted':1,'query_start0':0,'query_end0':100,
                    'target_start0':0,'target_end0':100,'score':200.,'score_scheme':'nt_blastn_v1','strand':'+',
                    'aligned_blocks':[(1,100,1,100)]}]})
        _apply_ordered_candidate_chains(rows,{r['occurrence_id']:r for r in occurrences},occurrences,paths)
        row=rows[1]
        self.assertEqual(row['membership_edge_eligible'],1)
        self.assertEqual(row['position_edge_eligible'],1)
        self.assertEqual(row['genomic_retained_candidate_ids'],'u')
        self.assertEqual(row['transcript_retained_candidate_ids'],'NA')
        self.assertEqual(row['correspondence_channels'],'DNA_only_no_observed_transcript_path')
        self.assertNotEqual(row['flanking_anchor_status'],'ordered_double_sided_homologous_flanks_same_path')

    def test_unknown_mid_exon_not_converted_to_jump_in_state_builder(self):
        from intraphy.observations.junctions import _junction_site_rows
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);inp=root/'in';out=root/'out';inp.mkdir();out.mkdir()
            occurrences=[];paths=[];elements=[]
            for rank,(name,start,end,role) in enumerate([('a',1,30,'exon'),('i1',31,60,'intron'),('b',61,90,'exon'),('i2',91,120,'intron'),('c',121,150,'exon')],1):
                occ={'occurrence_id':name,'family_id':'f','species':'A','gene_copy_id':'g','contig':'chr','start':str(start),'end':str(end),'strand':'+','role':role,'presence_status':'present','phase':'0','source_feature_id':name}
                occurrences.append(occ);paths.append({**occ,'transcript_id':'tx','path_rank':rank,'coding_role':'coding' if role=='exon' else 'intronic'})
                if name in ('a','c'):elements.append({'family_id':'f','element_id':name,'occurrence_id':name,'membership_call':'core_member','position_edge_eligible':'1','reference_occurrence_id':name,'element_class':'exon_like'})
            write_rows(inp/'transcript_paths.tsv',paths)
            rows=_junction_site_rows(inp,out,occurrences,elements,{'f'},{'f':{'A'}})
            self.assertEqual(rows,[])
            records=read_tsv(out/'splice_boundary_correspondence.tsv')
            self.assertTrue(any(r['unavailable_reason']=='unresolved_intervening_annotated_feature' for r in records))

    def test_ambiguous_membership_cannot_be_solid_or_ribbon(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp,prepared=build_example(tmp)
            occ=read_tsv(inp/'segment_occurrences.tsv');members=read_tsv(prepared/'element_correspondence.tsv')
            for row in members:
                if row['element_id']=='E3':row['position_edge_eligible']='0'
            target=FigureTarget('e3',FAMILY,'exon_presence','E3','E3','third','#0072B2')
            intervals=resolve_intervals(target,occ,members,[])
            self.assertTrue(intervals);self.assertTrue(all(i.coordinate_evidence=='ambiguous_membership_subinterval' for i in intervals))
            self.assertEqual(target_ribbons(target,intervals,occ,inp,prepared,{'A','B','C'}),[])

    def test_adding_context_features_does_not_add_inference_characters(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp,prepared=build_example(tmp);matrix=read_tsv(prepared/'structural_site_matrix.tsv')
            before=len(matrix)
            rows=read_tsv(inp/'raw_gene_features.tsv');repeat=dict(rows[-1]);repeat['id']='annotation_subdivision_only';rows.append(repeat)
            write_rows(inp/'raw_gene_features.tsv',rows)
            prepared_matrix=prepare_observation_matrix(inp,Path(tmp)/'out',prepared/'structural_site_matrix.tsv')
            self.assertEqual(len(prepared_matrix.full_rows),before)

    def test_unavailable_presence_does_not_harden_role(self):
        rows=[observation('U','A',1),observation('U','B',1),observation('U','A',None,'exon_role'),observation('U','B',None,'exon_role')]
        full,_,_,summaries=assess_scope(rows)
        self.assertEqual([r['state'] for r in full if r['layer']=='exon_role'],['unknown','unknown'])

    def test_unstranded_context_annotation_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp,prepared=build_example(tmp)
            raw=read_tsv(inp/'raw_gene_features.tsv')
            for row in raw:
                if row['type']=='repeat_region':row['strand']='.'
            write_rows(inp/'raw_gene_features.tsv',raw)
            scene=build_scene(inp,FAMILY,read_tsv(inp/'segment_occurrences.tsv'))
            self.assertTrue(any(f.kind=='repeat_region' and f.strand=='.' for f in scene.features))
            self.assertEqual(scene.track_for('C','C_gene').strand,'-')

    def test_zero_callable_does_not_become_global_absence(self):
        rows=[observation('X',sp,None) for sp in SPECIES]
        _,_,_,summary=assess_scope(rows)
        self.assertEqual(summary[0]['state0_n'],0);self.assertEqual(summary[0]['callable_n'],0)
        self.assertEqual(summary[0]['called_mrca'],'NA')

    def test_budget_truncation_reason_is_preserved(self):
        rows=[observation('X',sp,None,observation_reason='candidate_search_truncated_budget') for sp in SPECIES]
        full,_,audit,summary=assess_scope(rows)
        self.assertEqual(summary[0]['unknown_n'],4)
        self.assertTrue(all(r['reason']=='candidate_search_truncated_budget' for r in audit))


if __name__ == '__main__':
    unittest.main()
