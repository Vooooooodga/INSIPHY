#!/usr/bin/env python3
"""Build a transparent synthetic *prepared-input* example, not an accuracy benchmark.

The homology blocks/character states below are supplied inputs to exercise scope,
fixed-tree inference and plotting. The aligner did not discover their truth.
Separate unit tests exercise raw extraction, candidate chains and short matches.
"""
from pathlib import Path
import argparse
import json
import shutil

from intraphy.storage.tabular import write_tsv
from intraphy.observations.schema import write_structural_site_matrix

FAMILY = 'SYNTHETIC_GENE'
SPECIES = ('A', 'B', 'C', 'D')


def write_rows(path, rows, fields=None):
    fields = fields or sorted({key for row in rows for key in row})
    write_tsv(path, rows, fields)


def observation(site, species, state, layer='exon_presence', applicability=None, **kwargs):
    labels = ('not_exonic','exonic') if layer=='exon_role' else ('absent','present')
    row={'family_id':FAMILY, 'layer':layer, 'site_id':site,'species':species,
         'state':labels[state] if state in (0,1) else 'unknown','state_0':labels[0], 'state_1':labels[1],
         'schema_version':'3', 'annotation_view':'view_independent' if layer=='exon_presence' else 'repertoire',
         'transcript_scope':'gene_locus' if layer=='exon_presence' else 'annotated_transcript_repertoire',
         'applicability':applicability or ('applicable' if state in (0,1) else 'undetermined'),
         'observation_mask':'observed' if state in (0,1) else 'missing',
         'observation_reason':'supplied_synthetic_state' if state in (0,1) else 'unresolved_synthetic_homology',
         'evidence':'synthetic_prepared_input_not_real_biological_validation',
         'discovery_rule':'observed_at_least_one','linked_group_id':'NA','site_kind':'homologous_interval'}
    row.update(kwargs)
    return row


def build_example(root):
    root=Path(root); inp=root/'input'; prepared=root/'prepared'; inp.mkdir(parents=True,exist_ok=True);prepared.mkdir(exist_ok=True)
    tree=[{'node_id':'root','parent_id':'','label':'root','branch_length':'0'},
          {'node_id':'AB','parent_id':'root','label':'AB','branch_length':'.2'},
          {'node_id':'CD','parent_id':'root','label':'CD','branch_length':'.2'}]
    tree += [{'node_id':sp,'parent_id':'AB' if sp in 'AB' else 'CD','label':sp,'branch_length':'.1'} for sp in SPECIES]
    write_rows(inp/'species_tree.tsv',tree,['node_id','parent_id','label','branch_length'])
    occurrences=[];paths=[];raw=[];loci=[];members=[];sequences={}
    base=[('e1',101,250,'exon','E1'),('e2',401,550,'exon','E2'),
          ('q',561,680,'intron','Q'),('e3',701,820,'exon','E3'),('e4',1001,1180,'exon','E4')]
    for sp in SPECIES:
        strand='-' if sp=='C' else '+';copy=f'{sp}_gene';tx=f'{sp}_tx1';contig=f'chr{sp}'
        def pos(a,b):return (1401-b,1401-a) if strand=='-' else (a,b)
        loci.append({'species':sp,'gene_copy_id':copy,'contig':contig,'strand':strand,'annotation_start':101,'annotation_end':1300})
        raw.append({'family_id':FAMILY,'species':sp,'gene_copy_id':copy,'seqid':contig,'strand':strand,'type':'gene','id':copy,'parent':'NA','start':101,'end':1300,'source':'synthetic_native_annotation'})
        items=list(base)
        if sp=='B':items=[r for r in items if r[0]!='e3']+[('e3a',701,750,'exon','E3'),('e3b',766,820,'exon','E3')]
        if sp in {'C','D'}:items=[r for r in items if r[0]!='e4']
        items.sort(key=lambda x:x[1])
        for rank,(name,a,b,role,element) in enumerate(items,1):
            x,y=pos(a,b);oid=f'{sp}_{name}'
            # Q has reliable coordinates but no transcript membership: DNA-only candidate.
            txid='' if name=='q' else tx
            occ={'occurrence_id':oid,'family_id':FAMILY,'species':sp,'gene_copy_id':copy,'transcript_id':txid,
                 'segment_id':name,'contig':contig,'start':x,'end':y,'strand':strand,'role':role,'role_set':role,
                 'presence_status':'present','phase':'0','source_feature_id':oid,'transcript_order':rank,
                 'evidence':'synthetic_prepared_coordinate','boundary_state':'known','coding_status':'coding' if role=='exon' else 'noncoding'}
            occurrences.append(occ); sequences[oid]='ACGT'*((b-a+4)//4)
            if txid:paths.append({**occ,'path_rank':rank,'path_status':'annotated_transcript_path','path_role':'exonic','coding_role':'coding','partial_start':'0','partial_end':'0'})
            raw.append({'family_id':FAMILY,'species':sp,'gene_copy_id':copy,'seqid':contig,'strand':strand,'type':role,'id':oid,'parent':txid or 'DNA context','start':x,'end':y,'source':'synthetic_native_annotation'})
            if role=='exon':
                cx,cy=pos(a+12,b-8)
                raw.append({'family_id':FAMILY,'species':sp,'gene_copy_id':copy,'seqid':contig,'strand':strand,'type':'CDS','id':oid+'_cds','parent':tx,'start':cx,'end':cy,'source':'synthetic_native_annotation'})
            if name=='e1':
                ux,uy=pos(a,a+11)
                raw.append({'family_id':FAMILY,'species':sp,'gene_copy_id':copy,'seqid':contig,'strand':strand,'type':'five_prime_UTR','id':oid+'_utr','parent':tx,'start':ux,'end':uy,'source':'synthetic_native_annotation'})
            members.append({'family_id':FAMILY,'element_id':element,'occurrence_id':oid,'membership_call':'core_member',
                'element_class':'exon_like' if role=='exon' else 'candidate_source','member_start':1,'member_end':b-a+1,
                'matched_blocks':f'1-{b-a+1}','position_edge_eligible':'1','membership_status':'confirmed',
                'reference_occurrence_id':'A_'+ ('e3' if element=='E3' else name)})
        # Antisense nested annotation and an RNA/TE tag remain gray in every view.
        for label,a,b,kind,anti in [('antisense',725,770,'exon',True),('snoRNA',900,930,'ncRNA',False),('repeat',320,355,'repeat_region',False)]:
            x,y=pos(a,b)
            raw.append({'family_id':FAMILY,'species':sp,'gene_copy_id':copy,'seqid':contig,
                        'strand':('+' if strand=='-' else '-') if anti else strand,'type':kind,'id':sp+'_'+label,
                        'parent':'other_gene' if anti else 'other annotations','start':x,'end':y,'source':'synthetic_context'})
        # A second annotated path with a shorter valid first exon.
        x,y=pos(155,250)
        raw.append({'family_id':FAMILY,'species':sp,'gene_copy_id':copy,'seqid':contig,'strand':strand,'type':'exon',
                    'id':sp+'_alternative_first','parent':sp+'_tx2','start':x,'end':y,'source':'synthetic_alternative_transcript'})
    write_rows(inp/'segment_occurrences.tsv',occurrences)
    write_rows(inp/'transcript_paths.tsv',paths)
    write_rows(inp/'raw_gene_features.tsv',raw)
    write_rows(inp/'gene_loci.tsv',loci)
    (inp/'segment_sequences.fasta').write_text(''.join(f'>{r["occurrence_id"]}\n{sequences[r["occurrence_id"]][:int(r["end"])-int(r["start"])+1]}\n' for r in occurrences))
    write_rows(prepared/'element_correspondence.tsv',members)
    matches=[]
    for a,b in [('A_e3','B_e3a'),('A_e3','B_e3b'),('B_e3a','C_e3'),('B_e3b','C_e3'),('C_e3','D_e3'),('A_e4','B_e4'),('A_q','B_q'),('B_q','C_q'),('C_q','D_q')]:
        if a=='A_e3' and b=='B_e3a': blocks='1-50:1-50'
        elif a=='A_e3' and b=='B_e3b':blocks='66-120:1-55'
        elif a=='B_e3a':blocks='1-50:1-50'
        elif a=='B_e3b':blocks='1-55:66-120'
        else:
            n=min(int(next(r for r in occurrences if r['occurrence_id']==oid)['end'])-int(next(r for r in occurrences if r['occurrence_id']==oid)['start'])+1 for oid in (a,b))
            blocks=f'1-{n}:1-{n}'
        matches.append({'match_id':a+'__'+b,'query_occurrence_id':a,'subject_occurrence_id':b,'match_status':'mapped',
                        'alignment_strand':'+','matched_blocks':blocks,'alignment_backend':'synthetic_supplied_blocks',
                        'correspondence_basis':'DNA','position_edge_eligible':'1','enumeration_complete':'1'})
    write_rows(inp/'segment_matches.tsv',matches)
    states=[]
    for site,values,layer in [('E1',[1,1,1,1],'exon_presence'),('E3',[1,1,1,None],'exon_presence'),
        ('E4',[1,1,0,0],'exon_presence'),('Q',[1,1,None,None],'exon_presence'),
        ('E3',[1,1,None,None],'exon_role'),('J_E2_E3',[1,1,1,None],'splice_junction')]:
        states.extend(observation(site,sp,value,layer=layer) for sp,value in zip(SPECIES,values))
    write_structural_site_matrix(prepared/'structural_site_matrix.tsv',states)
    boundaries=[]
    for sp in SPECIES[:3]:
        boundaries.append({'family_id':FAMILY,'species':sp,'gene_copy_id':f'{sp}_gene','transcript_id':f'{sp}_tx1','site_id':'J_E2_E3',
           'element_id':'E2;E3','left_occurrence_id':sp+'_e2','right_occurrence_id':sp+('_e3a' if sp=='B' else '_e3'),
           'position_edge_eligible':'1'})
    write_rows(prepared/'splice_boundary_correspondence.tsv',boundaries)
    target_rows=[{'target_id':'exon_3','family_id':FAMILY,'layer':'exon_presence','site_id':'E3','label':'Exon 3 homologous DNA'},
                 {'target_id':'exon_4','family_id':FAMILY,'layer':'exon_presence','site_id':'E4','label':'Exon 4 homologous DNA'},
                 {'target_id':'intron_segment','family_id':FAMILY,'layer':'exon_presence','site_id':'Q','label':'Unannotated intronic DNA segment'},
                 {'target_id':'splice_position','family_id':FAMILY,'layer':'splice_junction','site_id':'J_E2_E3','label':'Splice position between E2 and E3'}]
    write_rows(root/'targets.tsv',target_rows)
    (root/'README.md').write_text('''# Synthetic 0.17 prepared-input example\n\nThis is **not real organism data and not an alignment-accuracy benchmark**.\nA/B/C/D are fictitious species. C is on the negative strand; B has split E3 pieces;\nD has an unresolved E3 correspondence. C and D have explicit E4 absences;
E4 therefore has two equally parsimonious gain/loss histories on the supplied tree.\nA shorter annotated transcript, antisense exon, ncRNA and repeat are gray context.\nQ has DNA coordinates without transcript membership. Its low-callability character\nremains in the full catalogue and background, even when excluded from a high-coverage fit.\nThe provided homology blocks and observations test downstream contracts, not inferred truth.\n\nThe builder runs fixed-tree parsimony and the actual grouped renderer.\nOriginal coordinates and all supplied features must be identical in every group.\n''')
    return inp,prepared


def run_example(root):
    from intraphy.phylogeny import infer_phylogeny
    from intraphy.visualize import visualize_results
    root=Path(root);inp,prepared=build_example(root)
    for scope in ('all','high-coverage'):
        out=root/('results_'+scope);out.mkdir(exist_ok=True)
        for name in ('element_correspondence.tsv','splice_boundary_correspondence.tsv'):
            shutil.copyfile(prepared/name,out/name)
        infer_phylogeny(inp,out,model='parsimony',structural_site_matrix_path=prepared/'structural_site_matrix.tsv',analysis_range=scope)
        visualize_results(inp,out,root/('figures_'+scope),target_manifest=root/'targets.tsv')
    return root


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    run_example(args.output)
    print(args.output)
