"""One target per figure group, unchanged native background, target-only ribbons.

Consumes stored results; no phylogeny inference and no annotation modification.
"""
from __future__ import annotations
from collections import defaultdict
from html import escape
import json
from pathlib import Path
import re

from insiphy.coordinates import Interval0
from insiphy.reporting.target_model import (LAYER_LABELS, load_targets, resolve_intervals, target_ribbons)
from insiphy.reporting.target_scene import (WIDTH, TRACK_LEFT, TRACK_RIGHT, attr, text, build_scene, scene_data)
from insiphy.run_result import result_model
from insiphy.storage.tabular import read_tsv,write_tsv


def _safe_name(value):
    return re.sub(r'[^A-Za-z0-9_.-]+','_',value).strip('._')[:150] or 'target'


def _state(target,matrix):
    return {row['species']:row for row in matrix if
            (row.get('family_id'),row.get('layer'),row.get('site_id'))==
            (target.family_id,target.layer,target.site_id)}


def _events(result_dir,target, *, model=None, rows=None):
    model=model or result_model(result_dir)
    if model=='parsimony':
        filename='branch_structural_events.tsv'
    elif model in {'er-ard','foreground'}:
        filename='structural_changes.tsv'
    else:
        return []
    rows=rows if rows is not None else read_tsv(Path(result_dir)/filename,optional=True)
    rows=[row for row in rows if (row.get('family_id'),row.get('layer'),row.get('site_id'))==
          (target.family_id,target.layer,target.site_id)]
    if model=='parsimony':
        return [row for row in rows if row.get('placement_status') in {'required','possible'}]
    return [row for row in rows if str(row.get('posterior_available','')).lower() in {'true','1','yes'}]


def _interval_lane_y(scene,item,occurrences):
    track=scene.track_for(item.species,item.gene_copy_id)
    if not track:return None,None
    occ=next((row for row in occurrences if row['occurrence_id']==item.occurrence_id),{})
    tx=occ.get('transcript_id','')
    native_ids={occ.get('source_feature_id'),item.occurrence_id}-{None,''}
    matching_lanes=[f.lane for f in scene.features if (f.species,f.copy)==(item.species,item.gene_copy_id)
                    and f.feature_id in native_ids and f.strand==item.strand]
    lane=tx if tx in dict(track.lane_y) else (matching_lanes[0] if matching_lanes else 'genomic DNA')
    y=dict(track.lane_y).get(lane,track.lane_y[-1][1])
    return track,y


def _render(scene,target,intervals,matrix,pairs,events,mode,occurrences,scope):
    states=_state(target,matrix)
    body=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{scene.height:.0f}" '
          f'viewBox="0 0 {WIDTH} {scene.height:.0f}" data-target="{attr(target.target_id)}" '
          f'data-family="{attr(target.family_id)}" data-layer="{attr(target.layer)}" '
          f'data-site="{attr(target.site_id)}" data-target-color="{target.color}" data-view="{mode}">',
          '<defs><pattern id="native-gap" width="6" height="6" patternUnits="userSpaceOnUse">'
          '<path d="M0,6 L6,0" stroke="#AAAAAA" stroke-width="1"/></pattern></defs>',
          '<rect width="100%" height="100%" fill="white"/>',
          text(24,28,f'{target.family_id} | {target.label}',18,weight='bold'),
          text(24,50,f'{LAYER_LABELS[target.layer]} — {mode} view',12),
          text(24,70,'Full supplied native background is identical across targets. Gray does not mean unanalysed, absent or conserved.',10),
          text(24,88,'Tracks retain genomic gaps; each locus has its own linear 5′→3′ scale. Ribbons show only this target’s actual matched bases.',10),
          f'<g id="target-legend" data-target="{attr(target.target_id)}" data-color="{target.color}">',
          f'<rect x="24" y="102" width="18" height="11" fill="{target.color}"/>',
          text(49,112,f'Target: {target.label}',10),
          '<rect x="375" y="102" width="18" height="11" fill="#CCCCCC"/>',
          text(400,112,'Other supplied annotations',10),
          f'<rect x="650" y="102" width="18" height="11" fill="none" stroke="{target.color}" stroke-dasharray="3 2"/>',
          text(675,112,'Located range; state unresolved / parent range only',10),
          text(1130,112,'Colors identify objects, not event probabilities',10),'</g>',
          '<g id="target-ribbons">']
    if mode in {'ribbons','phylogeny'}:
        order={sp:i for i,(sp,_) in enumerate(scene.species_y)}
        for a,b,ia,ib,mid in pairs:
            # Only join adjacent *panel* species, never jump across an unmeasured
            # intermediate track and visually conceal its uncertainty.
            if abs(order[a['species']]-order[b['species']])!=1:continue
            ta=scene.track_for(a['species'],a['gene_copy_id']);tb=scene.track_for(b['species'],b['gene_copy_id'])
            if not ta or not tb:continue
            from insiphy.reporting.target_model import TargetInterval
            _,ya=_interval_lane_y(scene,TargetInterval(a['species'],a['gene_copy_id'],a['contig'],ia,a['strand'],a['occurrence_id']),occurrences)
            _,yb=_interval_lane_y(scene,TargetInterval(b['species'],b['gene_copy_id'],b['contig'],ib,b['strand'],b['occurrence_id']),occurrences)
            x1,x2=sorted((ta.x(ia.start0),ta.x(ia.end0)));u1,u2=sorted((tb.x(ib.start0),tb.x(ib.end0)))
            midy=(ya+yb)/2
            body.append(f'<path d="M{x1:.3f},{ya+8:.3f} C{x1:.3f},{midy:.3f} {u1:.3f},{midy:.3f} {u1:.3f},{yb-8:.3f} '
                        f'L{u2:.3f},{yb-8:.3f} C{u2:.3f},{midy:.3f} {x2:.3f},{midy:.3f} {x2:.3f},{ya+8:.3f} Z" '
                        f'fill="{target.color}" fill-opacity=".18" stroke="{target.color}" stroke-opacity=".4" stroke-width=".6" '
                        f'data-ribbon-target="{attr(target.target_id)}" data-match="{attr(mid)}" '
                        f'data-query-occurrence="{attr(a["occurrence_id"])}" data-target-occurrence="{attr(b["occurrence_id"])}" '
                        f'data-query-start0="{ia.start0}" data-query-end0="{ia.end0}" data-target-start0="{ib.start0}" data-target-end0="{ib.end0}"/>')
    body.append('</g>')
    body.append(scene.background_svg)  # gray annotations remain gray over transparent ribbons
    body.append('<g id="target-highlights">')
    for item in intervals:
        track,y=_interval_lane_y(scene,item,occurrences)
        if track is None:continue
        if item.contig!=track.contig or item.strand!=track.strand:
            raise ValueError('highlight coordinate space disagrees with native locus')
        if item.interval.start0<track.start-1 or item.interval.end0>track.end:
            raise ValueError(f'target outside fixed native background: {target.target_id}')
        row=states.get(item.species,{})
        state=row.get('state','unknown');mask=row.get('observation_mask','missing')
        if target.layer=='exon_presence' and state==row.get('state_0') and mask!='missing':
            continue
        known=state!='unknown' and mask!='missing' and row.get('applicability')=='applicable'
        located=item.coordinate_evidence in {'actual_membership_subinterval','explicit_annotated_junction','user_specified_native_interval'}
        dashed=not known or not located
        attrs=(f'data-highlight-target="{attr(target.target_id)}" data-species="{attr(item.species)}" '
               f'data-copy="{attr(item.gene_copy_id)}" data-start0="{item.interval.start0}" data-end0="{item.interval.end0}" '
               f'data-state="{attr(state)}" data-coordinate-evidence="{attr(item.coordinate_evidence)}"')
        x1,x2=sorted((track.x(item.interval.start0),track.x(item.interval.end0)))
        if item.kind in {'donor','acceptor'}:
            body.append(f'<line x1="{x1:.3f}" x2="{x1:.3f}" y1="{y-12:.3f}" y2="{y+12:.3f}" stroke="{target.color}" '
                        f'stroke-width="2.5" stroke-dasharray="{"3 2" if dashed else "none"}" {attrs}><title>{item.kind}</title></line>')
        else:
            body.append(f'<rect x="{x1:.3f}" y="{y-10:.3f}" width="{max(.8,x2-x1):.3f}" height="20" '
                        f'fill="{"none" if dashed else target.color}" fill-opacity=".60" stroke="{target.color}" '
                        f'stroke-width="1.7" stroke-dasharray="{"4 2" if dashed else "none"}" {attrs}>'
                        f'<title>{attr(item.coordinate_evidence)}; {attr(state)}</title></rect>')
    body.append('</g>')
    body.append(text(1395,136,'Target state / scope',10,weight='bold'))
    for sp,y in scene.species_y:
        row=states.get(sp,{})
        state=row.get('state','unmeasured')
        if row.get('applicability')=='inapplicable':state='inapplicable'
        elif row.get('observation_mask')=='missing':state='unknown'
        body.append(text(1395,y-2,state,11,fill=target.color if state not in {'unknown','unmeasured','inapplicable'} else '#777777'))
        body.append(text(1395,y+13,('excluded from this fit' if str(scope.get('selected_for_analysis','1'))=='0' else 'scope: '+str(scope.get('analysis_range','unrecorded'))),8))
    if mode=='phylogeny' and scene.tree:
        xy={node:(x,y) for node,x,y in scene.node_xy}
        label_to_node={scene.tree.label[n]:n for n in scene.tree.parent}
        by_branch=defaultdict(list)
        for event in events:
            branch=event.get('branch_scope','')
            by_branch[branch].append(event)
        for branch,records in sorted(by_branch.items()):
            if '->' not in branch:continue
            pa,ch=(part.strip() for part in branch.split('->',1))
            parent=pa if pa in xy else label_to_node.get(pa);child=ch if ch in xy else label_to_node.get(ch)
            if parent not in xy or child not in xy or scene.tree.parent[child]!=parent:continue
            x=(xy[parent][0]+xy[child][0])/2;y=xy[child][1]
            required=any(r.get('placement_status')=='required' for r in records)
            label='; '.join(sorted({r.get('event_type',r.get('structural_change_type','conditional change')) for r in records}))
            body.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="4.5" fill="{target.color if required else "white"}" '
                        f'stroke="{target.color}" stroke-width="1.5" data-event-target="{attr(target.target_id)}" '
                        f'data-branch="{attr(branch)}"><title>{attr(label)}; {"required in all parsimonious histories" if required else "possible / model-conditional"}</title></circle>')
        body.append(text(24,scene.height-56,'Filled branch marker: required parsimony event. Open: possible or model-conditional. No invented ancestor exon track.',10))
    if target.layer=='splice_junction':
        note='A splice connection is represented by its two cuts. No intron-body ribbon is implied without homologous DNA evidence.'
    elif not pairs:
        note='No target-specific resolved matched blocks are available; no ribbons are fabricated from shared element labels.'
    else:
        note='Ribbons are matched DNA columns, not causal arrows or reconstructed evolutionary events.'
    body.extend([text(24,scene.height-36,note,10),text(24,scene.height-17,scene.annotation_source,9),'</svg>'])
    return ''.join(body)


def visualize_target_groups(input_dir,result_dir,output_dir,*,selectors=(),target_manifest=None):
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    targets,manual,matrix,members=load_targets(result_dir,selectors,target_manifest)
    if not targets:
        raise ValueError('no figure target is defined; provide a target manifest or structural matrix')
    occurrences=read_tsv(Path(input_dir)/'segment_occurrences.tsv')
    # All targets consume the same immutable run snapshot; do not reread large
    # event/match files for every SVG panel.
    from insiphy.reporting.inputs import ribbon_match_evidence
    match_data = ribbon_match_evidence(input_dir,result_dir,occurrences)[0]
    model = result_model(result_dir)
    event_file = 'branch_structural_events.tsv' if model=='parsimony' else 'structural_changes.tsv'
    event_index = defaultdict(list)
    for row in read_tsv(Path(result_dir)/event_file,optional=True):
        event_index[(row.get('family_id'),row.get('layer'),row.get('site_id'))].append(row)
    matrix_index = defaultdict(list)
    for row in matrix:
        matrix_index[(row.get('family_id'),row.get('layer'),row.get('site_id'))].append(row)
    boundaries=read_tsv(Path(result_dir)/'splice_boundary_correspondence.tsv',optional=True)
    eligibility=read_tsv(Path(result_dir)/'structural_character_eligibility.tsv',optional=True)
    scopes={(row['family_id'],row['layer'],row['site_id']):row for row in eligibility}
    scenes={family:build_scene(input_dir,family,occurrences) for family in sorted({t.family_id for t in targets})}
    for i,(family,scene) in enumerate(scenes.items()):
        (output/f'native_background_{i+1:03d}.json').write_text(json.dumps(scene_data(scene),indent=2,ensure_ascii=False)+'\n')
    manifest=[];interval_table=[];html=[]
    for index,target in enumerate(targets,1):
        directory=output/f'{index:04d}_{_safe_name(target.target_id)}';directory.mkdir(exist_ok=True)
        scene=scenes[target.family_id]
        intervals=resolve_intervals(target,occurrences,members,boundaries,manual.get(target.target_id,()))
        target_matrix=matrix_index[(target.family_id,target.layer,target.site_id)]
        states=_state(target,target_matrix)
        callable_species={sp for sp,row in states.items() if row.get('state')!='unknown'
                          and row.get('observation_mask')!='missing' and row.get('applicability')=='applicable'}
        pairs=target_ribbons(target,intervals,occurrences,input_dir,result_dir,callable_species,match_data=match_data)
        events=_events(result_dir,target,model=model,rows=event_index[(target.family_id,target.layer,target.site_id)])
        scope=scopes.get((target.family_id,target.layer,target.site_id),{})
        for mode in ('structure','ribbons','phylogeny'):
            path=directory/f'{mode}.svg'
            path.write_text(_render(scene,target,intervals,target_matrix,pairs,events,mode,occurrences,scope),encoding='utf-8')
            manifest.append({'target_id':target.target_id,'family_id':target.family_id,'layer':target.layer,
                'site_id':target.site_id,'label':target.label,'color':target.color,'view':mode,
                'path':str(path.relative_to(output)),'type':'svg','intervals':len(intervals),
                'description':'One target; full fixed native background; target-specific ribbons only',
                'selected_for_analysis':scope.get('selected_for_analysis','unrecorded')})
        for item in intervals:
            interval_table.append({'target_id':target.target_id,'species':item.species,'gene_copy_id':item.gene_copy_id,
                'contig':item.contig,'start0':item.interval.start0,'end0':item.interval.end0,'strand':item.strand,
                'occurrence_id':item.occurrence_id or 'NA','kind':item.kind,'coordinate_evidence':item.coordinate_evidence})
        rel=directory.relative_to(output).as_posix()
        html.append(f'<section><h2 style="border-left:8px solid {target.color};padding-left:12px">{escape(target.label)}</h2>'
                    f'<p>{escape(target.family_id)} · {escape(LAYER_LABELS[target.layer])} · {escape(target.site_id)}</p>'
                    f'<p><a href="{rel}/structure.svg">Structure</a> · <a href="{rel}/ribbons.svg">Target ribbons</a> · '
                    f'<a href="{rel}/phylogeny.svg">Conditional history</a></p>'
                    f'<object type="image/svg+xml" data="{rel}/phylogeny.svg" style="width:100%"></object></section>')
    write_tsv(output/'visualization_manifest.tsv',manifest,
              ['target_id','family_id','layer','site_id','label','color','view','path','type','intervals','selected_for_analysis','description'])
    write_tsv(output/'target_intervals.tsv',interval_table,['target_id','species','gene_copy_id','contig','start0','end0','strand','occurrence_id','kind','coordinate_evidence'])
    (output/'index.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>IntraPhy 0.17 target gallery</title>'
        '<style>body{font-family:system-ui;margin:25px}section{margin:35px 0;border-bottom:1px solid #ccc}a{margin-right:10px}</style>'
        '<h1>IntraPhy 0.17 · one object per figure group</h1><p>Gray is the unchanged supplied annotation, not proof of absence or conservation. '
        'Colors identify different targets; all three views of the same target share its color. No new inference is run by the renderer.</p>'
        +''.join(html)+'</html>',encoding='utf-8')
    return manifest
