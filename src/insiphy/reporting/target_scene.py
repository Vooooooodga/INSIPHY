"""Immutable native-background geometry, reused verbatim by every target figure."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from html import escape
from pathlib import Path
import json

from insiphy.storage.tabular import read_tsv
from insiphy.topology import SpeciesTree

WIDTH=1680
TRACK_LEFT=510
TRACK_RIGHT=1350


def attr(value):
    return escape(str(value),quote=True)


def text(x,y,value,size=11,fill="#333333",anchor="start",weight="normal"):
    return (f'<text x="{x:.2f}" y="{y:.2f}" fill="{fill}" font-size="{size}" '
            f'font-family="sans-serif" text-anchor="{anchor}" font-weight="{weight}">{escape(str(value))}</text>')


@dataclass(frozen=True)
class NativeFeature:
    species: str
    copy: str
    lane: str
    contig: str
    start: int
    end: int
    strand: str
    kind: str
    feature_id: str
    source: str


@dataclass(frozen=True)
class Track:
    species: str
    copy: str
    contig: str
    strand: str
    start: int
    end: int
    top: float
    lane_y: tuple

    def x(self,coordinate0):
        fraction=(self.end-coordinate0)/(self.end-self.start+1) if self.strand=='-' else (coordinate0-(self.start-1))/(self.end-self.start+1)
        return TRACK_LEFT+fraction*(TRACK_RIGHT-TRACK_LEFT)

    @property
    def center(self):
        return sum(y for _,y in self.lane_y)/len(self.lane_y)


@dataclass(frozen=True)
class NativeScene:
    family: str
    features: tuple
    tracks: tuple
    species_y: tuple
    tree: object
    node_xy: tuple
    background_svg: str
    height: float
    annotation_source: str

    def track_for(self,species,copy):
        return next((track for track in self.tracks if track.species==species and track.copy==copy),None)


def _features(input_dir,family,occurrences):
    raw=read_tsv(Path(input_dir)/'raw_gene_features.tsv',optional=True)
    raw=[row for row in raw if row.get('family_id')==family]
    features=[]
    if raw:
        for row in raw:
            kind=row.get('type','unknown')
            if kind.lower() in {'gene','mrna','transcript','pseudogene'}:
                continue
            parents=str(row.get('parents') or row.get('parent') or 'locus annotations').replace(',',';').split(';')
            for parent in parents:
                try:
                    features.append(NativeFeature(row['species'],row['gene_copy_id'],parent or 'locus annotations',
                        row.get('seqid',row.get('contig','NA')),int(row['start']),int(row['end']),row['strand'],
                        kind,row.get('id','NA'),row.get('source','supplied_annotation')))
                except (KeyError,ValueError) as exc:
                    raise ValueError('invalid raw native feature coordinate') from exc
        source='raw_gene_features.tsv; all supplied native annotations, not inferred ancestor structures'
    else:
        paths=read_tsv(Path(input_dir)/'transcript_paths.tsv',optional=True)
        by_occ=defaultdict(list)
        for row in paths:
            if row.get('family_id')==family:
                by_occ[row.get('occurrence_id')].append(row)
        for occ in occurrences:
            if occ.get('family_id')!=family or occ.get('presence_status')=='absent':
                continue
            paths_for=by_occ.get(occ['occurrence_id']) or [occ]
            for row in paths_for:
                row={**occ,**row}
                try:
                    features.append(NativeFeature(row['species'],row['gene_copy_id'],row.get('transcript_id') or 'unassigned DNA',
                        row.get('contig','NA'),int(row['start']),int(row['end']),row['strand'],row.get('role','unknown'),
                        row.get('source_feature_id') or row['occurrence_id'],'prepared_native_record'))
                except (KeyError,ValueError) as exc:
                    raise ValueError('invalid prepared feature coordinate') from exc
        source='prepared occurrence/path records; full raw annotation unavailable'
    for feature in features:
        if feature.start<1 or feature.end<feature.start or feature.strand not in {'+','-','.','?'}:
            raise ValueError(f'invalid native feature: {feature}')
    return tuple(sorted(set(features),key=lambda x:(x.species,x.copy,x.lane,x.start,-x.end,x.kind,x.feature_id))),source


def build_scene(input_dir,family,occurrences):
    features,source=_features(input_dir,family,occurrences)
    rows=read_tsv(Path(input_dir)/'species_tree.tsv',optional=True)
    tree=SpeciesTree(rows) if rows else None
    species=([tree.label[node] for node in tree.preorder() if node in tree.leaves] if tree else [])
    species+=sorted({f.species for f in features}-set(species))
    loci=read_tsv(Path(input_dir)/'gene_loci.tsv',optional=True)
    locus_by={(r.get('species'),r.get('gene_copy_id')):r for r in loci}
    grouped=defaultdict(list)
    for feature in features:grouped[(feature.species,feature.copy)].append(feature)
    tracks=[];species_y=[];cursor=160.0
    for sp in species:
        sp_start=cursor
        for key,items in sorted(grouped.items()):
            if key[0]!=sp:continue
            contigs={f.contig for f in items}
            if len(contigs)!=1:
                raise ValueError(f'one plot locus spans conflicting contigs: {key}')
            start=min(f.start for f in items);end=max(f.end for f in items)
            locus=locus_by.get(key,{})
            focal_strands={row.get('strand') for row in occurrences
                           if (row.get('species'),row.get('gene_copy_id'))==key
                           and row.get('strand') in {'+','-'}}
            orientation=locus.get('strand')
            if orientation not in {'+','-'}:
                if len(focal_strands)==1: orientation=next(iter(focal_strands))
                elif len({f.strand for f in items})==1 and items[0].strand in {'+','-'}: orientation=items[0].strand
                else: raise ValueError(f'focal locus strand is unresolved: {key}')
            try:
                start=min(start,int(locus.get('annotation_start',start)))
                end=max(end,int(locus.get('annotation_end',end)))
            except (ValueError,TypeError): pass
            lanes=tuple((lane,cursor+28+i*25) for i,lane in enumerate([*sorted({f.lane for f in items}), 'genomic DNA']))
            tracks.append(Track(sp,key[1],next(iter(contigs)),orientation,start,end,cursor,lanes))
            cursor+=max(1,len(lanes))*25+60
        if cursor==sp_start:
            cursor+=65
        species_y.append((sp,(sp_start+cursor)/2-12))
        cursor+=38
    ys=dict(species_y);node_xy=[]
    bg=['<g id="native-background" data-annotation-layout="fixed-per-family" data-filtering="none">']
    if tree:
        depth={tree.root:0}
        for node in tree.preorder():
            for child in tree.children.get(node,[]):depth[child]=depth[node]+1
        maxdepth=max(depth.values()) or 1
        node_y={node:ys[tree.label[node]] for node in tree.leaves}
        for node in tree.postorder():
            if node not in node_y:node_y[node]=sum(node_y[c] for c in tree.children[node])/len(tree.children[node])
        for node in tree.preorder():
            x=28+132*depth[node]/maxdepth;y=node_y[node];node_xy.append((node,x,y))
            parent=tree.parent[node]
            if parent:
                px=28+132*depth[parent]/maxdepth;py=node_y[parent]
                bg.append(f'<path d="M {px:.2f},{py:.2f} V {y:.2f} H {x:.2f}" fill="none" stroke="#888888"/>')
        bg.append(text(24,137,'Fixed tree; topology layout',9))
    for sp,y in species_y:
        bg.append(text(178,y+4,sp,11,weight='bold'))
        if not any(t.species==sp for t in tracks):bg.append(text(TRACK_LEFT,y,'No supplied native annotation (not gene absence)',10))
    for track in tracks:
        bg.append(text(325,track.top+11,track.copy,9))
        bg.append(text(TRACK_LEFT,track.top+11,f'{track.contig}:{track.start:,}–{track.end:,} ({track.strand})',9))
        lane_y=dict(track.lane_y)
        for lane,y in track.lane_y:
            label=lane if len(lane)<=24 else lane[:21]+'…'
            bg.append(text(TRACK_LEFT-12,y+4,label,9,anchor='end'))
            # This is a coordinate axis, not an inferred splice or adjacency edge.
            bg.append(f'<line x1="{TRACK_LEFT}" x2="{TRACK_RIGHT}" y1="{y:.2f}" y2="{y:.2f}" stroke="#DDDDDD" data-axis="genomic-coordinate"/>')
        for feature in features:
            if (feature.species,feature.copy)!=(track.species,track.copy):continue
            x1,x2=sorted((track.x(feature.start-1),track.x(feature.end)))
            y=lane_y[feature.lane]
            kind=feature.kind.lower()
            h=17 if kind in {'exon','noncoding_exon'} else 9 if kind=='cds' else 7
            fill='#C4C4C4' if kind in {'exon','cds'} else '#E5E5E5'
            if kind=='intron':h=3;fill='#B9B9B9'
            if kind in {'gap','assembly_gap'}:fill='url(#native-gap)';h=18
            data=(f'data-native="true" data-species="{attr(feature.species)}" data-copy="{attr(feature.copy)}" '
                  f'data-feature="{attr(feature.feature_id)}" data-kind="{attr(feature.kind)}" data-strand="{attr(feature.strand)}" '
                  f'data-start="{feature.start}" data-end="{feature.end}" data-lane="{attr(feature.lane)}"')
            bg.append(f'<rect x="{x1:.3f}" y="{y-h/2:.3f}" width="{max(.1,x2-x1):.3f}" height="{h}" fill="{fill}" stroke="#8A8A8A" stroke-width=".7" {data}>'
                      f'<title>{attr(feature.feature_id)} | {attr(feature.kind)} | {feature.start}–{feature.end} | {attr(feature.source)}</title></rect>')
        last_y=track.lane_y[-1][1]+19
        for fraction in (0,.25,.5,.75,1):
            x=TRACK_LEFT+fraction*(TRACK_RIGHT-TRACK_LEFT)
            value=round(track.end-fraction*(track.end-track.start)) if track.strand=='-' else round(track.start+fraction*(track.end-track.start))
            bg.append(text(x,last_y,f'{value:,}',8,fill='#666666',anchor='middle'))
    bg.append('</g>')
    return NativeScene(family,features,tuple(tracks),tuple(species_y),tree,tuple(node_xy),''.join(bg),cursor+82,source)


def scene_data(scene):
    """One plain JSON geometry record, not a checksum or a second annotation."""
    from dataclasses import asdict
    return {'family_id':scene.family,'annotation_source':scene.annotation_source,
            'width':WIDTH,'height':scene.height,'coordinate_compression':False,
            'scale':'each locus linearly scaled 5prime to 3prime; no trim',
            'tracks':[asdict(t) for t in scene.tracks],
            'features':[asdict(f) for f in scene.features]}
