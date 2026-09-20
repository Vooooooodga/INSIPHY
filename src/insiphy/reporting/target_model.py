"""Explicit per-object figure targets. Selection never changes inference inputs."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import colorsys
from pathlib import Path
import re

from insiphy.coordinates import Interval0, ClosedInterval1, local_interval_to_genome
from insiphy.reporting.inputs import _parse_membership_blocks, ribbon_match_evidence
from insiphy.storage.tabular import read_tsv

LAYER_LABELS = {"exon_presence": "Homologous DNA presence", "exon_role": "Annotated exonic use",
                "splice_junction": "Splice position / connection"}
PALETTE = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#5B4FA1", "#A56720", "#009EB5", "#9C3372")


def tokens(value):
    return {part.strip() for part in str(value or "").split(";") if part.strip() not in {"", "NA"}}


@dataclass(frozen=True)
class FigureTarget:
    target_id: str
    family_id: str
    layer: str
    site_id: str
    element_id: str
    label: str
    color: str


@dataclass(frozen=True)
class TargetInterval:
    species: str
    gene_copy_id: str
    contig: str
    interval: Interval0
    strand: str
    occurrence_id: str = ""
    coordinate_evidence: str = "explicit_interval"
    kind: str = "segment"


def _color(index):
    if index < len(PALETTE):
        return PALETTE[index]
    hue = (index * 0.618033988749895) % 1
    r,g,b = colorsys.hls_to_rgb(hue, 0.38 + 0.10 * (index % 2), 0.64)
    return '#%02X%02X%02X' % (round(r*255), round(g*255), round(b*255))


def load_targets(result_dir, selectors=(), manifest=None):
    """Colors are assigned over the complete catalogue, before selection.

    A manifest row declares one target identity, optionally repeated with native
    intervals (1-based inclusive start/end). A target cannot mix multiple sites.
    """
    result = Path(result_dir)
    matrix = read_tsv(result / "structural_site_matrix.tsv", optional=True)
    members = read_tsv(result / "element_correspondence.tsv", optional=True)
    catalog = {}
    for row in matrix:
        key = (row["family_id"], row["layer"], row["site_id"])
        catalog[key] = {"element_id": row["site_id"] if row["layer"] != "splice_junction" else "",
                        "label": row["site_id"]}
    if not catalog:
        for row in members:
            key = (row.get("family_id", "NA"), "exon_presence", row["element_id"])
            catalog[key] = {"element_id": row["element_id"], "label": row["element_id"]}
    manifest_rows = read_tsv(manifest) if manifest else []
    specifications = {}
    manual = defaultdict(list)
    for row in manifest_rows:
        ident = row.get("target_id", "")
        if not ident or ident in {".",".."}:
            raise ValueError("target manifest requires target_id")
        family = row.get("family_id", "")
        layer = row.get("layer") or "exon_presence"
        site = row.get("site_id") or row.get("element_id")
        if not family or layer not in LAYER_LABELS or not site:
            raise ValueError(f"target {ident}: family_id, valid layer and site_id/element_id required")
        element = row.get("element_id") or (site if layer != "splice_junction" else "")
        spec = (family, layer, site, element, row.get("label") or site, row.get("color") or "")
        if ident in specifications and specifications[ident] != spec:
            raise ValueError(f"target {ident} mixes identities, labels or colors")
        specifications[ident] = spec
        if any(row.get(field) not in {None,"","NA"} for field in ("start","end")):
            if not all(row.get(field) not in {None,"","NA"} for field in ("species","gene_copy_id","contig","start","end","strand")):
                raise ValueError(f"target {ident}: explicit interval requires species/copy/contig/start/end/strand")
            if row['strand'] not in {'+','-'}:
                raise ValueError("target interval strand must be + or -")
            manual[ident].append(TargetInterval(row['species'], row['gene_copy_id'], row['contig'],
                ClosedInterval1(int(row['start']),int(row['end'])).to_interval0(), row['strand'],
                row.get('occurrence_id',''), 'user_specified_native_interval'))
    if not manifest:
        for family, layer, site in sorted(catalog):
            ident = '/'.join((family,layer,site))
            item = catalog[(family,layer,site)]
            specifications[ident] = (family,layer,site,item['element_id'],item['label'],'')
    all_targets = []
    used_colors = set()
    for index,(ident,spec) in enumerate(sorted(specifications.items())):
        family, layer, site, element, label, supplied_color = spec
        color = supplied_color or _color(index)
        if not re.fullmatch(r'#[0-9A-Fa-f]{6}',color):
            raise ValueError(f"target {ident}: color must be #RRGGBB")
        color = color.upper()
        if color in used_colors:
            if supplied_color:
                raise ValueError("different targets must not share a color within one gallery")
            n = index + len(specifications)
            while color in used_colors:
                color = _color(n).upper(); n += 1
        used_colors.add(color)
        all_targets.append(FigureTarget(ident, family, layer, site, element, label, color))
    if selectors:
        chosen = set()
        for query in selectors:
            matches = [t for t in all_targets if query == t.target_id or query == t.site_id]
            if not matches:
                raise ValueError(f"unknown plot target: {query}")
            if len(matches)>1 and all(query != t.target_id for t in matches):
                raise ValueError(f"ambiguous plot target {query}; use family/layer/site")
            chosen.update(t.target_id for t in matches)
        all_targets = [t for t in all_targets if t.target_id in chosen]
    return all_targets, manual, matrix, members


def resolve_intervals(target, occurrences, members, boundaries, manual=()):
    """Use actual matched subranges when available; never color a bounding span
    of several disjoint matches as if every base were aligned.
    """
    if manual:
        return tuple(sorted(set(manual), key=lambda x:(x.species,x.gene_copy_id,x.interval)))
    occ_by_id = {row['occurrence_id']: row for row in occurrences}
    intervals = []
    if target.layer == 'splice_junction':
        for row in boundaries:
            if (row.get('family_id'), row.get('site_id')) != (target.family_id,target.site_id):
                continue
            if row.get('position_edge_eligible') != '1':
                continue
            left, right = occ_by_id.get(row.get('left_occurrence_id')), occ_by_id.get(row.get('right_occurrence_id'))
            if not left or not right or left.get('contig') != right.get('contig') or left.get('strand') != right.get('strand'):
                continue
            # A junction is two cuts, not a ribbon covering all intervening DNA.
            for occ,side in ((left,'donor'),(right,'acceptor')):
                cut = (int(occ['end']) if side=='donor' else int(occ['start'])-1) if occ['strand']=='+' else (
                       int(occ['start'])-1 if side=='donor' else int(occ['end']))
                intervals.append(TargetInterval(occ['species'],occ['gene_copy_id'],occ['contig'],
                    Interval0(cut,cut),occ['strand'],occ['occurrence_id'],'explicit_annotated_junction',side))
    else:
        for row in members:
            if row.get('element_id') != target.element_id or row.get('family_id') != target.family_id:
                continue
            occ = occ_by_id.get(row.get('occurrence_id'))
            if not occ or occ.get('presence_status') in {'absent','unknown','missing','NA'}:
                continue
            locus = ClosedInterval1(int(occ['start']),int(occ['end'])).to_interval0()
            ranges = _parse_membership_blocks(row.get('matched_blocks',''),locus.length)
            resolved = (row.get('membership_call','core_member') in {'core_member','resolved_member','resolved'}
                        and row.get('position_edge_eligible','1') not in {'0',0,False}
                        and row.get('correspondence_status') not in {'ambiguous','unresolved','candidate'})
            coordinate_evidence = 'actual_membership_subinterval' if resolved else 'ambiguous_membership_subinterval'
            for _key,start,end in ranges:
                interval = local_interval_to_genome(Interval0(start-1,end),locus,occ['strand'])
                intervals.append(TargetInterval(occ['species'],occ['gene_copy_id'],occ['contig'],interval,
                    occ['strand'],occ['occurrence_id'],coordinate_evidence))
            if not ranges:
                intervals.append(TargetInterval(occ['species'],occ['gene_copy_id'],occ['contig'],locus,
                    occ['strand'],occ['occurrence_id'],'annotated_parent_range_only'))
    return tuple(sorted(set(intervals),key=lambda x:(x.species,x.gene_copy_id,x.interval,x.kind)))


def target_ribbons(target, intervals, occurrences, input_dir, result_dir, callable_species, *, match_data=None):
    """Only actual paired alignment columns within this target. No all-to-all
    links through shared membership, no ribbons across an unresolved interval.
    Junction targets use cut markers; they are not DNA spans and have no ribbons.
    """
    if target.layer == 'splice_junction':
        return []
    selected = defaultdict(list)
    occs = {row['occurrence_id']:row for row in occurrences}
    for item in intervals:
        if (item.occurrence_id in occs and item.species in callable_species
                and item.coordinate_evidence in {'actual_membership_subinterval','user_specified_native_interval'}):
            selected[item.occurrence_id].append(item)
    matches = match_data if match_data is not None else ribbon_match_evidence(input_dir,result_dir,occurrences)[0]
    pairs = []
    seen = set()
    for (aid,bid),match in sorted(matches.items()):
        if aid not in selected or bid not in selected:
            continue
        a,b = occs[aid],occs[bid]
        if a['species'] == b['species']:
            continue
        la = ClosedInterval1(int(a['start']),int(a['end'])).to_interval0()
        lb = ClosedInterval1(int(b['start']),int(b['end'])).to_interval0()
        for qs,qe,ts,te in match['blocks']:
            qa, qb = Interval0(qs-1,qe),Interval0(ts-1,te)
            # The legacy matched-block contract is equal-length ungapped columns.
            if qa.length != qb.length:
                continue
            ga=local_interval_to_genome(qa,la,a['strand']); gb=local_interval_to_genome(qb,lb,b['strand'])
            for ia in selected[aid]:
                for ib in selected[bid]:
                    def offsets(block,highlight,strand):
                        lo=max(block.start0,highlight.start0); hi=min(block.end0,highlight.end0)
                        if lo>=hi:return None
                        return (lo-block.start0,hi-block.start0) if strand=='+' else (block.end0-hi,block.end0-lo)
                    oa=offsets(ga,ia.interval,a['strand']); ob=offsets(gb,ib.interval,b['strand'])
                    if oa is None or ob is None:continue
                    lo=max(oa[0],ob[0]); hi=min(oa[1],ob[1])
                    if lo>=hi:continue
                    aa=local_interval_to_genome(Interval0(qa.start0+lo,qa.start0+hi),la,a['strand'])
                    bb=local_interval_to_genome(Interval0(qb.start0+lo,qb.start0+hi),lb,b['strand'])
                    key=(aid,bid,aa,bb)
                    if key not in seen:
                        seen.add(key); pairs.append((a,b,aa,bb,match['match_id']))
    return pairs
