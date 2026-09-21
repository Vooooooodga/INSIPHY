"""The single registry of elementary edits, opportunities and counting semantics."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import replace
from .types import Catalogue, ElementaryEdit, ExonConfiguration, ExonSpan
from .material import normalize_exons, retained_length, valid_configuration

EDIT_KINDS = ("split", "fusion", "acceptor_shift", "donor_shift",
              "exonization", "exon_inactivation", "dna_insertion", "dna_deletion")


def elementary_edits(c: Catalogue, s: ExonConfiguration) -> tuple[ElementaryEdit, ...]:
    """Generate unique edits, distributing one opportunity across its destinations.

    Rates are per eligible exon, boundary, neighboring pair or DNA tract. Motif
    and coding scores never become transition weights. Material insertion can
    change exon boundaries in the same step; its consequences are not recounted.
    """
    candidates: dict[tuple[str, str, str], ElementaryEdit] = {}
    spans = frozenset(c.spans)

    def add(exons, kind, opportunity, footprint, affected, *, material=None, consequences=(), mid=None):
        target = ExonConfiguration(tuple(exons), s.material if material is None else tuple(material))
        if target == s or not valid_configuration(c, target) or not set(target.exons) <= spans:
            return
        edit = ElementaryEdit(s, target, kind, opportunity, footprint, tuple(affected), tuple(consequences), mid)
        candidates[(kind, opportunity, target.key)] = edit

    for i, e in enumerate(s.exons):
        before, after = s.exons[:i], s.exons[i+1:]
        add(before + after, "exon_inactivation", f"exon:{e.start}:{e.end}", (e.start, e.end), (e,))
        for d, a in c.junctions:
            if e.start < d < a < e.end and retained_length(c, s.material, d, a):
                add(before + (ExonSpan(e.start, d), ExonSpan(a, e.end)) + after,
                    "split", f"exon:{e.start}:{e.end}", (d, a), (e,), consequences=("exon_count_plus_one",))
        for dest in spans:
            if dest == e or not dest.overlaps(e):
                continue
            if before and before[-1].end >= dest.start or after and dest.end >= after[0].start:
                continue
            if dest.end == e.end:
                add(before + (dest,) + after, "acceptor_shift", f"acceptor:{e.start}",
                    (e.start, dest.start), (e,))
            if dest.start == e.start:
                add(before + (dest,) + after, "donor_shift", f"donor:{e.end}",
                    (e.end, dest.end), (e,))
    for i, (a, b) in enumerate(zip(s.exons, s.exons[1:])):
        add(s.exons[:i] + (ExonSpan(a.start, b.end),) + s.exons[i+2:],
            "fusion", f"pair:{a.start}:{a.end}:{b.start}:{b.end}", (a.end, b.start), (a, b),
            consequences=("exon_count_minus_one", "intervening_DNA_retained"))
    # A new exon must be one of the observed/candidate exon identities, not any
    # arbitrary cross-product of donor and acceptor coordinates.
    for e in c.boundary_candidates:
        if all(not e.overlaps(old) and e.start != old.end and e.end != old.start for old in s.exons):
            add(tuple(sorted((*s.exons, e))), "exonization", f"exon:{e.start}:{e.end}", (e.start, e.end), (e,))
    for k, m in enumerate(c.material):
        if s.material[k] == 1:
            material = list(s.material)
            material[k] = 2
            exons = normalize_exons(c, tuple(material), s.exons)
            effects = []
            lost = sum(retained_length(c, tuple(material), e.start, e.end) == 0 for e in s.exons)
            if lost:
                effects.append(f"deleted_exons:{lost}")
            if len(exons) < len(s.exons) - lost:
                effects.append("exon_fusion")
            if exons != s.exons and not effects:
                effects.append("exon_truncation_or_internal_deletion")
            affected = tuple(e for e in s.exons if e.start <= m.end and m.start <= e.end)
            add(exons, "dna_deletion", m.id, (m.start, m.end), affected,
                material=material, consequences=effects, mid=m.id)
        elif s.material[k] == 0:
            material = list(s.material)
            material[k] = 1
            affected = tuple(e for e in s.exons if e.start <= m.end and m.start <= e.end)
            add(s.exons, "dna_insertion", m.id, (m.start, m.end), affected,
                material=material, consequences=("material_introduction",), mid=m.id)
            # A tract inserted within an exon may introduce a spliced separator.
            for i, e in enumerate(s.exons):
                if e.start < m.start < m.end < e.end:
                    for d, a in c.junctions:
                        if e.start < d <= m.start and m.end <= a < e.end:
                            exons = s.exons[:i] + (ExonSpan(e.start, d), ExonSpan(a, e.end)) + s.exons[i+1:]
                            add(exons, "dna_insertion", m.id, (m.start, m.end), (e,),
                                material=material, consequences=("exon_split",), mid=m.id)
            # Inserting a source exon at a nonexonic position is one interval edit.
            for e in c.boundary_candidates:
                if m.start <= e.start < e.end <= m.end and all(not e.overlaps(x) for x in s.exons):
                    add(tuple(sorted((*s.exons, e))), "dna_insertion", m.id, (m.start, m.end), (e,),
                        material=material, consequences=("introduced_exon",), mid=m.id)
            # Only an explicitly evidenced module can arrive as a whole. The
            # presence of several target exons alone is not source evidence.
            for payload in c.insertion_payloads:
                if payload.material_id != m.id:
                    continue
                joined = tuple(sorted((*s.exons, *payload.exons)))
                if any(a.end >= b.start for a, b in zip(joined, joined[1:])):
                    continue
                add(joined, "dna_insertion", m.id, (m.start, m.end), payload.exons,
                    material=material, consequences=("introduced_exon_module",
                    f"introduced_exons:{len(payload.exons)}"), mid=m.id)
    opportunity_counts: dict[tuple[str, str], int] = defaultdict(int)
    for edit in candidates.values():
        opportunity_counts[(edit.kind, edit.opportunity)] += 1
    return tuple(replace(e, weight=1.0/opportunity_counts[(e.kind, e.opportunity)])
                 for _, e in sorted(candidates.items()))
