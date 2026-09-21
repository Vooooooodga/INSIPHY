"""Sequence-qualified whole annotation alternatives, independent of tree history.

This applies the same rule to missing, shortened, shifted, split and fused
annotations. Candidates are whole configurations from actual annotated paths,
not arbitrary combinations of independently masked boundary intervals.
"""
from __future__ import annotations
from .types import ConfigurationAlternative, ExonConfiguration
from .alignment import paired_support
from .material import normalize_exons


def _changed_cuts(first, second):
    a = {p for e in first.exons for p in (e.start, e.end)}
    b = {p for e in second.exons for p in (e.start, e.end)}
    return a ^ b


def structure_alternatives(catalogue, alignment, species, natives, observed_native,
                           material_observations, *, minimum_identity, anchor_bases,
                           anchor_identity):
    """Return candidates and their qualification trace, never a corrected GFF.

    Retained material, paired exon sequence and flanking location are required.
    At every changed cut, the two bases on each side must be identical to the
    source. This is a conservative candidate rule, not a splice-function test or
    a universal canonical-site requirement. Sequence-divergent alternatives can
    remain unresolved; failure of this rule does not prove biological change.
    """
    row = alignment.rows[species]
    offset = catalogue.alignment_offset
    presence = material_observations.get(species, ())
    accepted, trace = [], []
    for source, paths in sorted(observed_native.items()):
        if source == species:
            continue
        other_presence = material_observations.get(source, ())
        if presence != other_presence or any(p is None for p in presence):
            continue  # Never restore absent DNA by changing annotation.
        other = alignment.rows[source]
        for tx, exons in paths:
            if tx == "unresolved" or not exons:
                continue
            material = tuple(1 if p == 1 else 0 for p in presence)
            candidate = ExonConfiguration(normalize_exons(catalogue, material, exons), material)
            for native in natives:
                if candidate.exons == native.exons:
                    continue
                compared = tuple(sorted(set((*candidate.exons, *native.exons))))
                if not compared:
                    continue
                left, right = min(e.start for e in compared)+offset, max(e.end for e in compared)+offset
                reason, identities = "qualified", []
                for e in compared:
                    start, end = e.start+offset, e.end+offset
                    paired, identity = paired_support(row, other, start, end)
                    possible = sum(b in "ACGT" for b in other[start:end])
                    if not possible or paired != possible or identity < minimum_identity:
                        reason = "exon_sequence_not_fully_supported"
                        break
                    identities.append(identity)
                if reason == "qualified":
                    for start, end in ((left-anchor_bases, left), (right, right+anchor_bases)):
                        paired, identity = paired_support(row, other, max(0, start), end)
                        if start < 0 or end > len(row) or paired != anchor_bases or identity < anchor_identity:
                            reason = "paired_location_flanks_unavailable"
                            break
                if reason == "qualified":
                    for cut in _changed_cuts(native, candidate):
                        position = cut+offset
                        a, b = row[position-2:position+2], other[position-2:position+2]
                        if position < 2 or len(a) != 4 or a != b or set(a) - set("ACGT"):
                            reason = "changed_boundary_context_not_preserved"
                            break
                identity = min(identities, default=0.)
                trace.append({"family_id": catalogue.family, "unit_id": catalogue.unit,
                    "species": species, "source_species": source, "source_transcript": tx,
                    "native_configuration": native.key, "candidate_configuration": candidate.key,
                    "candidate_exons": [[e.start, e.end] for e in candidate.exons],
                    "evidence": reason, "identity": identity,
                    "status": "predicted_only" if reason == "qualified" else "not_qualified",
                    "promoted_to_annotation": False})
                if reason == "qualified":
                    accepted.append(ConfigurationAlternative(candidate, native.key, source, tx,
                        "paired_sequence_flanks_and_preserved_boundary_context", identity))
    unique = {(a.replaces_key, a.configuration.key, a.source_species, a.source_transcript): a for a in accepted}
    return tuple(unique[k] for k in sorted(unique)), trace
