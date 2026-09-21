"""One biological-unit validation contract shared by every catalogue entry point."""
from __future__ import annotations
from collections import defaultdict
from .types import Catalogue


def _check_disjoint(records, label):
    last_end, owner = -1, None
    for start, end, key in sorted(set(records)):
        if start < last_end and key != owner:
            raise ValueError(f"Overlapping physical regions in {label}: {owner} and {key}; "
                             "merge dependent regions, do not rename or count them twice")
        if end > last_end:
            last_end, owner = end, key


def validate_collection(catalogues: tuple[Catalogue, ...], *, allow_empty=False):
    """Each family's offsets refer to one common alignment coordinate system.

    Native coordinate evidence adds a second check, including renamed families.
    No arbitrary biological independence is inferred merely from distinct labels.
    """
    if not catalogues and not allow_empty:
        raise ValueError("Configuration collection is empty")
    keys = [(c.family, c.unit) for c in catalogues]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate family/unit in configuration collection")
    axes, native = defaultdict(list), defaultdict(list)
    for c in catalogues:
        key = (c.family, c.unit)
        axes[c.family].append((c.alignment_offset, c.alignment_offset+c.length, key))
        for e in c.exon_instances:
            # Locus names are metadata, not proof of different physical regions.
            native[(e.species, e.contig, e.strand)].append((e.start, e.end, key))
    for family, intervals in axes.items():
        _check_disjoint(intervals, f"family {family}'s common coordinate axis")
    for position, intervals in native.items():
        _check_disjoint(intervals, f"native genome {position}")
    return catalogues
