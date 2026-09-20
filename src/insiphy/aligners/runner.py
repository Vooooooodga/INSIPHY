"""aligners / runner: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations
from typing import Optional

import shutil


def _backend_version(backend: str) -> Optional[str]:
    if backend != "internal":
        return None
    from Bio import __version__ as biopython_version

    return str(biopython_version)


def available_alignment_backends():
    rows = [{"aligner": "auto", "available": 1, "notes": "global uses MAFFT; local selects LASTZ/minimap2 when available"}]
    rows.append({"aligner": "internal", "available": 1, "notes": "Biopython PairwiseAligner, explicit small-pair backend"})
    rows.append({"aligner": "minimap2", "available": int(shutil.which("minimap2") is not None), "notes": "external nucleotide aligner"})
    rows.append({"aligner": "miniprot", "available": int(shutil.which("miniprot") is not None), "notes": "external protein-to-genome aligner"})
    rows.append({"aligner": "mafft", "available": int(shutil.which("mafft") is not None), "notes": "external global nucleotide aligner; local requests use terminal-overhang-trimmed overlap projection"})
    rows.append({"aligner": "lastz", "available": int(shutil.which("lastz") is not None), "notes": "external local nucleotide aligner"})
    return rows
