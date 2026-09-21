"""CDS consequences from actual transcript-specific sequence, never an ORF filter."""
from __future__ import annotations
from Bio.Seq import Seq
from Bio.Data import CodonTable
from .native import NativeLocus


def native_cds(locus: NativeLocus, transcript: str, genetic_code: int = 1):
    if genetic_code not in CodonTable.unambiguous_dna_by_id:
        raise ValueError("Unknown genetic code identifier")
    if locus.translation_exceptions:
        return {"status": "translation_exception_not_modelled", "cds": None, "protein": None,
                "annotation_exceptions": list(locus.translation_exceptions),
                "used_as_evolutionary_filter": False}
    if transcript in locus.coding_by_transcript:
        positions = locus.coding_by_transcript[transcript]
    else:
        by_id = {e.id: e for e in locus.exons}
        positions = tuple(p for identifier in locus.paths[transcript]
            for p in sorted(by_id[identifier].cds, key=lambda p: p[0], reverse=locus.strand == "-"))
    if not positions:
        return {"status": "noncoding_or_CDS_unavailable", "cds": None, "protein": None}
    pieces = []
    for a, b, phase in positions:
        if locus.strand == "+":
            left, right = a-(locus.search_start-1), b-(locus.search_start-1)
        else:
            left, right = locus.search_end-b, locus.search_end-a
        if not 0 <= left < right <= len(locus.sequence):
            raise ValueError("CDS interval is outside the genomic sequence")
        pieces.append(locus.sequence[left:right])
    first_phase = positions[0][2]
    if first_phase is None:
        return {"status": "initial_CDS_phase_unknown", "cds": "".join(pieces), "protein": None}
    # Internal phases describe split codons; those bases must not be discarded.
    sequence = "".join(pieces)[first_phase:]
    if set(sequence)-set("ACGTN"):
        return {"status": "ambiguous_coding_bases", "cds": sequence, "protein": None}
    length = len(sequence)-len(sequence) % 3
    protein = str(Seq(sequence[:length]).translate(table=genetic_code))
    return {"status": "native_annotation_conditioned", "cds": sequence, "protein": protein,
            "length": len(sequence), "partial_terminal_codon_bases": len(sequence) % 3,
            "internal_stop_count": protein[:-1].count("*"), "genetic_code": genetic_code,
            "NMD_inferred": False, "used_as_evolutionary_filter": False}
