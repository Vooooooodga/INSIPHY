"""evidence / reference queries: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.evidence.protein import _reference_protein_context
from intraphy.evidence.search_context import _safe_int


def _prepare_reference_queries(rows_by_element, sequences, aligner, transcript_paths, proteins, occ_by_id):
    representative_infos = []
    for (family, element), members in sorted(rows_by_element.items()):
        representative_row, representative = max(
            members,
            key=lambda pair: len(sequences.get(pair[1]["occurrence_id"], "")),
        )
        if aligner == "miniprot":
            protein_members = [
                pair for pair in members
                if _reference_protein_context(pair[1], transcript_paths, proteins, occ_by_id) is not None
            ]
            if protein_members:
                representative_row, representative = min(
                    protein_members,
                    key=lambda pair: (-_safe_int(pair[1].get("cds_length"), 0), pair[1]["occurrence_id"]),
                )
        query = sequences.get(representative["occurrence_id"], "")
        if not query:
            continue
        representative_infos.append(
            {
                "family": family,
                "element": element,
                "members": members,
                "representative_row": representative_row,
                "representative": representative,
                "query": query,
                "protein_context": _reference_protein_context(representative, transcript_paths, proteins, occ_by_id)
                if aligner == "miniprot"
                else None,
            }
        )

    protein_contexts_by_family = defaultdict(dict)
    if aligner == "miniprot":
        for info in representative_infos:
            context = info.get("protein_context")
            if context is not None:
                protein_contexts_by_family[info["family"]][context["protein_id"]] = context
    return representative_infos, protein_contexts_by_family
