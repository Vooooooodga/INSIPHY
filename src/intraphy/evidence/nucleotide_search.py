"""evidence / nucleotide search: explicit implementation ownership."""
from __future__ import annotations

from intraphy.aligners.pairwise import local_alignment_stats
from intraphy.aligners.types import AlignmentBackendError
from intraphy.evidence.alignment_provenance import _offset_alignment_target
from intraphy.evidence.alignment_provenance import _set_alignment_context


def _align_nucleotide_evidence(aligner, anchor_context, query, short_context_max_length, threads, representative, source_transcript_ids, locus, min_identity, min_coverage):
    alignment = None
    alignment_error = ""
    alignment_evidence_scope = "unavailable"
    identity = coverage = 0.0
    dna_identity = dna_coverage = 0.0
    nucleotide_aligner = "minimap2" if aligner == "miniprot" else aligner
    bounded_target = anchor_context.get("sequence", "")
    bounded_backend = (
        "internal"
        if len(query) <= int(short_context_max_length)
        else nucleotide_aligner
    )
    if anchor_context["valid_double_flank"] and bounded_target:
        try:
            alignment = local_alignment_stats(
                query,
                bounded_target,
                backend=bounded_backend,
                threads=threads,
                query_occurrence_id=representative.get("occurrence_id"),
                target_occurrence_id=None,
                query_transcript_id=(source_transcript_ids[0] if len(source_transcript_ids) == 1 else None),
                left_anchor_id=anchor_context["left_anchor_id"],
                right_anchor_id=anchor_context["right_anchor_id"],
                search_interval=anchor_context["search_interval"],
            )
            alignment_evidence_scope = (
                "anchor_bounded_short_local"
                if bounded_backend == "internal"
                else "anchor_bounded_external_local"
            )
        except AlignmentBackendError as exc:
            if bounded_backend == "internal" and nucleotide_aligner != "internal":
                try:
                    alignment = local_alignment_stats(
                        query,
                        bounded_target,
                        backend=nucleotide_aligner,
                        threads=threads,
                    )
                    alignment_evidence_scope = "anchor_bounded_external_local"
                except AlignmentBackendError as fallback_exc:
                    alignment_error = f"bounded_alignment_unavailable:{fallback_exc}"
            else:
                alignment_error = f"bounded_alignment_unavailable:{exc}"
        if alignment is not None:
            _set_alignment_context(
                alignment,
                representative.get("occurrence_id", "NA"),
                None,
                source_transcript_ids[0] if len(source_transcript_ids) == 1 else None,
                anchor_context["left_anchor_id"],
                anchor_context["right_anchor_id"],
                anchor_context["search_interval"],
            )
            _offset_alignment_target(alignment, anchor_context["target_offset0"])
    if alignment is None and not anchor_context["valid_double_flank"]:
        try:
            alignment = local_alignment_stats(
                query, locus, backend=nucleotide_aligner, threads=threads
            )
            alignment_evidence_scope = "whole_locus_descriptive_fallback"
        except AlignmentBackendError as exc:
            alignment_error = str(exc)
    if alignment is not None:
        dna_identity = alignment.identity
        dna_coverage = alignment.query_coverage or alignment.coverage
        identity = dna_identity
        coverage = dna_coverage
    supported = (
        alignment_evidence_scope.startswith("anchor_bounded_")
        and alignment is not None
        and (bool(alignment.aligned_blocks) or int(alignment.aligned_pairs or 0) > 0)
        and identity >= min_identity
        and coverage >= min_coverage
    )
    return supported, alignment, dna_identity, dna_coverage, alignment_evidence_scope, nucleotide_aligner, alignment_error, bounded_target
