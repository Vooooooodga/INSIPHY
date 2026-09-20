"""evidence / fields: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations




_PROTEIN_QUERY_INTERVAL_FIELDS = (
    "protein_cds_query_start", "protein_cds_query_end",
    "protein_overlap_query_start", "protein_overlap_query_end",
)


_EVIDENCE_PROVENANCE_FIELDS = (
    "homologous_dna_presence",
    "homologous_dna_evidence",
    "predicted_exonic_role",
    "supplied_annotation_role",
    "supplied_annotation_roles",
    "source_parent_occurrence_id",
    "source_parent_transcript_ids",
    "target_parent_occurrence_ids",
    "target_parent_transcript_ids",
    "dna_aligned_blocks",
    "predicted_role_blocks",
    "supplied_annotation_overlaps",
    "annotation_conflict_blocks",
    "candidate_resolution_status",
    "candidate_search_complete",
    "candidate_search_incomplete_reason",
    "left_anchor_id",
    "right_anchor_id",
    "anchor_interval_status",
    "search_interval",
    "alignment_evidence_scope",
    "alignment_sequence_kind",
    "alignment_backend_version",
    "alignment_score_scheme",
    "alignment_nt_identity",
    "alignment_aa_identity",
    "alignment_known_aligned_pairs",
    "alignment_unknown_aligned_pairs",
    "alignment_query_covered_bases",
    "alignment_target_covered_bases",
    "alignment_query_length",
    "alignment_target_length",
    "alignment_gap_blocks",
    "alignment_relative_strand",
    "alignment_candidate_id",
    "alignment_alternative_candidate_ids",
    "alignment_enumeration_complete",
    "alignment_incomplete_reason",
    "original_annotation_start",
    "original_annotation_end",
    "linked_feature_start",
    "linked_feature_end",
    "search_bound_start",
    "search_bound_end",
    "search_left_limit_status",
    "search_right_limit_status",
    "hit_search_limit_status",
)


_SEARCH_COMPAT_FIELDS = (
    "search_expansion_status",
    "search_original_start",
    "search_original_end",
    "search_expanded_start",
    "search_expanded_end",
    "range_status",
)
