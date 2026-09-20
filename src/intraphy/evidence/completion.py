"""evidence / completion: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import Counter
from intraphy.elements import EXON_LIKE_ROLES
from intraphy.evidence.fields import _EVIDENCE_PROVENANCE_FIELDS
from intraphy.evidence.fields import _PROTEIN_QUERY_INTERVAL_FIELDS
from intraphy.evidence.fields import _SEARCH_COMPAT_FIELDS
from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from intraphy.storage.values import to_float
from pathlib import Path


def support_score(row):
    if row.get("evidence_status") == "supports_absence":
        return None
    identity = to_float(row.get("sequence_score"))
    coverage = to_float(row.get("sequence_coverage"), 1.0)
    return min(identity, coverage)


def completion_call(row, score, threshold):
    status = row.get("evidence_status", "ambiguous")
    anchor_status = row.get("anchor_interval_status", "unknown")
    if anchor_status not in {
        "ordered_double_flank_bounded_interval",
        "ordered_double_flank_empty_interval",
    }:
        return "ambiguous_evidence"
    event = row.get("inferred_event", "")
    frame = row.get("frame_status", "")
    inferred_role = row.get("inferred_role", "unknown")
    if (
        status == "supports_hidden_segment"
        and event == "protein_cds_projection"
        and row.get("predicted_role") in EXON_LIKE_ROLES
        and score >= threshold
    ):
        return "predicted_exon_candidate"
    if (
        status == "conflicts_annotation"
        and row.get("annotation_status") == "protein_projection_boundary_conflict"
        and row.get("predicted_role") in EXON_LIKE_ROLES
        and score >= threshold
    ):
        return "boundary_conflict_candidate"
    if status == "supports_hidden_segment" and inferred_role in EXON_LIKE_ROLES and score >= threshold:
        if event == "shifted_splice_site":
            return "shifted_splice_site_candidate"
        if event == "intron_deletion_joined_exon":
            return "joined_exon_candidate"
        if frame == "frameshift_or_stop_risk":
            return "hidden_segment_with_frame_disruption"
        return "hidden_segment_candidate"
    if status == "conflicts_annotation" and inferred_role in EXON_LIKE_ROLES and score >= threshold:
        if row.get("annotation_status") == "protein_projection_boundary_conflict":
            return "boundary_conflict_candidate"
        return "annotation_conflict_candidate"
    if status == "supports_annotation":
        return "supports_annotation"
    if status == "supports_absence":
        return "supports_true_absence"
    if status == "homologous_sequence_candidate":
        return "homologous_sequence_candidate"
    return "ambiguous_evidence"


def complete_annotation(input_dir, output_dir, threshold=0.55):
    generated = Path(output_dir) / "sequence_synteny_evidence.tsv"
    evidence_path = generated if generated.exists() else Path(input_dir) / "sequence_synteny_evidence.tsv"
    evidence = read_tsv(
        evidence_path,
        ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "annotation_status", "evidence_status"],
        optional=True,
    )
    rows = []
    summary = Counter()
    for row in evidence:
        score = support_score(row)
        status = row.get("evidence_status", "ambiguous")
        call = completion_call(row, score, threshold)
        summary[call] += 1
        rows.append(
            {
                "evidence_id": row["evidence_id"],
                "family_id": row["family_id"],
                "species": row["species"],
                "gene_copy_id": row["gene_copy_id"],
                "homology_id": row["homology_id"],
                "interval": f"{row.get('contig', 'NA')}:{row.get('start', 'NA')}-{row.get('end', 'NA')}:{row.get('strand', 'NA')}",
                "annotation_status": row.get("annotation_status", "unknown"),
                "inferred_role": row.get("inferred_role", "unknown"),
                "predicted_role": row.get("predicted_role", "unknown"),
                **{
                    field: row.get(
                        field,
                        "unknown"
                        if field in {
                            "homologous_dna_presence",
                            "homologous_dna_evidence",
                            "predicted_exonic_role",
                            "supplied_annotation_role",
                            "candidate_resolution_status",
                            "candidate_search_complete",
                        }
                        else "NA",
                    )
                    for field in _EVIDENCE_PROVENANCE_FIELDS
                },
                **{field: row.get(field, "NA") for field in _SEARCH_COMPAT_FIELDS},
                "primary_mapping_status": row.get("primary_mapping_status", "unknown"),
                # Legacy evidence lacks this field; only explicit unknown marks unresolved correspondence.
                "correspondence_status": row.get("correspondence_status", "unspecified"),
                "interval_scope": row.get("interval_scope", "unspecified"),
                **{field: row.get(field, "NA") for field in _PROTEIN_QUERY_INTERVAL_FIELDS},
                "reference_protein_id": row.get("reference_protein_id", "NA"),
                "reference_transcript_id": row.get("reference_transcript_id", "NA"),
                "reference_protein_query_start": row.get("reference_protein_query_start", "NA"),
                "reference_protein_query_end": row.get("reference_protein_query_end", "NA"),
                "interval_candidates": row.get("interval_candidates", "NA"),
                "support_score": "NA" if score is None else f"{score:.6g}",
                "completion_call": call,
                "evidence_status": status,
                "inferred_event": row.get("inferred_event", "NA"),
                "frame_status": row.get("frame_status", "NA"),
                "evidence_conclusion": row.get("evidence_conclusion", "unknown"),
                "confidence_flag": row.get("confidence_flag", "low"),
                "absence_evidence": row.get("absence_evidence", "NA"),
            }
        )
    write_tsv(
        f"{output_dir}/annotation_completion_candidates.tsv",
        rows,
        [
            "evidence_id", "family_id", "species", "gene_copy_id", "homology_id",
            "interval", "annotation_status", "inferred_role", "predicted_role", "support_score",
            *_EVIDENCE_PROVENANCE_FIELDS,
            *_SEARCH_COMPAT_FIELDS,
            "primary_mapping_status", "correspondence_status", "interval_scope", "interval_candidates",
            *_PROTEIN_QUERY_INTERVAL_FIELDS,
            "reference_protein_id", "reference_transcript_id",
            "reference_protein_query_start", "reference_protein_query_end",
            "completion_call", "evidence_status", "inferred_event", "frame_status",
            "evidence_conclusion", "confidence_flag", "absence_evidence",
        ],
    )
    write_tsv(
        f"{output_dir}/annotation_completion_summary.tsv",
        [{"completion_call": key, "count": value} for key, value in sorted(summary.items())],
        ["completion_call", "count"],
    )
    return rows
