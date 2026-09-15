"""Sequence-supported annotation completion."""

from collections import Counter

from .io import read_tsv, to_float, write_tsv


def support_score(row):
    seq = to_float(row.get("sequence_score"))
    left = to_float(row.get("left_synteny_score"))
    right = to_float(row.get("right_synteny_score"))
    motif = to_float(row.get("splice_motif_score"))
    phase_bonus = 0.05 if row.get("phase_compatibility") == "compatible" else 0.0
    frame_bonus = 0.05 if row.get("frame_status") in {"coding_frame_preserved", "coding_frame_annotated"} else 0.0
    return min(1.0, 0.36 * seq + 0.18 * left + 0.18 * right + 0.18 * motif + phase_bonus + frame_bonus)


def completion_call(row, score, threshold):
    status = row.get("evidence_status", "ambiguous")
    event = row.get("inferred_event", "")
    frame = row.get("frame_status", "")
    if status == "supports_hidden_segment" and score >= threshold:
        if event == "shifted_splice_site":
            return "shifted_splice_site_candidate"
        if event == "intron_deletion_joined_exon":
            return "joined_exon_candidate"
        if frame == "frameshift_or_stop_risk":
            return "hidden_segment_with_frame_disruption"
        return "hidden_segment_candidate"
    if status == "conflicts_annotation" and score >= threshold:
        return "annotation_conflict_candidate"
    if status == "supports_annotation":
        return "supports_annotation"
    if status == "supports_absence" and score >= threshold:
        return "supports_true_absence"
    return "ambiguous_evidence"


def complete_annotation(input_dir, output_dir, threshold=0.55):
    evidence = read_tsv(
        f"{input_dir}/sequence_synteny_evidence.tsv",
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
                "support_score": f"{score:.6g}",
                "completion_call": call,
                "evidence_status": status,
                "inferred_event": row.get("inferred_event", "NA"),
                "frame_status": row.get("frame_status", "NA"),
            }
        )
    write_tsv(
        f"{output_dir}/annotation_completion_candidates.tsv",
        rows,
        ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "interval", "annotation_status", "inferred_role", "support_score", "completion_call", "evidence_status", "inferred_event", "frame_status"],
    )
    write_tsv(
        f"{output_dir}/annotation_completion_summary.tsv",
        [{"completion_call": key, "count": value} for key, value in sorted(summary.items())],
        ["completion_call", "count"],
    )
    return rows
