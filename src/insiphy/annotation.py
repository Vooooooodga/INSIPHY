"""Sequence-supported annotation completion."""

from collections import Counter

from .io import read_tsv, to_float, write_tsv


def support_score(row):
    seq = to_float(row.get("sequence_score"))
    left = to_float(row.get("left_synteny_score"))
    right = to_float(row.get("right_synteny_score"))
    motif = to_float(row.get("splice_motif_score"))
    phase_bonus = 0.05 if row.get("phase_compatibility") == "compatible" else 0.0
    return min(1.0, 0.40 * seq + 0.20 * left + 0.20 * right + 0.15 * motif + phase_bonus)


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
        if status == "supports_hidden_segment" and score >= threshold:
            call = "hidden_segment_candidate"
        elif status == "conflicts_annotation" and score >= threshold:
            call = "annotation_conflict_candidate"
        elif status == "supports_annotation":
            call = "supports_annotation"
        else:
            call = "ambiguous_evidence"
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
            }
        )
    write_tsv(
        f"{output_dir}/annotation_completion_candidates.tsv",
        rows,
        ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "interval", "annotation_status", "inferred_role", "support_score", "completion_call", "evidence_status"],
    )
    write_tsv(
        f"{output_dir}/annotation_completion_summary.tsv",
        [{"completion_call": key, "count": value} for key, value in sorted(summary.items())],
        ["completion_call", "count"],
    )
    return rows
