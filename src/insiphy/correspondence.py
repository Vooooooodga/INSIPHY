"""Homologous segment group assignment and segment correspondence scoring."""

from collections import Counter, defaultdict

from .io import parse_fasta, read_tsv, to_float, uniq, write_tsv


def simple_identity(seq_a, seq_b):
    if not seq_a or not seq_b:
        return 0.0
    n = min(len(seq_a), len(seq_b))
    matches = sum(1 for a, b in zip(seq_a[:n], seq_b[:n]) if a == b and a not in "-N" and b not in "-N")
    return matches / max(1, n)


def match_total(row):
    if row.get("total_score") not in ("", "NA", None):
        return to_float(row.get("total_score"))
    return (
        0.25 * to_float(row.get("alignment_score"))
        + 0.15 * to_float(row.get("left_context_score"))
        + 0.15 * to_float(row.get("right_context_score"))
        + 0.20 * to_float(row.get("left_boundary_score"))
        + 0.20 * to_float(row.get("right_boundary_score"))
        + 0.05 * to_float(row.get("size_ratio"), 1.0)
    )


def infer_correspondence(input_dir, output_dir):
    homology = read_tsv(f"{input_dir}/segment_homology.tsv", ["homology_id", "occurrence_id", "support_type", "confidence"])
    occurrences = read_tsv(f"{input_dir}/segment_occurrences.tsv", ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"])
    matches = read_tsv(f"{input_dir}/segment_matches.tsv", ["match_id", "query_occurrence_id", "subject_occurrence_id", "match_status"], optional=True)
    seqs = parse_fasta(f"{input_dir}/segment_sequences.fasta")
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}

    hsg_rows = []
    by_hsg = defaultdict(list)
    for row in homology:
        occ = occ_by_id.get(row["occurrence_id"], {})
        by_hsg[row["homology_id"]].append(row["occurrence_id"])
        hsg_rows.append(
            {
                "homology_id": row["homology_id"],
                "occurrence_id": row["occurrence_id"],
                "family_id": occ.get("family_id", "NA"),
                "species": occ.get("species", "NA"),
                "gene_copy_id": occ.get("gene_copy_id", "NA"),
                "source_label": row.get("source_label", "NA"),
                "support_type": row.get("support_type", "NA"),
                "confidence": row.get("confidence", "NA"),
            }
        )

    conservation = []
    for hsg, occ_ids in sorted(by_hsg.items()):
        ids_with_seq = [occ_id for occ_id in occ_ids if occ_id in seqs]
        identities = []
        for i, left in enumerate(ids_with_seq):
            for right in ids_with_seq[i + 1 :]:
                identities.append(simple_identity(seqs[left], seqs[right]))
        mean_identity = sum(identities) / len(identities) if identities else 0.0
        roles = Counter(occ_by_id.get(occ_id, {}).get("role", "unknown") for occ_id in occ_ids)
        conservation.append(
            {
                "homology_id": hsg,
                "occurrence_count": len(occ_ids),
                "species_count": uniq(occ_by_id.get(occ_id, {}).get("species", "") for occ_id in occ_ids),
                "mean_pairwise_identity": f"{mean_identity:.6g}",
                "role_spectrum": ";".join(f"{key}:{roles[key]}" for key in sorted(roles)),
                "conservation_call": "conserved" if mean_identity >= 0.75 or len(occ_ids) >= 2 else "single_or_low_sequence_support",
            }
        )

    scored = []
    grouped = defaultdict(list)
    for row in matches:
        score = match_total(row)
        grouped[row["query_occurrence_id"]].append((score, row["match_id"]))
        scored.append(
            {
                "match_id": row["match_id"],
                "query_occurrence_id": row["query_occurrence_id"],
                "subject_occurrence_id": row["subject_occurrence_id"],
                "total_score": f"{score:.6g}",
                "match_status": row.get("match_status", "ambiguous"),
                "correspondence_call": "unclassified",
            }
        )
    best_status = {}
    for vals in grouped.values():
        best = max(score for score, _ in vals)
        best_ids = [match_id for score, match_id in vals if abs(score - best) < 1e-12]
        for match_id in best_ids:
            best_status[match_id] = "best_tie" if len(best_ids) > 1 else "best_unique"
    for row in scored:
        row["correspondence_call"] = row["match_status"] if row["match_status"] != "mapped" else best_status.get(row["match_id"], "not_best")

    write_tsv(f"{output_dir}/hsg_assignments.tsv", hsg_rows, ["homology_id", "occurrence_id", "family_id", "species", "gene_copy_id", "source_label", "support_type", "confidence"])
    write_tsv(f"{output_dir}/segment_conservation.tsv", conservation, ["homology_id", "occurrence_count", "species_count", "mean_pairwise_identity", "role_spectrum", "conservation_call"])
    write_tsv(f"{output_dir}/segment_correspondence.tsv", scored, ["match_id", "query_occurrence_id", "subject_occurrence_id", "total_score", "match_status", "correspondence_call"])
    return hsg_rows, conservation, scored
