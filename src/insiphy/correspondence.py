"""Homologous segment group assignment and segment correspondence scoring."""

from collections import Counter, defaultdict

from .alignment import global_alignment_stats
from .io import parse_fasta, read_tsv, to_float, uniq, write_tsv


def simple_identity(seq_a, seq_b):
    return global_alignment_stats(seq_a, seq_b).identity


def match_total(row):
    if row.get("total_score") not in ("", "NA", None):
        return to_float(row.get("total_score"))
    boundary = to_float(row.get("boundary_score"), None)
    if boundary is None:
        boundary = 0.5 * (to_float(row.get("left_boundary_score"), 0.5) + to_float(row.get("right_boundary_score"), 0.5))
    return (
        0.40 * to_float(row.get("alignment_score"))
        + 0.15 * to_float(row.get("coverage_score"), to_float(row.get("size_ratio"), 1.0))
        + 0.10 * to_float(row.get("left_context_score"), 0.5)
        + 0.10 * to_float(row.get("right_context_score"), 0.5)
        + 0.10 * boundary
        + 0.10 * to_float(row.get("phase_score"), 0.5)
        + 0.05 * to_float(row.get("order_score"), 0.5)
    )


def normalized_components(row):
    boundary = to_float(row.get("boundary_score"), None)
    if boundary is None:
        boundary = 0.5 * (to_float(row.get("left_boundary_score"), 0.5) + to_float(row.get("right_boundary_score"), 0.5))
    return {
        "alignment_score": to_float(row.get("alignment_score")),
        "coverage_score": to_float(row.get("coverage_score"), to_float(row.get("size_ratio"), 1.0)),
        "left_context_score": to_float(row.get("left_context_score"), 0.5),
        "right_context_score": to_float(row.get("right_context_score"), 0.5),
        "boundary_score": boundary,
        "phase_score": to_float(row.get("phase_score"), 0.5),
        "order_score": to_float(row.get("order_score"), 0.5),
        "strand_score": to_float(row.get("strand_score"), 0.5),
        "splice_score": to_float(row.get("splice_score"), 0.5),
    }


def reciprocal_status(scored):
    best_subject_by_query = defaultdict(list)
    best_query_by_subject = defaultdict(list)
    for row in scored:
        score = to_float(row["total_score"])
        q = row["query_occurrence_id"]
        s = row["subject_occurrence_id"]
        best_subject_by_query[q].append((score, s, row["match_id"]))
        best_query_by_subject[s].append((score, q, row["match_id"]))
        best_subject_by_query[s].append((score, q, row["match_id"]))
        best_query_by_subject[q].append((score, s, row["match_id"]))

    query_best = {}
    for query, vals in best_subject_by_query.items():
        best = max(score for score, _node, _mid in vals)
        query_best[query] = {_node for score, _node, _mid in vals if abs(score - best) < 1e-12}
    subject_best = {}
    for subject, vals in best_query_by_subject.items():
        best = max(score for score, _node, _mid in vals)
        subject_best[subject] = {_node for score, _node, _mid in vals if abs(score - best) < 1e-12}

    out = {}
    for row in scored:
        q = row["query_occurrence_id"]
        s = row["subject_occurrence_id"]
        if s in query_best.get(q, set()) and q in subject_best.get(s, set()):
            out[row["match_id"]] = "reciprocal_best"
        elif s in query_best.get(q, set()):
            out[row["match_id"]] = "query_best_only"
        else:
            out[row["match_id"]] = "not_reciprocal_best"
    return out


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
                "membership_score": row.get("confidence", "NA"),
                "membership_call": "core_member" if to_float(row.get("confidence"), 0.0) >= 0.7 else "ambiguous_member",
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
        components = normalized_components(row)
        grouped[row["query_occurrence_id"]].append((score, row["match_id"]))
        scored.append(
            {
                "match_id": row["match_id"],
                "query_occurrence_id": row["query_occurrence_id"],
                "subject_occurrence_id": row["subject_occurrence_id"],
                "alignment_score": f"{components['alignment_score']:.6g}",
                "coverage_score": f"{components['coverage_score']:.6g}",
                "left_context_score": f"{components['left_context_score']:.6g}",
                "right_context_score": f"{components['right_context_score']:.6g}",
                "boundary_score": f"{components['boundary_score']:.6g}",
                "phase_score": f"{components['phase_score']:.6g}",
                "order_score": f"{components['order_score']:.6g}",
                "strand_score": f"{components['strand_score']:.6g}",
                "splice_score": f"{components['splice_score']:.6g}",
                "total_score": f"{score:.6g}",
                "match_status": row.get("match_status", "ambiguous"),
                "distance_class": row.get("distance_class", "NA"),
                "threshold": row.get("threshold", "NA"),
                "alignment_cigar": row.get("alignment_cigar", "NA"),
                "reciprocal_status": "unclassified",
                "correspondence_call": "unclassified",
            }
        )
    best_status = {}
    for vals in grouped.values():
        best = max(score for score, _ in vals)
        best_ids = [match_id for score, match_id in vals if abs(score - best) < 1e-12]
        for match_id in best_ids:
            best_status[match_id] = "best_tie" if len(best_ids) > 1 else "best_unique"
    reciprocal = reciprocal_status(scored)
    graph_edges = []
    degree = Counter()
    for row in scored:
        row["reciprocal_status"] = reciprocal.get(row["match_id"], "not_reciprocal_best")
        if row["match_status"] == "mapped":
            degree[row["query_occurrence_id"]] += 1
            degree[row["subject_occurrence_id"]] += 1
            graph_edges.append(
                {
                    "edge_id": row["match_id"],
                    "query_occurrence_id": row["query_occurrence_id"],
                    "subject_occurrence_id": row["subject_occurrence_id"],
                    "total_score": row["total_score"],
                    "reciprocal_status": row["reciprocal_status"],
                    "edge_call": "high_confidence_correspondence" if row["reciprocal_status"] == "reciprocal_best" and to_float(row["total_score"]) >= 0.7 else "supporting_correspondence",
                }
            )
        if row["match_status"] != "mapped":
            row["correspondence_call"] = row["match_status"]
        elif row["reciprocal_status"] == "reciprocal_best":
            row["correspondence_call"] = "reciprocal_best"
        else:
            row["correspondence_call"] = best_status.get(row["match_id"], "not_best")

    for row in hsg_rows:
        deg = degree.get(row["occurrence_id"], 0)
        base = to_float(row.get("membership_score"), 0.5)
        adjusted = min(1.0, 0.75 * base + 0.25 * min(1.0, deg / 2))
        row["membership_score"] = f"{adjusted:.6g}"
        row["membership_call"] = "core_member" if adjusted >= 0.7 else "ambiguous_member"

    write_tsv(f"{output_dir}/hsg_assignments.tsv", hsg_rows, ["homology_id", "occurrence_id", "family_id", "species", "gene_copy_id", "source_label", "support_type", "confidence", "membership_score", "membership_call"])
    write_tsv(f"{output_dir}/segment_conservation.tsv", conservation, ["homology_id", "occurrence_count", "species_count", "mean_pairwise_identity", "role_spectrum", "conservation_call"])
    write_tsv(
        f"{output_dir}/segment_correspondence.tsv",
        scored,
        [
            "match_id",
            "query_occurrence_id",
            "subject_occurrence_id",
            "alignment_score",
            "coverage_score",
            "left_context_score",
            "right_context_score",
            "boundary_score",
            "phase_score",
            "order_score",
            "strand_score",
            "splice_score",
            "total_score",
            "distance_class",
            "threshold",
            "alignment_cigar",
            "match_status",
            "reciprocal_status",
            "correspondence_call",
        ],
    )
    write_tsv(f"{output_dir}/hsg_graph_edges.tsv", graph_edges, ["edge_id", "query_occurrence_id", "subject_occurrence_id", "total_score", "reciprocal_status", "edge_call"])
    return hsg_rows, conservation, scored
