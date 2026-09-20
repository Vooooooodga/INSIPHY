"""mapping / match context: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.mapping.fields import EXON_LIKE_ROLES
from intraphy.mapping.fields import _transcript_id_set
from intraphy.mapping.policies import occurrence_copy_key
from intraphy.preparation.transcripts import transcript_sort_key
from intraphy.storage.tabular import read_tsv
from intraphy.storage.values import to_float


def segment_length(row):
    return max(1, int(row.get("end", "0")) - int(row.get("start", "0")) + 1)


def phase_score(left, right):
    left_phase = left.get("phase", ".")
    right_phase = right.get("phase", ".")
    if left_phase in {".", "NA", ""} or right_phase in {".", "NA", ""}:
        return 0.6
    return 1.0 if left_phase == right_phase else 0.1


def role_boundary_score(left, right):
    if left.get("role") == right.get("role"):
        return 1.0
    exon_like = {"CDS", "exon", "UTR", "noncoding_exon"}
    noncoding_like = {"intron", "intergenic"}
    if left.get("role") in exon_like and right.get("role") in exon_like:
        return 0.75
    if {left.get("role"), right.get("role")} & exon_like and {left.get("role"), right.get("role")} & noncoding_like:
        return 0.45
    return 0.25


def _path_order_context(occurrences, transcript_paths=None):
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    by_tx = defaultdict(list)
    for path in transcript_paths or []:
        occ = occ_by_id.get(path.get("occurrence_id"))
        if not occ:
            continue
        key = (occ["family_id"], occ["species"], occ["gene_copy_id"], path.get("transcript_id", "NA"))
        by_tx[key].append(path)
    if not transcript_paths:
        for occ in occurrences:
            for transcript in _transcript_id_set(occ):
                key = (occ["family_id"], occ["species"], occ["gene_copy_id"], transcript)
                by_tx[key].append(occ)
    by_occ = {
        row["occurrence_id"]: {"indices": [], "left_roles": set(), "right_roles": set(), "copy_sizes": []}
        for row in occurrences
    }
    for _key, rows in by_tx.items():
        if transcript_paths:
            ordered = sorted(rows, key=lambda row: int(row.get("path_rank", "0") or 0))
        else:
            ordered = sorted(rows, key=lambda row: transcript_sort_key(row, row.get("strand", "+")))
        total = max(1, len(ordered) - 1)
        for idx, path in enumerate(ordered):
            occ_id = path.get("occurrence_id")
            if occ_id not in by_occ:
                continue
            rec = by_occ[occ_id]
            rec["indices"].append(idx / total)
            rec["left_roles"].add(ordered[idx - 1].get("role", "unknown") if idx > 0 else "terminal")
            rec["right_roles"].add(ordered[idx + 1].get("role", "unknown") if idx + 1 < len(ordered) else "terminal")
            rec["copy_sizes"].append(len(ordered))
    return by_occ


def copy_order_context(occurrences, transcript_paths=None):
    by_copy = defaultdict(list)
    for row in occurrences:
        by_copy[(row["family_id"], row["species"], row["gene_copy_id"])].append(row)
    context = {}
    path_context = _path_order_context(occurrences, transcript_paths)
    for key, rows in by_copy.items():
        by_occ = {row["occurrence_id"]: {"indices": [], "left_roles": set(), "right_roles": set(), "copy_sizes": []} for row in rows}
        for row in rows:
            if path_context and path_context.get(row["occurrence_id"], {}).get("indices"):
                by_occ[row["occurrence_id"]] = path_context[row["occurrence_id"]]
        if not transcript_paths and not any(_transcript_id_set(row) for row in rows):
            strand = rows[0].get("strand", "+") if rows else "+"
            ordered = sorted(
                rows,
                key=lambda row: (
                    -int(row.get("end", "0")) if strand == "-" else int(row.get("start", "0")),
                    -int(row.get("start", "0")) if strand == "-" else int(row.get("end", "0")),
                ),
            )
            total = max(1, len(ordered) - 1)
            for idx, row in enumerate(ordered):
                rec = by_occ[row["occurrence_id"]]
                rec["indices"].append(idx / total)
                rec["left_roles"].add(ordered[idx - 1]["role"] if idx > 0 else "terminal")
                rec["right_roles"].add(ordered[idx + 1]["role"] if idx + 1 < len(ordered) else "terminal")
                rec["copy_sizes"].append(len(ordered))
        for occ_id, rec in by_occ.items():
            left_roles = rec["left_roles"] or {"unknown"}
            right_roles = rec["right_roles"] or {"unknown"}
            context[occ_id] = {
                "index": "NA",
                "scaled_index": sum(rec["indices"]) / len(rec["indices"]) if rec["indices"] else 0.5,
                "left_role": next(iter(left_roles)) if len(left_roles) == 1 else "alternative",
                "right_role": next(iter(right_roles)) if len(right_roles) == 1 else "alternative",
                "copy_size": max(rec["copy_sizes"]) if rec["copy_sizes"] else len(rows),
            }
    return context


def context_score(left_ctx, right_ctx, side):
    key = f"{side}_role"
    if {left_ctx.get(key), right_ctx.get(key)} & {"alternative", "unknown", None}:
        return 0.5
    if left_ctx.get(key) == right_ctx.get(key):
        return 1.0
    if "terminal" in {left_ctx.get(key), right_ctx.get(key)}:
        return 0.5
    return 0.25


def load_distance_table(distance_table):
    if not distance_table:
        return {}
    rows = read_tsv(distance_table, optional=True)
    out = {}
    for row in rows:
        sp1 = row.get("species1") or row.get("sp1") or row.get("left_species") or row.get("query_species")
        sp2 = row.get("species2") or row.get("sp2") or row.get("right_species") or row.get("subject_species")
        dist = row.get("distance") or row.get("evolutionary_distance") or row.get("range") or row.get("distance_class")
        if sp1 and sp2 and dist:
            out[frozenset([sp1, sp2])] = dist
    return out


def pair_threshold(left, right, base_threshold, distance_lookup):
    category = distance_lookup.get(frozenset([left.get("species", ""), right.get("species", "")]), "medium")
    category = category.lower()
    if category == "short":
        return max(base_threshold, 0.75), category
    if category == "long":
        return max(0.45, base_threshold - 0.15), category
    return base_threshold, category


def cheap_match_evidence(left, right, context, alignment_backend="prefilter"):
    left_ctx = context.get(left["occurrence_id"], {})
    right_ctx = context.get(right["occurrence_id"], {})
    order = 1.0 - abs(to_float(left_ctx.get("scaled_index"), 0.5) - to_float(right_ctx.get("scaled_index"), 0.5))
    left_context = context_score(left_ctx, right_ctx, "left")
    right_context = context_score(left_ctx, right_ctx, "right")
    boundary = role_boundary_score(left, right)
    phase = phase_score(left, right)
    strand = 1.0
    splice = 1.0 - abs(to_float(left.get("splice_motif_score"), 0.5) - to_float(right.get("splice_motif_score"), 0.5))
    size_ratio = min(segment_length(left), segment_length(right)) / max(segment_length(left), segment_length(right))
    sequence_score = 0.0
    structural_context_score = 0.5 * left_context + 0.5 * right_context
    total = (
        0.34 * 0.0
        + 0.14 * 0.0
        + 0.10 * left_context
        + 0.10 * right_context
        + 0.10 * boundary
        + 0.08 * phase
        + 0.06 * order
        + 0.04 * strand
        + 0.04 * splice
    )
    return {
        "alignment_score": 0.0,
        "coverage_score": 0.0,
        "sequence_score": sequence_score,
        "structural_context_score": structural_context_score,
        "left_context_score": left_context,
        "right_context_score": right_context,
        "boundary_score": boundary,
        "phase_score": phase,
        "order_score": order,
        "strand_score": strand,
        "splice_score": splice,
        "size_ratio": size_ratio,
        "alignment_cigar": "NA",
        "alignment_backend": alignment_backend,
        "alignment_mode": "none",
        "alignment_meaning": "pair not aligned",
        "alignment_requested_backend": "NA",
        "alignment_strand": "NA",
        "projected_reference_blocks": "NA",
        "total_score": total,
    }


def should_align_pair(left, right, seqs, min_size_ratio=0.25):
    if left["family_id"] != right["family_id"]:
        return False, "different_family"
    left_role = left.get("role", "")
    right_role = right.get("role", "")
    if occurrence_copy_key(left) == occurrence_copy_key(right):
        return False, "same_copy_projection_only"
    if left_role == "intron" and right_role == "intron":
        return False, "intron_intron_skipped"
    size_ratio = min(segment_length(left), segment_length(right)) / max(segment_length(left), segment_length(right))
    has_exon_like = left_role in EXON_LIKE_ROLES or right_role in EXON_LIKE_ROLES
    if not has_exon_like and size_ratio < min_size_ratio:
        return False, "length_ratio_prefilter"
    if not seqs.get(left["occurrence_id"]) or not seqs.get(right["occurrence_id"]):
        return False, "missing_sequence"
    return True, "aligned_candidate"
