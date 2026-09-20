"""preparation / copy_context: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict


def make_adjacencies(occurrences):
    rows = []
    by_copy = defaultdict(list)
    for row in occurrences:
        if row.get("presence_status") == "present":
            by_copy[(row["family_id"], row["species"], row["gene_copy_id"])].append(row)
    for (family, species, copy), vals in sorted(by_copy.items()):
        vals = sorted(vals, key=lambda row: int(row.get("transcript_order", "0") or 0))
        for idx, left in enumerate(vals[:-1], start=1):
            right = vals[idx]
            rows.append(
                {
                    "adjacency_id": f"{species}_{copy}_adj_{idx:03d}",
                    "family_id": family,
                    "species": species,
                    "gene_copy_id": copy,
                    "left_occurrence_id": left["occurrence_id"],
                    "right_occurrence_id": right["occurrence_id"],
                    "adjacency_status": "present",
                }
            )
    return rows


def copy_spans(occurrences):
    spans = {}
    for row in occurrences:
        if row.get("presence_status") != "present":
            continue
        key = (row["family_id"], row["species"], row["gene_copy_id"])
        start = int(row.get("start", "0"))
        end = int(row.get("end", "0"))
        if key not in spans:
            spans[key] = {"family_id": row["family_id"], "species": row["species"], "gene_copy_id": row["gene_copy_id"], "contig": row.get("contig", "NA"), "start": start, "end": end, "roles": set()}
        spans[key]["start"] = min(spans[key]["start"], start)
        spans[key]["end"] = max(spans[key]["end"], end)
        spans[key]["roles"].add(row.get("role", "unknown"))
    return spans


def infer_copy_class(copies):
    if len(copies) == 1:
        return "single_copy"
    contigs = {copy["contig"] for copy in copies}
    if len(contigs) == 1:
        ordered = sorted(copies, key=lambda row: row["start"])
        max_gap = max(max(0, ordered[idx]["start"] - ordered[idx - 1]["end"]) for idx in range(1, len(ordered))) if len(ordered) > 1 else 0
        return "tandem_multi_copy" if max_gap <= 250_000 else "same_contig_multi_copy"
    return "dispersed_multi_copy"


def make_copy_context(occurrences):
    spans = copy_spans(occurrences)
    by_species = defaultdict(list)
    for span in spans.values():
        by_species[(span["family_id"], span["species"])].append(span)
    rows = []
    for (family, species), copies in sorted(by_species.items()):
        copy_class = infer_copy_class(copies)
        for copy in sorted(copies, key=lambda row: row["gene_copy_id"]):
            roles = copy["roles"]
            if "intron" not in roles and len(roles & {"CDS", "exon", "UTR"}) > 0 and len(copies) > 1:
                subtype = "processed_or_intronless_copy_candidate"
            else:
                subtype = copy_class
            rows.append({"family_id": family, "species": species, "gene_copy_id": copy["gene_copy_id"], "copy_class": copy_class, "copy_subclass": subtype, "copy_span": f"{copy['contig']}:{copy['start']}-{copy['end']}"})
    return rows


def make_copy_relationships(occurrences):
    spans = copy_spans(occurrences)
    by_species = defaultdict(list)
    for span in spans.values():
        by_species[(span["family_id"], span["species"])].append(span)
    rows = []
    for (family, species), copies in sorted(by_species.items()):
        copies = sorted(copies, key=lambda row: (row["contig"], row["start"], row["gene_copy_id"]))
        for i, left in enumerate(copies):
            for right in copies[i + 1 :]:
                same_contig = left["contig"] == right["contig"]
                gap = max(0, max(left["start"], right["start"]) - min(left["end"], right["end"]))
                if same_contig and gap <= 250_000:
                    rel = "tandem_duplication_candidate"
                    score = 0.9
                elif same_contig:
                    rel = "same_contig_duplication_candidate"
                    score = 0.65
                else:
                    rel = "dispersed_or_retrocopy_candidate"
                    score = 0.45
                rows.append(
                    {
                        "family_id": family,
                        "species": species,
                        "query_copy_id": left["gene_copy_id"],
                        "subject_copy_id": right["gene_copy_id"],
                        "relationship_class": rel,
                        "synteny_score": f"{score:.6g}",
                        "distance_bp": gap if same_contig else "NA",
                        "evidence": "copy_span_geometry",
                    }
                )
    return rows
