"""Genome annotation extraction and first-pass intragenic table builders."""

from collections import defaultdict
from pathlib import Path

from .correspondence import simple_identity
from .io import parse_fasta, read_tsv, to_float, write_tsv


SEGMENT_FIELDS = [
    "occurrence_id",
    "family_id",
    "species",
    "gene_copy_id",
    "role",
    "presence_status",
    "contig",
    "start",
    "end",
    "strand",
    "phase",
    "source_feature_id",
]


def parse_attributes(raw):
    attrs = {}
    for part in raw.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = value.strip().strip('"')
    return attrs


def read_annotation(path):
    rows = []
    with Path(path).open() as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            attrs = parse_attributes(parts[8])
            rows.append(
                {
                    "seqid": parts[0],
                    "source": parts[1],
                    "type": parts[2],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "score": parts[5],
                    "strand": parts[6],
                    "phase": parts[7],
                    "attrs": attrs,
                    "id": attrs.get("ID") or attrs.get("gene_id") or attrs.get("transcript_id") or "",
                    "parent": attrs.get("Parent") or attrs.get("gene_id") or "",
                    "name": attrs.get("Name") or attrs.get("gene_name") or "",
                }
            )
    return rows


def feature_matches_gene(feature, gene_ids):
    values = {feature.get("id", ""), feature.get("name", "")}
    values.update(v for v in feature.get("parent", "").replace(",", ";").split(";") if v)
    return bool(values & gene_ids)


def locate_gene(features, gene_id):
    gene_ids = {gene_id}
    genes = [row for row in features if row["type"].lower() == "gene" and feature_matches_gene(row, gene_ids)]
    if genes:
        return genes[0], gene_ids | {genes[0].get("id", ""), genes[0].get("name", "")}
    children = [row for row in features if feature_matches_gene(row, gene_ids)]
    if not children:
        raise SystemExit(f"gene_id not found in annotation: {gene_id}")
    seqids = {row["seqid"] for row in children}
    if len(seqids) != 1:
        raise SystemExit(f"gene_id maps to multiple contigs: {gene_id}")
    strand = children[0]["strand"]
    gene = {
        "seqid": children[0]["seqid"],
        "start": min(row["start"] for row in children),
        "end": max(row["end"] for row in children),
        "strand": strand,
        "id": gene_id,
        "name": gene_id,
    }
    return gene, gene_ids


def overlaps_gene(feature, gene):
    return feature["seqid"] == gene["seqid"] and feature["start"] <= gene["end"] and feature["end"] >= gene["start"]


def sequence_slice(seqs, contig, start, end, strand):
    seq = seqs.get(contig, "")
    if not seq:
        return ""
    sub = seq[start - 1 : end]
    if strand == "-":
        table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        sub = sub.translate(table)[::-1]
    return sub.upper()


def exon_features_for_gene(features, gene, gene_ids):
    rows = [
        row
        for row in features
        if row["type"].lower() in {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr"}
        and overlaps_gene(row, gene)
        and (feature_matches_gene(row, gene_ids) or gene["start"] <= row["start"] <= gene["end"])
    ]
    has_cds = any(row["type"].lower() == "cds" for row in rows)
    if has_cds:
        rows = [row for row in rows if row["type"].lower() in {"cds", "utr", "five_prime_utr", "three_prime_utr"}]
    if not rows:
        rows = [
            {
                "seqid": gene["seqid"],
                "type": "gene_body",
                "start": gene["start"],
                "end": gene["end"],
                "strand": gene["strand"],
                "phase": ".",
                "id": gene.get("id", "gene_body"),
            }
        ]
    return sorted(rows, key=lambda row: (row["start"], row["end"]))


def introns_from_exons(exons, gene):
    merged = []
    for exon in sorted(exons, key=lambda row: (row["start"], row["end"])):
        if not merged or exon["start"] > merged[-1]["end"] + 1:
            merged.append({"start": exon["start"], "end": exon["end"]})
        else:
            merged[-1]["end"] = max(merged[-1]["end"], exon["end"])
    introns = []
    for idx, left in enumerate(merged[:-1], start=1):
        right = merged[idx]
        start = left["end"] + 1
        end = right["start"] - 1
        if start <= end:
            introns.append(
                {
                    "seqid": gene["seqid"],
                    "type": "intron",
                    "start": start,
                    "end": end,
                    "strand": gene["strand"],
                    "phase": ".",
                    "id": f"intron_{idx}",
                }
            )
    return introns


def load_existing_segments(output_dir):
    path = Path(output_dir) / "segment_occurrences.tsv"
    if not path.exists():
        return []
    return read_tsv(path, SEGMENT_FIELDS)


def extract_gene(genome_fasta, annotation_path, gene_id, family_id, species, gene_copy_id, output_dir, append=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(genome_fasta)
    features = read_annotation(annotation_path)
    gene, gene_ids = locate_gene(features, gene_id)
    exons = exon_features_for_gene(features, gene, gene_ids)
    introns = introns_from_exons([row for row in exons if row["type"].lower() in {"exon", "cds", "gene_body"}], gene)
    ordered = sorted(exons + introns, key=lambda row: (row["start"], row["end"]))

    rows = load_existing_segments(output_dir) if append else []
    fasta_path = output_dir / "segment_sequences.fasta"
    fasta_mode = "a" if append and fasta_path.exists() else "w"
    with fasta_path.open(fasta_mode) as fasta:
        for idx, feat in enumerate(ordered, start=1):
            role = "CDS" if feat["type"].lower() == "cds" else "intron" if feat["type"].lower() == "intron" else "exon"
            occ_id = f"{species}_{gene_copy_id}_{idx:03d}_{role}"
            rows.append(
                {
                    "occurrence_id": occ_id,
                    "family_id": family_id,
                    "species": species,
                    "gene_copy_id": gene_copy_id,
                    "role": role,
                    "presence_status": "present",
                    "contig": feat["seqid"],
                    "start": feat["start"],
                    "end": feat["end"],
                    "strand": feat["strand"],
                    "phase": feat.get("phase", "."),
                    "source_feature_id": feat.get("id", "NA") or "NA",
                }
            )
            fasta.write(f">{occ_id}\n{sequence_slice(seqs, feat['seqid'], feat['start'], feat['end'], feat['strand'])}\n")
    write_tsv(output_dir / "segment_occurrences.tsv", rows, SEGMENT_FIELDS)
    return rows


def make_adjacencies(occurrences):
    rows = []
    by_copy = defaultdict(list)
    for row in occurrences:
        if row.get("presence_status") == "present":
            by_copy[(row["family_id"], row["species"], row["gene_copy_id"])].append(row)
    for (family, species, copy), vals in sorted(by_copy.items()):
        vals = sorted(vals, key=lambda row: (row.get("contig", ""), int(row.get("start", "0")), int(row.get("end", "0"))))
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


def make_copy_context(occurrences):
    rows = []
    by_species = defaultdict(set)
    for row in occurrences:
        if row.get("presence_status") == "present":
            by_species[(row["family_id"], row["species"])].add(row["gene_copy_id"])
    for (family, species), copies in sorted(by_species.items()):
        copy_class = "single_copy" if len(copies) == 1 else "tandem_multi_copy"
        for copy in sorted(copies):
            rows.append({"family_id": family, "species": species, "gene_copy_id": copy, "copy_class": copy_class})
    return rows


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
    if left.get("role") in exon_like and right.get("role") in exon_like:
        return 0.75
    return 0.25


def copy_order_context(occurrences):
    by_copy = defaultdict(list)
    for row in occurrences:
        by_copy[(row["family_id"], row["species"], row["gene_copy_id"])].append(row)
    context = {}
    for key, rows in by_copy.items():
        ordered = sorted(rows, key=lambda row: (row.get("contig", ""), int(row.get("start", "0")), int(row.get("end", "0"))))
        total = max(1, len(ordered) - 1)
        for idx, row in enumerate(ordered):
            context[row["occurrence_id"]] = {
                "index": idx,
                "scaled_index": idx / total,
                "left_role": ordered[idx - 1]["role"] if idx > 0 else "terminal",
                "right_role": ordered[idx + 1]["role"] if idx + 1 < len(ordered) else "terminal",
                "copy_size": len(ordered),
            }
    return context


def context_score(left_ctx, right_ctx, side):
    key = f"{side}_role"
    if left_ctx.get(key) == right_ctx.get(key):
        return 1.0
    if "terminal" in {left_ctx.get(key), right_ctx.get(key)}:
        return 0.5
    return 0.25


def match_evidence(left, right, seqs, context):
    left_seq = seqs.get(left["occurrence_id"], "")
    right_seq = seqs.get(right["occurrence_id"], "")
    identity = simple_identity(left_seq, right_seq)
    coverage = min(len(left_seq), len(right_seq)) / max(1, max(len(left_seq), len(right_seq)))
    left_ctx = context.get(left["occurrence_id"], {})
    right_ctx = context.get(right["occurrence_id"], {})
    order = 1.0 - abs(to_float(left_ctx.get("scaled_index"), 0.5) - to_float(right_ctx.get("scaled_index"), 0.5))
    left_context = context_score(left_ctx, right_ctx, "left")
    right_context = context_score(left_ctx, right_ctx, "right")
    boundary = 0.5 * role_boundary_score(left, right) + 0.5 * phase_score(left, right)
    phase = phase_score(left, right)
    total = 0.40 * identity + 0.15 * coverage + 0.10 * left_context + 0.10 * right_context + 0.10 * boundary + 0.10 * phase + 0.05 * order
    return {
        "alignment_score": identity,
        "coverage_score": coverage,
        "left_context_score": left_context,
        "right_context_score": right_context,
        "boundary_score": boundary,
        "phase_score": phase,
        "order_score": order,
        "size_ratio": min(segment_length(left), segment_length(right)) / max(segment_length(left), segment_length(right)),
        "total_score": total,
    }


def cluster_segments(occurrences, seqs, identity_threshold=0.7):
    parent = {row["occurrence_id"]: row["occurrence_id"] for row in occurrences}
    context = copy_order_context(occurrences)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    matches = []
    for i, left in enumerate(occurrences):
        for right in occurrences[i + 1 :]:
            if left["family_id"] != right["family_id"]:
                continue
            evidence = match_evidence(left, right, seqs, context)
            score = evidence["total_score"]
            same_role = left.get("role") == right.get("role")
            if score >= identity_threshold or (same_role and evidence["alignment_score"] >= identity_threshold - 0.1 and evidence["coverage_score"] >= 0.5):
                union(left["occurrence_id"], right["occurrence_id"])
                status = "mapped"
            else:
                status = "low_similarity"
            matches.append(
                {
                    "match_id": f"match_{len(matches) + 1:05d}",
                    "query_occurrence_id": left["occurrence_id"],
                    "subject_occurrence_id": right["occurrence_id"],
                    "alignment_score": f"{evidence['alignment_score']:.6g}",
                    "coverage_score": f"{evidence['coverage_score']:.6g}",
                    "left_context_score": f"{evidence['left_context_score']:.6g}",
                    "right_context_score": f"{evidence['right_context_score']:.6g}",
                    "boundary_score": f"{evidence['boundary_score']:.6g}",
                    "phase_score": f"{evidence['phase_score']:.6g}",
                    "order_score": f"{evidence['order_score']:.6g}",
                    "size_ratio": f"{evidence['size_ratio']:.6g}",
                    "total_score": f"{score:.6g}",
                    "match_status": status,
                }
            )
    clusters = defaultdict(list)
    for row in occurrences:
        clusters[find(row["occurrence_id"])].append(row["occurrence_id"])
    cluster_id = {root: f"HSG_{idx:04d}" for idx, root in enumerate(sorted(clusters), start=1)}
    homology = []
    for root, occ_ids in sorted(clusters.items()):
        for occ_id in sorted(occ_ids):
            homology.append(
                {
                    "homology_id": cluster_id[root],
                    "occurrence_id": occ_id,
                    "support_type": "sequence_context_cluster",
                    "confidence": "0.75",
                    "source_label": "unknown_source",
                }
            )
    return homology, matches


def derive_tables(input_dir, output_dir=None, identity_threshold=0.7):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir or input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", SEGMENT_FIELDS)
    seqs = parse_fasta(input_dir / "segment_sequences.fasta")
    homology, matches = cluster_segments(occurrences, seqs, identity_threshold)
    write_tsv(output_dir / "segment_homology.tsv", homology, ["homology_id", "occurrence_id", "support_type", "confidence", "source_label"])
    write_tsv(
        output_dir / "segment_matches.tsv",
        matches,
        ["match_id", "query_occurrence_id", "subject_occurrence_id", "alignment_score", "coverage_score", "left_context_score", "right_context_score", "boundary_score", "phase_score", "order_score", "size_ratio", "total_score", "match_status"],
    )
    write_tsv(output_dir / "physical_adjacencies.tsv", make_adjacencies(occurrences), ["adjacency_id", "family_id", "species", "gene_copy_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status"])
    write_tsv(output_dir / "copy_context.tsv", make_copy_context(occurrences), ["family_id", "species", "gene_copy_id", "copy_class"])
    if not (output_dir / "sequence_synteny_evidence.tsv").exists():
        evidence = []
        for row in occurrences:
            evidence.append(
                {
                    "evidence_id": f"ev_{row['occurrence_id']}",
                    "family_id": row["family_id"],
                    "species": row["species"],
                    "gene_copy_id": row["gene_copy_id"],
                    "homology_id": "NA",
                    "annotation_status": "annotated",
                    "evidence_status": "supports_annotation",
                    "inferred_role": row["role"],
                    "contig": row.get("contig", "NA"),
                    "start": row.get("start", "NA"),
                    "end": row.get("end", "NA"),
                    "strand": row.get("strand", "NA"),
                    "sequence_score": "1.0",
                    "left_synteny_score": "1.0",
                    "right_synteny_score": "1.0",
                    "splice_motif_score": "0.5",
                    "phase_compatibility": "compatible",
                }
            )
        write_tsv(
            output_dir / "sequence_synteny_evidence.tsv",
            evidence,
            ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "annotation_status", "evidence_status", "inferred_role", "contig", "start", "end", "strand", "sequence_score", "left_synteny_score", "right_synteny_score", "splice_motif_score", "phase_compatibility"],
        )
    return homology, matches
