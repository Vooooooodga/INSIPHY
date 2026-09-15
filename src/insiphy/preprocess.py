"""Genome annotation extraction and first-pass intragenic table builders."""

from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .alignment import available_alignment_backends, global_alignment_stats, phase_compatibility, splice_motif_score
from .io import open_text, parse_fasta, read_tsv, to_float, write_tsv


SEGMENT_FIELDS = [
    "occurrence_id",
    "family_id",
    "species",
    "gene_copy_id",
    "transcript_id",
    "role",
    "role_set",
    "presence_status",
    "contig",
    "start",
    "end",
    "strand",
    "phase",
    "source_feature_id",
    "boundary_class",
    "splice_motif_score",
    "splice_donor",
    "splice_acceptor",
    "frame_status",
]

SEGMENT_OUTPUT_FIELDS = SEGMENT_FIELDS + ["source_label", "copy_role"]

TRANSCRIPT_PATH_FIELDS = [
    "path_id",
    "family_id",
    "species",
    "gene_copy_id",
    "transcript_id",
    "path_rank",
    "occurrence_id",
    "role",
    "contig",
    "start",
    "end",
    "strand",
    "phase",
    "path_status",
]

INTRON_SITE_FIELDS = [
    "intron_id",
    "family_id",
    "species",
    "gene_copy_id",
    "transcript_id",
    "contig",
    "start",
    "end",
    "strand",
    "left_feature_id",
    "right_feature_id",
    "left_phase",
    "right_phase",
    "phase_compatibility",
    "splice_donor",
    "splice_acceptor",
    "splice_motif_score",
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


def split_ids(value):
    out = []
    for part in str(value or "").replace("|", ",").replace(";", ",").split(","):
        part = part.strip().strip('"')
        if part:
            out.append(part)
    return out


def read_annotation(path):
    rows = []
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            parts = raw.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            attrs = parse_attributes(parts[8])
            ftype = parts[2].lower()
            transcript_types = {"mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna"}
            if ftype == "gene":
                feature_id = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name") or ""
                parent = attrs.get("Parent") or ""
            elif ftype in transcript_types:
                feature_id = attrs.get("ID") or attrs.get("transcript_id") or attrs.get("Name") or ""
                parent = attrs.get("Parent") or attrs.get("gene_id") or ""
            else:
                feature_id = attrs.get("ID") or attrs.get("exon_id") or attrs.get("protein_id") or attrs.get("transcript_id") or attrs.get("gene_id") or ""
                parent = attrs.get("Parent") or attrs.get("transcript_id") or attrs.get("gene_id") or ""
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
                    "id": feature_id,
                    "parent": parent,
                    "parents": split_ids(parent),
                    "name": attrs.get("Name") or attrs.get("gene_name") or "",
                }
            )
    return rows


def feature_tokens(feature):
    attrs = feature.get("attrs", {})
    values = [feature.get("id", ""), feature.get("name", ""), feature.get("parent", "")]
    for key in ["ID", "Name", "Alias", "gene_id", "gene_name", "transcript_id", "Parent", "Dbxref"]:
        values.extend(split_ids(attrs.get(key, "")))
    tokens = set()
    for value in values:
        for token in split_ids(value):
            tokens.add(token)
            if ":" in token:
                tokens.add(token.split(":")[-1])
    return {token for token in tokens if token}


def feature_matches_gene(feature, gene_ids):
    return bool(feature_tokens(feature) & set(gene_ids))


def locate_gene(features, gene_id):
    gene_ids = {gene_id}
    genes = [row for row in features if row["type"].lower() == "gene" and feature_matches_gene(row, gene_ids)]
    if genes:
        gene = genes[0]
        return gene, gene_ids | feature_tokens(gene)
    children = [row for row in features if feature_matches_gene(row, gene_ids)]
    if not children:
        raise SystemExit(f"gene_id not found in annotation: {gene_id}")
    seqids = {row["seqid"] for row in children}
    if len(seqids) != 1:
        raise SystemExit(f"gene_id maps to multiple contigs: {gene_id}")
    strand = children[0]["strand"]
    gene = {
        "seqid": children[0]["seqid"],
        "type": "gene",
        "start": min(row["start"] for row in children),
        "end": max(row["end"] for row in children),
        "strand": strand,
        "phase": ".",
        "attrs": {"ID": gene_id, "Name": gene_id},
        "id": gene_id,
        "parent": "",
        "parents": [],
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


def transcript_features(features, gene, gene_ids):
    transcript_types = {"mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna"}
    transcripts = []
    for row in features:
        if row["type"].lower() not in transcript_types:
            continue
        if not overlaps_gene(row, gene):
            continue
        if set(row.get("parents", [])) & gene_ids or feature_matches_gene(row, gene_ids):
            transcripts.append(row)
    if not transcripts:
        transcripts = [
            {
                "seqid": gene["seqid"],
                "source": gene.get("source", "INSIPHY"),
                "type": "synthetic_transcript",
                "start": gene["start"],
                "end": gene["end"],
                "score": ".",
                "strand": gene["strand"],
                "phase": ".",
                "attrs": {"ID": f"{gene.get('id', 'gene')}.synthetic_tx"},
                "id": f"{gene.get('id', 'gene')}.synthetic_tx",
                "parent": gene.get("id", ""),
                "parents": [gene.get("id", "")],
                "name": f"{gene.get('id', 'gene')}.synthetic_tx",
            }
        ]
    return sorted(transcripts, key=lambda row: (row["start"], row["end"], row.get("id", "")))


def feature_role(feature):
    ftype = feature["type"].lower()
    if ftype == "cds":
        return "CDS"
    if ftype in {"utr", "five_prime_utr", "three_prime_utr"}:
        return "UTR"
    if ftype == "intron":
        return "intron"
    if ftype == "gene_body":
        return "exon"
    return "exon"


def child_features_for_transcript(features, gene, gene_ids, transcript):
    tx_id = transcript.get("id", "")
    selected = []
    allowed = {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr"}
    synthetic = transcript["type"] == "synthetic_transcript"
    for row in features:
        if row["type"].lower() not in allowed or not overlaps_gene(row, gene):
            continue
        parents = set(row.get("parents", []))
        if (tx_id and tx_id in parents) or (synthetic and (parents & gene_ids or feature_matches_gene(row, gene_ids))):
            selected.append(row)
    if not selected:
        selected = [
            {
                "seqid": gene["seqid"],
                "type": "gene_body",
                "start": gene["start"],
                "end": gene["end"],
                "strand": gene["strand"],
                "phase": ".",
                "id": gene.get("id", "gene_body"),
                "parent": tx_id,
                "parents": [tx_id],
                "attrs": {"ID": gene.get("id", "gene_body")},
            }
        ]
    has_cds = any(row["type"].lower() == "cds" for row in selected)
    if has_cds:
        selected = [row for row in selected if row["type"].lower() in {"cds", "utr", "five_prime_utr", "three_prime_utr"}]
    return sorted(selected, key=lambda row: (row["start"], row["end"], feature_role(row)))


def transcript_sort_key(row, strand):
    key = (row["start"], row["end"])
    if strand == "-":
        key = (-row["end"], -row["start"])
    return key


def transcript_cds_length(tx_features):
    return sum(max(0, row["end"] - row["start"] + 1) for row in tx_features if row["type"].lower() == "cds")


def select_transcripts(transcripts, features_by_tx, transcript_policy="canonical", canonical_rule="longest_cds"):
    if transcript_policy == "all":
        return transcripts
    if not transcripts:
        return []
    if canonical_rule == "longest_span":
        return [max(transcripts, key=lambda tx: (tx["end"] - tx["start"] + 1, tx.get("id", "")))]
    return [max(transcripts, key=lambda tx: (transcript_cds_length(features_by_tx[tx["id"]]), tx["end"] - tx["start"] + 1, tx.get("id", "")))]


def introns_from_path(path_features, gene, transcript_id, seqs=None):
    exonic = [row for row in path_features if feature_role(row) in {"CDS", "UTR", "exon"}]
    merged = []
    for exon in sorted(exonic, key=lambda row: (row["start"], row["end"])):
        if not merged or exon["start"] > merged[-1]["end"] + 1:
            merged.append({"start": exon["start"], "end": exon["end"], "left": exon, "right": exon})
        else:
            merged[-1]["end"] = max(merged[-1]["end"], exon["end"])
            merged[-1]["right"] = exon
    introns = []
    for idx, left in enumerate(merged[:-1], start=1):
        right = merged[idx]
        start = left["end"] + 1
        end = right["start"] - 1
        if start <= end:
            intron_seq = sequence_slice(seqs or {}, gene["seqid"], start, end, gene["strand"])
            motif_score, donor, acceptor = splice_motif_score(intron_seq)
            introns.append(
                {
                    "seqid": gene["seqid"],
                    "type": "intron",
                    "start": start,
                    "end": end,
                    "strand": gene["strand"],
                    "phase": ".",
                    "id": f"{transcript_id}.intron_{idx}",
                    "parent": transcript_id,
                    "parents": [transcript_id],
                    "attrs": {"ID": f"{transcript_id}.intron_{idx}"},
                    "left_feature_id": left["right"].get("id", "NA") or "NA",
                    "right_feature_id": right["left"].get("id", "NA") or "NA",
                    "left_phase": left["right"].get("phase", "."),
                    "right_phase": right["left"].get("phase", "."),
                    "splice_motif_score": motif_score,
                    "splice_donor": donor,
                    "splice_acceptor": acceptor,
                }
            )
    return introns


def exon_features_for_gene(features, gene, gene_ids):
    transcripts = transcript_features(features, gene, gene_ids)
    features_by_tx = {tx["id"]: child_features_for_transcript(features, gene, gene_ids, tx) for tx in transcripts}
    selected = select_transcripts(transcripts, features_by_tx)
    rows = []
    for tx in selected:
        rows.extend(features_by_tx[tx["id"]])
    return sorted(rows, key=lambda row: (row["start"], row["end"], feature_role(row)))


def load_existing_segments(output_dir):
    path = Path(output_dir) / "segment_occurrences.tsv"
    if not path.exists():
        return []
    return read_tsv(path, optional=True)


def _feature_key(feature):
    return (
        feature["seqid"],
        int(feature["start"]),
        int(feature["end"]),
        feature.get("strand", "."),
        feature_role(feature),
        feature.get("phase", "."),
    )


def _occurrence_id(species, gene_copy_id, index, role):
    safe_role = role.replace("/", "_")
    return f"{species}_{gene_copy_id}_{index:03d}_{safe_role}"


def extract_gene(
    genome_fasta,
    annotation_path,
    gene_id,
    family_id,
    species,
    gene_copy_id,
    output_dir,
    append=False,
    transcript_policy="canonical",
    canonical_rule="longest_cds",
    source_label="unknown_source",
    copy_role="candidate",
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(genome_fasta)
    features = read_annotation(annotation_path)
    gene, gene_ids = locate_gene(features, gene_id)
    gene_ids = set(gene_ids)
    transcripts = transcript_features(features, gene, gene_ids)
    features_by_tx = {tx["id"]: child_features_for_transcript(features, gene, gene_ids, tx) for tx in transcripts}
    selected = select_transcripts(transcripts, features_by_tx, transcript_policy, canonical_rule)

    rows = list(load_existing_segments(output_dir) if append else [])
    tx_path_rows = list(read_tsv(output_dir / "transcript_paths.tsv", optional=True) if append else [])
    intron_rows = list(read_tsv(output_dir / "intron_sites.tsv", optional=True) if append else [])

    unique = {}
    tx_paths = []
    intron_records = []
    for tx in selected:
        tx_id = tx.get("id", "") or f"{gene_id}.tx"
        child_features = features_by_tx[tx_id]
        introns = introns_from_path(child_features, gene, tx_id, seqs)
        path_features = sorted(child_features + introns, key=lambda row: transcript_sort_key(row, gene["strand"]))
        for rank, feat in enumerate(path_features, start=1):
            key = _feature_key(feat)
            unique.setdefault(key, {"feature": feat, "transcripts": set(), "source_ids": set()})
            unique[key]["transcripts"].add(tx_id)
            unique[key]["source_ids"].add(feat.get("id", "NA") or "NA")
            tx_paths.append((tx_id, rank, key, feat))
        for intron in introns:
            intron_records.append((tx_id, intron))

    start_index = len([row for row in rows if row.get("species") == species and row.get("gene_copy_id") == gene_copy_id])
    key_to_occ = {}
    fasta_path = output_dir / "segment_sequences.fasta"
    fasta_mode = "a" if append and fasta_path.exists() else "w"
    with fasta_path.open(fasta_mode) as fasta:
        for offset, key in enumerate(sorted(unique, key=lambda item: (item[1], item[2], item[4])), start=1):
            entry = unique[key]
            feat = entry["feature"]
            role = feature_role(feat)
            occ_id = _occurrence_id(species, gene_copy_id, start_index + offset, role)
            key_to_occ[key] = occ_id
            seq = sequence_slice(seqs, feat["seqid"], feat["start"], feat["end"], feat["strand"])
            motif_score = feat.get("splice_motif_score", "0.5" if role != "intron" else "0")
            donor = feat.get("splice_donor", "NA")
            acceptor = feat.get("splice_acceptor", "NA")
            rows.append(
                {
                    "occurrence_id": occ_id,
                    "family_id": family_id,
                    "species": species,
                    "gene_copy_id": gene_copy_id,
                    "transcript_id": ";".join(sorted(entry["transcripts"])),
                    "role": role,
                    "role_set": role,
                    "presence_status": "present",
                    "contig": feat["seqid"],
                    "start": feat["start"],
                    "end": feat["end"],
                    "strand": feat["strand"],
                    "phase": feat.get("phase", "."),
                    "source_feature_id": ";".join(sorted(entry["source_ids"])),
                    "boundary_class": "internal_intron" if role == "intron" else "annotated_segment",
                    "splice_motif_score": f"{to_float(motif_score):.6g}",
                    "splice_donor": donor,
                    "splice_acceptor": acceptor,
                    "frame_status": "coding_frame_annotated" if role == "CDS" and feat.get("phase", ".") not in {".", "NA"} else "not_coding_or_unknown",
                    "source_label": source_label or "unknown_source",
                    "copy_role": copy_role or "candidate",
                }
            )
            fasta.write(f">{occ_id}\n{seq}\n")

    for tx_id, rank, key, feat in tx_paths:
        tx_path_rows.append(
            {
                "path_id": f"{species}_{gene_copy_id}_{tx_id}_{rank:03d}",
                "family_id": family_id,
                "species": species,
                "gene_copy_id": gene_copy_id,
                "transcript_id": tx_id,
                "path_rank": rank,
                "occurrence_id": key_to_occ[key],
                "role": feature_role(feat),
                "contig": feat["seqid"],
                "start": feat["start"],
                "end": feat["end"],
                "strand": feat["strand"],
                "phase": feat.get("phase", "."),
                "path_status": "annotated_transcript_path" if transcript_policy == "all" else "canonical_transcript_path",
            }
        )
    for tx_id, intron in intron_records:
        intron_rows.append(
            {
                "intron_id": f"{species}_{gene_copy_id}_{intron['id']}",
                "family_id": family_id,
                "species": species,
                "gene_copy_id": gene_copy_id,
                "transcript_id": tx_id,
                "contig": intron["seqid"],
                "start": intron["start"],
                "end": intron["end"],
                "strand": intron["strand"],
                "left_feature_id": intron.get("left_feature_id", "NA"),
                "right_feature_id": intron.get("right_feature_id", "NA"),
                "left_phase": intron.get("left_phase", "."),
                "right_phase": intron.get("right_phase", "."),
                "phase_compatibility": phase_compatibility(intron.get("left_phase", "."), intron.get("right_phase", ".")),
                "splice_donor": intron.get("splice_donor", "NA"),
                "splice_acceptor": intron.get("splice_acceptor", "NA"),
                "splice_motif_score": f"{to_float(intron.get('splice_motif_score')):.6g}",
            }
        )

    write_tsv(output_dir / "segment_occurrences.tsv", rows, SEGMENT_OUTPUT_FIELDS)
    write_tsv(output_dir / "transcript_paths.tsv", tx_path_rows, TRANSCRIPT_PATH_FIELDS)
    write_tsv(output_dir / "intron_sites.tsv", intron_rows, INTRON_SITE_FIELDS)
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
    noncoding_like = {"intron", "regulatory", "intergenic", "noncoding", "intron_or_noncoding"}
    if left.get("role") in exon_like and right.get("role") in exon_like:
        return 0.75
    if {left.get("role"), right.get("role")} & exon_like and {left.get("role"), right.get("role")} & noncoding_like:
        return 0.45
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
    boundary = 0.5 * role_boundary_score(left, right) + 0.5 * phase_score(left, right)
    phase = phase_score(left, right)
    strand = 1.0 if left.get("strand") == right.get("strand") else 0.6
    splice = 1.0 - abs(to_float(left.get("splice_motif_score"), 0.5) - to_float(right.get("splice_motif_score"), 0.5))
    size_ratio = min(segment_length(left), segment_length(right)) / max(segment_length(left), segment_length(right))
    total = 0.10 * left_context + 0.10 * right_context + 0.10 * boundary + 0.08 * phase + 0.06 * order + 0.04 * strand + 0.04 * splice
    return {
        "alignment_score": 0.0,
        "coverage_score": 0.0,
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
        "total_score": total,
    }


def match_evidence(left, right, seqs, context, aligner="internal", threads=1):
    left_seq = seqs.get(left["occurrence_id"], "")
    right_seq = seqs.get(right["occurrence_id"], "")
    aln = global_alignment_stats(left_seq, right_seq, backend=aligner, threads=threads)
    identity = aln.identity
    coverage = aln.coverage
    left_ctx = context.get(left["occurrence_id"], {})
    right_ctx = context.get(right["occurrence_id"], {})
    order = 1.0 - abs(to_float(left_ctx.get("scaled_index"), 0.5) - to_float(right_ctx.get("scaled_index"), 0.5))
    left_context = context_score(left_ctx, right_ctx, "left")
    right_context = context_score(left_ctx, right_ctx, "right")
    boundary = 0.5 * role_boundary_score(left, right) + 0.5 * phase_score(left, right)
    phase = phase_score(left, right)
    strand = 1.0 if left.get("strand") == right.get("strand") else 0.6
    splice = 1.0 - abs(to_float(left.get("splice_motif_score"), 0.5) - to_float(right.get("splice_motif_score"), 0.5))
    total = 0.34 * identity + 0.14 * coverage + 0.10 * left_context + 0.10 * right_context + 0.10 * boundary + 0.08 * phase + 0.06 * order + 0.04 * strand + 0.04 * splice
    return {
        "alignment_score": identity,
        "coverage_score": coverage,
        "left_context_score": left_context,
        "right_context_score": right_context,
        "boundary_score": boundary,
        "phase_score": phase,
        "order_score": order,
        "strand_score": strand,
        "splice_score": splice,
        "size_ratio": min(segment_length(left), segment_length(right)) / max(segment_length(left), segment_length(right)),
        "alignment_cigar": aln.cigar,
        "alignment_backend": aln.backend,
        "total_score": total,
    }


EXON_LIKE_ROLES = {"CDS", "exon", "UTR", "noncoding_exon"}
NONCODING_ROLES = {"intron", "regulatory", "intergenic", "noncoding", "intron_or_noncoding"}
STRUCTURAL_ROLES = EXON_LIKE_ROLES | NONCODING_ROLES
UNKNOWN_SOURCE_LABELS = {"", "NA", "unknown", "unknown_source", "ambiguous", "unresolved"}


def occurrence_copy_key(row):
    return (row.get("family_id", ""), row.get("species", ""), row.get("gene_copy_id", ""))


def roles_compatible(left, right):
    left_role = left.get("role", "")
    right_role = right.get("role", "")
    if left_role == right_role:
        return True
    if left_role in EXON_LIKE_ROLES and right_role in EXON_LIKE_ROLES:
        return True
    return left_role in STRUCTURAL_ROLES and right_role in STRUCTURAL_ROLES


def sequence_supported_mapping(evidence, threshold):
    identity = evidence["alignment_score"]
    coverage = evidence["coverage_score"]
    size_ratio = evidence["size_ratio"]
    if size_ratio < 0.35:
        return False
    if identity >= threshold and coverage >= 0.45:
        return True
    if identity >= threshold + 0.15 and coverage >= 0.30:
        return True
    return False


def split_source_labels(value):
    labels = []
    for part in str(value or "").replace("|", ";").replace(",", ";").split(";"):
        label = part.strip()
        if label and label not in UNKNOWN_SOURCE_LABELS:
            labels.append(label)
    return sorted(set(labels))


def known_source_labels(row):
    if row.get("copy_role") == "derived":
        return []
    return split_source_labels(row.get("source_label"))


def inferred_source_label(row, support):
    own = split_source_labels(row.get("source_label"))
    if own and row.get("copy_role") != "derived":
        return ";".join(own)
    if not support:
        return "unknown_source"
    best = max(support.values())
    top = sorted(source for source, score in support.items() if score >= best * 0.90)
    return ";".join(top) if top else "unknown_source"


def graph_components(nodes, edges, occurrence_by_id=None):
    if occurrence_by_id:
        parent = {node: node for node in nodes}
        members = {
            node: {occurrence_copy_key(occurrence_by_id[node])}
            for node in nodes
            if node in occurrence_by_id
        }

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra == rb:
                return
            if members.get(ra, set()) & members.get(rb, set()):
                return
            parent[rb] = ra
            members[ra] = members.get(ra, set()) | members.get(rb, set())

        for left, right, _score in sorted(edges, key=lambda row: row[2], reverse=True):
            union(left, right)
        groups = defaultdict(list)
        for node in nodes:
            groups[find(node)].append(node)
        return [sorted(vals) for vals in groups.values()]

    try:
        import networkx as nx  # type: ignore

        graph = nx.Graph()
        graph.add_nodes_from(nodes)
        graph.add_weighted_edges_from(edges)
        return [sorted(comp) for comp in nx.connected_components(graph)]
    except Exception:
        parent = {node: node for node in nodes}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        for left, right, _score in edges:
            union(left, right)
        groups = defaultdict(list)
        for node in nodes:
            groups[find(node)].append(node)
        return [sorted(vals) for vals in groups.values()]


def should_align_pair(left, right, seqs, min_size_ratio=0.25):
    if left["family_id"] != right["family_id"]:
        return False, "different_family"
    if occurrence_copy_key(left) == occurrence_copy_key(right):
        return False, "same_copy_excluded"
    size_ratio = min(segment_length(left), segment_length(right)) / max(segment_length(left), segment_length(right))
    if size_ratio < min_size_ratio:
        return False, "length_ratio_prefilter"
    if not seqs.get(left["occurrence_id"]) or not seqs.get(right["occurrence_id"]):
        return False, "missing_sequence"
    return True, "aligned_candidate"


def cluster_segments(occurrences, seqs, identity_threshold=0.7, distance_table=None, aligner="internal", threads=1, min_size_ratio=0.25):
    context = copy_order_context(occurrences)
    distance_lookup = load_distance_table(distance_table)
    occurrence_by_id = {row["occurrence_id"]: row for row in occurrences}
    matches = []
    accepted_edges = []
    score_by_occ = defaultdict(list)
    source_support = defaultdict(lambda: defaultdict(float))
    raw_pairs = [(left, right) for i, left in enumerate(occurrences) for right in occurrences[i + 1 :] if left["family_id"] == right["family_id"]]
    pairs = [(idx, left, right) for idx, (left, right) in enumerate(raw_pairs)]

    def score_pair(item):
        idx, left, right = item
        should_align, prefilter_status = should_align_pair(left, right, seqs, min_size_ratio)
        if should_align:
            evidence = match_evidence(left, right, seqs, context, aligner=aligner, threads=1)
        else:
            evidence = cheap_match_evidence(left, right, context, alignment_backend=prefilter_status)
        return idx, left, right, evidence, prefilter_status

    scored_pairs = []
    worker_count = max(1, int(threads or 1))
    if worker_count > 1 and len(pairs) > 1:
        with ThreadPoolExecutor(max_workers=worker_count) as pool:
            futures = [pool.submit(score_pair, pair) for pair in pairs]
            for future in as_completed(futures):
                scored_pairs.append(future.result())
        scored_pairs.sort(key=lambda row: row[0])
    else:
        scored_pairs = [score_pair(pair) for pair in pairs]

    for _idx, left, right, evidence, prefilter_status in scored_pairs:
            threshold, distance_class = pair_threshold(left, right, identity_threshold, distance_lookup)
            score = evidence["total_score"]
            same_copy = occurrence_copy_key(left) == occurrence_copy_key(right)
            compatible = roles_compatible(left, right)
            sequence_ok = sequence_supported_mapping(evidence, threshold)
            mapped = (not same_copy) and compatible and sequence_ok and score >= threshold
            if mapped:
                accepted_edges.append((left["occurrence_id"], right["occurrence_id"], score))
                score_by_occ[left["occurrence_id"]].append(score)
                score_by_occ[right["occurrence_id"]].append(score)
                left_sources = known_source_labels(left)
                right_sources = known_source_labels(right)
                if left_sources and not right_sources:
                    for source in left_sources:
                        source_support[right["occurrence_id"]][source] += score
                if right_sources and not left_sources:
                    for source in right_sources:
                        source_support[left["occurrence_id"]][source] += score
                status = "mapped"
            else:
                status = prefilter_status if prefilter_status != "aligned_candidate" else "low_similarity"
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
                    "strand_score": f"{evidence['strand_score']:.6g}",
                    "splice_score": f"{evidence['splice_score']:.6g}",
                    "size_ratio": f"{evidence['size_ratio']:.6g}",
                    "total_score": f"{score:.6g}",
                    "distance_class": distance_class,
                    "threshold": f"{threshold:.6g}",
                    "alignment_cigar": evidence["alignment_cigar"],
                    "alignment_backend": evidence["alignment_backend"],
                    "match_status": status,
                }
            )
    nodes = [row["occurrence_id"] for row in occurrences]
    components = graph_components(nodes, accepted_edges, occurrence_by_id)
    homology = []
    for idx, occ_ids in enumerate(sorted(components, key=lambda vals: vals[0]), start=1):
        hsg = f"HSG_{idx:04d}"
        for occ_id in occ_ids:
            scores = score_by_occ.get(occ_id, [])
            confidence = sum(scores) / len(scores) if scores else 0.5
            homology.append(
                {
                    "homology_id": hsg,
                    "occurrence_id": occ_id,
                    "support_type": "sequence_boundary_context_graph",
                    "confidence": f"{confidence:.6g}",
                    "source_label": inferred_source_label(occurrence_by_id.get(occ_id, {}), source_support.get(occ_id, {})),
                }
            )
    return homology, matches


def derive_tables(input_dir, output_dir=None, identity_threshold=0.7, distance_table=None, aligner="internal", threads=1, min_size_ratio=0.25):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir or input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", SEGMENT_FIELDS)
    seqs = parse_fasta(input_dir / "segment_sequences.fasta")
    homology, matches = cluster_segments(occurrences, seqs, identity_threshold, distance_table, aligner=aligner, threads=threads, min_size_ratio=min_size_ratio)
    write_tsv(output_dir / "segment_homology.tsv", homology, ["homology_id", "occurrence_id", "support_type", "confidence", "source_label"])
    match_fields = [
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
        "size_ratio",
        "total_score",
        "distance_class",
        "threshold",
        "alignment_cigar",
        "alignment_backend",
        "match_status",
    ]
    write_tsv(output_dir / "segment_matches.tsv", matches, match_fields)
    backend_rows = []
    for row in available_alignment_backends():
        backend_rows.append(
            {
                **row,
                "selected": int(row["aligner"] == aligner),
                "threads": threads,
                "min_size_ratio": f"{min_size_ratio:.6g}",
            }
        )
    write_tsv(output_dir / "alignment_backend_report.tsv", backend_rows, ["aligner", "available", "selected", "threads", "min_size_ratio", "notes"])
    write_tsv(output_dir / "physical_adjacencies.tsv", make_adjacencies(occurrences), ["adjacency_id", "family_id", "species", "gene_copy_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status"])
    write_tsv(output_dir / "copy_context.tsv", make_copy_context(occurrences), ["family_id", "species", "gene_copy_id", "copy_class", "copy_subclass", "copy_span"])
    write_tsv(output_dir / "copy_relationships.tsv", make_copy_relationships(occurrences), ["family_id", "species", "query_copy_id", "subject_copy_id", "relationship_class", "synteny_score", "distance_bp", "evidence"])
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
                    "splice_motif_score": row.get("splice_motif_score", "0.5"),
                    "phase_compatibility": "compatible" if row.get("phase") not in {".", "NA", ""} else "unknown",
                    "inferred_event": "annotated_segment",
                    "frame_status": row.get("frame_status", "unknown"),
                }
            )
        write_tsv(
            output_dir / "sequence_synteny_evidence.tsv",
            evidence,
            [
                "evidence_id",
                "family_id",
                "species",
                "gene_copy_id",
                "homology_id",
                "annotation_status",
                "evidence_status",
                "inferred_role",
                "contig",
                "start",
                "end",
                "strand",
                "sequence_score",
                "left_synteny_score",
                "right_synteny_score",
                "splice_motif_score",
                "phase_compatibility",
                "inferred_event",
                "frame_status",
            ],
        )
    return homology, matches
