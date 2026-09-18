"""Genome annotation extraction and first-pass intragenic table builders."""

from __future__ import annotations

import math
import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from itertools import islice
from pathlib import Path

from .alignment import available_alignment_backends, local_alignment_stats, overlap_alignment_stats, phase_compatibility, splice_motif_score
from .coding_correspondence import CodingProjectionIndex
from .io import fasta_record_length, open_text, parse_fasta, read_fasta_interval, read_tsv, to_float, write_tsv


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
SEGMENT_OUTPUT_FIELDS += [
    "transcript_order",
    "coding_status",
    "cds_start",
    "cds_end",
    "cds_length",
    "cds_phase",
    "cds_intervals",
    "utr_intervals",
    "utr_status",
    "source_feature_type",
    "source_parent",
    "source_parents",
    "feature_ownership",
]

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
    "transcript_order",
    "coding_status",
    "cds_length",
    "cds_phase",
    "cds_intervals",
    "utr_intervals",
    "path_status",
]

RAW_FEATURE_FIELDS = [
    "family_id",
    "species",
    "gene_copy_id",
    "gene_id",
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "strand",
    "phase",
    "id",
    "name",
    "parent",
    "parents",
    "ownership",
    "attrs",
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
    "left_cds_length",
    "expected_right_phase",
    "phase_compatibility",
    "splice_donor",
    "splice_acceptor",
    "splice_motif_score",
]

GENE_LOCUS_FIELDS = [
    "species",
    "gene_copy_id",
    "contig",
    "strand",
    "annotation_start",
    "annotation_end",
    "linked_start",
    "linked_end",
    "search_start",
    "search_end",
    "contig_length",
    "genome_fasta",
    "max_extension",
    "range_status",
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
    for part in str(value or "").replace(";", ",").split(","):
        part = part.strip().strip('"')
        if part:
            out.append(part)
    return out


def iter_annotation(path):
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
            yield {
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


def read_annotation(path):
    return list(iter_annotation(path))


def _annotation_row_key(row):
    return (
        row.get("seqid", ""),
        row.get("type", ""),
        int(row.get("start", 0)),
        int(row.get("end", 0)),
        row.get("strand", "."),
        row.get("id", ""),
        row.get("parent", ""),
    )


def _copy_annotation_row(row):
    copied = dict(row)
    copied["attrs"] = dict(row.get("attrs", {}))
    copied["parents"] = list(row.get("parents", []))
    return copied


@lru_cache(maxsize=128)
def _read_annotation_for_gene_cached(path, gene_id):
    """Collect one gene hierarchy with bounded per-gene GFF/GTF caching."""
    requested = {gene_id}
    exact_genes = []
    alias_genes = []
    matching_children = []
    for row in iter_annotation(path):
        if not feature_matches_gene(row, requested):
            continue
        if row["type"].lower() == "gene":
            attrs = row.get("attrs", {})
            exact_ids = {row.get("id", ""), attrs.get("ID", ""), attrs.get("gene_id", "")}
            if gene_id in exact_ids:
                exact_genes.append(row)
            else:
                alias_genes.append(row)
        else:
            matching_children.append(row)
    if len(exact_genes) > 1:
        raise SystemExit(f"gene_id maps to multiple exact gene loci: {gene_id}")
    if exact_genes:
        gene = exact_genes[0]
    else:
        if len(alias_genes) > 1:
            loci = ",".join(
                f"{row.get('seqid')}:{row.get('start')}-{row.get('end')}:{row.get('id') or row.get('name')}"
                for row in alias_genes
            )
            raise SystemExit(f"gene_id ambiguously matches multiple gene loci: {gene_id} ({loci})")
        gene = alias_genes[0] if alias_genes else None
    if gene is None:
        if not matching_children:
            raise SystemExit(f"gene_id not found in annotation: {gene_id}")
        gene, _gene_ids = locate_gene(matching_children, gene_id)
    else:
        matching_children = []
    # Aliases select the seed; only exact identifiers extend its Parent hierarchy.
    gene_ids = {
        value for value in (
            gene.get("id"), gene.get("attrs", {}).get("ID"), gene.get("attrs", {}).get("gene_id")
        ) if value
    }

    linked = {}
    link_ids = set(gene_ids)
    if gene:
        linked[_annotation_row_key(gene)] = gene
    for row in matching_children:
        if row["seqid"] == gene["seqid"]:
            linked[_annotation_row_key(row)] = row
            if row.get("id"):
                link_ids.add(row["id"])

    changed = True
    while changed:
        changed = False
        for row in iter_annotation(path):
            if row["seqid"] != gene["seqid"] or row["type"].lower() == "gene":
                continue
            parents = set(row.get("parents", []))
            if not parents & link_ids:
                continue
            key = _annotation_row_key(row)
            if key not in linked:
                linked[key] = row
                changed = True
            before = len(link_ids)
            if row.get("id"):
                link_ids.add(row["id"])
            if len(link_ids) != before:
                changed = True

    if not linked:
        linked[_annotation_row_key(gene)] = gene
    linked_start = min(row["start"] for row in linked.values())
    linked_end = max(row["end"] for row in linked.values())
    region_rows = [
        row for row in iter_annotation(path)
        if row["seqid"] == gene["seqid"] and row["start"] <= linked_end and row["end"] >= linked_start
    ]
    for row in linked.values():
        region_rows.append(row)
    deduped = {}
    for row in region_rows:
        deduped[_annotation_row_key(row)] = row
    bounds = {
        "annotation_start": gene["start"],
        "annotation_end": gene["end"],
        "linked_start": linked_start,
        "linked_end": linked_end,
    }
    return tuple(deduped.values()), gene, frozenset(gene_ids | link_ids), bounds


def read_annotation_for_gene(path, gene_id):
    rows, gene, gene_ids, bounds = _read_annotation_for_gene_cached(str(Path(path)), gene_id)
    return (
        [_copy_annotation_row(row) for row in rows],
        _copy_annotation_row(gene),
        set(gene_ids),
        dict(bounds),
    )


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
        exact_genes = [
            row for row in genes
            if gene_id in {row.get("id", ""), row.get("attrs", {}).get("ID", ""), row.get("attrs", {}).get("gene_id", "")}
        ]
        if len(exact_genes) > 1:
            raise SystemExit(f"gene_id maps to multiple exact gene loci: {gene_id}")
        if len(exact_genes) == 1:
            gene = exact_genes[0]
            return gene, gene_ids | feature_tokens(gene)
        if len(genes) > 1:
            raise SystemExit(f"gene_id ambiguously matches multiple gene loci: {gene_id}")
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
    bounds = (seqs.get("__bounds__") or {}).get(contig)
    if bounds:
        lower, upper = bounds
        if start < lower or end > upper:
            raise SystemExit(f"requested interval {contig}:{start}-{end} outside loaded FASTA window {lower}-{upper}")
        local_start = start - lower + 1
        local_end = end - lower + 1
    else:
        local_start = start
        local_end = end
    sub = seq[local_start - 1 : local_end]
    if strand == "-":
        table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        sub = sub.translate(table)[::-1]
    return sub.upper()


def translate_cds(sequence):
    sequence = (sequence or "").upper().replace("U", "T")
    try:
        from Bio.Seq import Seq
    except ImportError as exc:
        raise SystemExit("Biopython is required to translate reconstructed CDS sequences") from exc
    usable = sequence[: len(sequence) - (len(sequence) % 3)]
    return str(Seq(usable).translate()) if usable else ""


def transcript_features(features, gene, gene_ids):
    transcript_types = {"mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna"}
    transcripts = []
    for row in features:
        if row["type"].lower() not in transcript_types:
            continue
        if not overlaps_gene(row, gene):
            continue
        if set(row.get("parents", [])) & gene_ids or row.get("id") in gene_ids:
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


def _merge_exonic_intervals(features, transcript):
    """Reconstruct exon intervals when an annotation omits explicit exon rows."""
    intervals = sorted(features, key=lambda row: (row["start"], row["end"]))
    merged = []
    for feature in intervals:
        if not merged or feature["start"] > merged[-1]["end"] + 1:
            merged.append(
                {
                    **feature,
                    "type": "exon",
                    "id": f"{transcript.get('id', 'tx')}.reconstructed_exon_{len(merged) + 1}",
                    "attrs": {"ID": f"{transcript.get('id', 'tx')}.reconstructed_exon_{len(merged) + 1}"},
                    "components": [feature],
                    "reconstructed": True,
                }
            )
        else:
            merged[-1]["end"] = max(merged[-1]["end"], feature["end"])
            merged[-1]["components"].append(feature)
    return merged


def _annotate_exon(exon, child_rows):
    item = dict(exon)
    overlapping = [
        row for row in child_rows
        if row["start"] <= item["end"] and row["end"] >= item["start"]
    ]
    cds = sorted(
        [row for row in overlapping if row["type"].lower() == "cds"],
        key=lambda row: transcript_sort_key(row, item.get("strand", "+")),
    )
    utr = [row for row in overlapping if row["type"].lower() in {"utr", "five_prime_utr", "three_prime_utr"}]
    item["cds_intervals"] = [(row["start"], row["end"]) for row in cds]
    item["utr_intervals"] = [(row["start"], row["end"]) for row in utr]
    item["cds_length"] = sum(row["end"] - row["start"] + 1 for row in cds)
    item["cds_start"] = min((row["start"] for row in cds), default=None)
    item["cds_end"] = max((row["end"] for row in cds), default=None)
    item["cds_phase"] = cds[0].get("phase", ".") if cds else "."
    item["coding_status"] = "coding" if cds else "noncoding"
    item["utr_status"] = "contains_utr" if utr else "no_annotated_utr"
    item["phase"] = item["cds_phase"]
    return item


def _format_intervals(intervals):
    if not intervals:
        return "NA"
    return ";".join(f"{int(start)}-{int(end)}" for start, end in sorted(intervals))


def _feature_parents(feature):
    parents = set(feature.get("parents", []))
    parent = feature.get("parent")
    if parent:
        parents.add(parent)
    return sorted(parent for parent in parents if parent)


def _format_attrs(attrs):
    return ";".join(f"{key}={attrs[key]}" for key in sorted(attrs)) if attrs else "NA"


def _interval_union_length(intervals):
    cleaned = sorted((int(start), int(end)) for start, end in intervals if int(end) >= int(start))
    if not cleaned:
        return 0
    merged = []
    for start, end in cleaned:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start + 1 for start, end in merged)


def child_features_for_transcript(features, gene, gene_ids, transcript):
    tx_id = transcript.get("id", "")
    selected = []
    allowed = {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr"}
    synthetic = transcript["type"] == "synthetic_transcript"
    for row in features:
        if row["type"].lower() not in allowed or not overlaps_gene(row, gene):
            continue
        parents = set(row.get("parents", []))
        if (tx_id and tx_id in parents) or (synthetic and (parents & gene_ids or row.get("id") in gene_ids)):
            selected.append(row)
    if not selected:
        return []
    explicit_exons = [row for row in selected if row["type"].lower() == "exon"]
    components = [row for row in selected if row["type"].lower() != "exon"]
    exons = explicit_exons or _merge_exonic_intervals(components, transcript)
    return sorted(
        [_annotate_exon(exon, components) for exon in exons],
        key=lambda row: transcript_sort_key(row, gene["strand"]),
    )


def transcript_sort_key(row, strand):
    key = (int(row["start"]), int(row["end"]))
    if strand == "-":
        key = (-int(row["end"]), -int(row["start"]))
    return key


def transcript_cds_length(tx_features):
    intervals = []
    for row in tx_features:
        intervals.extend(row.get("cds_intervals") or [])
        if row.get("type", "").lower() == "cds":
            intervals.append((row["start"], row["end"]))
    return _interval_union_length(intervals)


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
            upstream, downstream = left["right"], right["left"]
            if gene["strand"] == "-":
                upstream, downstream = downstream, upstream
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
                    "left_feature_id": upstream.get("id", "NA") or "NA",
                    "right_feature_id": downstream.get("id", "NA") or "NA",
                    "left_phase": upstream.get("phase", "."),
                    "right_phase": downstream.get("phase", "."),
                    "left_cds_length": upstream.get("cds_length", 0),
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
    transcript_policy="all",
    canonical_rule="longest_cds",
    source_label="unknown_source",
    copy_role="candidate",
    flank=1000,
    max_extension=10000,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    features, gene, gene_ids, annotation_bounds = read_annotation_for_gene(annotation_path, gene_id)
    contig = gene["seqid"]
    contig_length = fasta_record_length(str(Path(genome_fasta)), contig)
    annotation_start = int(annotation_bounds["annotation_start"])
    annotation_end = int(annotation_bounds["annotation_end"])
    linked_start = int(annotation_bounds["linked_start"])
    linked_end = int(annotation_bounds["linked_end"])
    flank = int(flank or 0)
    max_extension = int(max_extension or 0)
    if flank < 0:
        raise SystemExit("flank must be non-negative")
    if max_extension < 0:
        raise SystemExit("max_extension must be non-negative")
    if flank > max_extension:
        raise SystemExit("flank must be less than or equal to max_extension")
    requested_start = min(linked_start, annotation_start - flank)
    requested_end = max(linked_end, annotation_end + flank)
    search_start = max(1, requested_start)
    search_end = min(contig_length, requested_end)
    if requested_start < 1 and requested_end > contig_length:
        range_status = "both_truncated"
    elif requested_start < 1:
        range_status = "left_truncated"
    elif requested_end > contig_length:
        range_status = "right_truncated"
    else:
        range_status = "complete"
    seqs = {
        contig: read_fasta_interval(genome_fasta, contig, search_start, search_end),
        "__bounds__": {contig: (search_start, search_end)},
    }
    annotation_gene = {**gene, "start": linked_start, "end": linked_end}
    gene_ids = set(gene_ids)
    transcripts = transcript_features(features, annotation_gene, gene_ids)
    features_by_tx = {tx["id"]: child_features_for_transcript(features, annotation_gene, gene_ids, tx) for tx in transcripts}
    selected = select_transcripts(transcripts, features_by_tx, transcript_policy, canonical_rule)
    if not selected or any(not features_by_tx.get(tx.get("id", "")) for tx in selected):
        raise SystemExit(
            f"gene {gene_id} has no resolvable exon/CDS/UTR structure; the gene span is not treated as an exon"
        )

    rows = list(load_existing_segments(output_dir) if append else [])
    tx_path_rows = list(read_tsv(output_dir / "transcript_paths.tsv", optional=True) if append else [])
    intron_rows = list(read_tsv(output_dir / "intron_sites.tsv", optional=True) if append else [])
    locus_rows = list(read_tsv(output_dir / "gene_loci.tsv", optional=True) if append else [])
    raw_feature_rows = list(read_tsv(output_dir / "raw_gene_features.tsv", optional=True) if append else [])
    raw_seen = {
        (
            row.get("species"),
            row.get("gene_copy_id"),
            row.get("seqid"),
            row.get("type"),
            row.get("start"),
            row.get("end"),
            row.get("strand"),
            row.get("id"),
            row.get("parent"),
        )
        for row in raw_feature_rows
    }
    for feat in sorted(features, key=lambda row: (row.get("seqid", ""), int(row.get("start", 0)), int(row.get("end", 0)), row.get("type", ""), row.get("id", ""))):
        if feat.get("seqid") != contig or int(feat.get("end", 0)) < search_start or int(feat.get("start", 0)) > search_end:
            continue
        ownership = (
            "target_gene_descendant"
            if (set(feat.get("parents", [])) & gene_ids or feat.get("id") in gene_ids)
            else "overlapping_context"
        )
        raw_row = {
            "family_id": family_id,
            "species": species,
            "gene_copy_id": gene_copy_id,
            "gene_id": gene_id,
            "seqid": feat.get("seqid", "NA"),
            "source": feat.get("source", "NA"),
            "type": feat.get("type", "NA"),
            "start": feat.get("start", "NA"),
            "end": feat.get("end", "NA"),
            "strand": feat.get("strand", "NA"),
            "phase": feat.get("phase", "."),
            "id": feat.get("id", "NA") or "NA",
            "name": feat.get("name", "NA") or "NA",
            "parent": feat.get("parent", "NA") or "NA",
            "parents": ";".join(_feature_parents(feat)) or "NA",
            "ownership": ownership,
            "attrs": _format_attrs(feat.get("attrs", {})),
        }
        key = (
            raw_row["species"],
            raw_row["gene_copy_id"],
            raw_row["seqid"],
            raw_row["type"],
            str(raw_row["start"]),
            str(raw_row["end"]),
            raw_row["strand"],
            raw_row["id"],
            raw_row["parent"],
        )
        if key not in raw_seen:
            raw_seen.add(key)
            raw_feature_rows.append(raw_row)

    unique = {}
    tx_paths = []
    intron_records = []
    for tx in selected:
        tx_id = tx.get("id", "") or f"{gene_id}.tx"
        child_features = features_by_tx[tx_id]
        introns = introns_from_path(child_features, annotation_gene, tx_id, seqs)
        path_features = sorted(child_features + introns, key=lambda row: transcript_sort_key(row, annotation_gene["strand"]))
        for rank, feat in enumerate(path_features, start=1):
            feat = dict(feat)
            feat["transcript_order"] = rank
            key = _feature_key(feat)
            unique.setdefault(
                key,
                {
                    "feature": feat,
                    "transcripts": set(),
                    "source_ids": set(),
                    "source_types": set(),
                    "source_parents": set(),
                    "cds_intervals": set(),
                    "utr_intervals": set(),
                },
            )
            unique[key]["transcripts"].add(tx_id)
            unique[key]["source_ids"].add(feat.get("id", "NA") or "NA")
            unique[key]["source_types"].add(feat.get("type", "NA") or "NA")
            unique[key]["source_parents"].update(_feature_parents(feat))
            unique[key]["cds_intervals"].update(tuple(interval) for interval in feat.get("cds_intervals", []))
            unique[key]["utr_intervals"].update(tuple(interval) for interval in feat.get("utr_intervals", []))
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
            cds_intervals = sorted(entry["cds_intervals"] or feat.get("cds_intervals", []))
            cds_length = _interval_union_length(cds_intervals)
            cds_start = min((start for start, _end in cds_intervals), default=None)
            cds_end = max((end for _start, end in cds_intervals), default=None)
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
                    "transcript_order": feat.get("transcript_order", "NA"),
                    "coding_status": feat.get("coding_status", "not_applicable"),
                    "cds_start": cds_start or "NA",
                    "cds_end": cds_end or "NA",
                    "cds_length": cds_length,
                    "cds_phase": feat.get("cds_phase", "."),
                    "cds_intervals": _format_intervals(cds_intervals),
                    "utr_intervals": _format_intervals(entry["utr_intervals"] or feat.get("utr_intervals", [])),
                    "utr_status": feat.get("utr_status", "not_applicable"),
                    "source_feature_type": ";".join(sorted(entry["source_types"])) or feat.get("type", "NA"),
                    "source_parent": feat.get("parent", "NA") or "NA",
                    "source_parents": ";".join(sorted(entry["source_parents"])) or "NA",
                    "feature_ownership": (
                        "shared_by_observed_transcripts"
                        if len(entry["transcripts"]) > 1
                        else "transcript_specific"
                    ),
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
                "transcript_order": rank,
                "coding_status": feat.get("coding_status", "not_applicable"),
                "cds_length": _interval_union_length(feat.get("cds_intervals", [])),
                "cds_phase": feat.get("cds_phase", "."),
                "cds_intervals": _format_intervals(feat.get("cds_intervals", [])),
                "utr_intervals": _format_intervals(feat.get("utr_intervals", [])),
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
                "left_cds_length": intron.get("left_cds_length", 0),
                "expected_right_phase": (
                    (3 - ((int(intron.get("left_cds_length", 0)) - int(intron.get("left_phase", 0))) % 3)) % 3
                    if str(intron.get("left_phase", ".")) in {"0", "1", "2"} and int(intron.get("left_cds_length", 0)) > 0
                    else "NA"
                ),
                "phase_compatibility": phase_compatibility(
                    intron.get("left_phase", "."), intron.get("right_phase", "."), intron.get("left_cds_length")
                ),
                "splice_donor": intron.get("splice_donor", "NA"),
                "splice_acceptor": intron.get("splice_acceptor", "NA"),
                "splice_motif_score": f"{to_float(intron.get('splice_motif_score')):.6g}",
            }
        )

    write_tsv(output_dir / "segment_occurrences.tsv", rows, SEGMENT_OUTPUT_FIELDS)
    write_tsv(output_dir / "transcript_paths.tsv", tx_path_rows, TRANSCRIPT_PATH_FIELDS)
    write_tsv(output_dir / "intron_sites.tsv", intron_rows, INTRON_SITE_FIELDS)
    write_tsv(output_dir / "raw_gene_features.tsv", raw_feature_rows, RAW_FEATURE_FIELDS)
    locus_rows.append(
        {
            "species": species,
            "gene_copy_id": gene_copy_id,
            "contig": contig,
            "strand": gene["strand"],
            "annotation_start": annotation_start,
            "annotation_end": annotation_end,
            "linked_start": linked_start,
            "linked_end": linked_end,
            "search_start": search_start,
            "search_end": search_end,
            "contig_length": contig_length,
            "genome_fasta": str(genome_fasta),
            "max_extension": max_extension,
            "range_status": range_status,
        }
    )
    write_tsv(output_dir / "gene_loci.tsv", locus_rows, GENE_LOCUS_FIELDS)
    locus_fasta = output_dir / "gene_loci.fasta"
    locus_mode = "a" if append and locus_fasta.exists() else "w"
    with locus_fasta.open(locus_mode) as handle:
        handle.write(
            f">{species}|{gene_copy_id}|{contig}:{search_start}-{search_end}:{gene['strand']}\n"
            f"{sequence_slice(seqs, contig, search_start, search_end, gene['strand'])}\n"
        )
    transcript_fasta = output_dir / "transcript_sequences.fasta"
    protein_fasta = output_dir / "protein_sequences.fasta"
    sequence_mode = "a" if append and transcript_fasta.exists() else "w"
    protein_mode = "a" if append and protein_fasta.exists() else "w"
    with transcript_fasta.open(sequence_mode) as tx_handle, protein_fasta.open(protein_mode) as protein_handle:
        for tx in selected:
            tx_id = tx.get("id", "") or f"{gene_id}.tx"
            exons = sorted(features_by_tx[tx_id], key=lambda row: transcript_sort_key(row, annotation_gene["strand"]))
            transcript_sequence = "".join(
                sequence_slice(seqs, exon["seqid"], exon["start"], exon["end"], exon["strand"])
                for exon in exons
            )
            cds_parts = []
            for exon in exons:
                intervals = exon.get("cds_intervals", [])
                intervals = sorted(intervals, reverse=annotation_gene["strand"] == "-")
                cds_parts.extend(
                    sequence_slice(seqs, exon["seqid"], start, end, exon["strand"])
                    for start, end in intervals
                )
            cds_sequence = "".join(cds_parts)
            phase = next((exon.get("cds_phase") for exon in exons if exon.get("cds_phase") in {"0", "1", "2"}), "0")
            protein = translate_cds(cds_sequence[int(phase) :]) if cds_sequence else ""
            record_id = f"{species}|{gene_copy_id}|{tx_id}"
            tx_handle.write(f">{record_id}\n{transcript_sequence}\n")
            if protein:
                protein_handle.write(f">{record_id}\n{protein}\n")
    return rows


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


def _alignment_blocks(aln):
    blocks = getattr(aln, "aligned_blocks", None)
    if blocks:
        return [(int(qs), int(qe), int(ts), int(te)) for qs, qe, ts, te in blocks]
    return []


def _format_alignment_blocks(blocks):
    if not blocks:
        return "NA"
    return ";".join(f"{qs}-{qe}:{ts}-{te}" for qs, qe, ts, te in blocks)


def _mapped_genomic_interval(row, rel_start, rel_end):
    if rel_start in {"NA", None, ""} or rel_end in {"NA", None, ""}:
        return "NA", "NA", "NA"
    rel_start = int(rel_start)
    rel_end = int(rel_end)
    if rel_end < rel_start:
        return "NA", "NA", "NA"
    start = int(row["start"])
    end = int(row["end"])
    if row.get("strand") == "-":
        genomic_start = end - rel_end + 1
        genomic_end = end - rel_start + 1
    else:
        genomic_start = start + rel_start - 1
        genomic_end = start + rel_end - 1
    genomic_start, genomic_end = min(genomic_start, genomic_end), max(genomic_start, genomic_end)
    return genomic_start, genomic_end, genomic_end - genomic_start + 1


def match_evidence(left, right, seqs, context, aligner="auto", threads=1, context_aligner="minimap2"):
    left_seq = seqs.get(left["occurrence_id"], "")
    right_seq = seqs.get(right["occurrence_id"], "")
    exon_pair = left.get("role") in EXON_LIKE_ROLES and right.get("role") in EXON_LIKE_ROLES
    requested_backend = aligner if exon_pair else context_aligner
    # Keep left/query and right/target order for blocks, coverage, and genomic projection.
    if exon_pair and aligner in {"mafft", "auto"}:
        aln = overlap_alignment_stats(left_seq, right_seq, backend="mafft", threads=threads)
    else:
        if not exon_pair and context_aligner not in {"internal", "minimap2", "lastz"}:
            raise ValueError("context_aligner must be internal, minimap2, or lastz")
        aln = local_alignment_stats(left_seq, right_seq, backend=requested_backend, threads=threads)
    identity = aln.identity
    aligned_pairs = int(getattr(aln, "aligned_pairs", 0) or 0)
    coverage = aligned_pairs / max(1, min(len(left_seq), len(right_seq)))
    blocks = _alignment_blocks(aln)
    qstart = min((block[0] for block in blocks), default="NA")
    qend = max((block[1] for block in blocks), default="NA")
    tstart = min((block[2] for block in blocks), default="NA")
    tend = max((block[3] for block in blocks), default="NA")
    q_gen_start, q_gen_end, q_gen_len = _mapped_genomic_interval(left, qstart, qend)
    t_gen_start, t_gen_end, t_gen_len = _mapped_genomic_interval(right, tstart, tend)
    left_ctx = context.get(left["occurrence_id"], {})
    right_ctx = context.get(right["occurrence_id"], {})
    order = 1.0 - abs(to_float(left_ctx.get("scaled_index"), 0.5) - to_float(right_ctx.get("scaled_index"), 0.5))
    left_context = context_score(left_ctx, right_ctx, "left")
    right_context = context_score(left_ctx, right_ctx, "right")
    boundary = role_boundary_score(left, right)
    phase = phase_score(left, right)
    alignment_strand = getattr(aln, "strand", "+")
    strand = 1.0 if alignment_strand == "+" else 0.0
    splice = 1.0 - abs(to_float(left.get("splice_motif_score"), 0.5) - to_float(right.get("splice_motif_score"), 0.5))
    sequence_score = 0.70 * identity + 0.30 * coverage
    structural_context_score = 0.5 * left_context + 0.5 * right_context
    total = (
        0.34 * identity
        + 0.14 * coverage
        + 0.10 * left_context
        + 0.10 * right_context
        + 0.10 * boundary
        + 0.08 * phase
        + 0.06 * order
        + 0.04 * strand
        + 0.04 * splice
    )
    return {
        "alignment_score": identity,
        "coverage_score": coverage,
        "sequence_score": sequence_score,
        "structural_context_score": structural_context_score,
        "query_coverage": aln.query_coverage,
        "target_coverage": aln.target_coverage,
        "aligned_pairs": aligned_pairs,
        "query_alignment_start": qstart,
        "query_alignment_end": qend,
        "target_alignment_start": tstart,
        "target_alignment_end": tend,
        "query_mapped_contig": left.get("contig", "NA"),
        "query_mapped_start": q_gen_start,
        "query_mapped_end": q_gen_end,
        "query_mapped_strand": left.get("strand", "NA"),
        "query_mapped_length": q_gen_len,
        "subject_mapped_contig": right.get("contig", "NA"),
        "subject_mapped_start": t_gen_start,
        "subject_mapped_end": t_gen_end,
        "subject_mapped_strand": right.get("strand", "NA"),
        "subject_mapped_length": t_gen_len,
        "alignment_strand": alignment_strand,
        "projected_reference_occurrence_id": right["occurrence_id"],
        "projected_reference_start": tstart,
        "projected_reference_end": tend,
        "projected_reference_blocks": _format_alignment_blocks(blocks),
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
        "alignment_mode": aln.alignment_mode,
        "alignment_meaning": aln.alignment_meaning,
        "alignment_requested_backend": requested_backend,
        "total_score": total,
    }


EXON_LIKE_ROLES = {"CDS", "exon", "UTR", "noncoding_exon"}
NONCODING_ROLES = {"intron", "intergenic"}
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


def graph_components(nodes, edges, occurrence_by_id=None, species_distances=None, same_copy_compatible=None, cross_copy_compatible=None):
    if occurrence_by_id is None:
        import networkx as nx
        graph = nx.Graph()
        graph.add_nodes_from(nodes)
        graph.add_weighted_edges_from(edges)
        return [sorted(component) for component in nx.connected_components(graph)]
    parent = {node: node for node in nodes}
    members = {node: {node} for node in nodes}
    direct = {frozenset((left, right)) for left, right, _score in edges}
    same_copy_compatible = same_copy_compatible or set()
    cross_copy_compatible = cross_copy_compatible or set()

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def compatible(left_root, right_root):
        for left_id in members[left_root]:
            for right_id in members[right_root]:
                left = occurrence_by_id.get(left_id, {})
                right = occurrence_by_id.get(right_id, {})
                if occurrence_copy_key(left) == occurrence_copy_key(right):
                    if frozenset((left_id, right_id)) not in same_copy_compatible:
                        return False
                elif (
                    frozenset((left_id, right_id)) not in direct
                    and frozenset((left_id, right_id)) not in cross_copy_compatible
                ):
                    return False
        return True

    species_distances = species_distances or {}
    def progressive_key(edge):
        left, right, score = edge
        left_species = occurrence_by_id.get(left, {}).get("species")
        right_species = occurrence_by_id.get(right, {}).get("species")
        return (species_distances.get((left_species, right_species), math.inf), -score, left, right)

    for left, right, _score in sorted(edges, key=progressive_key):
        left_root, right_root = find(left), find(right)
        if left_root == right_root or not compatible(left_root, right_root):
            continue
        if len(members[left_root]) < len(members[right_root]):
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        members[left_root].update(members.pop(right_root))
    return [sorted(component) for component in members.values()]


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


def _species_tree_distances(path):
    """Prioritize correspondence by supplied lengths, or topology if any are missing.

    Unit edges apply to this ordering only; CTMC branch lengths are unchanged.
    """
    if not Path(path).exists():
        return {}
    from .tree import SpeciesTree
    tree = SpeciesTree(read_tsv(path))
    use_topology = any(tree.length[node] is None for node in tree.parent if node != tree.root)
    root_dist = {tree.root: 0.0}
    ancestors = {}
    for node in tree.preorder():
        ancestors[node] = [node] + (ancestors.get(tree.parent.get(node), []))
        for child in tree.children.get(node, []):
            root_dist[child] = root_dist[node] + (1.0 if use_topology else tree.branch_length(child))
    result = {}
    for left_label, left_node in tree.leaf_by_label.items():
        for right_label, right_node in tree.leaf_by_label.items():
            right_ancestors = set(ancestors[right_node])
            lca = next(node for node in ancestors[left_node] if node in right_ancestors)
            result[(left_label, right_label)] = root_dist[left_node] + root_dist[right_node] - 2 * root_dist[lca]
    return result


MATCH_FIELDS = [
    "match_id",
    "query_occurrence_id",
    "subject_occurrence_id",
    "alignment_score",
    "coverage_score",
    "sequence_score",
    "structural_context_score",
    "query_coverage",
    "target_coverage",
    "aligned_pairs",
    "query_alignment_start",
    "query_alignment_end",
    "target_alignment_start",
    "target_alignment_end",
    "query_mapped_contig",
    "query_mapped_start",
    "query_mapped_end",
    "query_mapped_strand",
    "query_mapped_length",
    "subject_mapped_contig",
    "subject_mapped_start",
    "subject_mapped_end",
    "subject_mapped_strand",
    "subject_mapped_length",
    "alignment_strand",
    "projected_reference_occurrence_id",
    "projected_reference_start",
    "projected_reference_end",
    "projected_reference_blocks",
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
    "alignment_mode",
    "alignment_meaning",
    "alignment_requested_backend",
    "dna_match_status",
    "correspondence_basis",
    "correspondence_score",
    "protein_status",
    "protein_unavailable_reason",
    "protein_aa_identity",
    "protein_query_cds_coverage",
    "protein_target_cds_coverage",
    "protein_metrics_scope",
    "protein_best_query_transcript",
    "protein_best_target_transcript",
    "protein_supporting_transcripts",
    "protein_projected_blocks",
]


def _projection_interval(evidence, side):
    if evidence.get("correspondence_basis") == "annotated_CDS_protein":
        blocks = evidence.get("protein_projected_blocks", ())
        offset = 0 if side == "query" else 2
        if not blocks:
            return None
        return min(block[offset] for block in blocks), max(block[offset + 1] for block in blocks)
    if side == "query":
        start = evidence.get("query_alignment_start", "NA")
        end = evidence.get("query_alignment_end", "NA")
    else:
        start = evidence.get("target_alignment_start", "NA")
        end = evidence.get("target_alignment_end", "NA")
    if start in {"NA", None, ""} or end in {"NA", None, ""}:
        return None
    start, end = int(start), int(end)
    if end < start:
        return None
    return start, end


def _projection_record(evidence, side):
    interval = _projection_interval(evidence, side)
    if interval is None:
        return None
    return {
        "interval": interval,
        "strand": "+" if evidence.get("correspondence_basis") == "annotated_CDS_protein" else evidence.get("alignment_strand", "NA"),
    }


def _ordered_projection_compatible(left, right, left_ref_interval, right_ref_interval):
    if not left_ref_interval or not right_ref_interval:
        return False
    left_ref_start, left_ref_end = left_ref_interval
    right_ref_start, right_ref_end = right_ref_interval
    if not _disjoint_reference_intervals(left_ref_interval, right_ref_interval):
        return False
    if left.get("contig") != right.get("contig") or left.get("strand") != right.get("strand"):
        return False
    left_order = int(left["start"])
    right_order = int(right["start"])
    if left.get("strand") == "-":
        left_order, right_order = -left_order, -right_order
    copy_order = -1 if left_order < right_order else 1
    ref_order = -1 if left_ref_start < right_ref_start else 1
    return copy_order == ref_order


def _disjoint_reference_intervals(left_ref_interval, right_ref_interval):
    if not left_ref_interval or not right_ref_interval:
        return False
    left_ref_start, left_ref_end = left_ref_interval
    right_ref_start, right_ref_end = right_ref_interval
    return left_ref_end < right_ref_start or right_ref_end < left_ref_start


def _projection_compatibility(projection_by_occ_ref, occurrence_by_id):
    same_copy_compatible = set()
    cross_copy_compatible = set()
    by_reference = defaultdict(list)
    for (occ_id, ref_id), projection in projection_by_occ_ref.items():
        if projection:
            by_reference[ref_id].append((occ_id, projection))
    for projected in by_reference.values():
        for left_index, (left_id, left_projection) in enumerate(projected):
            left = occurrence_by_id.get(left_id, {})
            for right_id, right_projection in projected[left_index + 1 :]:
                right = occurrence_by_id.get(right_id, {})
                pair = frozenset((left_id, right_id))
                if occurrence_copy_key(left) == occurrence_copy_key(right):
                    if (
                        left_projection.get("strand") == "+"
                        and right_projection.get("strand") == "+"
                        and _ordered_projection_compatible(
                            left,
                            right,
                            left_projection["interval"],
                            right_projection["interval"],
                        )
                    ):
                        same_copy_compatible.add(pair)
                elif _disjoint_reference_intervals(left_projection["interval"], right_projection["interval"]):
                    cross_copy_compatible.add(pair)
    return same_copy_compatible, cross_copy_compatible


def _match_row(left, right, evidence, score, threshold, distance_class, status, match_index):
    return {
        "match_id": f"match_{match_index:05d}",
        "query_occurrence_id": left["occurrence_id"],
        "subject_occurrence_id": right["occurrence_id"],
        "alignment_score": f"{evidence['alignment_score']:.6g}",
        "coverage_score": f"{evidence['coverage_score']:.6g}",
        "sequence_score": f"{evidence.get('sequence_score', 0.0):.6g}",
        "structural_context_score": f"{evidence.get('structural_context_score', 0.5):.6g}",
        "query_coverage": f"{evidence.get('query_coverage', 0.0):.6g}",
        "target_coverage": f"{evidence.get('target_coverage', 0.0):.6g}",
        "aligned_pairs": evidence.get("aligned_pairs", 0),
        "query_alignment_start": evidence.get("query_alignment_start", "NA"),
        "query_alignment_end": evidence.get("query_alignment_end", "NA"),
        "target_alignment_start": evidence.get("target_alignment_start", "NA"),
        "target_alignment_end": evidence.get("target_alignment_end", "NA"),
        "query_mapped_contig": evidence.get("query_mapped_contig", "NA"),
        "query_mapped_start": evidence.get("query_mapped_start", "NA"),
        "query_mapped_end": evidence.get("query_mapped_end", "NA"),
        "query_mapped_strand": evidence.get("query_mapped_strand", "NA"),
        "query_mapped_length": evidence.get("query_mapped_length", "NA"),
        "subject_mapped_contig": evidence.get("subject_mapped_contig", "NA"),
        "subject_mapped_start": evidence.get("subject_mapped_start", "NA"),
        "subject_mapped_end": evidence.get("subject_mapped_end", "NA"),
        "subject_mapped_strand": evidence.get("subject_mapped_strand", "NA"),
        "subject_mapped_length": evidence.get("subject_mapped_length", "NA"),
        "alignment_strand": evidence.get("alignment_strand", "NA"),
        "projected_reference_occurrence_id": evidence.get("projected_reference_occurrence_id", "NA"),
        "projected_reference_start": evidence.get("projected_reference_start", "NA"),
        "projected_reference_end": evidence.get("projected_reference_end", "NA"),
        "projected_reference_blocks": evidence.get("projected_reference_blocks", "NA"),
        "left_context_score": f"{evidence['left_context_score']:.6g}",
        "right_context_score": f"{evidence['right_context_score']:.6g}",
        "boundary_score": f"{evidence['boundary_score']:.6g}",
        "phase_score": f"{evidence['phase_score']:.6g}",
        "order_score": f"{evidence['order_score']:.6g}",
        "strand_score": f"{evidence['strand_score']:.6g}",
        "splice_score": f"{evidence['splice_score']:.6g}",
        "size_ratio": f"{evidence['size_ratio']:.6g}",
        "total_score": f"{evidence['total_score']:.6g}",
        "distance_class": distance_class,
        "threshold": f"{threshold:.6g}",
        "alignment_cigar": evidence["alignment_cigar"],
        "alignment_backend": evidence["alignment_backend"],
        "match_status": status,
        "alignment_mode": evidence["alignment_mode"],
        "alignment_meaning": evidence["alignment_meaning"],
        "alignment_requested_backend": evidence["alignment_requested_backend"],
        "dna_match_status": evidence.get("dna_match_status", status),
        "correspondence_basis": evidence.get("correspondence_basis", "DNA"),
        "correspondence_score": f"{score:.6g}",
        "protein_status": evidence.get("protein_status", "not_requested"),
        "protein_unavailable_reason": evidence.get("protein_unavailable_reason", "NA"),
        "protein_metrics_scope": evidence.get("protein_metrics_scope", "NA"),
        "protein_best_query_transcript": evidence.get("protein_best_query_transcript", "NA"),
        "protein_best_target_transcript": evidence.get("protein_best_target_transcript", "NA"),
        "protein_supporting_transcripts": evidence.get("protein_supporting_transcripts", "NA"),
        "protein_projected_blocks": _format_alignment_blocks(evidence.get("protein_projected_blocks", ())),
        **{
            field: f"{evidence[field]:.6g}" if field in evidence else "NA"
            for field in ("protein_aa_identity", "protein_query_cds_coverage", "protein_target_cds_coverage")
        },
    }


def _transcript_id_set(row):
    return {token for token in str(row.get("transcript_id", "")).split(";") if token and token != "NA"}


def _relative_overlap_interval(row, overlap_start, overlap_end):
    start = int(row["start"])
    end = int(row["end"])
    if row.get("strand") == "-":
        rel_start = end - int(overlap_end) + 1
        rel_end = end - int(overlap_start) + 1
    else:
        rel_start = int(overlap_start) - start + 1
        rel_end = int(overlap_end) - start + 1
    return min(rel_start, rel_end), max(rel_start, rel_end)


def alternative_overlap_evidence(left, right, context):
    if occurrence_copy_key(left) != occurrence_copy_key(right):
        return None
    if left.get("role") not in EXON_LIKE_ROLES or right.get("role") not in EXON_LIKE_ROLES:
        return None
    if left.get("contig") != right.get("contig") or left.get("strand") != right.get("strand"):
        return None
    left_tx = _transcript_id_set(left)
    right_tx = _transcript_id_set(right)
    if not left_tx or not right_tx or left_tx & right_tx:
        return None
    try:
        overlap_start = max(int(left["start"]), int(right["start"]))
        overlap_end = min(int(left["end"]), int(right["end"]))
    except (KeyError, TypeError, ValueError):
        return None
    if overlap_start > overlap_end:
        return None
    query_start, query_end = _relative_overlap_interval(left, overlap_start, overlap_end)
    target_start, target_end = _relative_overlap_interval(right, overlap_start, overlap_end)
    overlap_len = overlap_end - overlap_start + 1
    left_len = segment_length(left)
    right_len = segment_length(right)
    coverage = overlap_len / max(1, min(left_len, right_len))
    left_ctx = context.get(left["occurrence_id"], {})
    right_ctx = context.get(right["occurrence_id"], {})
    order = 1.0 - abs(to_float(left_ctx.get("scaled_index"), 0.5) - to_float(right_ctx.get("scaled_index"), 0.5))
    left_context = context_score(left_ctx, right_ctx, "left")
    right_context = context_score(left_ctx, right_ctx, "right")
    boundary = role_boundary_score(left, right)
    phase = phase_score(left, right)
    splice = 1.0 - abs(to_float(left.get("splice_motif_score"), 0.5) - to_float(right.get("splice_motif_score"), 0.5))
    total = (
        0.34 * 1.0
        + 0.14 * coverage
        + 0.10 * left_context
        + 0.10 * right_context
        + 0.10 * boundary
        + 0.08 * phase
        + 0.06 * order
        + 0.04 * 1.0
        + 0.04 * splice
    )
    return {
        "alignment_score": 1.0,
        "coverage_score": coverage,
        "sequence_score": 0.70 * 1.0 + 0.30 * coverage,
        "structural_context_score": 0.5 * left_context + 0.5 * right_context,
        "query_coverage": overlap_len / max(1, left_len),
        "target_coverage": overlap_len / max(1, right_len),
        "aligned_pairs": overlap_len,
        "query_alignment_start": query_start,
        "query_alignment_end": query_end,
        "target_alignment_start": target_start,
        "target_alignment_end": target_end,
        "query_mapped_contig": left.get("contig", "NA"),
        "query_mapped_start": overlap_start,
        "query_mapped_end": overlap_end,
        "query_mapped_strand": left.get("strand", "NA"),
        "query_mapped_length": overlap_len,
        "subject_mapped_contig": right.get("contig", "NA"),
        "subject_mapped_start": overlap_start,
        "subject_mapped_end": overlap_end,
        "subject_mapped_strand": right.get("strand", "NA"),
        "subject_mapped_length": overlap_len,
        "alignment_strand": "+",
        "projected_reference_occurrence_id": right["occurrence_id"],
        "projected_reference_start": target_start,
        "projected_reference_end": target_end,
        "projected_reference_blocks": f"{query_start}-{query_end}:{target_start}-{target_end}",
        "left_context_score": left_context,
        "right_context_score": right_context,
        "boundary_score": boundary,
        "phase_score": phase,
        "order_score": order,
        "strand_score": 1.0,
        "splice_score": splice,
        "size_ratio": min(left_len, right_len) / max(left_len, right_len),
        "alignment_cigar": f"{overlap_len}=",
        "alignment_backend": "genomic_overlap",
        "alignment_mode": "coordinate_overlap",
        "alignment_meaning": "same-copy alternative isoform shared genomic interval",
        "alignment_requested_backend": "genomic_overlap",
        "total_score": total,
    }


def cluster_segments(occurrences, seqs, identity_threshold=0.7, distance_table=None, aligner="mafft", threads=1, min_size_ratio=0.25, species_distances=None, match_writer=None, context_aligner="minimap2", transcript_paths=None, raw_features=None):
    context = copy_order_context(occurrences, transcript_paths)
    distance_lookup = load_distance_table(distance_table)
    occurrence_by_id = {row["occurrence_id"]: row for row in occurrences}
    matches = []
    accepted_edges = []
    projection_by_occ_ref = {}
    score_by_occ = defaultdict(list)
    source_support = defaultdict(lambda: defaultdict(float))
    genomic_overlap_compatible = set()
    protein_index = None
    if aligner in {"mafft", "auto"} and transcript_paths:
        parsed_features = [
            {**row, "attrs": parse_attributes(row["attrs"]) if isinstance(row.get("attrs"), str) else row.get("attrs", {})}
            for row in raw_features or []
        ]
        protein_index = CodingProjectionIndex(occurrences, seqs, transcript_paths, parsed_features)

    def emit_match(row):
        if match_writer is not None:
            match_writer(row)
            if row.get("match_status") == "mapped":
                matches.append(row)
        else:
            matches.append(row)

    def iter_pairs():
        index = 0
        by_family = defaultdict(list)
        for occurrence in occurrences:
            by_family[occurrence["family_id"]].append(occurrence)
        for family_rows in by_family.values():
            for left_index, left in enumerate(family_rows):
                for right in family_rows[left_index + 1 :]:
                    if occurrence_copy_key(left) == occurrence_copy_key(right):
                        continue
                    if left.get("role") == "intron" and right.get("role") == "intron":
                        continue
                    yield index, left, right
                    index += 1

    match_count = 0
    by_copy = defaultdict(list)
    for occurrence in occurrences:
        by_copy[occurrence_copy_key(occurrence)].append(occurrence)
    for copy_rows in by_copy.values():
        for left_index, left in enumerate(copy_rows):
            for right in copy_rows[left_index + 1 :]:
                evidence = alternative_overlap_evidence(left, right, context)
                if evidence is None:
                    continue
                score = evidence["total_score"]
                accepted_edges.append((left["occurrence_id"], right["occurrence_id"], score))
                pair = frozenset((left["occurrence_id"], right["occurrence_id"]))
                genomic_overlap_compatible.add(pair)
                projection_by_occ_ref[(left["occurrence_id"], right["occurrence_id"])] = _projection_record(evidence, "target")
                projection_by_occ_ref[(right["occurrence_id"], left["occurrence_id"])] = _projection_record(evidence, "query")
                score_by_occ[left["occurrence_id"]].append(score)
                score_by_occ[right["occurrence_id"]].append(score)
                match_count += 1
                emit_match(_match_row(left, right, evidence, score, 1.0, "same_copy_alternative_overlap", "mapped", match_count))

    def score_pair(item):
        idx, left, right = item
        should_align, prefilter_status = should_align_pair(left, right, seqs, min_size_ratio)
        if should_align:
            evidence = match_evidence(left, right, seqs, context, aligner=aligner, threads=1, context_aligner=context_aligner)
        else:
            evidence = cheap_match_evidence(left, right, context, alignment_backend=prefilter_status)
        return idx, left, right, evidence, prefilter_status

    worker_count = max(1, int(threads or 1))
    if worker_count > 1:
        pool = ThreadPoolExecutor(max_workers=worker_count)
        pair_iterator = iter(iter_pairs())

        def bounded_scores():
            batch_size = max(8, worker_count * 4)
            while True:
                batch = list(islice(pair_iterator, batch_size))
                if not batch:
                    break
                yield from pool.map(score_pair, batch)

        scored_pairs = bounded_scores()
    else:
        pool = None
        scored_pairs = map(score_pair, iter_pairs())

    try:
        for _idx, left, right, evidence, prefilter_status in scored_pairs:
            threshold, distance_class = pair_threshold(left, right, identity_threshold, distance_lookup)
            score = evidence["total_score"]
            compatible = roles_compatible(left, right)
            sequence_ok = sequence_supported_mapping(evidence, threshold)
            mapped = compatible and sequence_ok and score >= threshold
            evidence["dna_match_status"] = "mapped" if mapped else prefilter_status if prefilter_status != "aligned_candidate" else "low_similarity"
            exon_pair = left.get("role") in EXON_LIKE_ROLES and right.get("role") in EXON_LIKE_ROLES
            if not mapped and exon_pair:
                if protein_index is not None:
                    evidence.update(protein_index.evidence(left["occurrence_id"], right["occurrence_id"]))
                    if evidence.get("protein_status") == "supported":
                        # Keep the original weights and raw DNA statistics. These
                        # sequence components use annotated CDS evidence only.
                        coding_score = (
                            score - 0.34 * evidence["alignment_score"] - 0.14 * evidence["coverage_score"]
                            + 0.34 * evidence["protein_aa_identity"]
                            + 0.14 * max(evidence["protein_query_cds_coverage"], evidence["protein_target_cds_coverage"])
                        )
                        if compatible and coding_score >= threshold:
                            mapped = True
                            score = coding_score
                            evidence["correspondence_basis"] = "annotated_CDS_protein"
                else:
                    evidence["protein_status"] = "unavailable" if aligner in {"mafft", "auto"} else "disabled_backend"
                    evidence["protein_unavailable_reason"] = "no_CDS_transcript_path" if aligner in {"mafft", "auto"} else "MAFFT_exon_backend_required"
            if mapped:
                accepted_edges.append((left["occurrence_id"], right["occurrence_id"], score))
                projection_by_occ_ref[(left["occurrence_id"], right["occurrence_id"])] = _projection_record(evidence, "target")
                projection_by_occ_ref[(right["occurrence_id"], left["occurrence_id"])] = _projection_record(evidence, "query")
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
            match_count += 1
            emit_match(_match_row(left, right, evidence, score, threshold, distance_class, status, match_count))
    finally:
        if pool is not None:
            pool.shutdown(wait=True)
    nodes = [row["occurrence_id"] for row in occurrences]
    same_copy_compatible, cross_copy_compatible = _projection_compatibility(projection_by_occ_ref, occurrence_by_id)
    same_copy_compatible |= genomic_overlap_compatible
    components = graph_components(
        nodes,
        accepted_edges,
        occurrence_by_id,
        species_distances,
        same_copy_compatible,
        cross_copy_compatible,
    )
    homology = []
    for idx, occ_ids in enumerate(sorted(components, key=lambda vals: vals[0]), start=1):
        component_id = f"HC_{idx:04d}"
        for occ_id in occ_ids:
            scores = score_by_occ.get(occ_id, [])
            confidence = sum(scores) / len(scores) if scores else 0.5
            homology.append(
                {
                    "homology_id": component_id,
                    "occurrence_id": occ_id,
                    "support_type": "sequence_boundary_context_graph",
                    "confidence": f"{confidence:.6g}",
                    "source_label": inferred_source_label(occurrence_by_id.get(occ_id, {}), source_support.get(occ_id, {})),
                }
            )
    return homology, matches


def derive_tables(input_dir, output_dir=None, identity_threshold=0.7, distance_table=None, aligner="mafft", threads=1, min_size_ratio=0.25, context_aligner="minimap2"):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir or input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", SEGMENT_FIELDS)
    transcript_paths = read_tsv(input_dir / "transcript_paths.tsv", optional=True)
    raw_features = read_tsv(input_dir / "raw_gene_features.tsv", optional=True)
    seqs = parse_fasta(input_dir / "segment_sequences.fasta")
    species_distances = _species_tree_distances(input_dir / "species_tree.tsv")
    match_path = output_dir / "segment_matches.tsv"
    with match_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, MATCH_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        homology, matches = cluster_segments(
            occurrences, seqs, identity_threshold, distance_table, aligner=aligner,
            threads=threads, min_size_ratio=min_size_ratio, species_distances=species_distances,
            match_writer=lambda row: writer.writerow({field: row.get(field, "NA") for field in MATCH_FIELDS}),
            context_aligner=context_aligner,
            transcript_paths=transcript_paths,
            raw_features=raw_features,
        )
    write_tsv(output_dir / "segment_homology.tsv", homology, ["homology_id", "occurrence_id", "support_type", "confidence", "source_label"])
    backend_rows = []
    for row in available_alignment_backends():
        selected_exon = row["aligner"] == aligner
        selected_context = row["aligner"] == context_aligner
        backend_rows.append(
            {
                **row,
                "selected": int(selected_exon or selected_context),
                "selected_exon": int(selected_exon),
                "selected_context": int(selected_context),
                "alignment_mode": (
                    "overlap_projection" if selected_exon and aligner in {"mafft", "auto"}
                    else "local" if selected_exon or selected_context else "NA"
                ),
                "threads": threads,
                "min_size_ratio": f"{min_size_ratio:.6g}",
                "notes": (
                    f"{row['notes']}; min_size_ratio is restricted to non-exon-like prefiltering"
                    if selected_exon or selected_context
                    else row["notes"]
                ),
            }
        )
    write_tsv(output_dir / "alignment_backend_report.tsv", backend_rows, ["aligner", "available", "selected", "threads", "min_size_ratio", "notes", "selected_exon", "selected_context", "alignment_mode"])
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
