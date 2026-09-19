"""Genome annotation extraction and first-pass intragenic table builders."""

from __future__ import annotations

import math
import csv
import json
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import lru_cache
from itertools import islice
from pathlib import Path

from .alignment import (
    AlignmentBackendError,
    AlignmentStats,
    NT_BLASTN_V1_GAP_EXTEND,
    NT_BLASTN_V1_GAP_OPEN,
    NT_BLASTN_V1_MATCH,
    NT_BLASTN_V1_MISMATCH,
    anchored_short_alignment,
    available_alignment_backends,
    local_alignment_stats,
    overlap_alignment_stats,
    phase_compatibility,
    splice_motif_score,
)
from .candidate_chain import (
    DEFAULT_CHAIN_CONFIGURATION,
    ChainCandidate,
    ChainPathMembership,
    ordered_candidate_chain,
)
from .coding_correspondence import CodingProjectionIndex
from .coordinates import (
    ClosedInterval1,
    CoordinateBlock,
    Interval0,
    format_legacy_blocks,
    genome_interval_to_local,
    local_interval_to_genome,
    parse_legacy_blocks,
)
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


@dataclass(frozen=True)
class CorrespondenceCriteria:
    regular_min_coverage: float = 0.45
    high_identity_offset: float = 0.15
    high_identity_min_coverage: float = 0.30
    short_min_identity: float = 0.70
    short_min_query_coverage: float = 0.80
    short_min_aligned_pairs: int = 12


DEFAULT_CORRESPONDENCE_CRITERIA = CorrespondenceCriteria()

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
    "path_roles",
    "coding_roles",
    "position_roles",
    "annotation_source",
    "original_attributes",
    "partial_start",
    "partial_end",
    "path_role_records",
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
    "canonical_selection_rule",
    "source_feature_id",
    "source_feature_type",
    "path_role",
    "coding_role",
    "position_role",
    "annotation_source",
    "original_attributes",
    "partial_start",
    "partial_end",
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
    item["utr_types"] = [row["type"] for row in utr]
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


def _partial_boundary(feature, boundary):
    attrs = feature.get("attrs", {}) or {}
    keys = (
        ("partial_start", "start_partial", "partial5")
        if boundary == "start"
        else ("partial_end", "end_partial", "partial3")
    )
    for key in keys:
        if key in attrs:
            return attrs[key]
    range_value = attrs.get(f"{boundary}_range")
    if range_value not in {None, "", "."}:
        return range_value
    if str(attrs.get("partial", "")).lower() in {"1", "true", "yes"}:
        return "1"
    return "0"


def _coding_role(feature):
    if feature_role(feature) == "intron":
        return "unknown"
    cds = feature.get("cds_intervals", ())
    utr_types = {str(value).lower() for value in feature.get("utr_types", ())}
    if cds and utr_types:
        return "mixed"
    if cds:
        return "CDS"
    if utr_types == {"five_prime_utr"}:
        return "five_prime_UTR"
    if utr_types == {"three_prime_utr"}:
        return "three_prime_UTR"
    if utr_types:
        return "mixed" if len(utr_types) > 1 else "unknown"
    return "noncoding_exon"


def _position_role(path_features, index):
    feature = path_features[index]
    if feature_role(feature) == "intron":
        return "internal"
    exonic_indices = [
        item_index
        for item_index, item in enumerate(path_features)
        if feature_role(item) != "intron"
    ]
    if len(exonic_indices) == 1:
        return "single"
    if index == exonic_indices[0]:
        return "first"
    if index == exonic_indices[-1]:
        return "last"
    return "internal"


def _path_role_record(transcript_id, rank, path_features, index, feature, transcript=None):
    position_role = _position_role(path_features, index)
    partial_start = _partial_boundary(feature, "start")
    partial_end = _partial_boundary(feature, "end")
    if transcript is not None and position_role in {"first", "single"} and partial_start == "0":
        partial_start = _partial_boundary(transcript, "start")
    if transcript is not None and position_role in {"last", "single"} and partial_end == "0":
        partial_end = _partial_boundary(transcript, "end")
    # A transcript may be wholly contained in its linked gene while its
    # terminal CDS/exon has no explicit partial attribute.  The uncovered
    # gene boundary is still evidence of a partial path.
    if transcript is not None:
        gene_start = transcript.get("_gene_start")
        gene_end = transcript.get("_gene_end")
        strand = transcript.get("strand", "+")
        if position_role in {"first", "single"} and gene_start is not None and gene_end is not None:
            if (strand == "+" and int(feature["start"]) > int(gene_start)) or (
                strand == "-" and int(feature["end"]) < int(gene_end)
            ):
                partial_start = "1"
        if position_role in {"last", "single"} and gene_start is not None and gene_end is not None:
            if (strand == "+" and int(feature["end"]) < int(gene_end)) or (
                strand == "-" and int(feature["start"]) > int(gene_start)
            ):
                partial_end = "1"
    return {
        "transcript_id": transcript_id,
        "path_rank": rank,
        "source_feature_id": feature.get("id", "NA") or "NA",
        "source_feature_type": feature.get("type", "NA") or "NA",
        "path_role": "intronic" if feature_role(feature) == "intron" else "exonic",
        "coding_role": _coding_role(feature),
        "position_role": position_role,
        "annotation_source": feature.get("source", (transcript or {}).get("source", "NA")) or "NA",
        "original_attributes": _format_attrs(feature.get("attrs", {})),
        "partial_start": partial_start,
        "partial_end": partial_end,
    }


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
    canonical_ids = {
        tx["id"]
        for tx in select_transcripts(
            transcripts,
            features_by_tx,
            transcript_policy="canonical",
            canonical_rule=canonical_rule,
        )
    }
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
            path_record = _path_role_record(
                tx_id, rank, path_features, rank - 1, feat,
                transcript={**tx, "_gene_start": annotation_gene["start"], "_gene_end": annotation_gene["end"]},
            )
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
                    "path_records": [],
                },
            )
            unique[key]["transcripts"].add(tx_id)
            unique[key]["source_ids"].add(feat.get("id", "NA") or "NA")
            unique[key]["source_types"].add(feat.get("type", "NA") or "NA")
            unique[key]["source_parents"].update(_feature_parents(feat))
            unique[key]["cds_intervals"].update(tuple(interval) for interval in feat.get("cds_intervals", []))
            unique[key]["utr_intervals"].update(tuple(interval) for interval in feat.get("utr_intervals", []))
            unique[key]["path_records"].append(path_record)
            tx_paths.append((tx_id, rank, key, feat, path_record))
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
            path_records = sorted(
                entry["path_records"],
                key=lambda item: (item["transcript_id"], int(item["path_rank"])),
            )
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
                    "role_set": ";".join(sorted({role, *(item["coding_role"] for item in path_records)})),
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
                    "path_roles": ";".join(sorted({item["path_role"] for item in path_records})) or "unknown",
                    "coding_roles": ";".join(sorted({item["coding_role"] for item in path_records})) or "unknown",
                    "position_roles": ";".join(sorted({item["position_role"] for item in path_records})) or "unknown",
                    "annotation_source": ";".join(sorted({item["annotation_source"] for item in path_records})) or "NA",
                    "original_attributes": json.dumps(
                        [item["original_attributes"] for item in path_records],
                        separators=(",", ":"),
                    ),
                    "partial_start": ";".join(sorted({str(item["partial_start"]) for item in path_records})),
                    "partial_end": ";".join(sorted({str(item["partial_end"]) for item in path_records})),
                    "path_role_records": json.dumps(path_records, sort_keys=True, separators=(",", ":")),
                }
            )
            fasta.write(f">{occ_id}\n{seq}\n")

    for tx_id, rank, key, feat, path_record in tx_paths:
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
                "path_status": (
                    "canonical_transcript_path"
                    if tx_id in canonical_ids
                    else "annotated_transcript_path"
                ),
                "canonical_selection_rule": canonical_rule,
                **path_record,
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


def _coordinate_block0(block):
    if isinstance(block, CoordinateBlock):
        return block
    if isinstance(block, dict):
        return CoordinateBlock(
            Interval0(int(block["query_start0"]), int(block["query_end0"])),
            Interval0(int(block["target_start0"]), int(block["target_end0"])),
        )
    query_start, query_end, target_start, target_end = block
    return CoordinateBlock(
        ClosedInterval1(int(query_start), int(query_end)).to_interval0(),
        ClosedInterval1(int(target_start), int(target_end)).to_interval0(),
    )


def _alignment_blocks0(aln):
    blocks = getattr(aln, "aligned_blocks", None)
    if blocks:
        return tuple(_coordinate_block0(block) for block in blocks)
    return tuple()


def _alignment_blocks(aln):
    """Compatibility view for callers that still consume public 1-based blocks."""

    public = []
    for block in _alignment_blocks0(aln):
        query = ClosedInterval1.from_interval0(block.query)
        target = ClosedInterval1.from_interval0(block.target)
        public.append((query.start, query.end, target.start, target.end))
    return public


def _format_alignment_blocks(blocks):
    if not blocks:
        return "NA"
    return format_legacy_blocks(tuple(_coordinate_block0(block) for block in blocks))


def _block_signature(block):
    block = _coordinate_block0(block)
    return (
        block.query.start0, block.query.end0,
        block.target.start0, block.target.end0,
    )


def _valid_local_boundary_range(row, sequence):
    try:
        start = int(row["start"])
        end = int(row["end"])
    except (KeyError, TypeError, ValueError):
        return False
    return (
        row.get("contig") not in {None, "", "NA"}
        and row.get("strand") in {"+", "-"}
        and start <= end
        and len(sequence) == end - start + 1
    )


def _explicit_bounded_target(row, sequence):
    explicit = row.get("target_interval_bounded", row.get("search_interval_bounded"))
    if explicit not in {None, "", "NA"}:
        return explicit in {True, 1, "1", "true", "True", "yes"} and _valid_local_boundary_range(row, sequence)
    if str(row.get("boundary_class", "")).lower() in {
        "whole_locus",
        "whole_locus_search_interval",
        "unbounded_locus_search",
    }:
        return False
    return _valid_local_boundary_range(row, sequence) and (
        row.get("role") in STRUCTURAL_ROLES
        or row.get("source_feature_id") not in {None, "", "NA"}
        or row.get("boundary_class") not in {None, "", "NA"}
    )


def _genomic_matched_blocks(row, blocks, side):
    mapped = []
    for block in blocks:
        block = _coordinate_block0(block)
        local = block.query if side == "query" else block.target
        try:
            locus = ClosedInterval1(int(row["start"]), int(row["end"])).to_interval0()
            genome = local_interval_to_genome(local, locus, row.get("strand"))
            public = ClosedInterval1.from_interval0(genome)
        except (KeyError, TypeError, ValueError):
            continue
        mapped.append(
            f"{row.get('contig', 'NA')}:{public.start}-{public.end}:{row.get('strand', 'NA')}"
        )
    return ";".join(mapped) if mapped else "NA"


def _genomic_blocks0(row, blocks, side):
    mapped = []
    try:
        locus = ClosedInterval1(int(row["start"]), int(row["end"])).to_interval0()
    except (KeyError, TypeError, ValueError):
        return tuple()
    for block in blocks:
        block = _coordinate_block0(block)
        local = block.query if side == "query" else block.target
        try:
            mapped.append(local_interval_to_genome(local, locus, row.get("strand")))
        except ValueError:
            continue
    return tuple(mapped)


def _format_genomic_blocks(contig, strand, intervals):
    blocks = []
    for interval in intervals or ():
        try:
            public = ClosedInterval1.from_interval0(interval)
        except ValueError:
            continue
        blocks.append(f"{contig}:{public.start}-{public.end}:{strand}")
    return ";".join(blocks) if blocks else "NA"


def _candidate_value(candidate, key, default="NA"):
    if isinstance(candidate, dict):
        return candidate.get(key, default)
    return getattr(candidate, key, default)


def _alignment_gap_blocks(candidate):
    gaps = []
    for gap in _candidate_value(candidate, "gap_blocks", ()) or ():
        if isinstance(gap, dict):
            if gap.get("gap_in") == "query":
                gaps.append({
                    "gap_in": "query",
                    "query_cut0": int(gap.get("query_cut0", gap.get("query_start0", 0))),
                    "target_start0": int(gap["target_start0"]),
                    "target_end0": int(gap["target_end0"]),
                })
            elif gap.get("gap_in") == "target":
                gaps.append({
                    "gap_in": "target",
                    "query_start0": int(gap["query_start0"]),
                    "query_end0": int(gap["query_end0"]),
                    "target_cut0": int(gap.get("target_cut0", gap.get("target_start0", 0))),
                })
            continue
        query = getattr(gap, "query", None)
        target = getattr(gap, "target", None)
        if query is None or target is None:
            continue
        gaps.append(
            ({
                "gap_in": "query",
                "query_cut0": int(query.start0),
                "target_start0": int(target.start0),
                "target_end0": int(target.end0),
            } if query.length == 0 else {
                "gap_in": "target",
                "query_start0": int(query.start0),
                "query_end0": int(query.end0),
                "target_cut0": int(target.start0),
            })
        )
    return gaps


def _covered_bases(blocks, side):
    intervals = sorted(
        (
            (_coordinate_block0(block).query if side in {"query", 0} else _coordinate_block0(block).target)
        )
        for block in blocks
    )
    if not intervals:
        return 0
    covered = 0
    start, end = intervals[0].start0, intervals[0].end0
    for interval in intervals[1:]:
        if interval.start0 <= end:
            end = max(end, interval.end0)
        else:
            covered += end - start
            start, end = interval.start0, interval.end0
    return covered + end - start


def _unknown_pair_count(candidate):
    explicit = _candidate_value(candidate, "unknown_aligned_pairs", None)
    if explicit not in {None, "", "NA"}:
        return int(explicit)
    return sum(
        int(length)
        for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", str(_candidate_value(candidate, "cigar", "")))
        if operation == "M"
    )


def _candidate_record(candidate, rank, fallback_backend, fallback_scheme):
    if isinstance(candidate, dict):
        record = dict(candidate)
        blocks = tuple(_coordinate_block0(block) for block in record.get("aligned_blocks", ()))
        record["aligned_blocks"] = blocks
        if blocks:
            record["query_start0"] = min(block.query.start0 for block in blocks)
            record["query_end0"] = max(block.query.end0 for block in blocks)
            record["target_start0"] = min(block.target.start0 for block in blocks)
            record["target_end0"] = max(block.target.end0 for block in blocks)
        elif record.get("query_start") not in {None, "", "NA"}:
            query = ClosedInterval1(
                int(record["query_start"]), int(record["query_end"]),
            ).to_interval0()
            target = ClosedInterval1(
                int(record["target_start"]), int(record["target_end"]),
            ).to_interval0()
            record["query_start0"], record["query_end0"] = query.start0, query.end0
            record["target_start0"], record["target_end0"] = target.start0, target.end0
        for field in ("query_start", "query_end", "target_start", "target_end"):
            record.pop(field, None)
        record.setdefault("rank", rank)
        record.setdefault("backend", fallback_backend)
        record.setdefault("query_coverage", record.get("coverage", "NA"))
        record.setdefault("target_coverage", record.get("coverage", "NA"))
        record.setdefault("aligned_pairs", sum(block.query.length for block in blocks))
        record["gap_blocks"] = _alignment_gap_blocks(record)
        if record.get("score_scheme") in {None, "", "unspecified"}:
            record["score_scheme"] = fallback_scheme
        sequence_kind = record.setdefault("sequence_kind", "nucleotide")
        record.setdefault("backend_version", "NA")
        record.setdefault("raw_score", record.get("score", "NA"))
        record.setdefault("nt_identity", record.get("identity", "NA") if sequence_kind == "nucleotide" else "NA")
        record.setdefault("aa_identity", record.get("identity", "NA") if sequence_kind == "amino_acid" else "NA")
        unknown = _unknown_pair_count(record)
        record.setdefault("unknown_aligned_pairs", unknown)
        record.setdefault("known_aligned_pairs", int(record.get("aligned_pairs", 0) or 0))
        record.setdefault("query_covered_bases", _covered_bases(blocks, "query"))
        record.setdefault("target_covered_bases", _covered_bases(blocks, "target"))
        record.setdefault("query_length", "NA")
        record.setdefault("target_length", "NA")
        record.setdefault("relative_strand", record.get("strand", "+"))
        record.setdefault("mapq", record.get("mapping_quality", "NA"))
        if record.get("mapq") is None:
            record["mapq"] = "NA"
        record.setdefault("search_interval_side", "target")
        return record
    adapter_fields = dict(vars(candidate)) if hasattr(candidate, "__dict__") else {}
    for field in (
        "aligned_blocks", "gap_blocks", "alternative_hits",
        "query_interval", "target_interval",
    ):
        adapter_fields.pop(field, None)
    blocks = _alignment_blocks0(candidate)
    aligned_pairs = int(_candidate_value(candidate, "aligned_pairs", 0) or 0)
    unknown_pairs = _unknown_pair_count(candidate)
    sequence_kind = _candidate_value(candidate, "sequence_kind", "nucleotide")
    mapq = _candidate_value(
        candidate, "mapq", _candidate_value(candidate, "mapping_quality", "NA"),
    )
    if mapq is None:
        mapq = "NA"
    return {
        **adapter_fields,
        "candidate_id": _candidate_value(candidate, "candidate_id", "NA"),
        "rank": rank,
        "identity": _candidate_value(candidate, "identity", "NA"),
        "coverage": _candidate_value(candidate, "coverage", "NA"),
        "query_coverage": _candidate_value(candidate, "query_coverage", "NA"),
        "target_coverage": _candidate_value(candidate, "target_coverage", "NA"),
        "aligned_pairs": aligned_pairs,
        "query_start0": min((block.query.start0 for block in blocks), default="NA"),
        "query_end0": max((block.query.end0 for block in blocks), default="NA"),
        "target_start0": min((block.target.start0 for block in blocks), default="NA"),
        "target_end0": max((block.target.end0 for block in blocks), default="NA"),
        "strand": _candidate_value(candidate, "strand", "+"),
        "mapping_quality": _candidate_value(candidate, "mapping_quality", "NA"),
        "is_secondary": int(bool(_candidate_value(candidate, "is_secondary", rank > 1))),
        "score": _candidate_value(candidate, "score", "NA"),
        "cigar": _candidate_value(candidate, "cigar", "NA"),
        "aligned_blocks": blocks,
        "gap_blocks": _alignment_gap_blocks(candidate),
        "sequence_kind": sequence_kind,
        "backend": _candidate_value(candidate, "backend", fallback_backend),
        "backend_version": _candidate_value(candidate, "backend_version", "NA"),
        "score_scheme": _candidate_value(candidate, "score_scheme", fallback_scheme),
        "raw_score": _candidate_value(candidate, "raw_score", _candidate_value(candidate, "score", "NA")),
        "nt_identity": _candidate_value(
            candidate,
            "nt_identity",
            _candidate_value(candidate, "identity", "NA") if sequence_kind == "nucleotide" else "NA",
        ),
        "aa_identity": _candidate_value(
            candidate,
            "aa_identity",
            _candidate_value(candidate, "identity", "NA") if sequence_kind == "amino_acid" else "NA",
        ),
        "known_aligned_pairs": _candidate_value(
            candidate, "known_aligned_pairs", aligned_pairs,
        ),
        "unknown_aligned_pairs": unknown_pairs,
        "query_covered_bases": _candidate_value(
            candidate, "query_covered_bases", _covered_bases(blocks, "query"),
        ),
        "target_covered_bases": _candidate_value(
            candidate, "target_covered_bases", _covered_bases(blocks, "target"),
        ),
        "query_length": _candidate_value(candidate, "query_length", "NA"),
        "target_length": _candidate_value(candidate, "target_length", "NA"),
        "relative_strand": _candidate_value(
            candidate, "relative_strand", _candidate_value(candidate, "strand", "+"),
        ),
        "mapq": mapq,
        "left_anchor_id": _candidate_value(candidate, "left_anchor_id", "NA"),
        "right_anchor_id": _candidate_value(candidate, "right_anchor_id", "NA"),
        "search_interval": _candidate_value(candidate, "search_interval", "NA"),
        "search_interval_side": "target",
        "enumeration_complete": _candidate_value(
            candidate, "enumeration_complete", "NA",
        ),
        "incomplete_reason": _candidate_value(candidate, "incomplete_reason", "NA"),
    }


def _public_interval(interval):
    """The sole adapter for nonempty half-open intervals written to TSV JSON."""

    if not isinstance(interval, dict) or interval.get("start0") in {None, "", "NA"}:
        return interval
    internal = Interval0(int(interval["start0"]), int(interval["end0"]))
    public = ClosedInterval1.from_interval0(internal)
    return {
        key: value
        for key, value in interval.items()
        if key not in {"coordinate_system", "start0", "end0"}
    } | {
        "coordinate_system": "1-based-closed",
        "start": public.start,
        "end": public.end,
    }


def _public_gap_blocks(gaps):
    public = []
    for gap in gaps or ():
        if gap.get("gap_in") == "query":
            interval = ClosedInterval1.from_interval0(Interval0(
                int(gap["target_start0"]), int(gap["target_end0"]),
            ))
            public.append({
                "gap_in": "query",
                "query_cut0": int(gap["query_cut0"]),
                "target_start": interval.start,
                "target_end": interval.end,
            })
        elif gap.get("gap_in") == "target":
            interval = ClosedInterval1.from_interval0(Interval0(
                int(gap["query_start0"]), int(gap["query_end0"]),
            ))
            public.append({
                "gap_in": "target",
                "query_start": interval.start,
                "query_end": interval.end,
                "target_cut0": int(gap["target_cut0"]),
            })
    return public


def _public_candidate_record(record):
    public = dict(record)
    blocks = tuple(_coordinate_block0(block) for block in record.get("aligned_blocks", ()))
    public["aligned_blocks"] = []
    for block in blocks:
        query = ClosedInterval1.from_interval0(block.query)
        target = ClosedInterval1.from_interval0(block.target)
        public["aligned_blocks"].append(
            (query.start, query.end, target.start, target.end)
        )
    for side in ("query", "target"):
        start0 = record.get(f"{side}_start0")
        end0 = record.get(f"{side}_end0")
        if start0 not in {None, "", "NA"} and int(end0) > int(start0):
            interval = ClosedInterval1.from_interval0(Interval0(int(start0), int(end0)))
            public[f"{side}_start"] = interval.start
            public[f"{side}_end"] = interval.end
        public.pop(f"{side}_start0", None)
        public.pop(f"{side}_end0", None)
    public["gap_blocks"] = _public_gap_blocks(record.get("gap_blocks", ()))
    public["search_interval"] = _public_interval(record.get("search_interval"))
    for side in ("query", "target"):
        genomic = record.get(f"{side}_genomic_blocks0")
        if genomic:
            public[f"{side}_genomic_blocks"] = [
                _public_interval({"start0": interval.start0, "end0": interval.end0})
                for interval in genomic
            ]
        public.pop(f"{side}_genomic_blocks0", None)
    return public


def _nt_column_score(aln):
    score = (
        float(getattr(aln, "matches", 0) or 0) * NT_BLASTN_V1_MATCH
        + float(getattr(aln, "mismatches", 0) or 0) * NT_BLASTN_V1_MISMATCH
    )
    for length_text, operation in re.findall(r"(\d+)([ID])", str(getattr(aln, "cigar", ""))):
        length = int(length_text)
        score += NT_BLASTN_V1_GAP_OPEN
        score += max(0, length - 1) * NT_BLASTN_V1_GAP_EXTEND
    return score


def _set_explicit_alignment_score(aln):
    if getattr(aln, "backend", "") == "mafft_overlap":
        aln.score = _nt_column_score(aln)
        aln.score_scheme = "nt_blastn_v1/sum_of_column_scores"


def _alignment_candidate_records(aln, candidate_set=None):
    backend = getattr(aln, "backend", "internal")
    score_scheme = getattr(aln, "score_scheme", "unspecified")
    if candidate_set is not None:
        return [
            _candidate_record(candidate, rank, backend, candidate_set.score_scheme)
            for rank, candidate in enumerate(candidate_set.candidates, start=1)
        ]
    records = [_candidate_record(aln, 1, backend, score_scheme)] if _alignment_blocks(aln) else []
    records.extend(
        _candidate_record(candidate, rank, backend, score_scheme)
        for rank, candidate in enumerate(getattr(aln, "alternative_hits", ()) or (), start=2)
    )
    return records


def _transpose_cigar(cigar):
    return str(cigar).translate(str.maketrans({"I": "D", "D": "I"}))


def _transpose_gap_blocks(gaps):
    transposed = []
    for gap in gaps or ():
        if gap.get("gap_in") == "query":
            transposed.append({
                "gap_in": "target",
                "query_start0": int(gap["target_start0"]),
                "query_end0": int(gap["target_end0"]),
                "target_cut0": int(gap["query_cut0"]),
            })
        elif gap.get("gap_in") == "target":
            transposed.append({
                "gap_in": "query",
                "query_cut0": int(gap["target_cut0"]),
                "target_start0": int(gap["query_start0"]),
                "target_end0": int(gap["query_end0"]),
            })
    return transposed


def _transpose_candidate_record(record):
    transposed = dict(record)
    transposed["query_start0"], transposed["target_start0"] = (
        record.get("target_start0", "NA"), record.get("query_start0", "NA"),
    )
    transposed["query_end0"], transposed["target_end0"] = (
        record.get("target_end0", "NA"), record.get("query_end0", "NA"),
    )
    transposed["query_coverage"], transposed["target_coverage"] = (
        record.get("target_coverage", "NA"), record.get("query_coverage", "NA"),
    )
    transposed["query_covered_bases"], transposed["target_covered_bases"] = (
        record.get("target_covered_bases", "NA"), record.get("query_covered_bases", "NA"),
    )
    transposed["query_length"], transposed["target_length"] = (
        record.get("target_length", "NA"), record.get("query_length", "NA"),
    )
    transposed["query_occurrence_id"], transposed["target_occurrence_id"] = (
        record.get("target_occurrence_id", "NA"),
        record.get("query_occurrence_id", "NA"),
    )
    transposed["query_transcript_id"], transposed["target_transcript_id"] = (
        record.get("target_transcript_id", "NA"),
        record.get("query_transcript_id", "NA"),
    )
    transposed["aligned_blocks"] = tuple(
        CoordinateBlock(
            _coordinate_block0(block).target,
            _coordinate_block0(block).query,
        )
        for block in record.get("aligned_blocks", ())
    )
    transposed["gap_blocks"] = _transpose_gap_blocks(record.get("gap_blocks", ()))
    transposed["cigar"] = _transpose_cigar(record.get("cigar", "NA"))
    transposed["short_sequence_coverage"] = record.get("query_coverage", "NA")
    transposed["alignment_input_transposed"] = 1
    transposed["search_interval_side"] = "query"
    return transposed


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


def _mapped_genomic_interval0(row, local):
    if local is None:
        return "NA", "NA", "NA"
    try:
        locus = ClosedInterval1(int(row["start"]), int(row["end"])).to_interval0()
        genome = local_interval_to_genome(local, locus, row.get("strand"))
        public = ClosedInterval1.from_interval0(genome)
    except (KeyError, TypeError, ValueError):
        return "NA", "NA", "NA"
    return public.start, public.end, public.length


def match_evidence(left, right, seqs, context, aligner="auto", threads=1, context_aligner="minimap2", short_context_max_length=300):
    left_seq = seqs.get(left["occurrence_id"], "")
    right_seq = seqs.get(right["occurrence_id"], "")
    exon_pair = left.get("role") in EXON_LIKE_ROLES and right.get("role") in EXON_LIKE_ROLES
    requested_backend = aligner if exon_pair else context_aligner
    candidate_set = None
    left_is_short_query = (
        len(left_seq) <= int(short_context_max_length)
        and _valid_local_boundary_range(left, left_seq)
        and _explicit_bounded_target(right, right_seq)
    )
    right_is_short_query = (
        not left_is_short_query
        and len(right_seq) <= int(short_context_max_length)
        and _valid_local_boundary_range(right, right_seq)
        and _explicit_bounded_target(left, left_seq)
    )
    bounded_short_context = left_is_short_query or right_is_short_query
    alignment_input_transposed = bool(right_is_short_query)
    bounded_row = left if alignment_input_transposed else right
    # Keep left/query and right/target order for blocks, coverage, and genomic projection.
    if bounded_short_context:
        adapter_query_row, adapter_target_row = (
            (right, left) if alignment_input_transposed else (left, right)
        )
        adapter_query, adapter_target = (
            (right_seq, left_seq) if alignment_input_transposed else (left_seq, right_seq)
        )
        bounded_interval = {
            "coordinate_system": "0-based-half-open",
            "contig": bounded_row.get("contig", "NA"),
            "start0": int(bounded_row["start"]) - 1,
            "end0": int(bounded_row["end"]),
            "strand": bounded_row.get("strand", "NA"),
        }
        candidate_set = anchored_short_alignment(
            adapter_query,
            adapter_target,
            mode="local",
            query_occurrence_id=adapter_query_row.get("occurrence_id"),
            target_occurrence_id=adapter_target_row.get("occurrence_id"),
            query_transcript_id=adapter_query_row.get("transcript_id"),
            target_transcript_id=adapter_target_row.get("transcript_id"),
            search_interval=bounded_interval,
        )
        aln = candidate_set.primary
        requested_backend = "anchored_short_alignment"
        if aln is None:
            aln = AlignmentStats(
                0.0,
                0.0,
                0.0,
                backend="internal",
                alignment_mode="local",
                alignment_meaning="bounded short nucleotide local alignment",
                score_scheme=candidate_set.score_scheme,
                enumeration_complete=candidate_set.enumeration_complete,
                incomplete_reason=candidate_set.incomplete_reason,
            )
    elif exon_pair and aligner in {"mafft", "auto"}:
        aln = overlap_alignment_stats(left_seq, right_seq, backend="mafft", threads=threads)
    else:
        if not exon_pair and context_aligner not in {"internal", "minimap2", "lastz"}:
            raise ValueError("context_aligner must be internal, minimap2, or lastz")
        aln = local_alignment_stats(left_seq, right_seq, backend=requested_backend, threads=threads)
    _set_explicit_alignment_score(aln)
    identity = aln.identity
    aligned_pairs = int(getattr(aln, "aligned_pairs", 0) or 0)
    coverage = aligned_pairs / max(1, min(len(left_seq), len(right_seq)))
    blocks = _alignment_blocks0(aln)
    if alignment_input_transposed:
        blocks = tuple(CoordinateBlock(block.target, block.query) for block in blocks)
    query_interval = (
        Interval0(
            min(block.query.start0 for block in blocks),
            max(block.query.end0 for block in blocks),
        ) if blocks else None
    )
    target_interval = (
        Interval0(
            min(block.target.start0 for block in blocks),
            max(block.target.end0 for block in blocks),
        ) if blocks else None
    )
    public_query = (
        ClosedInterval1.from_interval0(query_interval)
        if query_interval is not None else None
    )
    public_target = (
        ClosedInterval1.from_interval0(target_interval)
        if target_interval is not None else None
    )
    qstart = public_query.start if public_query is not None else "NA"
    qend = public_query.end if public_query is not None else "NA"
    tstart = public_target.start if public_target is not None else "NA"
    tend = public_target.end if public_target is not None else "NA"
    q_gen_start, q_gen_end, q_gen_len = _mapped_genomic_interval0(left, query_interval)
    t_gen_start, t_gen_end, t_gen_len = _mapped_genomic_interval0(right, target_interval)
    score_scheme = getattr(aln, "score_scheme", "unspecified")
    if score_scheme in {None, "", "unspecified"}:
        score_scheme = f"{getattr(aln, 'backend', 'unknown')}_{getattr(aln, 'alignment_mode', 'alignment')}_raw_score"
        aln.score_scheme = score_scheme
    candidate_records = _alignment_candidate_records(aln, candidate_set)
    if alignment_input_transposed:
        candidate_records = [
            _transpose_candidate_record(record) for record in candidate_records
        ]
    for record in candidate_records:
        record["coverage"] = (
            int(record.get("aligned_pairs", 0) or 0)
            / max(1, min(len(left_seq), len(right_seq)))
        )
        if record.get("short_sequence_coverage") in {None, "", "NA"}:
            record["short_sequence_coverage"] = (
                record.get("query_coverage", "NA") if bounded_short_context else "NA"
            )
        record.setdefault("alignment_input_transposed", int(alignment_input_transposed))
        if record.get("query_length") in {None, "", "NA", 0, "0"}:
            record["query_length"] = len(left_seq)
        if record.get("target_length") in {None, "", "NA", 0, "0"}:
            record["target_length"] = len(right_seq)
        if record.get("backend_version") in {None, "", "NA"}:
            record["backend_version"] = getattr(aln, "backend_version", None) or "NA"
        if record.get("raw_score") in {None, "", "NA"}:
            record["raw_score"] = getattr(aln, "raw_score", None)
            if record["raw_score"] is None:
                record["raw_score"] = getattr(aln, "score", "NA")
        if (
            record.get("sequence_kind", "nucleotide") == "nucleotide"
            and record.get("nt_identity") in {None, "", "NA"}
        ):
            record["nt_identity"] = record.get("identity", identity)
        if (
            record.get("sequence_kind") == "amino_acid"
            and record.get("aa_identity") in {None, "", "NA"}
        ):
            record["aa_identity"] = record.get("identity", "NA")
        if record.get("query_covered_bases") in {None, "", "NA", 0, "0"}:
            record["query_covered_bases"] = _covered_bases(record.get("aligned_blocks", ()), "query")
        if record.get("target_covered_bases") in {None, "", "NA", 0, "0"}:
            record["target_covered_bases"] = _covered_bases(record.get("aligned_blocks", ()), "target")
        search_interval = record.get("search_interval")
        if not bounded_short_context:
            record["search_interval"] = "NA"
        elif search_interval is None or search_interval == "" or search_interval == "NA":
            record["search_interval"] = bounded_interval
        record.setdefault(
            "search_interval_side", "query" if alignment_input_transposed else "target",
        )
        if record.get("left_anchor_id") in {None, "", "NA"}:
            record["left_anchor_id"] = "NA"
        if record.get("right_anchor_id") in {None, "", "NA"}:
            record["right_anchor_id"] = "NA"
    primary_record = candidate_records[0] if candidate_records else {}
    unknown_pairs = int(
        getattr(aln, "unknown_aligned_pairs", getattr(aln, "unknown_bases", 0)) or 0
    )
    if candidate_set is not None or getattr(aln, "backend", "") == "internal":
        enumeration_complete = (
            candidate_set.enumeration_complete
            if candidate_set is not None
            else bool(getattr(aln, "enumeration_complete", True))
        )
        enumeration_status = "complete" if enumeration_complete else "incomplete"
        incomplete_reason = (
            candidate_set.incomplete_reason
            if candidate_set is not None
            else getattr(aln, "incomplete_reason", "")
        )
    else:
        enumeration_complete = False
        enumeration_status = "unassessed"
        incomplete_reason = "external_backend_candidate_enumeration_unassessed"
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
        "query_coverage": aln.target_coverage if alignment_input_transposed else aln.query_coverage,
        "target_coverage": aln.query_coverage if alignment_input_transposed else aln.target_coverage,
        "aligned_pairs": aligned_pairs,
        "sequence_kind": primary_record.get("sequence_kind", "nucleotide"),
        "backend": primary_record.get("backend", getattr(aln, "backend", "NA")),
        "backend_version": primary_record.get("backend_version", getattr(aln, "backend_version", "NA")),
        "raw_score": primary_record.get("raw_score", getattr(aln, "score", "NA")),
        "nt_identity": primary_record.get("nt_identity", identity),
        "aa_identity": primary_record.get("aa_identity", "NA"),
        "known_aligned_pairs": primary_record.get(
            "known_aligned_pairs", getattr(aln, "known_aligned_pairs", aligned_pairs),
        ),
        "unknown_aligned_pairs": primary_record.get("unknown_aligned_pairs", unknown_pairs),
        "query_covered_bases": primary_record.get("query_covered_bases", _covered_bases(blocks, "query")),
        "target_covered_bases": primary_record.get("target_covered_bases", _covered_bases(blocks, "target")),
        "query_length": len(left_seq),
        "target_length": len(right_seq),
        "relative_strand": primary_record.get("relative_strand", getattr(aln, "strand", "+")),
        "gap_blocks": primary_record.get("gap_blocks", []),
        "search_interval": bounded_interval if bounded_short_context else "NA",
        "search_interval_side": (
            "query" if alignment_input_transposed else "target"
        ) if bounded_short_context else "NA",
        "alignment_input_transposed": int(alignment_input_transposed),
        "short_sequence_coverage": primary_record.get(
            "short_sequence_coverage",
            aln.query_coverage if bounded_short_context else "NA",
        ),
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
        "matched_blocks": _format_alignment_blocks(blocks),
        "query_genomic_matched_blocks": _genomic_matched_blocks(left, blocks, "query"),
        "subject_genomic_matched_blocks": _genomic_matched_blocks(right, blocks, "subject"),
        "query_parent_feature_ids": left.get("source_feature_id", "NA"),
        "subject_parent_feature_ids": right.get("source_feature_id", "NA"),
        "query_transcript_ids": left.get("transcript_id", "NA"),
        "subject_transcript_ids": right.get("transcript_id", "NA"),
        "left_context_score": left_context,
        "right_context_score": right_context,
        "boundary_score": boundary,
        "phase_score": phase,
        "order_score": order,
        "strand_score": strand,
        "splice_score": splice,
        "size_ratio": min(segment_length(left), segment_length(right)) / max(segment_length(left), segment_length(right)),
        "alignment_cigar": _transpose_cigar(aln.cigar) if alignment_input_transposed else aln.cigar,
        "alignment_backend": aln.backend,
        "alignment_mode": getattr(aln, "alignment_mode", "local"),
        "alignment_meaning": getattr(aln, "alignment_meaning", "bounded short nucleotide local alignment"),
        "alignment_requested_backend": requested_backend,
        "mapping_quality": (
            getattr(aln, "mapping_quality", "NA")
            if getattr(aln, "backend", "internal") not in {"internal", "mafft"}
            else "NA"
        ),
        "hit_count": (
            len(candidate_set.candidates)
            if candidate_set is not None
            else getattr(aln, "hit_count", 0)
        ),
        "ambiguous_hit_count": (
            max(0, len(candidate_set.candidates) - 1)
            if candidate_set is not None
            else getattr(aln, "ambiguous_hit_count", 0)
        ),
        "alternative_hits": candidate_records[1:],
        "score_scheme": getattr(aln, "score_scheme", "unspecified"),
        "raw_alignment_score": getattr(aln, "score", "NA"),
        "enumeration_complete": enumeration_complete,
        "candidate_enumeration_status": enumeration_status,
        "incomplete_reason": incomplete_reason or "NA",
        "short_context_route": "bounded_local" if bounded_short_context else "not_used",
        "local_boundary_range": bounded_interval if bounded_short_context else "not_evaluated",
        "flanking_anchor_status": "not_evaluated_pre_chain",
        "true_absence_eligible": 0,
        "true_absence_evidence_status": "insufficient_evidence",
        "true_absence_reason": "ordered_double_flanks_not_evaluated;anchor_interval_sequence_not_extracted;assembly_continuity_unassessed;ambiguous_base_status_unassessed;query_only_deletion_gap_unassessed;alternative_alignment_concordance_unassessed",
        "candidate_records": candidate_records,
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


def sequence_supported_mapping(
    evidence,
    threshold,
    criteria=DEFAULT_CORRESPONDENCE_CRITERIA,
):
    identity = evidence["alignment_score"]
    coverage = evidence["coverage_score"]
    if identity >= threshold and coverage >= criteria.regular_min_coverage:
        return True
    if (
        identity >= threshold + criteria.high_identity_offset
        and coverage >= criteria.high_identity_min_coverage
    ):
        return True
    return False


def _candidate_sequence_accepted(
    record,
    threshold,
    short_context=False,
    criteria=DEFAULT_CORRESPONDENCE_CRITERIA,
):
    identity = to_float(record.get("identity"), 0.0)
    coverage = to_float(record.get("coverage"), 0.0)
    accepted = (
        (identity >= threshold and coverage >= criteria.regular_min_coverage)
        or (
            identity >= threshold + criteria.high_identity_offset
            and coverage >= criteria.high_identity_min_coverage
        )
    )
    if short_context:
        accepted = bool(
            accepted
            and identity >= max(criteria.short_min_identity, threshold)
            and to_float(
                record.get("short_sequence_coverage", record.get("query_coverage")),
                0.0,
            )
            >= criteria.short_min_query_coverage
            and int(to_float(
                record.get("known_aligned_pairs", record.get("aligned_pairs", 0)),
                0.0,
            ))
            >= criteria.short_min_aligned_pairs
        )
    return accepted


def assess_short_candidate_thresholds(
    segment_matches,
    output_path,
    short_identity_thresholds=(0.60, 0.70, 0.80),
    short_query_coverage_thresholds=(0.60, 0.80),
):
    identities = sorted({float(value) for value in short_identity_thresholds})
    coverages = sorted({float(value) for value in short_query_coverage_thresholds})
    if not identities or not coverages:
        raise ValueError("at least one identity and coverage threshold is required")
    if any(value < 0.0 or value > 1.0 for value in identities + coverages):
        raise ValueError("identity and coverage thresholds must be within [0, 1]")

    rows = []
    for match in read_tsv(segment_matches):
        if match.get("short_context_route") not in {
            "feature_bounded_candidate", "bounded_local", "anchor_bounded_local",
        }:
            continue
        encoded = match.get("dna_candidate_assessments")
        if encoded in {None, "", "NA"}:
            continue
        candidates = json.loads(encoded)
        if not isinstance(candidates, list):
            raise ValueError("dna_candidate_assessments must encode a JSON list")
        for candidate in candidates:
            source = candidate.get("source", "")
            score_scheme = candidate.get("score_scheme", "")
            if source != "nucleotide_alignment" or not str(score_scheme).startswith("nt_"):
                raise ValueError(
                    "dna_candidate_assessments contains a non-nucleotide candidate: "
                    f"source={source!r}, score_scheme={score_scheme!r}"
                )
            saved_pair_threshold = to_float(candidate.get("acceptance_threshold"), None)
            if saved_pair_threshold is None:
                raise ValueError(
                    "DNA candidate is missing its saved pair-specific acceptance_threshold"
                )
            for identity in identities:
                for coverage in coverages:
                    criteria = CorrespondenceCriteria(
                        short_min_identity=identity,
                        short_min_query_coverage=coverage,
                    )
                    effective_identity_cutoff = max(saved_pair_threshold, identity)
                    rows.append({
                        "match_id": match.get("match_id", "NA"),
                        "query_occurrence_id": match.get("query_occurrence_id", "NA"),
                        "subject_occurrence_id": match.get("subject_occurrence_id", "NA"),
                        "candidate_id": candidate.get("candidate_id", "NA"),
                        "saved_pair_threshold": f"{saved_pair_threshold:.6g}",
                        "short_identity_threshold": f"{identity:.6g}",
                        "effective_identity_cutoff": f"{effective_identity_cutoff:.6g}",
                        "query_coverage_threshold": f"{coverage:.6g}",
                        "candidate_identity": candidate.get("identity", "NA"),
                        "candidate_query_coverage": candidate.get(
                            "short_sequence_coverage", candidate.get("query_coverage", "NA")
                        ),
                        "candidate_aligned_pairs": candidate.get("aligned_pairs", "NA"),
                        "acceptance": int(_candidate_sequence_accepted(
                            candidate,
                            saved_pair_threshold,
                            short_context=True,
                            criteria=criteria,
                        )),
                        "interpretation": "short_DNA_acceptance_rule_sensitivity_only",
                    })
    write_tsv(
        output_path,
        rows,
        [
            "match_id", "query_occurrence_id", "subject_occurrence_id",
            "candidate_id", "saved_pair_threshold", "short_identity_threshold",
            "effective_identity_cutoff", "query_coverage_threshold",
            "candidate_identity", "candidate_query_coverage",
            "candidate_aligned_pairs", "acceptance", "interpretation",
        ],
    )
    return rows


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
    "matched_blocks",
    "query_genomic_matched_blocks",
    "subject_genomic_matched_blocks",
    "query_parent_feature_ids",
    "subject_parent_feature_ids",
    "query_transcript_ids",
    "subject_transcript_ids",
    "candidate_id",
    "alternative_candidate_ids",
    "aligned_blocks",
    "gap_blocks",
    "sequence_kind",
    "backend",
    "backend_version",
    "raw_score",
    "nt_identity",
    "aa_identity",
    "known_aligned_pairs",
    "unknown_aligned_pairs",
    "query_covered_bases",
    "target_covered_bases",
    "query_length",
    "target_length",
    "relative_strand",
    "is_secondary",
    "left_anchor_id",
    "right_anchor_id",
    "search_interval",
    "search_interval_side",
    "short_sequence_coverage",
    "alignment_input_transposed",
    "mapq",
    "mapping_quality",
    "hit_count",
    "ambiguous_hit_count",
    "alternative_hits",
    "candidate_assessments",
    "dna_candidate_assessments",
    "candidate_ids",
    "score_scheme",
    "raw_alignment_score",
    "enumeration_complete",
    "candidate_enumeration_status",
    "incomplete_reason",
    "short_context_route",
    "local_boundary_range",
    "candidate_resolution",
    "retained_candidate_ids",
    "best_path_candidate_ids",
    "chain_best_score",
    "chain_score_delta",
    "chain_configuration",
    "chain_delta_rule",
    "chain_local_mode",
    "chain_status",
    "chain_ambiguity",
    "chain_start_anchor_ids",
    "chain_end_anchor_ids",
    "chain_retained_edges",
    "chain_best_path_count_capped",
    "chain_near_optimal_path_count_capped",
    "flanking_anchor_status",
    "membership_edge_eligible",
    "membership_edge_reason",
    "position_edge_eligible",
    "position_edge_reason",
    "true_absence_eligible",
    "true_absence_evidence_status",
    "true_absence_reason",
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
    "protein_mapping_status",
    "protein_candidate_evidence_available",
    "protein_membership_eligible",
    "protein_position_eligible",
    "protein_hard_observation_eligible",
    "protein_known_aa_pairs",
    "protein_blosum62_score",
    "protein_gap_fraction",
    "protein_left_anchor_pairs",
    "protein_right_anchor_pairs",
    "protein_left_anchor_score",
    "protein_right_anchor_score",
    "protein_left_anchor_supported",
    "protein_right_anchor_supported",
    "protein_terminal_side",
    "protein_msa_column_start0",
    "protein_msa_column_end0",
    "protein_msa_column_interval",
    "protein_msa_mode",
    "protein_candidate_mapping_count",
    "protein_candidate_coordinate_consensus",
    "protein_candidate_details",
    "protein_competing_occurrences",
    "protein_query_source_features",
    "protein_target_source_features",
]


def _projection_interval(evidence, side):
    if side not in {"query", "target"}:
        raise ValueError(f"invalid projection side: {side}")
    blocks = ()
    if (
        "annotated_CDS_protein" in str(evidence.get("correspondence_basis", ""))
        and evidence.get("protein_hard_observation_eligible") in {1, "1", True}
    ):
        blocks = tuple(
            _coordinate_block0(block)
            for block in evidence.get("protein_projected_blocks", ())
        )
    if not blocks:
        accepted = [
            record
            for record in evidence.get("candidate_records", ())
            if record.get("accepted", 1) in {1, "1", True}
        ]
        if accepted:
            blocks = tuple(
                _coordinate_block0(block)
                for block in accepted[0].get("aligned_blocks", ())
            )
    if blocks:
        intervals = [getattr(block, side) for block in blocks]
        return Interval0(
            min(interval.start0 for interval in intervals),
            max(interval.end0 for interval in intervals),
        )
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
    return ClosedInterval1(start, end).to_interval0()


def _projection_record(evidence, side):
    interval = _projection_interval(evidence, side)
    if interval is None:
        return None
    return {
        "interval": interval,
        "strand": (
            "+"
            if "annotated_CDS_protein" in str(evidence.get("correspondence_basis", ""))
            and evidence.get("protein_position_eligible") in {1, "1", True}
            else evidence.get("alignment_strand", "NA")
        ),
    }


def _ordered_projection_compatible(left, right, left_ref_interval, right_ref_interval):
    if not left_ref_interval or not right_ref_interval:
        return False
    if not _disjoint_reference_intervals(left_ref_interval, right_ref_interval):
        return False
    if left.get("contig") != right.get("contig") or left.get("strand") != right.get("strand"):
        return False
    left_order = int(left["start"])
    right_order = int(right["start"])
    if left.get("strand") == "-":
        left_order, right_order = -left_order, -right_order
    copy_order = -1 if left_order < right_order else 1
    ref_order = -1 if left_ref_interval.start0 < right_ref_interval.start0 else 1
    return copy_order == ref_order


def _disjoint_reference_intervals(left_ref_interval, right_ref_interval):
    if not left_ref_interval or not right_ref_interval:
        return False
    return not left_ref_interval.overlaps(right_ref_interval)


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


def _candidate_records_for_match(evidence, match_id):
    records = [dict(record) for record in evidence.get("candidate_records", ())]
    if not records and evidence.get("projected_reference_blocks") not in {None, "", "NA"}:
        sequence_kind = (
            "amino_acid"
            if str(evidence.get("score_scheme", "")).startswith("blosum")
            else evidence.get("sequence_kind", "nucleotide")
        )
        blocks = list(parse_legacy_blocks(evidence["projected_reference_blocks"]))
        records.append(
            {
                "rank": 1,
                "identity": evidence.get("alignment_score", "NA"),
                "coverage": evidence.get("coverage_score", "NA"),
                "query_start0": min((block.query.start0 for block in blocks), default="NA"),
                "query_end0": max((block.query.end0 for block in blocks), default="NA"),
                "target_start0": min((block.target.start0 for block in blocks), default="NA"),
                "target_end0": max((block.target.end0 for block in blocks), default="NA"),
                "strand": evidence.get("alignment_strand", "+"),
                "mapping_quality": evidence.get("mapping_quality", "NA"),
                "is_secondary": 0,
                "score": evidence.get("raw_alignment_score", evidence.get("alignment_score", "NA")),
                "cigar": evidence.get("alignment_cigar", "NA"),
                "aligned_blocks": blocks,
                "query_transcript_id": evidence.get(
                    "protein_best_query_transcript",
                    evidence.get("query_transcript_ids", "NA"),
                ),
                "target_transcript_id": evidence.get(
                    "protein_best_target_transcript",
                    evidence.get("subject_transcript_ids", "NA"),
                ),
                "backend": evidence.get("alignment_backend", "NA"),
                "score_scheme": evidence.get("score_scheme", "unspecified"),
                "source": (
                    "protein_msa_projection"
                    if str(evidence.get("score_scheme", "")).startswith("blosum")
                    else "nucleotide_alignment"
                ),
                "query_coverage": evidence.get("query_coverage", "NA"),
                "target_coverage": evidence.get("target_coverage", "NA"),
                "aligned_pairs": evidence.get("aligned_pairs", 0),
                "known_aligned_pairs": evidence.get(
                    "protein_known_aa_pairs" if sequence_kind == "amino_acid" else "known_aligned_pairs",
                    evidence.get("aligned_pairs", 0),
                ),
                "unknown_aligned_pairs": evidence.get("unknown_aligned_pairs", 0),
                "query_covered_bases": _covered_bases(blocks, "query"),
                "target_covered_bases": _covered_bases(blocks, "target"),
                "query_length": evidence.get("query_length", "NA"),
                "target_length": evidence.get("target_length", "NA"),
                "gap_blocks": evidence.get("gap_blocks", []),
                "sequence_kind": sequence_kind,
                "backend_version": evidence.get("backend_version", "NA"),
                "raw_score": evidence.get("raw_score", evidence.get("raw_alignment_score", "NA")),
                "nt_identity": (
                    evidence.get("nt_identity", evidence.get("alignment_score", "NA"))
                    if sequence_kind == "nucleotide" else "NA"
                ),
                "aa_identity": (
                    evidence.get("aa_identity", evidence.get("protein_aa_identity", "NA"))
                    if sequence_kind == "amino_acid" else "NA"
                ),
                "relative_strand": evidence.get("relative_strand", evidence.get("alignment_strand", "+")),
                "mapq": evidence.get("mapq", evidence.get("mapping_quality", "NA")),
                "search_interval": evidence.get("search_interval", evidence.get("local_boundary_range", "NA")),
                "search_interval_side": evidence.get("search_interval_side", "target"),
                "accepted": int(bool(evidence.get("candidate_accepted", True))),
            }
        )
    for index, record in enumerate(records, start=1):
        if record.get("candidate_id") in {None, "", "NA"}:
            record["candidate_id"] = f"{match_id}.candidate_{index:03d}"
    return records


def _format_optional_number(value):
    if value in {None, "", "NA"}:
        return "NA"
    return f"{float(value):.6g}"


def _format_contract_value(value):
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def _match_row(left, right, evidence, score, threshold, distance_class, status, match_index):
    match_id = f"match_{match_index:05d}"
    candidate_records = _candidate_records_for_match(evidence, match_id)
    dna_candidate_records = [
        dict(record)
        for record in evidence.get("dna_candidate_assessments", evidence.get("candidate_records", ()))
    ]
    for index, record in enumerate(dna_candidate_records, start=1):
        if record.get("candidate_id") in {None, "", "NA"}:
            record["candidate_id"] = f"{match_id}.dna_candidate_{index:03d}"
    for record in candidate_records + dna_candidate_records:
        blocks = record.get("aligned_blocks", ())
        record.setdefault("query_genomic_blocks0", _genomic_blocks0(left, blocks, "query"))
        record.setdefault("target_genomic_blocks0", _genomic_blocks0(right, blocks, "target"))
        record["query_genomic_matched_blocks"] = _format_genomic_blocks(
            left.get("contig", "NA"), left.get("strand", "NA"),
            record["query_genomic_blocks0"],
        )
        record["subject_genomic_matched_blocks"] = _format_genomic_blocks(
            right.get("contig", "NA"), right.get("strand", "NA"),
            record["target_genomic_blocks0"],
        )
        record["query_parent_feature_ids"] = evidence.get(
            "query_parent_feature_ids", left.get("source_feature_id", "NA"),
        )
        record["subject_parent_feature_ids"] = evidence.get(
            "subject_parent_feature_ids", right.get("source_feature_id", "NA"),
        )
        record["query_transcript_ids"] = evidence.get(
            "query_transcript_ids", left.get("transcript_id", "NA"),
        )
        record["subject_transcript_ids"] = evidence.get(
            "subject_transcript_ids", right.get("transcript_id", "NA"),
        )
    primary = candidate_records[0] if candidate_records else {}
    public_candidates = [_public_candidate_record(record) for record in candidate_records]
    public_dna_candidates = [
        _public_candidate_record(record) for record in dna_candidate_records
    ]
    public_primary = public_candidates[0] if public_candidates else {}
    row = {
        "match_id": match_id,
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
        "matched_blocks": evidence.get("matched_blocks", evidence.get("projected_reference_blocks", "NA")),
        "query_genomic_matched_blocks": evidence.get("query_genomic_matched_blocks", "NA"),
        "subject_genomic_matched_blocks": evidence.get("subject_genomic_matched_blocks", "NA"),
        "query_parent_feature_ids": evidence.get("query_parent_feature_ids", left.get("source_feature_id", "NA")),
        "subject_parent_feature_ids": evidence.get("subject_parent_feature_ids", right.get("source_feature_id", "NA")),
        "query_transcript_ids": evidence.get("query_transcript_ids", left.get("transcript_id", "NA")),
        "subject_transcript_ids": evidence.get("subject_transcript_ids", right.get("transcript_id", "NA")),
        "candidate_id": primary.get("candidate_id", "NA"),
        "alternative_candidate_ids": ";".join(
            record["candidate_id"] for record in public_candidates[1:]
        ) or "NA",
        "aligned_blocks": _format_alignment_blocks(primary.get("aligned_blocks", ())),
        "gap_blocks": json.dumps(public_primary.get("gap_blocks", []), sort_keys=True, separators=(",", ":")),
        "sequence_kind": primary.get("sequence_kind", evidence.get("sequence_kind", "nucleotide")),
        "backend": primary.get("backend", evidence.get("backend", evidence.get("alignment_backend", "NA"))),
        "backend_version": primary.get("backend_version", evidence.get("backend_version", "NA")),
        "raw_score": _format_optional_number(primary.get("raw_score", primary.get("score", evidence.get("raw_score")))),
        "nt_identity": _format_optional_number(primary.get("nt_identity", evidence.get("nt_identity"))),
        "aa_identity": _format_optional_number(primary.get("aa_identity", evidence.get("aa_identity"))),
        "known_aligned_pairs": primary.get("known_aligned_pairs", evidence.get("known_aligned_pairs", "NA")),
        "unknown_aligned_pairs": primary.get("unknown_aligned_pairs", evidence.get("unknown_aligned_pairs", "NA")),
        "query_covered_bases": primary.get("query_covered_bases", evidence.get("query_covered_bases", "NA")),
        "target_covered_bases": primary.get("target_covered_bases", evidence.get("target_covered_bases", "NA")),
        "query_length": primary.get("query_length", evidence.get("query_length", "NA")),
        "target_length": primary.get("target_length", evidence.get("target_length", "NA")),
        "relative_strand": primary.get("relative_strand", evidence.get("relative_strand", evidence.get("alignment_strand", "NA"))),
        "is_secondary": primary.get("is_secondary", 0),
        "left_anchor_id": primary.get("left_anchor_id", "NA"),
        "right_anchor_id": primary.get("right_anchor_id", "NA"),
        "search_interval": _format_contract_value(
            _public_interval(primary.get(
                "search_interval", evidence.get("search_interval", evidence.get("local_boundary_range", "NA")),
            ))
        ),
        "search_interval_side": primary.get(
            "search_interval_side", evidence.get("search_interval_side", "NA"),
        ),
        "short_sequence_coverage": primary.get("short_sequence_coverage", evidence.get("short_sequence_coverage", "NA")),
        "alignment_input_transposed": primary.get("alignment_input_transposed", evidence.get("alignment_input_transposed", 0)),
        "mapq": primary.get("mapq", evidence.get("mapping_quality", "NA")),
        "mapping_quality": evidence.get("mapping_quality", "NA"),
        "hit_count": evidence.get("hit_count", len(candidate_records)),
        "ambiguous_hit_count": evidence.get("ambiguous_hit_count", max(0, len(candidate_records) - 1)),
        "alternative_hits": json.dumps(public_candidates[1:], sort_keys=True, separators=(",", ":")),
        "candidate_assessments": json.dumps(public_candidates, sort_keys=True, separators=(",", ":")),
        "dna_candidate_assessments": json.dumps(
            public_dna_candidates,
            sort_keys=True,
            separators=(",", ":"),
        ),
        "candidate_ids": ";".join(record["candidate_id"] for record in candidate_records) or "NA",
        "score_scheme": evidence.get("score_scheme", "unspecified"),
        "raw_alignment_score": _format_optional_number(evidence.get("raw_alignment_score")),
        "enumeration_complete": int(bool(evidence.get("enumeration_complete", False))),
        "candidate_enumeration_status": evidence.get("candidate_enumeration_status", "unassessed"),
        "incomplete_reason": evidence.get("incomplete_reason", "NA") or "NA",
        "short_context_route": evidence.get("short_context_route", "not_used"),
        "local_boundary_range": _format_contract_value(_public_interval(
            evidence.get("local_boundary_range", "not_evaluated")
        )),
        "candidate_resolution": "unassessed",
        "retained_candidate_ids": "NA",
        "best_path_candidate_ids": "NA",
        "chain_best_score": "NA",
        "chain_score_delta": "NA",
        "chain_configuration": DEFAULT_CHAIN_CONFIGURATION.name,
        "chain_delta_rule": DEFAULT_CHAIN_CONFIGURATION.delta_rule,
        "chain_local_mode": "NA",
        "chain_status": "unassessed",
        "chain_ambiguity": "unassessed",
        "chain_start_anchor_ids": "NA",
        "chain_end_anchor_ids": "NA",
        "chain_retained_edges": "NA",
        "chain_best_path_count_capped": 0,
        "chain_near_optimal_path_count_capped": 0,
        "flanking_anchor_status": evidence.get("flanking_anchor_status", "not_evaluated_pre_chain"),
        "membership_edge_eligible": 0,
        "membership_edge_reason": "candidate_chain_not_evaluated",
        "position_edge_eligible": 0,
        "position_edge_reason": "candidate_chain_not_evaluated",
        "true_absence_eligible": 0,
        "true_absence_evidence_status": evidence.get("true_absence_evidence_status", "insufficient_evidence"),
        "true_absence_reason": evidence.get("true_absence_reason", "ordered_double_flanks_not_evaluated;anchor_interval_sequence_not_extracted;assembly_continuity_unassessed;ambiguous_base_status_unassessed;query_only_deletion_gap_unassessed;alternative_alignment_concordance_unassessed"),
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
        "protein_mapping_status": evidence.get("protein_mapping_status", "uncovered"),
        "protein_candidate_evidence_available": int(bool(
            evidence.get("protein_candidate_evidence_available", False)
        )),
        "protein_membership_eligible": int(bool(evidence.get("protein_membership_eligible", False))),
        "protein_position_eligible": int(bool(evidence.get("protein_position_eligible", False))),
        "protein_hard_observation_eligible": int(bool(evidence.get("protein_hard_observation_eligible", False))),
        "protein_terminal_side": evidence.get("protein_terminal_side", "NA"),
        "protein_msa_column_start0": evidence.get("protein_msa_column_start0", "NA"),
        "protein_msa_column_end0": evidence.get("protein_msa_column_end0", "NA"),
        "protein_msa_column_interval": evidence.get("protein_msa_column_interval", "NA"),
        "protein_msa_mode": evidence.get("protein_msa_mode", "NA"),
        "protein_competing_occurrences": evidence.get("protein_competing_occurrences", "NA"),
        "protein_candidate_coordinate_consensus": int(bool(evidence.get("protein_candidate_coordinate_consensus", False))),
        "protein_candidate_details": evidence.get("protein_candidate_details", "NA"),
        "protein_query_source_features": evidence.get("protein_query_source_features", "NA"),
        "protein_target_source_features": evidence.get("protein_target_source_features", "NA"),
        **{
            field: f"{evidence[field]:.6g}" if field in evidence else "NA"
            for field in (
                "protein_aa_identity",
                "protein_query_cds_coverage",
                "protein_target_cds_coverage",
                "protein_blosum62_score",
                "protein_gap_fraction",
                "protein_left_anchor_score",
                "protein_right_anchor_score",
            )
        },
        **{
            field: evidence.get(field, "NA")
            for field in (
                "protein_known_aa_pairs",
                "protein_left_anchor_pairs",
                "protein_right_anchor_pairs",
                "protein_left_anchor_supported",
                "protein_right_anchor_supported",
                "protein_candidate_mapping_count",
            )
        },
    }
    row["_candidate_records"] = candidate_records
    return row


def _copy_transcription_bounds(occurrences):
    bounds = {}
    grouped = defaultdict(list)
    for occurrence in occurrences:
        grouped[occurrence_copy_key(occurrence)].append(occurrence)
    for key, rows in grouped.items():
        bounds[key] = (
            min(int(row["start"]) for row in rows),
            max(int(row["end"]) for row in rows),
        )
    return bounds


def _candidate_copy_interval(occurrence, record, side, copy_bounds):
    genomic_blocks = tuple(record.get(f"{side}_genomic_blocks0", ()) or ())
    if genomic_blocks:
        genomic = Interval0(
            min(interval.start0 for interval in genomic_blocks),
            max(interval.end0 for interval in genomic_blocks),
        )
    else:
        start0 = record.get(f"{side}_start0")
        end0 = record.get(f"{side}_end0")
        if start0 in {None, "", "NA"} or end0 in {None, "", "NA"}:
            return None
        try:
            locus = ClosedInterval1(
                int(occurrence["start"]), int(occurrence["end"]),
            ).to_interval0()
            genomic = local_interval_to_genome(
                Interval0(int(start0), int(end0)), locus, occurrence.get("strand"),
            )
        except (KeyError, TypeError, ValueError):
            return None
    if genomic.length == 0:
        return None
    copy_start, copy_end = copy_bounds[occurrence_copy_key(occurrence)]
    if occurrence.get("strand") == "-":
        return Interval0(copy_end - genomic.end0, copy_end - genomic.start0)
    return Interval0(genomic.start0 - (copy_start - 1), genomic.end0 - (copy_start - 1))


def _chain_precedes(left, right):
    return left.query.end0 <= right.query.start0 and left.target.end0 <= right.target.start0


def _candidate_path_memberships(
    query,
    subject,
    transcript_paths,
    reverse=False,
    query_transcript_ids=None,
    subject_transcript_ids=None,
):
    def transcript_tokens(value):
        return {
            token
            for token in str(value or "").replace(",", ";").split(";")
            if token and token != "NA"
        }

    allowed_query_transcripts = transcript_tokens(query_transcript_ids)
    allowed_subject_transcripts = transcript_tokens(subject_transcript_ids)
    by_occurrence = defaultdict(list)
    for path in transcript_paths or ():
        by_occurrence[path.get("occurrence_id")].append(path)
    memberships = []
    for query_path in by_occurrence.get(query.get("occurrence_id"), ()):
        if (
            allowed_query_transcripts
            and query_path.get("transcript_id") not in allowed_query_transcripts
        ):
            continue
        for subject_path in by_occurrence.get(subject.get("occurrence_id"), ()):
            if (
                allowed_subject_transcripts
                and subject_path.get("transcript_id") not in allowed_subject_transcripts
            ):
                continue
            if query_path.get("transcript_id") in {None, "", "NA"}:
                continue
            if subject_path.get("transcript_id") in {None, "", "NA"}:
                continue
            try:
                query_order = int(query_path.get("path_rank", query_path.get("transcript_order")))
                subject_order = int(subject_path.get("path_rank", subject_path.get("transcript_order")))
            except (TypeError, ValueError):
                continue
            query_membership = {
                "path_id": "|".join(
                    str(query_path.get(field, query.get(field, "NA")))
                    for field in ("family_id", "species", "gene_copy_id", "transcript_id")
                ),
                "order": query_order,
                "contig": query_path.get("contig", query.get("contig", "NA")),
                "strand": query_path.get("strand", query.get("strand", "NA")),
            }
            subject_membership = {
                "path_id": "|".join(
                    str(subject_path.get(field, subject.get(field, "NA")))
                    for field in ("family_id", "species", "gene_copy_id", "transcript_id")
                ),
                "order": subject_order,
                "contig": subject_path.get("contig", subject.get("contig", "NA")),
                "strand": subject_path.get("strand", subject.get("strand", "NA")),
            }
            if reverse:
                query_membership, subject_membership = subject_membership, query_membership
            if (
                query_membership["strand"] not in {"+", "-"}
                or subject_membership["strand"] not in {"+", "-"}
                or query_membership["contig"] in {None, "", "NA"}
                or subject_membership["contig"] in {None, "", "NA"}
            ):
                continue
            memberships.append(
                ChainPathMembership(
                    query_path_id=query_membership["path_id"],
                    target_path_id=subject_membership["path_id"],
                    query_order=query_membership["order"],
                    target_order=subject_membership["order"],
                    query_contig=query_membership["contig"],
                    target_contig=subject_membership["contig"],
                    query_strand=query_membership["strand"],
                    target_strand=subject_membership["strand"],
                )
            )
    return tuple(sorted(set(memberships), key=lambda item: (item.context, item.query_order, item.target_order)))


def _apply_ordered_candidate_chains(rows, occurrence_by_id, occurrences, transcript_paths=None):
    for row in rows:
        row["membership_edge_eligible"] = 0
        row["position_edge_eligible"] = 0
        row["_membership_edge_eligible"] = False
        row["_position_edge_eligible"] = False
        row["retained_candidate_ids"] = "NA"
        row["best_path_candidate_ids"] = "NA"
        row["chain_best_path_count_capped"] = 0
        row["chain_near_optimal_path_count_capped"] = 0
        row["chain_status"] = "unassessed"
        row["chain_ambiguity"] = "unassessed"
        row["chain_start_anchor_ids"] = "NA"
        row["chain_end_anchor_ids"] = "NA"
        row["chain_retained_edges"] = "NA"
        row["left_anchor_id"] = "NA"
        row["right_anchor_id"] = "NA"
    copy_bounds = _copy_transcription_bounds(occurrences)
    groups = defaultdict(list)
    owner_by_candidate = {}
    candidate_by_id = {}
    record_by_candidate = {}
    original_intervals = {}
    for row in rows:
        if row.get("match_status") != "mapped":
            continue
        query = occurrence_by_id.get(row.get("query_occurrence_id"))
        subject = occurrence_by_id.get(row.get("subject_occurrence_id"))
        if not query or not subject:
            continue
        query_key = occurrence_copy_key(query)
        subject_key = occurrence_copy_key(subject)
        reverse = subject_key < query_key
        group_copies = tuple(sorted((query_key, subject_key)))
        for record in row.get("_candidate_records", ()):
            if record.get("accepted") not in {1, "1", True}:
                continue
            query_interval = _candidate_copy_interval(
                query, record, "query", copy_bounds,
            )
            target_interval = _candidate_copy_interval(
                subject, record, "target", copy_bounds,
            )
            if query_interval is None or target_interval is None:
                continue
            original_query_interval = query_interval
            original_target_interval = target_interval
            if reverse:
                query_interval, target_interval = target_interval, query_interval
            try:
                raw_score = float(record.get("score"))
            except (TypeError, ValueError):
                continue
            score_scheme = str(record.get("score_scheme") or row.get("score_scheme") or "unspecified")
            alignment_strand = str(record.get("strand") or row.get("alignment_strand") or "+")
            candidate = ChainCandidate(
                candidate_id=record["candidate_id"],
                query=query_interval,
                target=target_interval,
                score=raw_score,
                score_scheme=score_scheme,
                relative_strand=alignment_strand if alignment_strand in {"+", "-"} else "+",
                path_memberships=_candidate_path_memberships(
                    query,
                    subject,
                    transcript_paths,
                    reverse=reverse,
                    query_transcript_ids=record.get(
                        "query_transcript_id", row.get("query_transcript_ids"),
                    ),
                    subject_transcript_ids=record.get(
                        "target_transcript_id", row.get("subject_transcript_ids"),
                    ),
                ),
            )
            owner_by_candidate[candidate.candidate_id] = row
            candidate_by_id[candidate.candidate_id] = candidate
            record_by_candidate[candidate.candidate_id] = record
            original_intervals[candidate.candidate_id] = (
                original_query_interval, original_target_interval,
            )
            groups[(*group_copies, score_scheme)].append(candidate)

    retained_ids = set()
    best_ids = set()
    chain_metadata = {}
    anchor_status_by_candidate = {}
    left_flank_ids_by_candidate = defaultdict(set)
    right_flank_ids_by_candidate = defaultdict(set)

    def owner_has_resolved_anchor_position(owner):
        if owner.get("enumeration_complete") not in {1, "1", True}:
            return False
        if owner.get("short_context_route") in {
            "feature_bounded_candidate", "anchor_bounded_unavailable",
        }:
            return False
        if (
            "annotated_CDS_protein" in str(owner.get("correspondence_basis", ""))
            and owner.get("protein_hard_observation_eligible") not in {1, "1", True}
        ):
            return False
        signatures = {
            (
                original_intervals[candidate_id][0],
                original_intervals[candidate_id][1],
                tuple(
                    _block_signature(block)
                    for block in record_by_candidate[candidate_id].get("aligned_blocks", ())
                ),
            )
            for record in owner.get("_candidate_records", ())
            for candidate_id in (record.get("candidate_id"),)
            if record.get("accepted") in {1, "1", True}
            and candidate_id in original_intervals
        }
        return len(signatures) == 1

    def fixed_anchor_ids(group):
        anchors = {}
        contexts = {
            membership.context
            for candidate in group
            for membership in candidate.path_memberships
        }
        for context in contexts:
            start_ids, end_ids = set(), set()
            path_candidates = []
            for candidate in group:
                membership = next(
                    (item for item in candidate.path_memberships if item.context == context),
                    None,
                )
                if membership is None or not owner_has_resolved_anchor_position(
                    owner_by_candidate[candidate.candidate_id]
                ):
                    continue
                path_candidates.append((candidate, membership))
            if len({
                (
                    owner_by_candidate[candidate.candidate_id]["query_occurrence_id"],
                    owner_by_candidate[candidate.candidate_id]["subject_occurrence_id"],
                )
                for candidate, _membership in path_candidates
            }) < 2:
                anchors[context] = (start_ids, end_ids)
                continue
            ordered = sorted(
                path_candidates,
                key=lambda item: (
                    item[1].query_order,
                    item[1].target_order,
                    item[0].query.start0,
                    item[0].target.start0,
                    item[0].candidate_id,
                ),
            )
            first_order = (ordered[0][1].query_order, ordered[0][1].target_order)
            last_order = (ordered[-1][1].query_order, ordered[-1][1].target_order)
            if first_order == last_order:
                anchors[context] = (start_ids, end_ids)
                continue
            start_ids.update(
                candidate.candidate_id
                for candidate, membership in ordered
                if (membership.query_order, membership.target_order) == first_order
            )
            end_ids.update(
                candidate.candidate_id
                for candidate, membership in ordered
                if (membership.query_order, membership.target_order) == last_order
            )
            anchors[context] = (start_ids, end_ids)
        return anchors

    for group in groups.values():
        collinear = [
            candidate
            for candidate in group
            if candidate.relative_strand == "+"
            and (not transcript_paths or candidate.path_memberships)
        ]
        for candidate in group:
            if candidate.relative_strand == "-":
                chain_metadata[candidate.candidate_id] = {
                    "status": "noncollinear_candidate",
                    "best_score": "NA",
                    "score_delta": "NA",
                    "local_mode": "NA",
                }
            elif transcript_paths and not candidate.path_memberships:
                chain_metadata[candidate.candidate_id] = {
                    "status": "incompatible_transcript_paths",
                    "best_score": "NA",
                    "score_delta": "NA",
                    "local_mode": "NA",
                }
        if not collinear:
            continue
        context_anchor_ids = fixed_anchor_ids(collinear)
        exact = ordered_candidate_chain(
            collinear,
            0.0,
            context_anchor_ids=context_anchor_ids,
            configuration_name=DEFAULT_CHAIN_CONFIGURATION.name,
        )
        scheme = collinear[0].score_scheme
        delta = DEFAULT_CHAIN_CONFIGURATION.score_delta(
            exact.best_score, scheme,
        )
        result = ordered_candidate_chain(
            collinear,
            delta,
            context_anchor_ids=context_anchor_ids,
            configuration_name=DEFAULT_CHAIN_CONFIGURATION.name,
        )
        retained_ids.update(result.retained_ids)
        best_ids.update(result.best_path_member_ids)
        for candidate in collinear:
            status = "outside_near_optimal_chain"
            if candidate.candidate_id in result.retained_ids:
                status = "retained_near_optimal"
            if candidate.candidate_id in result.best_path_member_ids:
                status = "best_path_member"
            chain_metadata[candidate.candidate_id] = {
                "status": status,
                "best_score": f"{result.best_score:.6g}",
                "score_delta": f"{result.score_delta:.6g}",
                "local_mode": int(result.local_mode),
                "configuration": result.configuration_name,
                "ambiguity": result.ambiguity_status,
                "candidate_ambiguous": candidate.candidate_id in result.ambiguous_ids,
                "start_anchor_ids": ";".join(sorted(result.start_anchor_ids)) or "NA",
                "end_anchor_ids": ";".join(sorted(result.end_anchor_ids)) or "NA",
                "retained_edges": ";".join(
                    f"{left}>{right}" for left, right in sorted(result.retained_edges)
                ) or "NA",
                "best_path_count_capped": result.best_path_count_capped,
                "near_optimal_path_count_capped": result.near_optimal_path_count_capped,
            }
            focus_owner = owner_by_candidate[candidate.candidate_id]

            def independent_structure_unit(peer):
                peer_owner = owner_by_candidate[peer.candidate_id]
                return (
                    peer_owner["query_occurrence_id"] != focus_owner["query_occurrence_id"]
                    and peer_owner["subject_occurrence_id"] != focus_owner["subject_occurrence_id"]
                    and owner_has_resolved_anchor_position(peer_owner)
                    and peer.candidate_id not in result.ambiguous_ids
                )

            saw_one_sided = False
            saw_double_sided = False
            for summary in result.context_summaries:
                if candidate.candidate_id not in summary.retained_ids:
                    continue
                context_candidates = [
                    candidate_by_id[candidate_id]
                    for candidate_id in summary.retained_ids
                ]

                def path_reachable(source_id, target_id):
                    pending = [source_id]
                    seen = set()
                    while pending:
                        current = pending.pop()
                        if current == target_id:
                            return True
                        if current in seen:
                            continue
                        seen.add(current)
                        pending.extend(
                            right for left, right in summary.retained_edges
                            if left == current
                        )
                    return False

                left_anchors = [
                    peer for peer in context_candidates
                    if peer.candidate_id != candidate.candidate_id
                    and independent_structure_unit(peer)
                    and (
                        path_reachable(peer.candidate_id, candidate.candidate_id)
                        if summary.retained_edges else _chain_precedes(peer, candidate)
                    )
                ]
                right_anchors = [
                    peer for peer in context_candidates
                    if peer.candidate_id != candidate.candidate_id
                    and independent_structure_unit(peer)
                    and (
                        path_reachable(candidate.candidate_id, peer.candidate_id)
                        if summary.retained_edges else _chain_precedes(candidate, peer)
                    )
                ]
                if left_anchors:
                    nearest_left = max(
                        left_anchors,
                        key=lambda item: (
                            item.query.end0, item.target.end0, item.candidate_id,
                        ),
                    )
                    left_flank_ids_by_candidate[candidate.candidate_id].add(
                        nearest_left.candidate_id
                    )
                if right_anchors:
                    nearest_right = min(
                        right_anchors,
                        key=lambda item: (
                            item.query.start0, item.target.start0, item.candidate_id,
                        ),
                    )
                    right_flank_ids_by_candidate[candidate.candidate_id].add(
                        nearest_right.candidate_id
                    )
                saw_one_sided |= bool(left_anchors or right_anchors)
                saw_double_sided |= bool(left_anchors and right_anchors)
            if saw_double_sided:
                anchor_status_by_candidate[candidate.candidate_id] = (
                    "ordered_double_sided_homologous_flanks_same_path"
                )
            elif saw_one_sided:
                anchor_status_by_candidate[candidate.candidate_id] = (
                    "ordered_one_sided_independent_homologous_flank"
                )
            else:
                anchor_status_by_candidate[candidate.candidate_id] = (
                    "no_independent_homologous_flanks_on_same_path"
                )

    partners_by_query = defaultdict(set)
    partners_by_subject = defaultdict(set)
    for row in rows:
        own_ids = {record["candidate_id"] for record in row.get("_candidate_records", ())}
        retained = own_ids & retained_ids
        best = own_ids & best_ids
        row["retained_candidate_ids"] = ";".join(sorted(retained)) or "NA"
        row["best_path_candidate_ids"] = ";".join(sorted(best)) or "NA"
        metadata = [chain_metadata[candidate_id] for candidate_id in retained if candidate_id in chain_metadata]
        if metadata:
            row["chain_best_score"] = max(
                metadata, key=lambda item: to_float(item["best_score"], float("-inf"))
            )["best_score"]
            row["chain_score_delta"] = max(
                metadata, key=lambda item: to_float(item["score_delta"], float("-inf"))
            )["score_delta"]
            row["chain_configuration"] = metadata[0]["configuration"]
            row["chain_delta_rule"] = DEFAULT_CHAIN_CONFIGURATION.delta_rule
            row["chain_local_mode"] = int(any(item["local_mode"] == 1 for item in metadata))
            row["chain_ambiguity"] = (
                "multiple_near_optimal_chains"
                if any(item["candidate_ambiguous"] for item in metadata)
                else "unique_within_reported_candidates"
            )
            row["chain_best_path_count_capped"] = max(
                item["best_path_count_capped"] for item in metadata
            )
            row["chain_near_optimal_path_count_capped"] = max(
                item["near_optimal_path_count_capped"] for item in metadata
            )
            row["chain_start_anchor_ids"] = ";".join(sorted({
                token
                for item in metadata
                for token in item["start_anchor_ids"].split(";")
                if token != "NA"
            })) or "NA"
            row["chain_end_anchor_ids"] = ";".join(sorted({
                token
                for item in metadata
                for token in item["end_anchor_ids"].split(";")
                if token != "NA"
            })) or "NA"
            row["left_anchor_id"] = ";".join(sorted({
                anchor_id
                for candidate_id in retained
                for anchor_id in left_flank_ids_by_candidate[candidate_id]
            })) or "NA"
            row["right_anchor_id"] = ";".join(sorted({
                anchor_id
                for candidate_id in retained
                for anchor_id in right_flank_ids_by_candidate[candidate_id]
            })) or "NA"
            row["chain_retained_edges"] = ";".join(sorted({
                token
                for item in metadata
                for token in item["retained_edges"].split(";")
                if token != "NA"
            })) or "NA"
            row["chain_status"] = (
                "best_path_member" if best else "retained_near_optimal"
            )
        elif own_ids & set(chain_metadata):
            own_statuses = {
                chain_metadata[candidate_id]["status"]
                for candidate_id in own_ids
                if candidate_id in chain_metadata
            }
            row["chain_status"] = (
                "incompatible_transcript_paths"
                if "incompatible_transcript_paths" in own_statuses
                else "noncollinear_candidate"
                if "noncollinear_candidate" in own_statuses
                else "outside_near_optimal_chain"
            )
        anchor_states = {
            anchor_status_by_candidate[candidate_id]
            for candidate_id in retained
            if candidate_id in anchor_status_by_candidate
        }
        if "ordered_double_sided_homologous_flanks_same_path" in anchor_states:
            row["flanking_anchor_status"] = "ordered_double_sided_homologous_flanks_same_path"
        elif "ordered_one_sided_independent_homologous_flank" in anchor_states:
            row["flanking_anchor_status"] = "ordered_one_sided_independent_homologous_flank"
        else:
            row["flanking_anchor_status"] = "no_independent_homologous_flanks_on_same_path"
        for record in row.get("_candidate_records", ()):
            if record.get("candidate_id") not in retained:
                continue
            candidate_id = record["candidate_id"]
            record["left_anchor_id"] = ";".join(sorted(
                left_flank_ids_by_candidate[candidate_id]
            )) or "NA"
            record["right_anchor_id"] = ";".join(sorted(
                right_flank_ids_by_candidate[candidate_id]
            )) or "NA"
            record["chain_configuration"] = row.get("chain_configuration", DEFAULT_CHAIN_CONFIGURATION.name)
            record["chain_score_delta"] = row.get("chain_score_delta", "NA")
        public_candidates = [
            _public_candidate_record(record)
            for record in row.get("_candidate_records", ())
        ]
        row["candidate_assessments"] = json.dumps(
            public_candidates, sort_keys=True, separators=(",", ":"),
        )
        row["alternative_hits"] = json.dumps(
            public_candidates[1:], sort_keys=True, separators=(",", ":"),
        )
        if retained and row.get("match_status") == "mapped":
            query_id = row["query_occurrence_id"]
            subject_id = row["subject_occurrence_id"]
            query_copy = occurrence_copy_key(occurrence_by_id[query_id])
            subject_copy = occurrence_copy_key(occurrence_by_id[subject_id])
            partners_by_query[(query_id, subject_copy)].add(subject_id)
            partners_by_subject[(subject_id, query_copy)].add(query_id)

    def compatible_partner_projections(shared_id, other_copy, shared_side):
        intervals = []
        for candidate_id in retained_ids:
            owner = owner_by_candidate[candidate_id]
            if candidate_id not in original_intervals:
                continue
            query_id = owner["query_occurrence_id"]
            subject_id = owner["subject_occurrence_id"]
            if shared_side == "query":
                if query_id != shared_id or occurrence_copy_key(occurrence_by_id[subject_id]) != other_copy:
                    continue
                intervals.append(original_intervals[candidate_id][0])
            else:
                if subject_id != shared_id or occurrence_copy_key(occurrence_by_id[query_id]) != other_copy:
                    continue
                intervals.append(original_intervals[candidate_id][1])
        intervals.sort()
        return all(left.end0 <= right.start0 for left, right in zip(intervals, intervals[1:]))

    for row in rows:
        retained = set(str(row.get("retained_candidate_ids", "NA")).split(";")) - {"NA", ""}
        # Direct callers and legacy rows may provide accepted candidate
        # records without the chain annotation pass.  Their membership still
        # carries evidence; coordinate eligibility is decided below from the
        # number and identity of placements.
        if not retained:
            retained = {
                record.get("candidate_id")
                for record in row.get("_candidate_records", ())
                if record.get("candidate_id") and record.get("accepted") in {1, "1", True}
            }
        query_id = row["query_occurrence_id"]
        subject_id = row["subject_occurrence_id"]
        query_copy = occurrence_copy_key(occurrence_by_id[query_id])
        subject_copy = occurrence_copy_key(occurrence_by_id[subject_id])
        query_partners = partners_by_query[(query_id, subject_copy)]
        subject_partners = partners_by_subject[(subject_id, query_copy)]
        query_partner_compatible = (
            len(query_partners) <= 1
            or compatible_partner_projections(query_id, subject_copy, "query")
        )
        subject_partner_compatible = (
            len(subject_partners) <= 1
            or compatible_partner_projections(subject_id, query_copy, "subject")
        )
        placement_signatures = {
            (
                original_intervals[candidate_id][0],
                original_intervals[candidate_id][1],
                tuple(
                    _block_signature(block)
                    for block in record_by_candidate[candidate_id].get("aligned_blocks", ())
                ),
            )
            for candidate_id in retained
            if candidate_id in record_by_candidate and candidate_id in original_intervals
        }
        unique_position = len(placement_signatures) == 1
        enumeration_complete = row.get("enumeration_complete") in {1, "1", True}
        chain_unambiguous = (
            row.get("chain_ambiguity", "unique_within_reported_candidates")
            == "unique_within_reported_candidates"
        )
        double_flanks = (
            row.get("flanking_anchor_status")
            == "ordered_double_sided_homologous_flanks_same_path"
        )
        protein_hard = row.get("protein_hard_observation_eligible") in {1, "1", True}
        protein_position = row.get("protein_position_eligible") in {1, "1", True}
        short_route = row.get("short_context_route")
        short_context_supported = short_route not in {
            "feature_bounded_candidate", "bounded_local", "anchor_bounded_local",
            "anchor_bounded_unavailable",
        } or (
            short_route == "anchor_bounded_local" and double_flanks
        )
        chain_membership = bool(
            retained
            and query_partner_compatible
            and subject_partner_compatible
            and short_context_supported
        )
        protein_basis = (
            "annotated_CDS_protein"
            in str(row.get("correspondence_basis", ""))
        )
        membership_eligible = bool(
            row.get("match_status") == "mapped"
            and (chain_membership or protein_hard)
            and (not protein_basis or protein_hard)
        )
        if row.get("match_status") != "mapped":
            membership_reason = "sequence_correspondence_not_accepted"
        elif protein_basis and not protein_hard:
            membership_reason = "protein_candidate_without_hard_coordinates"
        elif protein_hard and not chain_membership:
            membership_reason = "resolved_annotated_CDS_protein_membership"
        elif not retained:
            membership_reason = "no_retained_accepted_candidate"
        elif not unique_position:
            membership_reason = "multiple_accepted_retained_coordinate_placements"
        elif not chain_unambiguous:
            membership_reason = "multiple_near_optimal_candidate_chains"
        elif not query_partner_compatible or not subject_partner_compatible:
            membership_reason = "overlapping_partner_projections"
        elif not short_context_supported:
            membership_reason = "short_context_without_same_path_independent_double_flanks"
        else:
            membership_reason = "retained_sequence_membership"

        position_eligible = bool(
            membership_eligible
            and enumeration_complete
            and unique_position
            and (len(retained) >= 1 or (protein_basis and protein_position))
            and (
                not protein_basis
                or protein_position
            )
        )
        if not membership_eligible:
            position_reason = membership_reason
        elif not enumeration_complete:
            position_reason = "candidate_enumeration_unassessed_or_incomplete"
        elif not unique_position:
            position_reason = "multiple_accepted_retained_coordinate_placements"
        elif (
            protein_basis
            and not protein_position
        ):
            position_reason = "protein_membership_without_resolved_coordinates"
        else:
            position_reason = "unique_resolved_actual_coordinates"

        row["membership_edge_eligible"] = int(membership_eligible)
        row["membership_edge_reason"] = membership_reason
        row["position_edge_eligible"] = int(position_eligible)
        row["position_edge_reason"] = position_reason
        if row.get("match_status") != "mapped":
            row["candidate_resolution"] = "candidate"
        elif protein_hard and not retained:
            row["candidate_resolution"] = "ambiguous"
        elif not retained:
            row["candidate_resolution"] = "excluded"
        elif not position_eligible:
            row["candidate_resolution"] = "ambiguous"
        else:
            row["candidate_resolution"] = "resolved"
        row["_membership_edge_eligible"] = membership_eligible
        row["_position_edge_eligible"] = position_eligible
        if row.get("match_status") == "mapped" and not membership_eligible:
            if not retained:
                row["match_status"] = "candidate_chain_excluded"
            elif not short_context_supported:
                row["match_status"] = "candidate_unanchored"
            else:
                row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
        elif row.get("match_status") == "mapped" and not enumeration_complete:
            row["match_status"] = "candidate_search_incomplete"
            row["candidate_resolution"] = "ambiguous"

        row["true_absence_eligible"] = 0
        if row.get("alignment_backend") == "genomic_overlap":
            row["true_absence_evidence_status"] = "not_applicable"
            row["true_absence_reason"] = "same_locus_annotation_overlap_is_not_deletion_evidence"
            continue
        absence_reasons = []
        if not double_flanks:
            absence_reasons.append("same_path_independent_double_flanks_not_established")
        if not enumeration_complete:
            absence_reasons.append("acceptable_alternative_alignment_set_unassessed_or_incomplete")
        if not unique_position:
            absence_reasons.append("acceptable_alternatives_do_not_define_one_position")
        absence_reasons.extend(
            [
                "anchor_interval_sequence_not_extracted",
                "assembly_continuity_unassessed",
                "ambiguous_base_status_unassessed",
                "query_only_deletion_gap_unassessed",
                "alternative_alignment_concordance_unassessed",
            ]
        )
        row["true_absence_evidence_status"] = (
            "evidence_candidate" if double_flanks else "insufficient_evidence"
        )
        row["true_absence_reason"] = ";".join(absence_reasons)


def _gene_locus_records(input_dir):
    input_dir = Path(input_dir)
    sequences = parse_fasta(input_dir / "gene_loci.fasta")
    metadata = {
        (row.get("species"), row.get("gene_copy_id")): row
        for row in read_tsv(input_dir / "gene_loci.tsv", optional=True)
    }
    records = {}
    for header, sequence in sequences.items():
        try:
            species, gene_copy_id, geometry = header.split("|", 2)
            contig, bounds, strand = geometry.rsplit(":", 2)
            start, end = (int(value) for value in bounds.split("-", 1))
            interval = ClosedInterval1(start, end).to_interval0()
        except (TypeError, ValueError):
            continue
        if strand not in {"+", "-"} or interval.length != len(sequence):
            continue
        row = metadata.get((species, gene_copy_id), {})
        records[(species, gene_copy_id)] = {
            "header": header,
            "sequence": sequence,
            "contig": contig,
            "strand": strand,
            "interval": interval,
            "range_status": row.get("range_status", "NA"),
        }
    return records


def _value_tokens(value):
    return {
        token
        for token in str(value or "").replace(",", ";").split(";")
        if token and token != "NA"
    }


def _candidate_blocks_for_copy(record, owner, copy_key, occurrence_by_id):
    query = occurrence_by_id.get(owner.get("query_occurrence_id"), {})
    subject = occurrence_by_id.get(owner.get("subject_occurrence_id"), {})
    if occurrence_copy_key(query) == copy_key:
        return tuple(record.get("query_genomic_blocks0", ()) or ())
    if occurrence_copy_key(subject) == copy_key:
        return tuple(record.get("target_genomic_blocks0", ()) or ())
    return tuple()


def _owner_occurrence_for_copy(owner, copy_key, occurrence_by_id):
    query = occurrence_by_id.get(owner.get("query_occurrence_id"), {})
    subject = occurrence_by_id.get(owner.get("subject_occurrence_id"), {})
    if occurrence_copy_key(query) == copy_key:
        return query
    if occurrence_copy_key(subject) == copy_key:
        return subject
    return None


def _interval_between_transcript_flanks(left_blocks, right_blocks, strand):
    if not left_blocks or not right_blocks:
        return None
    left = Interval0(
        min(block.start0 for block in left_blocks),
        max(block.end0 for block in left_blocks),
    )
    right = Interval0(
        min(block.start0 for block in right_blocks),
        max(block.end0 for block in right_blocks),
    )
    if strand == "+" and left.end0 <= right.start0:
        return Interval0(left.end0, right.start0)
    if strand == "-" and right.end0 <= left.start0:
        return Interval0(right.end0, left.start0)
    return None


def _cut0_between_loci(cut0, inner_locus, outer_locus, strand):
    cut0 = int(cut0)
    if cut0 < 0 or cut0 > inner_locus.length:
        raise ValueError("cut lies outside the inner locus")
    genome_cut0 = (
        inner_locus.start0 + cut0
        if strand == "+"
        else inner_locus.end0 - cut0
    )
    if genome_cut0 < outer_locus.start0 or genome_cut0 > outer_locus.end0:
        raise ValueError("cut lies outside the parent feature")
    return (
        genome_cut0 - outer_locus.start0
        if strand == "+"
        else outer_locus.end0 - genome_cut0
    )


def _anchor_bounded_gap_blocks(gaps, search_interval, parent_interval, strand):
    projected = []
    for gap in gaps or ():
        if gap.get("gap_in") == "query":
            target = genome_interval_to_local(
                local_interval_to_genome(
                    Interval0(
                        int(gap["target_start0"]), int(gap["target_end0"]),
                    ),
                    search_interval,
                    strand,
                ),
                parent_interval,
                strand,
            )
            projected.append({
                "gap_in": "query",
                "query_cut0": int(gap["query_cut0"]),
                "target_start0": target.start0,
                "target_end0": target.end0,
            })
        elif gap.get("gap_in") == "target":
            projected.append({
                "gap_in": "target",
                "query_start0": int(gap["query_start0"]),
                "query_end0": int(gap["query_end0"]),
                "target_cut0": _cut0_between_loci(
                    gap["target_cut0"], search_interval, parent_interval, strand,
                ),
            })
    return projected


def _context_flank_pairs(
    focus_row,
    source,
    bounded,
    candidate_owner,
    candidate_record,
    occurrence_by_id,
    transcript_paths,
):
    if not transcript_paths:
        return set()
    if (
        source.get("contig") in {None, "", "NA"}
        or bounded.get("contig") in {None, "", "NA"}
        or source.get("strand") not in {"+", "-"}
        or bounded.get("strand") not in {"+", "-"}
    ):
        return set()
    focal_memberships = _candidate_path_memberships(
        source, bounded, transcript_paths,
    )
    if not focal_memberships:
        return set()
    source_copy = occurrence_copy_key(source)
    bounded_copy = occurrence_copy_key(bounded)

    def distinct_placement_ids(candidate_ids):
        by_placement = defaultdict(list)
        for candidate_id in candidate_ids:
            owner = candidate_owner[candidate_id]
            record = candidate_record[candidate_id]
            source_blocks = _candidate_blocks_for_copy(
                record, owner, source_copy, occurrence_by_id,
            )
            bounded_blocks = _candidate_blocks_for_copy(
                record, owner, bounded_copy, occurrence_by_id,
            )
            signature = (
                tuple((block.start0, block.end0) for block in source_blocks),
                tuple((block.start0, block.end0) for block in bounded_blocks),
            )
            by_placement[signature].append(candidate_id)
        return tuple(
            min(candidate_ids)
            for _signature, candidate_ids in sorted(by_placement.items())
        )

    anchors_by_context = defaultdict(list)
    for candidate_id, record in candidate_record.items():
        owner = candidate_owner[candidate_id]
        if owner is focus_row or not owner.get("_position_edge_eligible"):
            continue
        if candidate_id not in _value_tokens(owner.get("retained_candidate_ids")):
            continue
        query = occurrence_by_id.get(owner.get("query_occurrence_id"), {})
        subject = occurrence_by_id.get(owner.get("subject_occurrence_id"), {})
        query_copy = occurrence_copy_key(query)
        subject_copy = occurrence_copy_key(subject)
        if query_copy == source_copy and subject_copy == bounded_copy:
            anchor_source, anchor_bounded = query, subject
            source_transcript = record.get("query_transcript_id")
            bounded_transcript = record.get("target_transcript_id")
        elif subject_copy == source_copy and query_copy == bounded_copy:
            anchor_source, anchor_bounded = subject, query
            source_transcript = record.get("target_transcript_id")
            bounded_transcript = record.get("query_transcript_id")
        else:
            continue
        if (
            anchor_source.get("occurrence_id") == source.get("occurrence_id")
            or anchor_bounded.get("occurrence_id") == bounded.get("occurrence_id")
        ):
            continue
        if (
            anchor_source.get("contig") != source.get("contig")
            or anchor_source.get("strand") != source.get("strand")
            or anchor_bounded.get("contig") != bounded.get("contig")
            or anchor_bounded.get("strand") != bounded.get("strand")
        ):
            continue
        for membership in _candidate_path_memberships(
            anchor_source,
            anchor_bounded,
            transcript_paths,
            query_transcript_ids=source_transcript,
            subject_transcript_ids=bounded_transcript,
        ):
            anchors_by_context[membership.context].append((candidate_id, membership))

    flank_pairs = set()
    for focal in focal_memberships:
        contextual = anchors_by_context.get(focal.context, ())
        left = [
            (candidate_id, membership)
            for candidate_id, membership in contextual
            if membership.query_order < focal.query_order
            and membership.target_order < focal.target_order
        ]
        right = [
            (candidate_id, membership)
            for candidate_id, membership in contextual
            if membership.query_order > focal.query_order
            and membership.target_order > focal.target_order
        ]
        nearest_left = {
            candidate_id
            for candidate_id, membership in left
            if not any(
                (other.query_order >= membership.query_order)
                and (other.target_order >= membership.target_order)
                and (
                    other.query_order > membership.query_order
                    or other.target_order > membership.target_order
                )
                for _other_id, other in left
            )
        }
        nearest_right = {
            candidate_id
            for candidate_id, membership in right
            if not any(
                (other.query_order <= membership.query_order)
                and (other.target_order <= membership.target_order)
                and (
                    other.query_order < membership.query_order
                    or other.target_order < membership.target_order
                )
                for _other_id, other in right
            )
        }
        if nearest_left and nearest_right:
            flank_pairs.add((
                distinct_placement_ids(nearest_left),
                distinct_placement_ids(nearest_right),
            ))
    return flank_pairs


def _rerun_anchor_bounded_short_candidates(
    rows,
    occurrence_by_id,
    seqs,
    gene_loci,
    transcript_paths=None,
):
    candidate_owner = {}
    candidate_record = {}
    for owner in rows:
        for record in owner.get("_candidate_records", ()):
            candidate_id = record.get("candidate_id")
            if candidate_id not in {None, "", "NA"}:
                candidate_owner[candidate_id] = owner
                candidate_record[candidate_id] = record

    changed = False
    for row in rows:
        if row.get("short_context_route") != "feature_bounded_candidate":
            continue
        row["membership_edge_eligible"] = 0
        row["position_edge_eligible"] = 0
        row["_membership_edge_eligible"] = False
        row["_position_edge_eligible"] = False
        input_transposed = row.get("alignment_input_transposed") in {1, "1", True}
        bounded_side = "query" if input_transposed else "target"
        source_side = "target" if input_transposed else "query"
        bounded_id = (
            row["query_occurrence_id"] if bounded_side == "query"
            else row["subject_occurrence_id"]
        )
        source_id = (
            row["query_occurrence_id"] if source_side == "query"
            else row["subject_occurrence_id"]
        )
        bounded = occurrence_by_id.get(bounded_id, {})
        source = occurrence_by_id.get(source_id, {})
        retained = [
            record for record in row.get("_candidate_records", ())
            if record.get("candidate_id") in _value_tokens(row.get("retained_candidate_ids"))
        ]
        flank_pairs = {
            (
                tuple(sorted(_value_tokens(record.get("left_anchor_id")))),
                tuple(sorted(_value_tokens(record.get("right_anchor_id")))),
            )
            for record in retained
        }
        flank_pairs.discard((tuple(), tuple()))
        inferred_pairs = _context_flank_pairs(
            row,
            source=source,
            bounded=bounded,
            candidate_owner=candidate_owner,
            candidate_record=candidate_record,
            occurrence_by_id=occurrence_by_id,
            transcript_paths=transcript_paths,
        )
        if transcript_paths:
            flank_pairs = inferred_pairs
        if len(flank_pairs) != 1:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "unique_same_context_flank_pair_unavailable"
            continue
        left_ids, right_ids = next(iter(flank_pairs))
        if len(left_ids) != 1 or len(right_ids) != 1:
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "competing_same_context_flank_pairs"
            continue

        bounded_copy = occurrence_copy_key(bounded)
        locus = gene_loci.get((bounded.get("species"), bounded.get("gene_copy_id")))
        if (
            locus is None
            or locus.get("contig") != bounded.get("contig")
            or locus.get("strand") != bounded.get("strand")
        ):
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "target_gene_locus_sequence_unavailable"
            continue

        left_id, right_id = left_ids[0], right_ids[0]
        left_record = candidate_record.get(left_id)
        right_record = candidate_record.get(right_id)
        left_owner = candidate_owner.get(left_id)
        right_owner = candidate_owner.get(right_id)
        if not left_record or not right_record or not left_owner or not right_owner:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flank_candidate_coordinates_unavailable"
            continue
        flank_occurrences = (
            _owner_occurrence_for_copy(
                left_owner, occurrence_copy_key(source), occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                right_owner, occurrence_copy_key(source), occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                left_owner, bounded_copy, occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                right_owner, bounded_copy, occurrence_by_id,
            ),
        )
        expected_geometry = (
            (source.get("contig"), source.get("strand")),
            (source.get("contig"), source.get("strand")),
            (bounded.get("contig"), bounded.get("strand")),
            (bounded.get("contig"), bounded.get("strand")),
        )
        if any(
            occurrence is None
            or (occurrence.get("contig"), occurrence.get("strand")) != expected
            for occurrence, expected in zip(flank_occurrences, expected_geometry)
        ):
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flanks_are_not_on_matching_contigs_and_gene_strands"
            continue
        left_blocks = _candidate_blocks_for_copy(
            left_record, left_owner, bounded_copy, occurrence_by_id,
        )
        right_blocks = _candidate_blocks_for_copy(
            right_record, right_owner, bounded_copy, occurrence_by_id,
        )
        search_interval = _interval_between_transcript_flanks(
            left_blocks, right_blocks, bounded.get("strand"),
        )
        if (
            search_interval is None
            or search_interval.length == 0
            or search_interval.start0 < locus["interval"].start0
            or search_interval.end0 > locus["interval"].end0
        ):
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flanks_do_not_define_one_contained_genome_interval"
            continue
        local_search = genome_interval_to_local(
            search_interval, locus["interval"], locus["strand"],
        )
        target_sequence = locus["sequence"][local_search.start0:local_search.end0]
        query_sequence = seqs.get(source_id, "")
        search_metadata = {
            "coordinate_system": "0-based-half-open",
            "contig": locus["contig"],
            "start0": search_interval.start0,
            "end0": search_interval.end0,
            "strand": locus["strand"],
        }
        try:
            candidate_set = anchored_short_alignment(
                query_sequence,
                target_sequence,
                mode="local",
                query_occurrence_id=source_id,
                target_occurrence_id=bounded_id,
                query_transcript_id=source.get("transcript_id"),
                target_transcript_id=bounded.get("transcript_id"),
                left_anchor_id=left_id,
                right_anchor_id=right_id,
                search_interval=search_metadata,
            )
        except AlignmentBackendError as error:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["short_context_route"] = "anchor_bounded_unavailable"
            row["incomplete_reason"] = str(error)
            continue

        records = []
        try:
            bounded_locus = ClosedInterval1(
                int(bounded["start"]), int(bounded["end"]),
            ).to_interval0()
        except (KeyError, TypeError, ValueError):
            bounded_locus = None
        for rank, candidate in enumerate(candidate_set.candidates, start=1):
            record = _candidate_record(
                candidate, rank, candidate.backend, candidate.score_scheme,
            )
            source_genomic = _genomic_blocks0(
                source, record["aligned_blocks"], "query",
            )
            bounded_genomic = tuple(
                local_interval_to_genome(
                    block.target, search_interval, locus["strand"],
                )
                for block in record["aligned_blocks"]
            )
            try:
                if bounded_locus is None:
                    raise ValueError("bounded parent feature is unavailable")
                bounded_local = tuple(
                    genome_interval_to_local(block, bounded_locus, bounded.get("strand"))
                    for block in bounded_genomic
                )
                parent_gap_blocks = _anchor_bounded_gap_blocks(
                    record.get("gap_blocks", ()),
                    search_interval,
                    bounded_locus,
                    bounded.get("strand"),
                )
            except (KeyError, TypeError, ValueError):
                bounded_local = tuple()
                parent_gap_blocks = []
            if (
                len(source_genomic) == len(record["aligned_blocks"])
                and len(bounded_local) == len(record["aligned_blocks"])
            ):
                row_blocks = tuple(
                    CoordinateBlock(source_block.query, target_local)
                    for source_block, target_local in zip(
                        record["aligned_blocks"], bounded_local,
                    )
                )
                if input_transposed:
                    row_blocks = tuple(
                        CoordinateBlock(block.target, block.query) for block in row_blocks
                    )
                    parent_gap_blocks = _transpose_gap_blocks(parent_gap_blocks)
            else:
                row_blocks = tuple()
                parent_gap_blocks = []
            if input_transposed:
                record = _transpose_candidate_record(record)
                record["query_genomic_blocks0"] = bounded_genomic
                record["target_genomic_blocks0"] = source_genomic
            else:
                record["query_genomic_blocks0"] = source_genomic
                record["target_genomic_blocks0"] = bounded_genomic
            record["aligned_blocks"] = row_blocks
            record["gap_blocks"] = parent_gap_blocks
            if row_blocks:
                record["query_start0"] = min(block.query.start0 for block in row_blocks)
                record["query_end0"] = max(block.query.end0 for block in row_blocks)
                record["target_start0"] = min(block.target.start0 for block in row_blocks)
                record["target_end0"] = max(block.target.end0 for block in row_blocks)
            record["candidate_id"] = (
                f"{row['match_id']}.anchor_bounded_candidate_{rank:03d}"
            )
            record["left_anchor_id"] = left_id
            record["right_anchor_id"] = right_id
            record["search_interval"] = search_metadata
            record["search_interval_side"] = bounded_side
            record["source"] = "nucleotide_alignment"
            record["short_sequence_coverage"] = candidate.query_coverage
            record["accepted"] = int(bool(
                row_blocks
                and _candidate_sequence_accepted(
                    record, to_float(row.get("threshold"), 0.0), short_context=True,
                )
            ))
            record["acceptance_threshold"] = row.get("threshold", "NA")
            records.append(record)

        final_candidate_ids = [record["candidate_id"] for record in records]
        for record in records:
            record["hit_count"] = len(records)
            record["alternative_candidate_ids"] = tuple(
                candidate_id for candidate_id in final_candidate_ids
                if candidate_id != record["candidate_id"]
            )

        row["_candidate_records"] = records
        row["short_context_route"] = "anchor_bounded_local"
        row["flanking_anchor_status"] = "ordered_double_sided_homologous_flanks_same_path"
        row["left_anchor_id"] = left_id
        row["right_anchor_id"] = right_id
        row["search_interval"] = _format_contract_value(_public_interval(search_metadata))
        row["search_interval_side"] = bounded_side
        row["local_boundary_range"] = row["search_interval"]
        row["enumeration_complete"] = int(candidate_set.enumeration_complete)
        row["candidate_enumeration_status"] = (
            "complete" if candidate_set.enumeration_complete else "incomplete"
        )
        row["incomplete_reason"] = candidate_set.incomplete_reason or "NA"
        row["hit_count"] = len(records)
        row["ambiguous_hit_count"] = max(0, len(records) - 1)
        accepted = [record for record in records if record.get("accepted") == 1]
        row["match_status"] = "mapped" if accepted else "candidate_low_similarity"
        row["candidate_resolution"] = "unassessed"
        primary = records[0] if records else None
        if primary is not None:
            row["alignment_score"] = f"{to_float(primary.get('identity'), 0.0):.6g}"
            row["coverage_score"] = f"{to_float(primary.get('coverage'), 0.0):.6g}"
            row["sequence_score"] = f"{(0.70 * to_float(primary.get('identity'), 0.0) + 0.30 * to_float(primary.get('coverage'), 0.0)):.6g}"
            row["correspondence_score"] = row["sequence_score"]
            row["candidate_id"] = primary["candidate_id"]
            row["candidate_ids"] = ";".join(record["candidate_id"] for record in records)
            row["alternative_candidate_ids"] = ";".join(
                record["candidate_id"] for record in records[1:]
            ) or "NA"
            row["matched_blocks"] = _format_alignment_blocks(primary["aligned_blocks"])
            row["aligned_blocks"] = row["matched_blocks"]
            row["projected_reference_blocks"] = row["matched_blocks"]
            query_blocks = tuple(
                _coordinate_block0(block).query for block in primary["aligned_blocks"]
            )
            target_blocks = tuple(
                _coordinate_block0(block).target for block in primary["aligned_blocks"]
            )
            if query_blocks:
                public_query = ClosedInterval1.from_interval0(Interval0(
                    min(block.start0 for block in query_blocks),
                    max(block.end0 for block in query_blocks),
                ))
                public_target = ClosedInterval1.from_interval0(Interval0(
                    min(block.start0 for block in target_blocks),
                    max(block.end0 for block in target_blocks),
                ))
                row["query_alignment_start"] = public_query.start
                row["query_alignment_end"] = public_query.end
                row["target_alignment_start"] = public_target.start
                row["target_alignment_end"] = public_target.end
            row["query_genomic_matched_blocks"] = _format_genomic_blocks(
                occurrence_by_id[row["query_occurrence_id"]].get("contig", "NA"),
                occurrence_by_id[row["query_occurrence_id"]].get("strand", "NA"),
                primary.get("query_genomic_blocks0", ()),
            )
            row["subject_genomic_matched_blocks"] = _format_genomic_blocks(
                occurrence_by_id[row["subject_occurrence_id"]].get("contig", "NA"),
                occurrence_by_id[row["subject_occurrence_id"]].get("strand", "NA"),
                primary.get("target_genomic_blocks0", ()),
            )
            row["gap_blocks"] = json.dumps(
                _public_gap_blocks(primary.get("gap_blocks", ())),
                sort_keys=True, separators=(",", ":"),
            )
            row["alignment_cigar"] = primary.get("cigar", "NA")
            row["raw_alignment_score"] = _format_optional_number(primary.get("score"))
            row["raw_score"] = row["raw_alignment_score"]
            row["short_sequence_coverage"] = primary.get("short_sequence_coverage", "NA")
        public_records = [_public_candidate_record(record) for record in records]
        row["candidate_assessments"] = json.dumps(
            public_records, sort_keys=True, separators=(",", ":"),
        )
        row["dna_candidate_assessments"] = row["candidate_assessments"]
        row["alternative_hits"] = json.dumps(
            public_records[1:], sort_keys=True, separators=(",", ":"),
        )
        changed = True
    return changed


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
        "matched_blocks": f"{query_start}-{query_end}:{target_start}-{target_end}",
        "query_genomic_matched_blocks": f"{left.get('contig', 'NA')}:{overlap_start}-{overlap_end}:{left.get('strand', 'NA')}",
        "subject_genomic_matched_blocks": f"{right.get('contig', 'NA')}:{overlap_start}-{overlap_end}:{right.get('strand', 'NA')}",
        "query_parent_feature_ids": left.get("source_feature_id", "NA"),
        "subject_parent_feature_ids": right.get("source_feature_id", "NA"),
        "query_transcript_ids": left.get("transcript_id", "NA"),
        "subject_transcript_ids": right.get("transcript_id", "NA"),
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
        "mapping_quality": "NA",
        "hit_count": 1,
        "ambiguous_hit_count": 0,
        "alternative_hits": [],
        "score_scheme": "genomic_coordinate_overlap",
        "raw_alignment_score": overlap_len,
        "enumeration_complete": True,
        "candidate_enumeration_status": "complete",
        "incomplete_reason": "NA",
        "short_context_route": "not_used",
        "local_boundary_range": "available",
        "flanking_anchor_status": "same_locus_annotation_overlap",
        "true_absence_eligible": 0,
        "true_absence_evidence_status": "not_applicable",
        "true_absence_reason": "same_locus_annotation_overlap_is_not_deletion_evidence",
        "total_score": total,
    }


def cluster_segments(occurrences, seqs, identity_threshold=0.7, distance_table=None, aligner="mafft", threads=1, min_size_ratio=0.25, species_distances=None, match_writer=None, context_aligner="minimap2", transcript_paths=None, raw_features=None, coding_msa_mode="linsi", short_context_max_length=300, gene_loci=None):
    context = copy_order_context(occurrences, transcript_paths)
    distance_lookup = load_distance_table(distance_table)
    occurrence_by_id = {row["occurrence_id"]: row for row in occurrences}
    matches = []
    pending_match_rows = []
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
        protein_index = CodingProjectionIndex(
            occurrences,
            seqs,
            transcript_paths,
            parsed_features,
            threads=threads,
            msa_mode=coding_msa_mode,
        )

    def emit_match(row):
        pending_match_rows.append(row)

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
                score = evidence["sequence_score"]
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
            evidence = match_evidence(
                left,
                right,
                seqs,
                context,
                aligner=aligner,
                threads=1,
                context_aligner=context_aligner,
                short_context_max_length=short_context_max_length,
            )
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
            score = evidence["sequence_score"]
            compatible = roles_compatible(left, right)
            candidate_records = evidence.get("candidate_records", ())
            short_context = evidence.get("short_context_route") in {
                "feature_bounded_candidate", "anchor_bounded_local",
            }
            for record in candidate_records:
                if record.get("source") not in {None, "", "NA"}:
                    record.setdefault("alignment_source", record["source"])
                record["source"] = "nucleotide_alignment"
                record.setdefault("score_scheme", evidence.get("score_scheme", "unspecified"))
                record["accepted"] = int(
                    _candidate_sequence_accepted(record, threshold, short_context)
                )
                record["acceptance_threshold"] = f"{threshold:.6g}"
            sequence_ok = any(record.get("accepted") == 1 for record in candidate_records)
            if not candidate_records:
                sequence_ok = sequence_supported_mapping(evidence, threshold)
            mapped = compatible and sequence_ok
            evidence["dna_match_status"] = "mapped" if mapped else prefilter_status if prefilter_status != "aligned_candidate" else "low_similarity"
            exon_pair = left.get("role") in EXON_LIKE_ROLES and right.get("role") in EXON_LIKE_ROLES
            # Sensitivity analysis evaluates the original nucleotide evidence,
            # even when a protein projection later supplies the final position.
            evidence["dna_candidate_assessments"] = [
                dict(record) for record in candidate_records
            ]
            if exon_pair and protein_index is not None:
                evidence.update(protein_index.evidence(left["occurrence_id"], right["occurrence_id"]))
                protein_hard = bool(evidence.get("protein_hard_observation_eligible"))
                protein_position = bool(evidence.get("protein_position_eligible"))
                if mapped and protein_hard:
                    evidence["correspondence_basis"] = "DNA_and_annotated_CDS_protein"
                elif compatible and protein_hard:
                    mapped = True
                    score = (
                        0.70 * float(evidence["protein_aa_identity"])
                        + 0.30 * min(
                            float(evidence["protein_query_cds_coverage"]),
                            float(evidence["protein_target_cds_coverage"]),
                        )
                    )
                    evidence["correspondence_basis"] = "annotated_CDS_protein"
                if mapped and protein_position:
                    protein_blocks = tuple(
                        _coordinate_block0(block)
                        for block in evidence.get("protein_projected_blocks", ())
                    )
                    evidence["projected_reference_blocks"] = _format_alignment_blocks(protein_blocks)
                    evidence["matched_blocks"] = evidence["projected_reference_blocks"]
                    if protein_blocks:
                        query_interval = Interval0(
                            min(block.query.start0 for block in protein_blocks),
                            max(block.query.end0 for block in protein_blocks),
                        )
                        target_interval = Interval0(
                            min(block.target.start0 for block in protein_blocks),
                            max(block.target.end0 for block in protein_blocks),
                        )
                        public_query = ClosedInterval1.from_interval0(query_interval)
                        public_target = ClosedInterval1.from_interval0(target_interval)
                        evidence["query_alignment_start"] = public_query.start
                        evidence["query_alignment_end"] = public_query.end
                        evidence["target_alignment_start"] = public_target.start
                        evidence["target_alignment_end"] = public_target.end
                        evidence["query_genomic_matched_blocks"] = _genomic_matched_blocks(left, protein_blocks, "query")
                        evidence["subject_genomic_matched_blocks"] = _genomic_matched_blocks(right, protein_blocks, "subject")
                    evidence["score_scheme"] = "blosum62_cds_projection"
                    evidence["raw_alignment_score"] = evidence.get("protein_blosum62_score", "NA")
                    evidence["alignment_backend"] = "family_protein_msa_projection"
                    evidence["alignment_mode"] = "coding_projection"
                    evidence["alignment_meaning"] = (
                        "family protein MSA projected through transcript codons to CDS bases"
                    )
                    evidence["alignment_strand"] = "+"
                    evidence["candidate_records"] = []
                    evidence["candidate_accepted"] = True
                    evidence["enumeration_complete"] = (
                        evidence.get("protein_candidate_details") not in {None, "", "NA"}
                    )
                    evidence["candidate_enumeration_status"] = (
                        "complete" if evidence["enumeration_complete"] else "incomplete"
                    )
                    evidence["incomplete_reason"] = (
                        "NA" if evidence["enumeration_complete"]
                        else "protein_candidate_set_unavailable"
                    )
                elif evidence.get("protein_candidate_evidence_available"):
                    evidence["incomplete_reason"] = "protein_candidate_without_hard_coordinates"
            elif not mapped and exon_pair:
                if protein_index is None:
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

    # A split transcript can map to the two terminal portions of one short,
    # unsplit coding occurrence.  Each member then lacks enough independent
    # amino-acid columns for the ordinary per-pair anchor threshold, while the
    # complementary left/right terminal partition resolves the mapping jointly.
    terminal_groups = defaultdict(list)
    for row in pending_match_rows:
        if (
            row.get("protein_mapping_status") == "supported_unanchored"
            and row.get("protein_candidate_coordinate_consensus") in {1, "1", True}
            and row.get("protein_terminal_side") in {"left", "right"}
            and row.get("protein_projected_blocks") not in {None, "", "NA"}
        ):
            query = occurrence_by_id.get(row.get("query_occurrence_id"), {})
            subject = occurrence_by_id.get(row.get("subject_occurrence_id"), {})
            key = (
                tuple(sorted((occurrence_copy_key(query), occurrence_copy_key(subject)))),
                row.get("protein_best_query_transcript", "NA"),
                row.get("protein_best_target_transcript", "NA"),
            )
            terminal_groups[key].append(row)
    for group in terminal_groups.values():
        if {row.get("protein_terminal_side") for row in group} != {"left", "right"}:
            continue
        shared_occurrences = set.intersection(*(
            {row.get("query_occurrence_id"), row.get("subject_occurrence_id")}
            for row in group
        ))
        if len(shared_occurrences) != 1:
            continue
        reference_id = next(iter(shared_occurrences))
        reference_intervals = []
        parsed_by_row = {}
        for row in group:
            try:
                blocks = tuple(parse_legacy_blocks(row["protein_projected_blocks"]))
            except (TypeError, ValueError):
                blocks = tuple()
            if not blocks:
                break
            parsed_by_row[id(row)] = blocks
            reference_intervals.extend(
                block.query if row.get("query_occurrence_id") == reference_id else block.target
                for block in blocks
            )
        else:
            ordered = sorted(reference_intervals)
            if any(left.overlaps(right) for left, right in zip(ordered, ordered[1:])):
                continue
            for row in group:
                row["protein_mapping_status"] = "resolved_joint_terminal_partition"
                row["protein_membership_eligible"] = 1
                row["protein_position_eligible"] = 1
                row["protein_hard_observation_eligible"] = 1
                row["match_status"] = "mapped"
                row["correspondence_basis"] = "annotated_CDS_protein"
                blocks = parsed_by_row[id(row)]
                evidence = dict(row, protein_projected_blocks=blocks)
                left_id = row["query_occurrence_id"]
                right_id = row["subject_occurrence_id"]
                projection_by_occ_ref[(left_id, right_id)] = _projection_record(evidence, "target")
                projection_by_occ_ref[(right_id, left_id)] = _projection_record(evidence, "query")
    _apply_ordered_candidate_chains(
        pending_match_rows,
        occurrence_by_id,
        occurrences,
        transcript_paths=transcript_paths,
    )
    if gene_loci and _rerun_anchor_bounded_short_candidates(
        pending_match_rows,
        occurrence_by_id,
        seqs,
        gene_loci,
        transcript_paths=transcript_paths,
    ):
        _apply_ordered_candidate_chains(
            pending_match_rows,
            occurrence_by_id,
            occurrences,
            transcript_paths=transcript_paths,
        )
    position_pairs = {
        frozenset((row["query_occurrence_id"], row["subject_occurrence_id"]))
        for row in pending_match_rows
        if row.get("_position_edge_eligible") or (
            "annotated_CDS_protein" in str(row.get("correspondence_basis", ""))
            and row.get("match_status") == "mapped"
            and row.get("protein_mapping_status") in {
                "resolved_local", "resolved_joint_terminal_partition",
            }
        )
    }
    accepted_edges = []
    score_by_occ = defaultdict(list)
    source_support = defaultdict(lambda: defaultdict(float))
    for row in pending_match_rows:
        protein_hard = (
            "annotated_CDS_protein" in str(row.get("correspondence_basis", ""))
            and row.get("match_status") == "mapped"
            and row.get("protein_mapping_status") in {
                "resolved_local", "resolved_joint_terminal_partition",
            }
        )
        if protein_hard:
            row["match_status"] = "mapped"
            row["membership_edge_eligible"] = 1
            row["membership_edge_reason"] = "resolved_annotated_CDS_protein_membership"
            row["_membership_edge_eligible"] = True
        if not row.get("_membership_edge_eligible"):
            continue
        left_id = row["query_occurrence_id"]
        right_id = row["subject_occurrence_id"]
        edge_score = to_float(row.get("correspondence_score"), 0.0)
        accepted_edges.append((left_id, right_id, edge_score))
        score_by_occ[left_id].append(edge_score)
        score_by_occ[right_id].append(edge_score)
        left_sources = known_source_labels(occurrence_by_id.get(left_id, {}))
        right_sources = known_source_labels(occurrence_by_id.get(right_id, {}))
        if left_sources and not right_sources:
            for source in left_sources:
                source_support[right_id][source] += edge_score
        if right_sources and not left_sources:
            for source in right_sources:
                source_support[left_id][source] += edge_score
    projection_by_occ_ref = {
        key: value
        for key, value in projection_by_occ_ref.items()
        if frozenset(key) in position_pairs
    }
    for row in pending_match_rows:
        public_row = {key: value for key, value in row.items() if not key.startswith("_")}
        if match_writer is not None:
            match_writer(public_row)
        if public_row.get("match_status") == "mapped":
            matches.append(public_row)
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
                    "support_type": "ordered_sequence_correspondence_graph",
                    "confidence": f"{confidence:.6g}",
                    "source_label": inferred_source_label(occurrence_by_id.get(occ_id, {}), source_support.get(occ_id, {})),
                }
            )
    return homology, matches


def derive_tables(input_dir, output_dir=None, identity_threshold=0.7, distance_table=None, aligner="mafft", threads=1, min_size_ratio=0.25, context_aligner="minimap2", coding_msa_mode="linsi", short_context_max_length=300):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir or input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", SEGMENT_FIELDS)
    transcript_paths = read_tsv(input_dir / "transcript_paths.tsv", optional=True)
    raw_features = read_tsv(input_dir / "raw_gene_features.tsv", optional=True)
    seqs = parse_fasta(input_dir / "segment_sequences.fasta")
    gene_loci = _gene_locus_records(input_dir)
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
            coding_msa_mode=coding_msa_mode,
            short_context_max_length=short_context_max_length,
            gene_loci=gene_loci,
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
                "short_context_max_length": int(short_context_max_length),
                "coding_msa_mode": coding_msa_mode,
                "notes": (
                    f"{row['notes']}; min_size_ratio is restricted to non-exon-like prefiltering"
                    if selected_exon or selected_context
                    else row["notes"]
                ),
            }
        )
    write_tsv(
        output_dir / "alignment_backend_report.tsv",
        backend_rows,
        [
            "aligner",
            "available",
            "selected",
            "threads",
            "min_size_ratio",
            "short_context_max_length",
            "coding_msa_mode",
            "notes",
            "selected_exon",
            "selected_context",
            "alignment_mode",
        ],
    )
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
