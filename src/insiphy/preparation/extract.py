"""preparation / extract: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from insiphy.aligners.scoring import phase_compatibility
from insiphy.mapping.fields import GENE_LOCUS_FIELDS
from insiphy.mapping.fields import INTRON_SITE_FIELDS
from insiphy.mapping.fields import RAW_FEATURE_FIELDS
from insiphy.mapping.fields import SEGMENT_OUTPUT_FIELDS
from insiphy.mapping.fields import TRANSCRIPT_PATH_FIELDS
from insiphy.preparation.annotation_index import read_annotation_for_gene
from insiphy.preparation.transcripts import _feature_key
from insiphy.preparation.transcripts import _feature_parents
from insiphy.preparation.transcripts import _format_attrs
from insiphy.preparation.transcripts import _format_intervals
from insiphy.preparation.transcripts import _interval_union_length
from insiphy.preparation.transcripts import _occurrence_id
from insiphy.preparation.transcripts import _path_role_record
from insiphy.preparation.transcripts import child_features_for_transcript
from insiphy.preparation.transcripts import feature_role
from insiphy.preparation.transcripts import introns_from_path
from insiphy.preparation.transcripts import load_existing_segments
from insiphy.preparation.transcripts import select_transcripts
from insiphy.preparation.transcripts import sequence_slice
from insiphy.preparation.transcripts import transcript_features
from insiphy.preparation.transcripts import transcript_sort_key
from insiphy.preparation.transcripts import translate_cds
from insiphy.storage.fasta import fasta_record_length
from insiphy.storage.fasta import read_fasta_interval
from insiphy.storage.tabular import read_tsv
from insiphy.storage.tabular import write_tsv
from insiphy.storage.values import to_float
from pathlib import Path
import json


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
