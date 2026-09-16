"""Sequence-supported annotation completion."""

from collections import Counter, defaultdict
from pathlib import Path

from .alignment import AlignmentBackendError, local_alignment_stats
from .io import parse_fasta, read_tsv, to_float, write_tsv


def support_score(row):
    identity = to_float(row.get("sequence_score"))
    coverage = to_float(row.get("sequence_coverage"), 1.0)
    return min(identity, coverage)


def completion_call(row, score, threshold):
    status = row.get("evidence_status", "ambiguous")
    event = row.get("inferred_event", "")
    frame = row.get("frame_status", "")
    if status == "supports_hidden_segment" and score >= threshold:
        if event == "shifted_splice_site":
            return "shifted_splice_site_candidate"
        if event == "intron_deletion_joined_exon":
            return "joined_exon_candidate"
        if frame == "frameshift_or_stop_risk":
            return "hidden_segment_with_frame_disruption"
        return "hidden_segment_candidate"
    if status == "conflicts_annotation" and score >= threshold:
        return "annotation_conflict_candidate"
    if status == "supports_annotation":
        return "supports_annotation"
    if status == "supports_absence" and score >= threshold:
        return "supports_true_absence"
    return "ambiguous_evidence"


def _locus_key(header):
    parts = header.split("|", 2)
    return tuple(parts[:2]) if len(parts) >= 2 else ("NA", header)


def _overlaps_annotated_exon(start, end, locus_header, occurrences):
    locus = locus_header.split("|", 2)[2]
    contig, interval = locus.rsplit(":", 1)[0].rsplit(":", 1)
    lower, upper = (int(value) for value in interval.split("-", 1))
    if locus.endswith(":-"):
        hit_start, hit_end = upper - end + 1, upper - start + 1
    else:
        hit_start, hit_end = lower + start - 1, lower + end - 1
    return any(
        row.get("role") in {"exon", "CDS", "UTR", "noncoding_exon"}
        and row.get("contig") == contig
        and int(row["start"]) <= hit_end
        and int(row["end"]) >= hit_start
        for row in occurrences
    )


def generate_sequence_evidence(input_dir, result_dir, min_identity=0.70, min_coverage=0.60):
    """Search missing homologous exon sequences inside supplied homologous gene loci."""
    input_dir = Path(input_dir)
    result_dir = Path(result_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv")
    elements = read_tsv(result_dir / "element_correspondence.tsv", optional=True)
    sequences = parse_fasta(input_dir / "segment_sequences.fasta")
    locus_records = parse_fasta(input_dir / "gene_loci.fasta")
    loci = {_locus_key(name): (name, sequence) for name, sequence in locus_records.items()}
    if not elements or not loci:
        return []

    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    element_by_occ = {row["occurrence_id"]: row["element_id"] for row in elements}
    rows_by_element = defaultdict(list)
    copies_by_family_species = defaultdict(set)
    occurrences_by_copy = defaultdict(list)
    present_elements = defaultdict(set)
    order_by_copy = defaultdict(list)
    for occurrence in occurrences:
        key = (occurrence["family_id"], occurrence["species"])
        copies_by_family_species[key].add(occurrence["gene_copy_id"])
        occurrences_by_copy[(occurrence["family_id"], occurrence["species"], occurrence["gene_copy_id"])].append(occurrence)
        element = element_by_occ.get(occurrence["occurrence_id"])
        if element and occurrence.get("role") != "intron":
            present_elements[(occurrence["family_id"], occurrence["species"], occurrence["gene_copy_id"])].add(element)
            order_by_copy[(occurrence["family_id"], occurrence["species"], occurrence["gene_copy_id"])].append(
                (int(occurrence.get("transcript_order", 0) or 0), element)
            )
    for row in elements:
        occurrence = occ_by_id.get(row.get("occurrence_id", ""))
        if occurrence and row.get("element_class") == "exon_like":
            rows_by_element[(occurrence["family_id"], row["element_id"])].append((row, occurrence))

    evidence = []
    for (family, element), members in sorted(rows_by_element.items()):
        representative_row, representative = max(
            members,
            key=lambda pair: len(sequences.get(pair[1]["occurrence_id"], "")),
        )
        query = sequences.get(representative["occurrence_id"], "")
        if not query:
            continue
        source_order = [element_id for _rank, element_id in sorted(order_by_copy[(family, representative["species"], representative["gene_copy_id"])])]
        source_index = source_order.index(element) if element in source_order else -1
        left_element = source_order[source_index - 1] if source_index > 0 else None
        right_element = source_order[source_index + 1] if 0 <= source_index < len(source_order) - 1 else None
        for (candidate_family, species), copies in sorted(copies_by_family_species.items()):
            if candidate_family != family or len(copies) != 1:
                continue
            gene_copy = next(iter(copies))
            copy_key = (family, species, gene_copy)
            if element in present_elements[copy_key]:
                continue
            locus_header, locus = loci.get((species, gene_copy), ("", ""))
            if not locus:
                continue
            backend = "internal" if len(query) * len(locus) <= 250_000 else "minimap2"
            try:
                alignment = local_alignment_stats(query, locus, backend=backend)
                identity = alignment.identity
                coverage = alignment.query_coverage or alignment.coverage
            except AlignmentBackendError:
                alignment = None
                identity = coverage = 0.0
            anchored = bool(
                left_element
                and right_element
                and left_element in present_elements[copy_key]
                and right_element in present_elements[copy_key]
            )
            assembly_complete = locus.count("N") / max(1, len(locus)) <= 0.05
            supported = alignment is not None and identity >= min_identity and coverage >= min_coverage
            overlaps_exon = supported and _overlaps_annotated_exon(
                alignment.target_start, alignment.target_end,
                locus_header, occurrences_by_copy[copy_key],
            )
            if overlaps_exon:
                status = "ambiguous"
                annotation_status = "alignment_overlaps_annotated_exon"
                start, end = alignment.target_start, alignment.target_end
            elif supported:
                status = "supports_hidden_segment"
                annotation_status = "unannotated_homologous_sequence"
                start, end = alignment.target_start, alignment.target_end
            elif anchored and assembly_complete:
                status = "supports_absence"
                annotation_status = "sequence_absence_between_flanking_homologs"
                start = end = "NA"
            else:
                status = "ambiguous"
                annotation_status = "insufficient_sequence_or_flank_evidence"
                start = end = "NA"
            evidence.append(
                {
                    "evidence_id": f"completion_{family}_{element}_{species}",
                    "family_id": family,
                    "species": species,
                    "gene_copy_id": gene_copy,
                    "homology_id": representative_row.get("homology_id", "NA"),
                    "annotation_status": annotation_status,
                    "evidence_status": status,
                    "inferred_role": "exon",
                    "contig": "supplied_gene_locus",
                    "start": start,
                    "end": end,
                    "strand": "+",
                    "sequence_score": f"{identity:.6g}",
                    "sequence_coverage": f"{coverage:.6g}",
                    "left_synteny_score": "1" if left_element in present_elements[copy_key] else "0",
                    "right_synteny_score": "1" if right_element in present_elements[copy_key] else "0",
                    "splice_motif_score": "NA",
                    "phase_compatibility": "unknown",
                    "inferred_event": "homologous_exon_sequence_search",
                    "frame_status": "unknown",
                }
            )
    fields = [
        "evidence_id", "family_id", "species", "gene_copy_id", "homology_id",
        "annotation_status", "evidence_status", "inferred_role", "contig", "start", "end",
        "strand", "sequence_score", "sequence_coverage", "left_synteny_score",
        "right_synteny_score", "splice_motif_score", "phase_compatibility",
        "inferred_event", "frame_status",
    ]
    write_tsv(result_dir / "sequence_synteny_evidence.tsv", evidence, fields)
    return evidence


def complete_annotation(input_dir, output_dir, threshold=0.55):
    generated = Path(output_dir) / "sequence_synteny_evidence.tsv"
    evidence_path = generated if generated.exists() else Path(input_dir) / "sequence_synteny_evidence.tsv"
    evidence = read_tsv(
        evidence_path,
        ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "annotation_status", "evidence_status"],
        optional=True,
    )
    rows = []
    summary = Counter()
    for row in evidence:
        score = support_score(row)
        status = row.get("evidence_status", "ambiguous")
        call = completion_call(row, score, threshold)
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
                "inferred_event": row.get("inferred_event", "NA"),
                "frame_status": row.get("frame_status", "NA"),
            }
        )
    write_tsv(
        f"{output_dir}/annotation_completion_candidates.tsv",
        rows,
        ["evidence_id", "family_id", "species", "gene_copy_id", "homology_id", "interval", "annotation_status", "inferred_role", "support_score", "completion_call", "evidence_status", "inferred_event", "frame_status"],
    )
    write_tsv(
        f"{output_dir}/annotation_completion_summary.tsv",
        [{"completion_call": key, "count": value} for key, value in sorted(summary.items())],
        ["completion_call", "count"],
    )
    return rows
