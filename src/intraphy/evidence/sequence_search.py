"""evidence / sequence search: explicit implementation ownership."""
from __future__ import annotations

from collections import Counter
from collections import defaultdict
from intraphy.elements import EXON_LIKE_ROLES
from intraphy.evidence.alignment_provenance import _alignment_evidence_provenance
from intraphy.evidence.alignment_provenance import _ambiguous_repeated_mapping
from intraphy.evidence.alignment_provenance import _format_alternative_hits
from intraphy.evidence.alignment_provenance import _nucleotide_interval_candidates
from intraphy.evidence.alignment_provenance import _projection_within_anchor_interval
from intraphy.evidence.classification import _classify_sequence_evidence
from intraphy.evidence.deletion_support import _ordered_anchor_deletion_support
from intraphy.evidence.fields import _EVIDENCE_PROVENANCE_FIELDS
from intraphy.evidence.fields import _PROTEIN_QUERY_INTERVAL_FIELDS
from intraphy.evidence.nucleotide_search import _align_nucleotide_evidence
from intraphy.evidence.projection import _alignment_genome_blocks
from intraphy.evidence.projection import _annotation_overlaps_for_blocks
from intraphy.evidence.projection import _candidate_span_blocks
from intraphy.evidence.projection import _json_records
from intraphy.evidence.projection import _locus_key
from intraphy.evidence.projection import _protein_projection_blocks
from intraphy.evidence.projection import _supplied_annotation_role
from intraphy.evidence.projection import _transcripts_for_occurrence
from intraphy.evidence.protein import _cached_protein_projection
from intraphy.evidence.protein import _protein_context_provenance
from intraphy.evidence.reference_queries import _prepare_reference_queries
from intraphy.evidence.search_context import _extended_locus
from intraphy.evidence.search_context import _ordered_anchor_search_context
from intraphy.evidence.search_context import _read_locus_metadata
from intraphy.evidence.search_context import _search_provenance
from intraphy.evidence.search_context import _source_occurrence_path
from intraphy.storage.fasta import parse_fasta
from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from pathlib import Path
import json


def generate_sequence_evidence(
    input_dir,
    result_dir,
    min_identity=0.70,
    min_coverage=0.60,
    threads=1,
    aligner="internal",
    short_context_max_length=300,
):
    """Search missing homologous exon sequences inside supplied homologous gene loci."""
    input_dir = Path(input_dir)
    result_dir = Path(result_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv")
    elements = read_tsv(result_dir / "element_correspondence.tsv", optional=True)
    sequences = parse_fasta(input_dir / "segment_sequences.fasta")
    locus_records = parse_fasta(input_dir / "gene_loci.fasta")
    locus_metadata = _read_locus_metadata(input_dir)
    transcript_paths = read_tsv(input_dir / "transcript_paths.tsv", optional=True)
    transcript_ids_by_occurrence = defaultdict(set)
    for path_row in transcript_paths:
        if path_row.get("occurrence_id") and path_row.get("transcript_id"):
            transcript_ids_by_occurrence[path_row["occurrence_id"]].add(path_row["transcript_id"])
    proteins = parse_fasta(input_dir / "protein_sequences.fasta") if aligner == "miniprot" else {}
    loci = {_locus_key(name): (name, sequence) for name, sequence in locus_records.items()}
    if not elements or not loci:
        return []

    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    element_by_occ = {row["occurrence_id"]: row["element_id"] for row in elements}
    rows_by_element = defaultdict(list)
    copies_by_family_species = defaultdict(set)
    occurrences_by_copy = defaultdict(list)
    present_elements = defaultdict(set)
    core_exonic_occurrences = set()
    annotated_exon_counts = Counter()
    for occurrence in occurrences:
        key = (occurrence["family_id"], occurrence["species"])
        copies_by_family_species[key].add(occurrence["gene_copy_id"])
        occurrences_by_copy[(occurrence["family_id"], occurrence["species"], occurrence["gene_copy_id"])].append(occurrence)
        element = element_by_occ.get(occurrence["occurrence_id"])
        if element and occurrence.get("role") in EXON_LIKE_ROLES:
            present_elements[(occurrence["family_id"], occurrence["species"], occurrence["gene_copy_id"])].add(element)
    for row in elements:
        occurrence = occ_by_id.get(row.get("occurrence_id", ""))
        if (
            occurrence
            and row.get("element_class") == "exon_like"
            and occurrence.get("role") in EXON_LIKE_ROLES
        ):
            annotated_exon_counts[(occurrence["family_id"], row["element_id"])] += 1
    for row in elements:
        occurrence = occ_by_id.get(row.get("occurrence_id", ""))
        if (
            occurrence
            and row.get("element_class") == "exon_like"
            and occurrence.get("role") in EXON_LIKE_ROLES
        ):
            key = (occurrence["family_id"], row["element_id"])
            if row.get("membership_call", "core_member") == "core_member" or annotated_exon_counts[key] == 1:
                rows_by_element[key].append((row, occurrence))
            if row.get("membership_call", "core_member") == "core_member":
                core_exonic_occurrences.add(row.get("occurrence_id"))

    representative_infos, protein_contexts_by_family = _prepare_reference_queries(rows_by_element, sequences, aligner, transcript_paths, proteins, occ_by_id)

    protein_projection_cache = {}
    evidence = []
    for info in representative_infos:
        family = info["family"]
        element = info["element"]
        representative_row = info["representative_row"]
        representative = info["representative"]
        query = info["query"]
        protein_context = info.get("protein_context")
        protein_context_provenance = _protein_context_provenance(protein_context)
        source_transcript_ids = _transcripts_for_occurrence(
            representative, transcript_ids_by_occurrence
        )
        source_occurrences, source_path_status = _source_occurrence_path(
            family,
            representative,
            transcript_paths,
            occurrences_by_copy,
            occ_by_id,
            element_by_occ,
            core_exonic_occurrences,
        )
        source_elements = [element_by_occ.get(row.get("occurrence_id")) for row in source_occurrences]
        source_index = next(
            (
                index for index, row in enumerate(source_occurrences)
                if row.get("occurrence_id") == representative.get("occurrence_id")
            ),
            -1,
        )
        left_element = source_elements[source_index - 1] if source_index > 0 else None
        right_element = source_elements[source_index + 1] if 0 <= source_index < len(source_elements) - 1 else None
        terminal_need = None
        if source_elements and source_index == 0:
            terminal_need = "left"
        elif source_index >= 0 and source_index == len(source_elements) - 1:
            terminal_need = "right"
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
            metadata = locus_metadata.get((species, gene_copy), {})
            locus_header, locus, expansion = _extended_locus(
                input_dir, locus_header, locus, metadata, terminal_need
            )
            anchor_context = _ordered_anchor_search_context(
                source_occurrences,
                occurrences_by_copy[copy_key],
                element_by_occ,
                representative,
                left_element,
                right_element,
                locus_header,
                locus,
            )
            supported, alignment, dna_identity, dna_coverage, alignment_evidence_scope, nucleotide_aligner, alignment_error, bounded_target = _align_nucleotide_evidence(aligner, anchor_context, query, short_context_max_length, threads, representative, source_transcript_ids, locus, min_identity, min_coverage)
            identity, coverage = dna_identity, dna_coverage
            projected = None
            projected_candidate = None
            protein_candidates = []
            projection_reason = "not_requested"
            protein_identity = protein_coverage = 0.0
            if aligner == "miniprot":
                projected_candidate, projection_reason = _cached_protein_projection(
                    protein_context,
                    family,
                    copy_key,
                    locus_header,
                    locus,
                    protein_contexts_by_family,
                    protein_projection_cache,
                    min_identity,
                    min_coverage,
                    threads,
                    candidates=protein_candidates,
                )
                projected = (
                    projected_candidate
                    if _projection_within_anchor_interval(projected_candidate, anchor_context)
                    else None
                )
                if projected:
                    protein_identity = projected["identity"]
                    protein_coverage = projected["coverage"]
                    identity = protein_identity
                    coverage = protein_coverage
            if not source_occurrences:
                deletion_supported = False
                deletion_reason = source_path_status
                deletion_provenance = {"backend": "NA", "cigar": "NA"}
            else:
                deletion_supported, deletion_reason, deletion_provenance = _ordered_anchor_deletion_support(
                    source_occurrences,
                    occurrences_by_copy[copy_key],
                    element_by_occ,
                    representative,
                    left_element,
                    right_element,
                    loci.get((representative["species"], representative["gene_copy_id"]), ("", ""))[0],
                    loci.get((representative["species"], representative["gene_copy_id"]), ("", ""))[1],
                    locus_header,
                    locus,
                    min_identity,
                    min_coverage,
                    threads=threads,
                )
            bounded_protein_candidates = [
                candidate for candidate in protein_candidates
                if _projection_within_anchor_interval(candidate, anchor_context)
            ]
            if projected is None and len(bounded_protein_candidates) == 1:
                projected = bounded_protein_candidates[0]
                projection_reason = "protein_projection_resolved_by_ordered_anchor_interval"
                protein_identity = projected["identity"]
                protein_coverage = projected["coverage"]
                identity = protein_identity
                coverage = protein_coverage
            elif projected_candidate is not None and projected is None:
                projection_reason = "protein_projection_outside_ordered_anchor_interval"
            ambiguous_dna = supported and _ambiguous_repeated_mapping(alignment)
            ambiguous_protein = len(bounded_protein_candidates) > 1
            interval_candidates = _nucleotide_interval_candidates(alignment, locus_header) + protein_candidates
            dna_blocks = _alignment_genome_blocks(
                alignment, locus_header, representative, source_transcript_ids
            )
            dna_annotation_overlaps = _annotation_overlaps_for_blocks(
                dna_blocks,
                occurrences_by_copy[copy_key],
                transcript_ids_by_occurrence,
            )
            protein_blocks = _protein_projection_blocks(
                bounded_protein_candidates, representative, source_transcript_ids
            )
            projection_annotation_overlaps = _annotation_overlaps_for_blocks(
                protein_blocks,
                occurrences_by_copy[copy_key],
                transcript_ids_by_occurrence,
            )
            status, annotation_status, inferred_role, predicted_role, conflict_blocks, contig, start, end, strand, primary_score, primary_coverage, conclusion, confidence, primary_backend, primary_mapping_status, interval_scope, primary_cigar, correspondence_status = _classify_sequence_evidence(anchor_context, alignment, dna_identity, min_identity, dna_coverage, min_coverage, protein_candidates, alignment_evidence_scope, nucleotide_aligner, metadata, locus_header, ambiguous_dna, ambiguous_protein, bounded_protein_candidates, projected, protein_identity, protein_coverage, projection_annotation_overlaps, identity, coverage, supported, dna_annotation_overlaps, deletion_supported, deletion_provenance, alignment_error, aligner, projection_reason, deletion_reason)

            if protein_blocks and (projected or ambiguous_protein):
                observation_blocks = protein_blocks
                annotation_overlaps = projection_annotation_overlaps
            elif ambiguous_dna:
                observation_blocks = _candidate_span_blocks(
                    interval_candidates,
                    representative,
                    source_transcript_ids,
                    backend=alignment.backend,
                )
                annotation_overlaps = _annotation_overlaps_for_blocks(
                    observation_blocks,
                    occurrences_by_copy[copy_key],
                    transcript_ids_by_occurrence,
                )
            elif dna_blocks:
                observation_blocks = dna_blocks
                annotation_overlaps = dna_annotation_overlaps
            else:
                observation_blocks = []
                annotation_overlaps = []

            supplied_role = (
                "ambiguous_candidate_roles"
                if ambiguous_dna or ambiguous_protein
                else _supplied_annotation_role(annotation_overlaps)
            )
            supplied_roles = sorted({row["role"] for row in annotation_overlaps})
            target_parent_occurrence_ids = sorted(
                {row["occurrence_id"] for row in annotation_overlaps}
            )
            target_parent_transcript_ids = sorted(
                {
                    transcript_id
                    for row in annotation_overlaps
                    for transcript_id in row.get("transcript_ids", [])
                }
            )

            if status == "supports_absence":
                homologous_dna_presence = "absent"
                homologous_dna_evidence = "ordered_flank_deletion"
            elif supported or projected or bounded_protein_candidates:
                homologous_dna_presence = "present"
                dna_evidence_types = []
                if supported:
                    dna_evidence_types.append("anchor_bounded_nucleotide_alignment")
                if projected or bounded_protein_candidates:
                    dna_evidence_types.append("anchor_bounded_protein_coding_projection")
                homologous_dna_evidence = ";".join(dna_evidence_types)
            else:
                homologous_dna_presence = "unknown"
                homologous_dna_evidence = "no_resolved_sequence_evidence"

            search_provenance = _search_provenance(
                metadata, locus_header, terminal_need, interval_candidates
            )
            incomplete_reasons = []
            if alignment is not None and not alignment.enumeration_complete:
                incomplete_reasons.append(
                    alignment.incomplete_reason or "nucleotide_candidate_enumeration_truncated"
                )
            if search_provenance["hit_search_limit_status"].startswith("candidate_touches_"):
                incomplete_reasons.append(
                    search_provenance["hit_search_limit_status"]
                )
            if search_provenance["hit_search_limit_status"].startswith("no_hit_with_"):
                incomplete_reasons.append(
                    search_provenance["hit_search_limit_status"]
                )
            if incomplete_reasons:
                candidate_resolution_status = "candidate_search_incomplete"
                candidate_search_complete = "false"
                correspondence_status = "unknown"
            elif ambiguous_dna or ambiguous_protein:
                candidate_resolution_status = "ambiguous"
                candidate_search_complete = (
                    "true"
                    if not bounded_protein_candidates
                    and alignment is not None
                    and alignment.enumeration_complete
                    else "unknown"
                )
            elif correspondence_status == "resolved" and homologous_dna_presence in {"present", "absent"}:
                candidate_resolution_status = "resolved"
                candidate_search_complete = (
                    "true"
                    if supported
                    and not bounded_protein_candidates
                    and alignment is not None
                    and alignment.enumeration_complete
                    else "unknown"
                )
            else:
                candidate_resolution_status = "unresolved"
                candidate_search_complete = "unknown"

            evidence.append(
                {
                    "evidence_id": f"completion_{family}_{element}_{species}",
                    "family_id": family,
                    "species": species,
                    "gene_copy_id": gene_copy,
                    "homology_id": representative_row.get("homology_id", "NA"),
                    "annotation_status": annotation_status,
                    "evidence_status": status,
                    "inferred_role": inferred_role,
                    "predicted_role": predicted_role,
                    "homologous_dna_presence": homologous_dna_presence,
                    "homologous_dna_evidence": homologous_dna_evidence,
                    "predicted_exonic_role": predicted_role,
                    "supplied_annotation_role": supplied_role,
                    "supplied_annotation_roles": ";".join(supplied_roles) or "NA",
                    "source_parent_occurrence_id": representative.get("occurrence_id", "NA"),
                    "source_parent_transcript_ids": ";".join(source_transcript_ids) or "NA",
                    "target_parent_occurrence_ids": ";".join(target_parent_occurrence_ids) or "NA",
                    "target_parent_transcript_ids": ";".join(target_parent_transcript_ids) or "NA",
                    "dna_aligned_blocks": _json_records(dna_blocks),
                    "predicted_role_blocks": _json_records(protein_blocks if predicted_role in EXON_LIKE_ROLES else []),
                    "supplied_annotation_overlaps": _json_records(annotation_overlaps),
                    "annotation_conflict_blocks": _json_records(conflict_blocks),
                    "candidate_resolution_status": candidate_resolution_status,
                    "candidate_search_complete": candidate_search_complete,
                    "candidate_search_incomplete_reason": ";".join(sorted(set(incomplete_reasons))) or "NA",
                    "left_anchor_id": anchor_context["left_anchor_id"],
                    "right_anchor_id": anchor_context["right_anchor_id"],
                    "anchor_interval_status": anchor_context["status"],
                    "search_interval": (
                        json.dumps(anchor_context["search_interval"], separators=(",", ":"))
                        if anchor_context["search_interval"] is not None
                        else "NA"
                    ),
                    "alignment_evidence_scope": alignment_evidence_scope,
                    **_alignment_evidence_provenance(
                        alignment,
                        query_length=len(query),
                        target_length=(
                            len(bounded_target)
                            if alignment_evidence_scope.startswith("anchor_bounded_")
                            else len(locus)
                        ),
                    ),
                    **search_provenance,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "sequence_score": "NA" if primary_score is None else f"{primary_score:.6g}",
                    "sequence_coverage": "NA" if primary_coverage is None else f"{primary_coverage:.6g}",
                    "left_synteny_score": "1" if left_element in present_elements[copy_key] else "0",
                    "right_synteny_score": "1" if right_element in present_elements[copy_key] else "0",
                    "splice_motif_score": "NA",
                    "phase_compatibility": "unknown",
                    "inferred_event": "protein_cds_projection" if projected else "homologous_exon_sequence_search",
                    "frame_status": "unknown",
                    "evidence_conclusion": conclusion,
                    "confidence_flag": confidence,
                    "absence_evidence": deletion_reason,
                    "alignment_backend": primary_backend,
                    "primary_mapping_status": primary_mapping_status,
                    "correspondence_status": correspondence_status,
                    "interval_scope": interval_scope,
                    "interval_candidates": json.dumps(interval_candidates, separators=(",", ":")) if interval_candidates else "NA",
                    "evidence_aligner": aligner,
                    "alignment_error": alignment_error or "NA",
                    "alignment_cigar": primary_cigar,
                    "dna_sequence_score": f"{dna_identity:.6g}" if alignment is not None else "NA",
                    "dna_sequence_coverage": f"{dna_coverage:.6g}" if alignment is not None else "NA",
                    "dna_alignment_backend": alignment.backend if alignment is not None else nucleotide_aligner,
                    "dna_alignment_cigar": alignment.cigar if alignment is not None else "NA",
                    "protein_sequence_score": f"{protein_identity:.6g}" if projected else "NA",
                    "protein_sequence_coverage": f"{protein_coverage:.6g}" if projected else "NA",
                    "protein_alignment_backend": "miniprot" if projected else "NA",
                    "protein_alignment_cigar": projected.get("cigar", "NA") if projected else "NA",
                    "protein_projection_parent_id": projected.get("parent_id", "NA") if projected else "NA",
                    "protein_projection_reason": projection_reason,
                    **{field: projected.get(field, "NA") if projected else "NA" for field in _PROTEIN_QUERY_INTERVAL_FIELDS},
                    "protein_identity_scope": projected.get("identity_scope", "unknown") if projected else "NA",
                    "protein_cigar_scope": "parent_alignment" if projected and projected.get("cigar", "NA") != "NA" else "NA",
                    "protein_cds_phase": projected.get("phase", "NA") if projected else "NA",
                    **protein_context_provenance,
                    "alignment_hit_count": len(protein_candidates) if primary_backend == "miniprot" else alignment.hit_count if alignment is not None and primary_mapping_status != "ordered_flank_deletion" else "NA",
                    "alignment_ambiguous_hit_count": max(0, len(protein_candidates) - 1) if primary_backend == "miniprot" else alignment.ambiguous_hit_count if alignment is not None and primary_mapping_status != "ordered_flank_deletion" else "NA",
                    "alignment_mapping_quality": alignment.mapping_quality if alignment is not None and alignment.mapping_quality is not None and primary_backend != "miniprot" and primary_mapping_status != "ordered_flank_deletion" else "NA",
                    "alignment_alternative_hits": _format_alternative_hits(alignment) if primary_backend != "miniprot" and primary_mapping_status != "ordered_flank_deletion" else "NA",
                    "dna_alignment_hit_count": alignment.hit_count if alignment is not None else 0,
                    "dna_alignment_ambiguous_hit_count": alignment.ambiguous_hit_count if alignment is not None else 0,
                    "dna_alignment_mapping_quality": alignment.mapping_quality if alignment is not None and alignment.mapping_quality is not None else "NA",
                    "dna_alignment_alternative_hits": _format_alternative_hits(alignment),
                    "deletion_alignment_backend": deletion_provenance.get("backend", "NA"),
                    "deletion_alignment_cigar": deletion_provenance.get("cigar", "NA"),
                    "left_flank_identity": deletion_provenance.get("left_flank_identity", "NA"),
                    "left_flank_coverage": deletion_provenance.get("left_flank_coverage", "NA"),
                    "left_flank_paired_bases": deletion_provenance.get("left_flank_paired_bases", "NA"),
                    "left_flank_matches": deletion_provenance.get("left_flank_matches", "NA"),
                    "left_flank_mismatches": deletion_provenance.get("left_flank_mismatches", "NA"),
                    "right_flank_identity": deletion_provenance.get("right_flank_identity", "NA"),
                    "right_flank_coverage": deletion_provenance.get("right_flank_coverage", "NA"),
                    "right_flank_paired_bases": deletion_provenance.get("right_flank_paired_bases", "NA"),
                    "right_flank_matches": deletion_provenance.get("right_flank_matches", "NA"),
                    "right_flank_mismatches": deletion_provenance.get("right_flank_mismatches", "NA"),
                    "search_expansion_status": expansion.get("search_expansion_status", "not_requested"),
                    "search_original_start": expansion.get("search_original_start", metadata.get("search_start", "NA")),
                    "search_original_end": expansion.get("search_original_end", metadata.get("search_end", "NA")),
                    "search_expanded_start": expansion.get("search_expanded_start", "NA"),
                    "search_expanded_end": expansion.get("search_expanded_end", "NA"),
                    "range_status": expansion.get("range_status", metadata.get("range_status", "NA")),
                }
            )
    fields = [
        "evidence_id", "family_id", "species", "gene_copy_id", "homology_id",
        "annotation_status", "evidence_status", "inferred_role", "predicted_role", "contig", "start", "end",
        *_EVIDENCE_PROVENANCE_FIELDS,
        "strand", "sequence_score", "sequence_coverage", "left_synteny_score",
        "right_synteny_score", "splice_motif_score", "phase_compatibility",
        "inferred_event", "frame_status", "evidence_conclusion", "confidence_flag",
        "absence_evidence", "alignment_backend", "primary_mapping_status", "evidence_aligner", "alignment_error", "alignment_cigar",
        "correspondence_status", "interval_scope", "interval_candidates",
        "dna_sequence_score", "dna_sequence_coverage", "dna_alignment_backend", "dna_alignment_cigar",
        "protein_sequence_score", "protein_sequence_coverage", "protein_alignment_backend", "protein_alignment_cigar",
        "protein_projection_parent_id", "protein_projection_reason",
        *_PROTEIN_QUERY_INTERVAL_FIELDS,
        "protein_identity_scope", "protein_cigar_scope", "protein_cds_phase",
        "reference_protein_id", "reference_transcript_id", "reference_protein_query_start",
        "reference_protein_query_end", "reference_cds_length", "reference_cds_phase",
        "alignment_hit_count", "alignment_ambiguous_hit_count", "alignment_mapping_quality", "alignment_alternative_hits",
        "dna_alignment_hit_count", "dna_alignment_ambiguous_hit_count", "dna_alignment_mapping_quality", "dna_alignment_alternative_hits",
        "deletion_alignment_backend", "deletion_alignment_cigar",
        "left_flank_identity", "left_flank_coverage", "left_flank_paired_bases",
        "left_flank_matches", "left_flank_mismatches",
        "right_flank_identity", "right_flank_coverage", "right_flank_paired_bases",
        "right_flank_matches", "right_flank_mismatches",
        "search_expansion_status", "search_original_start", "search_original_end",
        "search_expanded_start", "search_expanded_end", "range_status",
    ]
    write_tsv(result_dir / "sequence_synteny_evidence.tsv", evidence, fields)
    return evidence
