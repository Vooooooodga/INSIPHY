"""annotation: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import Counter
from collections import defaultdict
from insiphy.aligners.pairwise import local_alignment_stats
from insiphy.aligners.types import AlignmentBackendError
from insiphy.elements import EXON_LIKE_ROLES
from insiphy.evidence.deletion import _alignment_evidence_provenance
from insiphy.evidence.deletion import _ambiguous_repeated_mapping
from insiphy.evidence.deletion import _format_alternative_hits
from insiphy.evidence.deletion import _nucleotide_interval_candidates
from insiphy.evidence.deletion import _offset_alignment_target
from insiphy.evidence.deletion import _ordered_anchor_deletion_support
from insiphy.evidence.deletion import _projection_within_anchor_interval
from insiphy.evidence.deletion import _set_alignment_context
from insiphy.evidence.fields import _EVIDENCE_PROVENANCE_FIELDS
from insiphy.evidence.fields import _PROTEIN_QUERY_INTERVAL_FIELDS
from insiphy.evidence.projection import _alignment_genome_blocks
from insiphy.evidence.projection import _annotation_overlaps_for_blocks
from insiphy.evidence.projection import _candidate_span_blocks
from insiphy.evidence.projection import _json_records
from insiphy.evidence.projection import _locus_geometry
from insiphy.evidence.projection import _locus_key
from insiphy.evidence.projection import _project_locus_interval
from insiphy.evidence.projection import _protein_projection_blocks
from insiphy.evidence.projection import _supplied_annotation_role
from insiphy.evidence.projection import _transcripts_for_occurrence
from insiphy.evidence.protein import _cached_protein_projection
from insiphy.evidence.protein import _protein_context_provenance
from insiphy.evidence.protein import _reference_protein_context
from insiphy.evidence.search_context import _extended_locus
from insiphy.evidence.search_context import _ordered_anchor_search_context
from insiphy.evidence.search_context import _read_locus_metadata
from insiphy.evidence.search_context import _safe_int
from insiphy.evidence.search_context import _search_provenance
from insiphy.evidence.search_context import _source_occurrence_path
from insiphy.storage.fasta import parse_fasta
from insiphy.storage.tabular import read_tsv
from insiphy.storage.tabular import write_tsv
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

    representative_infos = []
    for (family, element), members in sorted(rows_by_element.items()):
        representative_row, representative = max(
            members,
            key=lambda pair: len(sequences.get(pair[1]["occurrence_id"], "")),
        )
        if aligner == "miniprot":
            protein_members = [
                pair for pair in members
                if _reference_protein_context(pair[1], transcript_paths, proteins, occ_by_id) is not None
            ]
            if protein_members:
                representative_row, representative = min(
                    protein_members,
                    key=lambda pair: (-_safe_int(pair[1].get("cds_length"), 0), pair[1]["occurrence_id"]),
                )
        query = sequences.get(representative["occurrence_id"], "")
        if not query:
            continue
        representative_infos.append(
            {
                "family": family,
                "element": element,
                "members": members,
                "representative_row": representative_row,
                "representative": representative,
                "query": query,
                "protein_context": _reference_protein_context(representative, transcript_paths, proteins, occ_by_id)
                if aligner == "miniprot"
                else None,
            }
        )

    protein_contexts_by_family = defaultdict(dict)
    if aligner == "miniprot":
        for info in representative_infos:
            context = info.get("protein_context")
            if context is not None:
                protein_contexts_by_family[info["family"]][context["protein_id"]] = context

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
            alignment = None
            alignment_error = ""
            alignment_evidence_scope = "unavailable"
            identity = coverage = 0.0
            dna_identity = dna_coverage = 0.0
            nucleotide_aligner = "minimap2" if aligner == "miniprot" else aligner
            bounded_target = anchor_context.get("sequence", "")
            bounded_backend = (
                "internal"
                if len(query) <= int(short_context_max_length)
                else nucleotide_aligner
            )
            if anchor_context["valid_double_flank"] and bounded_target:
                try:
                    alignment = local_alignment_stats(
                        query,
                        bounded_target,
                        backend=bounded_backend,
                        threads=threads,
                        query_occurrence_id=representative.get("occurrence_id"),
                        target_occurrence_id=None,
                        query_transcript_id=(source_transcript_ids[0] if len(source_transcript_ids) == 1 else None),
                        left_anchor_id=anchor_context["left_anchor_id"],
                        right_anchor_id=anchor_context["right_anchor_id"],
                        search_interval=anchor_context["search_interval"],
                    )
                    alignment_evidence_scope = (
                        "anchor_bounded_short_local"
                        if bounded_backend == "internal"
                        else "anchor_bounded_external_local"
                    )
                except AlignmentBackendError as exc:
                    if bounded_backend == "internal" and nucleotide_aligner != "internal":
                        try:
                            alignment = local_alignment_stats(
                                query,
                                bounded_target,
                                backend=nucleotide_aligner,
                                threads=threads,
                            )
                            alignment_evidence_scope = "anchor_bounded_external_local"
                        except AlignmentBackendError as fallback_exc:
                            alignment_error = f"bounded_alignment_unavailable:{fallback_exc}"
                    else:
                        alignment_error = f"bounded_alignment_unavailable:{exc}"
                if alignment is not None:
                    _set_alignment_context(
                        alignment,
                        representative.get("occurrence_id", "NA"),
                        None,
                        source_transcript_ids[0] if len(source_transcript_ids) == 1 else None,
                        anchor_context["left_anchor_id"],
                        anchor_context["right_anchor_id"],
                        anchor_context["search_interval"],
                    )
                    _offset_alignment_target(alignment, anchor_context["target_offset0"])
            if alignment is None and not anchor_context["valid_double_flank"]:
                try:
                    alignment = local_alignment_stats(
                        query, locus, backend=nucleotide_aligner, threads=threads
                    )
                    alignment_evidence_scope = "whole_locus_descriptive_fallback"
                except AlignmentBackendError as exc:
                    alignment_error = str(exc)
            if alignment is not None:
                dna_identity = alignment.identity
                dna_coverage = alignment.query_coverage or alignment.coverage
                identity = dna_identity
                coverage = dna_coverage
            supported = (
                alignment_evidence_scope.startswith("anchor_bounded_")
                and alignment is not None
                and (bool(alignment.aligned_blocks) or int(alignment.aligned_pairs or 0) > 0)
                and identity >= min_identity
                and coverage >= min_coverage
            )
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
            conflict_blocks = []
            correspondence_status = "resolved"
            if not anchor_context["valid_double_flank"]:
                has_descriptive_candidate = (
                    (
                        alignment is not None
                        and (bool(alignment.aligned_blocks) or int(alignment.aligned_pairs or 0) > 0)
                        and dna_identity >= min_identity
                        and dna_coverage >= min_coverage
                    )
                    or bool(protein_candidates)
                )
                primary_mapping_status = (
                    "whole_locus_descriptive_candidate" if has_descriptive_candidate else "unresolved"
                )
                interval_scope = alignment_evidence_scope
                correspondence_status = "unknown"
                primary_backend = (
                    "miniprot"
                    if protein_candidates
                    else alignment.backend if alignment is not None else nucleotide_aligner
                )
                primary_cigar = alignment.cigar if alignment is not None else "NA"
                primary_score = dna_identity if alignment is not None else 0.0
                primary_coverage = dna_coverage if alignment is not None else 0.0
                predicted_role = "unknown"
                inferred_role = "unknown"
                status = "ambiguous"
                annotation_status = anchor_context["status"]
                conclusion = "unknown"
                contig = metadata.get("contig") or _locus_geometry(locus_header)[0]
                start = end = "NA"
                strand = metadata.get("strand") or _locus_geometry(locus_header)[3]
                confidence = "low"
            elif ambiguous_dna or ambiguous_protein:
                primary_mapping_status = "ambiguous_repeated_mapping"
                interval_scope = "ambiguous_candidates"
                correspondence_status = "unknown"
                primary_backend = "miniprot" if bounded_protein_candidates else alignment.backend
                primary_cigar = "NA"
                if bounded_protein_candidates:
                    primary_candidate = max(bounded_protein_candidates, key=lambda row: (row["coverage"], row["identity"]))
                    primary_score = primary_candidate["identity"]
                    primary_coverage = primary_candidate["coverage"]
                else:
                    primary_score, primary_coverage = dna_identity, dna_coverage
                predicted_role = "CDS" if bounded_protein_candidates else "unknown"
                inferred_role = "unknown"
                status = "homologous_sequence_candidate"
                annotation_status = "alignment_ambiguous_repeated_mapping"
                conclusion = "sequence_present_repeated_mapping_ambiguous"
                contig = _locus_geometry(locus_header)[0]
                start = end = strand = "NA"
                confidence = "low"
            elif projected:
                primary_mapping_status = "protein_projection"
                interval_scope = "projected_cds_container"
                primary_backend = "miniprot"
                primary_cigar = projected.get("cigar", "NA")
                primary_score = protein_identity
                primary_coverage = protein_coverage
                predicted_role = "CDS"
                inferred_role = "predicted_CDS"
                contig = projected["contig"]
                start = projected["start"]
                end = projected["end"]
                strand = projected["strand"]
                overlapping_exons = [
                    row for row in projection_annotation_overlaps
                    if row.get("role") in EXON_LIKE_ROLES
                    and row.get("strand_relation") == "sense"
                ]
                contained_roles = {
                    row["role"] for row in overlapping_exons if row.get("contains_block")
                }
                if strand != _locus_geometry(locus_header)[3]:
                    status = "homologous_sequence_candidate"
                    annotation_status = "protein_projection_antisense_to_gene_locus"
                    inferred_role = "unknown"
                    predicted_role = "CDS"
                    conclusion = "sequence_present_role_unknown"
                elif contained_roles:
                    status = "supports_annotation"
                    annotation_status = "protein_projection_contained_in_annotated_exon"
                    inferred_role = "CDS" if "CDS" in contained_roles else sorted(contained_roles)[0]
                    conclusion = "annotated_exon_sequence_present"
                elif overlapping_exons:
                    status = "conflicts_annotation"
                    annotation_status = "protein_projection_boundary_conflict"
                    inferred_role = "predicted_CDS"
                    conclusion = "predicted_cds_boundary_conflict"
                    conflict_blocks = overlapping_exons
                else:
                    status = "supports_hidden_segment"
                    annotation_status = "protein_projection_supports_missing_cds"
                    conclusion = "predicted_exon_candidate"
                confidence = "high" if identity >= min_identity and coverage >= min_coverage else "low"
            elif supported:
                interval_scope = "aligned_sequence"
                primary_backend = alignment.backend
                primary_cigar = alignment.cigar
                primary_score = dna_identity
                primary_coverage = dna_coverage
                predicted_role = "unknown"
                contig, start, end, strand = _project_locus_interval(
                    alignment.target_start, alignment.target_end, locus_header
                )
                if alignment.strand in {"+", "-"} and strand in {"+", "-"}:
                    strand = "+" if alignment.strand == strand else "-"
                primary_mapping_status = "unique_nucleotide_mapping"
                sense_roles = {
                    row["role"] for row in dna_annotation_overlaps
                    if row.get("strand_relation") == "sense"
                }
                antisense_exonic = any(
                    row.get("role") in EXON_LIKE_ROLES
                    and row.get("strand_relation") == "antisense"
                    for row in dna_annotation_overlaps
                )
                exonic_roles = sense_roles & EXON_LIKE_ROLES
                if exonic_roles:
                    status = "supports_annotation"
                    annotation_status = "alignment_overlaps_annotated_exon"
                    inferred_role = "CDS" if "CDS" in exonic_roles else sorted(exonic_roles)[0]
                    conclusion = "annotated_exon_sequence_present"
                else:
                    status = "homologous_sequence_candidate"
                    if antisense_exonic:
                        annotation_status = "homologous_sequence_overlaps_antisense_exon_annotation"
                    elif sense_roles:
                        annotation_status = "homologous_sequence_overlaps_non_exonic_annotation"
                    else:
                        annotation_status = "unannotated_homologous_sequence_candidate"
                    inferred_role = "unknown"
                    conclusion = "sequence_present_role_unknown"
                confidence = "medium"
            elif deletion_supported:
                primary_mapping_status = "ordered_flank_deletion"
                interval_scope = "not_applicable"
                primary_backend = deletion_provenance.get("backend", "NA")
                primary_cigar = deletion_provenance.get("cigar", "NA")
                primary_score = None
                primary_coverage = None
                predicted_role = "unknown"
                status = "supports_absence"
                annotation_status = "sequence_absence_between_ordered_flanking_homologs"
                inferred_role = "unknown"
                contig = metadata.get("contig") or "supplied_gene_locus"
                start = end = "NA"
                strand = metadata.get("strand") or "+"
                conclusion = "targeted_deletion_supported"
                confidence = "medium"
            else:
                primary_mapping_status = "unresolved"
                interval_scope = "unresolved"
                correspondence_status = "unknown"
                primary_backend = nucleotide_aligner
                primary_cigar = "NA"
                primary_score = 0.0
                primary_coverage = 0.0
                predicted_role = "unknown"
                status = "ambiguous"
                annotation_status = (
                    "alignment_backend_error" if alignment_error else
                    projection_reason if aligner == "miniprot" and projection_reason != "not_requested" else
                    deletion_reason
                )
                inferred_role = "unknown"
                contig = metadata.get("contig") or "supplied_gene_locus"
                start = end = "NA"
                strand = metadata.get("strand") or "+"
                conclusion = "unknown"
                confidence = "low"

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


# Backward-compatible symbol exports; no alternate implementations.
from insiphy.evidence.completion import (
    support_score,
    completion_call,
    complete_annotation,
)
from insiphy.evidence.deletion import (
    _offset_alignment_target,
    _set_alignment_context,
    _projection_within_anchor_interval,
    _cigar_ops,
    _supported_alignment_blocks,
    _query_only_gaps,
    _alternative_supports_expected_deletion,
    _projected_interval_coverage,
    _flank_support_from_spanning_blocks,
    _format_alternative_hits,
    _alignment_evidence_provenance,
    _ambiguous_repeated_mapping,
    _nucleotide_interval_candidates,
    _ordered_anchor_deletion_support,
)
from insiphy.evidence.fields import (
    _PROTEIN_QUERY_INTERVAL_FIELDS,
    _EVIDENCE_PROVENANCE_FIELDS,
    _SEARCH_COMPAT_FIELDS,
)
from insiphy.evidence.projection import (
    _locus_key,
    _locus_geometry,
    _project_locus_interval,
    _project_occurrence_interval,
    _json_records,
    _transcripts_for_occurrence,
    _alignment_genome_blocks,
    _protein_projection_blocks,
    _candidate_span_blocks,
    _annotation_overlaps_for_blocks,
    _supplied_annotation_role,
    _overlapping_annotation_role,
    _split_transcript_ids,
    _oriented_locus_slice,
    _relative_interval_within_span,
)
from insiphy.evidence.protein import (
    _protein_context_provenance,
    _coding_length_and_phase,
    _reference_protein_context,
    _protein_projection,
    _protein_projection_from_rows,
    _cached_protein_projection,
)
from insiphy.evidence.search_context import (
    _read_locus_metadata,
    _search_provenance,
    _extended_locus,
    _occurrence_for_element,
    _occurrences_for_element,
    _safe_int,
    _path_occurrences,
    _disjoint_genomic_occurrences,
    _source_occurrence_path,
    _ordered_occurrence_triplet,
    _ordered_occurrence_pair,
    _ordered_anchor_search_context,
)
import json
