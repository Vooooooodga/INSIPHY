"""evidence / classification: explicit implementation ownership."""
from __future__ import annotations

from intraphy.elements import EXON_LIKE_ROLES
from intraphy.evidence.projection import _locus_geometry
from intraphy.evidence.projection import _project_locus_interval


def _classify_sequence_evidence(anchor_context, alignment, dna_identity, min_identity, dna_coverage, min_coverage, protein_candidates, alignment_evidence_scope, nucleotide_aligner, metadata, locus_header, ambiguous_dna, ambiguous_protein, bounded_protein_candidates, projected, protein_identity, protein_coverage, projection_annotation_overlaps, identity, coverage, supported, dna_annotation_overlaps, deletion_supported, deletion_provenance, alignment_error, aligner, projection_reason, deletion_reason):
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
    return status, annotation_status, inferred_role, predicted_role, conflict_blocks, contig, start, end, strand, primary_score, primary_coverage, conclusion, confidence, primary_backend, primary_mapping_status, interval_scope, primary_cigar, correspondence_status
