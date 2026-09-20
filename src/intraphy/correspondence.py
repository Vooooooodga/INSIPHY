"""Correspondence workflow and backwards-compatible public helpers."""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from intraphy.coordinates import ClosedInterval1, Interval0, parse_legacy_blocks
from intraphy.elements import collect_element_profiles, element_class_for_occurrence, element_role_from_occurrences, membership_call
from intraphy.alignment import global_alignment_stats
from intraphy.io import norm_state, parse_fasta, read_tsv, to_float, uniq, write_tsv
from intraphy.topology import SpeciesTree
from intraphy.mapping.membership import CORRESPONDENCE_EVIDENCE_FIELDS
from intraphy.mapping.membership import _hard_membership_eligible
from intraphy.mapping.membership import _membership_match_details
from intraphy.mapping.profiles import match_total
from intraphy.mapping.profiles import normalized_components
from intraphy.mapping.profiles import progressive_tier
from intraphy.mapping.profiles import reciprocal_status
from intraphy.mapping.profiles import tree_distances
from intraphy.reporting.summaries import summarize_intragenic_paths
from intraphy.reporting.summaries import summarize_observed_element_tree_coverage
from intraphy.reporting.summaries import summarize_progressive_elements


def simple_identity(seq_a, seq_b):
    return global_alignment_stats(seq_a, seq_b, backend="auto").identity


def infer_correspondence(input_dir, output_dir, *, annotation_rows=None, annotation_completion_path=None):
    homology = read_tsv(f"{input_dir}/segment_homology.tsv", ["homology_id", "occurrence_id", "support_type", "confidence"])
    occurrences = read_tsv(f"{input_dir}/segment_occurrences.tsv", ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"])
    matches = read_tsv(f"{input_dir}/segment_matches.tsv", ["match_id", "query_occurrence_id", "subject_occurrence_id", "match_status"], optional=True)
    evidence_rows = read_tsv(f"{input_dir}/sequence_synteny_evidence.tsv", optional=True)
    if annotation_rows is not None and annotation_completion_path is not None:
        raise ValueError("provide annotation_rows or annotation_completion_path, not both")
    if annotation_completion_path is not None:
        annotation_rows = read_tsv(annotation_completion_path)
    annotation_rows = list(annotation_rows or ())
    seqs = parse_fasta(f"{input_dir}/segment_sequences.fasta")
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    distances = tree_distances(input_dir)
    profiles, element_by_homology = collect_element_profiles(homology, occurrences, evidence_rows + annotation_rows)

    component_rows = []
    element_rows = []
    by_component = defaultdict(list)
    for row in homology:
        occ = occ_by_id.get(row["occurrence_id"], {})
        by_component[row["homology_id"]].append(row["occurrence_id"])
        score, call = membership_call(row.get("confidence"), 0)
        component_rows.append(
            {
                "homology_id": row["homology_id"],
                "occurrence_id": row["occurrence_id"],
                "family_id": occ.get("family_id", "NA"),
                "species": occ.get("species", "NA"),
                "gene_copy_id": occ.get("gene_copy_id", "NA"),
                "source_label": row.get("source_label", "NA"),
                "support_type": row.get("support_type", "NA"),
                "confidence": row.get("confidence", "NA"),
                "membership_score": f"{score:.6g}",
                "membership_call": call,
            }
        )
        element_id = element_by_homology.get(row["homology_id"])
        if element_id:
            element_class = element_class_for_occurrence(occ, row["homology_id"], profiles, element_by_homology)
            element_rows.append(
                {
                    "element_id": element_id,
                    "family_id": occ.get("family_id", "NA"),
                    "homology_id": row["homology_id"],
                    "occurrence_id": row["occurrence_id"],
                    "species": occ.get("species", "NA"),
                    "gene_copy_id": occ.get("gene_copy_id", "NA"),
                    "element_class": element_class,
                    "display_role": element_role_from_occurrences([occ]),
                    "source_label": row.get("source_label", "NA"),
                    "support_type": row.get("support_type", "NA"),
                    "confidence": row.get("confidence", "NA"),
                    "membership_score": f"{score:.6g}",
                    "membership_call": call,
                }
            )

    conservation = []
    for component_id, occ_ids in sorted(by_component.items()):
        ids_with_seq = [occ_id for occ_id in occ_ids if occ_id in seqs]
        identities = []
        for i, left in enumerate(ids_with_seq):
            for right in ids_with_seq[i + 1 :]:
                identities.append(simple_identity(seqs[left], seqs[right]))
        mean_identity = sum(identities) / len(identities) if identities else 0.0
        roles = Counter(occ_by_id.get(occ_id, {}).get("role", "unknown") for occ_id in occ_ids)
        observed_species = {
            occ_by_id[occ_id]["species"] for occ_id in occ_ids
            if norm_state(occ_by_id.get(occ_id, {}).get("presence_status")) == "present"
            and norm_state(occ_by_id.get(occ_id, {}).get("species")) != "unknown"
        }
        observation_call = (
            "observed_in_multiple_species" if len(observed_species) > 1
            else "observed_single_species" if observed_species
            else "no_species_presence_observed"
        )
        conservation.append(
            {
                "homology_id": component_id,
                "occurrence_count": len(occ_ids),
                "species_count": len(observed_species),
                "mean_pairwise_identity": f"{mean_identity:.6g}",
                "role_spectrum": ";".join(f"{key}:{roles[key]}" for key in sorted(roles)),
                "conservation_call": observation_call,
            }
        )

    scored = []
    grouped = defaultdict(list)
    for row in matches:
        score = match_total(row)
        components = normalized_components(row)
        if components["sequence_score"] is None:
            components["sequence_score"] = 0.70 * components["alignment_score"] + 0.30 * components["coverage_score"]
        if components["structural_context_score"] is None:
            components["structural_context_score"] = 0.5 * components["left_context_score"] + 0.5 * components["right_context_score"]
        grouped[row["query_occurrence_id"]].append((score, row["match_id"]))
        scored.append(
            {
                "match_id": row["match_id"],
                "query_occurrence_id": row["query_occurrence_id"],
                "subject_occurrence_id": row["subject_occurrence_id"],
                "alignment_score": f"{components['alignment_score']:.6g}",
                "coverage_score": f"{components['coverage_score']:.6g}",
                "sequence_score": f"{components['sequence_score']:.6g}",
                "structural_context_score": f"{components['structural_context_score']:.6g}",
                "left_context_score": f"{components['left_context_score']:.6g}",
                "right_context_score": f"{components['right_context_score']:.6g}",
                "boundary_score": f"{components['boundary_score']:.6g}",
                "phase_score": f"{components['phase_score']:.6g}",
                "order_score": f"{components['order_score']:.6g}",
                "strand_score": f"{components['strand_score']:.6g}",
                "splice_score": f"{components['splice_score']:.6g}",
                "total_score": f"{score:.6g}",
                "dna_total_score": row.get("total_score", "NA"),
                "correspondence_basis": row.get("correspondence_basis", "DNA"),
                "protein_status": row.get("protein_status", "not_requested"),
                "protein_aa_identity": row.get("protein_aa_identity", "NA"),
                "protein_query_cds_coverage": row.get("protein_query_cds_coverage", "NA"),
                "protein_target_cds_coverage": row.get("protein_target_cds_coverage", "NA"),
                "protein_metrics_scope": row.get("protein_metrics_scope", "NA"),
                "protein_best_query_transcript": row.get("protein_best_query_transcript", "NA"),
                "protein_best_target_transcript": row.get("protein_best_target_transcript", "NA"),
                "protein_supporting_transcripts": row.get("protein_supporting_transcripts", "NA"),
                "match_status": row.get("match_status", "ambiguous"),
                "distance_class": row.get("distance_class", "NA"),
                "threshold": row.get("threshold", "NA"),
                "alignment_cigar": row.get("alignment_cigar", "NA"),
                "reciprocal_status": "unclassified",
                "correspondence_call": "unclassified",
                **{
                    field: row.get(field, "NA")
                    for field in CORRESPONDENCE_EVIDENCE_FIELDS
                },
            }
        )
    best_status = {}
    for vals in grouped.values():
        best = max(score for score, _ in vals)
        best_ids = [match_id for score, match_id in vals if abs(score - best) < 1e-12]
        for match_id in best_ids:
            best_status[match_id] = "best_tie" if len(best_ids) > 1 else "best_unique"
    reciprocal = reciprocal_status(scored)
    graph_edges = []
    degree = Counter()
    for row in scored:
        row["reciprocal_status"] = reciprocal.get(row["match_id"], "not_reciprocal_best")
        if row["match_status"] == "mapped" and _hard_membership_eligible(row):
            degree[row["query_occurrence_id"]] += 1
            degree[row["subject_occurrence_id"]] += 1
            graph_edges.append(
                {
                    "edge_id": row["match_id"],
                    "query_occurrence_id": row["query_occurrence_id"],
                    "subject_occurrence_id": row["subject_occurrence_id"],
                    "total_score": row["total_score"],
                    "reciprocal_status": row["reciprocal_status"],
                    "edge_call": (
                        "high_confidence_correspondence"
                        if row.get("candidate_resolution") == "resolved"
                        and row["reciprocal_status"] == "reciprocal_best"
                        and to_float(row["total_score"]) >= 0.7
                        else "supporting_correspondence"
                    ),
                }
            )
        if row["match_status"] != "mapped":
            row["correspondence_call"] = row["match_status"]
        elif row.get("candidate_resolution") == "ambiguous":
            row["correspondence_call"] = "ambiguous_candidate_mapping"
        elif row.get("candidate_resolution") == "excluded":
            row["correspondence_call"] = "outside_near_optimal_chain"
        elif row["reciprocal_status"] == "reciprocal_best":
            row["correspondence_call"] = "reciprocal_best"
        else:
            row["correspondence_call"] = best_status.get(row["match_id"], "not_best")

    _membership_match_details(element_rows, occurrences, scored)

    progressive_rows = []
    for row in scored:
        left = occ_by_id.get(row["query_occurrence_id"], {})
        right = occ_by_id.get(row["subject_occurrence_id"], {})
        pair_distance = distances.get((left.get("species"), right.get("species")), "NA")
        pair_distance_text = f"{pair_distance:.6g}" if pair_distance != "NA" else "NA"
        progressive_rows.append(
            {
                "match_id": row["match_id"],
                "query_species": left.get("species", "NA"),
                "subject_species": right.get("species", "NA"),
                "query_gene_copy_id": left.get("gene_copy_id", "NA"),
                "subject_gene_copy_id": right.get("gene_copy_id", "NA"),
                "tree_distance": pair_distance_text,
                "progressive_tier": progressive_tier(pair_distance_text),
                "total_score": row["total_score"],
                "correspondence_call": row["correspondence_call"],
            }
        )
    progressive_element_rows = summarize_progressive_elements(input_dir, occurrences, element_rows, scored)
    observed_coverage_rows = summarize_observed_element_tree_coverage(input_dir, element_rows)
    path_rows = summarize_intragenic_paths(occurrences, element_rows)

    for row in component_rows:
        deg = degree.get(row["occurrence_id"], 0)
        adjusted, call = membership_call(row.get("confidence"), deg)
        row["membership_score"] = f"{adjusted:.6g}"
        row["membership_call"] = call
    for row in element_rows:
        occ_id = row["occurrence_id"]
        adjusted, call = membership_call(row.get("confidence"), degree.get(occ_id, 0))
        row["membership_score"] = f"{adjusted:.6g}"
        row["membership_call"] = call

    write_tsv(f"{output_dir}/internal_homology_assignments.tsv", component_rows, ["homology_id", "occurrence_id", "family_id", "species", "gene_copy_id", "source_label", "support_type", "confidence", "membership_score", "membership_call"])
    write_tsv(
        f"{output_dir}/element_correspondence.tsv",
        element_rows,
        [
            "element_id",
            "parent_element_id",
            "member_interval_id",
            "family_id",
            "homology_id",
            "occurrence_id",
            "species",
            "gene_copy_id",
            "element_class",
            "display_role",
            "source_label",
            "support_type",
            "confidence",
            "membership_score",
            "membership_call",
            "matched_blocks",
            "genomic_matched_blocks",
            "actual_matched_blocks",
            "parent_feature_ids",
            "transcript_ids",
            "path_roles",
            "coding_roles",
            "position_roles",
            "annotation_source",
            "original_attributes",
            "partial_start",
            "partial_end",
            "path_role_records",
            "repeat_instance_id",
            "correspondence_status",
            "candidate_ids",
            "membership_edge_eligible",
            "membership_edge_reason",
            "position_edge_eligible",
            "position_edge_reason",
            "reference_occurrence_id",
            "reference_length",
            "reference_length_source",
            "reference_coverage_relation",
            "reference_covered_bases",
            "reference_overlap_bases",
            "reference_uncovered_bases",
        ],
    )
    write_tsv(f"{output_dir}/segment_conservation.tsv", conservation, ["homology_id", "occurrence_count", "species_count", "mean_pairwise_identity", "role_spectrum", "conservation_call"])
    write_tsv(
        f"{output_dir}/segment_correspondence.tsv",
        scored,
        [
            "match_id",
            "query_occurrence_id",
            "subject_occurrence_id",
            "alignment_score",
            "coverage_score",
            "sequence_score",
            "structural_context_score",
            "left_context_score",
            "right_context_score",
            "boundary_score",
            "phase_score",
            "order_score",
            "strand_score",
            "splice_score",
            "total_score",
            "dna_total_score",
            "correspondence_basis",
            "protein_status",
            "protein_aa_identity",
            "protein_query_cds_coverage",
            "protein_target_cds_coverage",
            "protein_metrics_scope",
            "protein_best_query_transcript",
            "protein_best_target_transcript",
            "protein_supporting_transcripts",
            "distance_class",
            "threshold",
            "alignment_cigar",
            "match_status",
            "reciprocal_status",
            "correspondence_call",
            *CORRESPONDENCE_EVIDENCE_FIELDS,
        ],
    )
    write_tsv(f"{output_dir}/internal_homology_graph_edges.tsv", graph_edges, ["edge_id", "query_occurrence_id", "subject_occurrence_id", "total_score", "reciprocal_status", "edge_call"])
    write_tsv(f"{output_dir}/progressive_correspondence.tsv", progressive_rows, ["match_id", "query_species", "subject_species", "query_gene_copy_id", "subject_gene_copy_id", "tree_distance", "progressive_tier", "total_score", "correspondence_call"])
    write_tsv(f"{output_dir}/progressive_element_correspondence.tsv", progressive_element_rows, ["element_id", "family_id", "homology_ids", "member_count", "species_count", "copy_count", "exon_like_members", "candidate_source_members", "nearest_species_support", "within_clade_support", "deep_tree_support", "best_pair_score", "mean_pair_score", "progressive_call"])
    write_tsv(f"{output_dir}/observed_element_tree_coverage.tsv", observed_coverage_rows, ["family_id", "element_id", "mrca_label", "present_species_count", "tree_tip_count", "coverage_class", "present_species", "present_copies"])
    write_tsv(f"{output_dir}/observed_intragenic_paths.tsv", path_rows, ["family_id", "path_scope", "node_label", "species", "gene_copy_id", "path_type", "element_path", "context_count", "path_support"])
    return component_rows, conservation, scored

# Public helper compatibility; each implementation has one owner.
from intraphy.reporting.summaries import _ancestor_chain
from intraphy.mapping.membership import _block_union_length
from intraphy.mapping.membership import _blocks_overlap
from intraphy.mapping.membership import _candidate_block_evidence
from intraphy.reporting.summaries import _coverage_class
from intraphy.mapping.membership import _edge_eligible
from intraphy.mapping.membership import _genomic_block_records
from intraphy.mapping.membership import _hard_position_eligible
from intraphy.reporting.summaries import _mrca_label
from intraphy.mapping.membership import _parse_genomic_blocks
from intraphy.mapping.membership import _tagged_membership_blocks
from intraphy.mapping.membership import _tokens
from intraphy.reporting.summaries import _tree_depths
