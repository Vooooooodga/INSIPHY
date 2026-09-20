"""preprocess: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from insiphy.aligners.runner import available_alignment_backends
from insiphy.coding_correspondence import CodingProjectionIndex
from insiphy.coordinates import ClosedInterval1
from insiphy.coordinates import Interval0
from insiphy.coordinates import parse_legacy_blocks
from insiphy.mapping.alternative import alternative_overlap_evidence
from insiphy.mapping.candidate_codec import _coordinate_block0
from insiphy.mapping.candidate_codec import _format_alignment_blocks
from insiphy.mapping.candidate_codec import _genomic_matched_blocks
from insiphy.mapping.chains import _apply_ordered_candidate_chains
from insiphy.mapping.fields import EXON_LIKE_ROLES
from insiphy.mapping.fields import MATCH_FIELDS
from insiphy.mapping.fields import SEGMENT_FIELDS
from insiphy.mapping.match_records import _match_row
from insiphy.mapping.match_records import _projection_compatibility
from insiphy.mapping.match_records import _projection_record
from insiphy.mapping.matches import cheap_match_evidence
from insiphy.mapping.matches import copy_order_context
from insiphy.mapping.matches import load_distance_table
from insiphy.mapping.matches import match_evidence
from insiphy.mapping.matches import pair_threshold
from insiphy.mapping.matches import should_align_pair
from insiphy.mapping.policies import _candidate_sequence_accepted
from insiphy.mapping.policies import _species_tree_distances
from insiphy.mapping.policies import graph_components
from insiphy.mapping.policies import inferred_source_label
from insiphy.mapping.policies import known_source_labels
from insiphy.mapping.policies import occurrence_copy_key
from insiphy.mapping.policies import roles_compatible
from insiphy.mapping.policies import sequence_supported_mapping
from insiphy.mapping.short_context import _gene_locus_records
from insiphy.mapping.short_context import _rerun_anchor_bounded_short_candidates
from insiphy.preparation.annotation_index import parse_attributes
from insiphy.preparation.copy_context import make_adjacencies
from insiphy.preparation.copy_context import make_copy_context
from insiphy.preparation.copy_context import make_copy_relationships
from insiphy.storage.fasta import parse_fasta
from insiphy.storage.tabular import read_tsv
from insiphy.storage.tabular import write_tsv
from insiphy.storage.values import to_float
from itertools import islice
from pathlib import Path
import csv


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
                row["correspondence_score"] = (
                    0.70 * float(row["protein_aa_identity"])
                    + 0.30 * min(float(row["protein_query_cds_coverage"]),
                                 float(row["protein_target_cds_coverage"]))
                )
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
            and row.get("protein_hard_observation_eligible") in {1, "1", True}
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
            and row.get("protein_hard_observation_eligible") in {1, "1", True}
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


# Backward-compatible symbol exports; no alternate implementations.
from insiphy.mapping.alternative import (
    _relative_overlap_interval,
    alternative_overlap_evidence,
)
from insiphy.mapping.candidate_codec import (
    _coordinate_block0,
    _alignment_blocks0,
    _alignment_blocks,
    _format_alignment_blocks,
    _block_signature,
    _valid_local_boundary_range,
    _explicit_bounded_target,
    _genomic_matched_blocks,
    _genomic_blocks0,
    _format_genomic_blocks,
    _candidate_value,
    _alignment_gap_blocks,
    _covered_bases,
    _unknown_pair_count,
    _candidate_record,
    _public_interval,
    _public_gap_blocks,
    _public_candidate_record,
    _nt_column_score,
    _set_explicit_alignment_score,
    _alignment_candidate_records,
    _transpose_cigar,
    _transpose_gap_blocks,
    _transpose_candidate_record,
    _mapped_genomic_interval,
    _mapped_genomic_interval0,
)
from insiphy.mapping.chains import (
    _copy_transcription_bounds,
    _candidate_copy_interval,
    _chain_precedes,
    _candidate_path_memberships,
    _apply_ordered_candidate_chains,
)
from insiphy.mapping.fields import (
    SEGMENT_FIELDS,
    CorrespondenceCriteria,
    DEFAULT_CORRESPONDENCE_CRITERIA,
    SEGMENT_OUTPUT_FIELDS,
    TRANSCRIPT_PATH_FIELDS,
    RAW_FEATURE_FIELDS,
    INTRON_SITE_FIELDS,
    GENE_LOCUS_FIELDS,
    EXON_LIKE_ROLES,
    NONCODING_ROLES,
    STRUCTURAL_ROLES,
    UNKNOWN_SOURCE_LABELS,
    MATCH_FIELDS,
    _transcript_id_set,
)
from insiphy.mapping.match_records import (
    _projection_interval,
    _projection_record,
    _ordered_projection_compatible,
    _disjoint_reference_intervals,
    _projection_compatibility,
    _candidate_records_for_match,
    _format_optional_number,
    _format_contract_value,
    _match_row,
)
from insiphy.mapping.matches import (
    segment_length,
    phase_score,
    role_boundary_score,
    _path_order_context,
    copy_order_context,
    context_score,
    load_distance_table,
    pair_threshold,
    cheap_match_evidence,
    match_evidence,
    should_align_pair,
)
from insiphy.mapping.policies import (
    occurrence_copy_key,
    roles_compatible,
    sequence_supported_mapping,
    _candidate_sequence_accepted,
    assess_short_candidate_thresholds,
    split_source_labels,
    known_source_labels,
    inferred_source_label,
    graph_components,
    _species_tree_distances,
)
from insiphy.mapping.short_context import (
    _gene_locus_records,
    _value_tokens,
    _candidate_blocks_for_copy,
    _owner_occurrence_for_copy,
    _interval_between_transcript_flanks,
    _cut0_between_loci,
    _anchor_bounded_gap_blocks,
    _context_flank_pairs,
    _rerun_anchor_bounded_short_candidates,
)
from insiphy.preparation.annotation_index import (
    parse_attributes,
    split_ids,
    iter_annotation,
    read_annotation,
    _annotation_row_key,
    _copy_annotation_row,
    _read_annotation_for_gene_cached,
    read_annotation_for_gene,
    feature_tokens,
    feature_matches_gene,
    locate_gene,
    overlaps_gene,
)
from insiphy.preparation.copy_context import (
    make_adjacencies,
    copy_spans,
    infer_copy_class,
    make_copy_context,
    make_copy_relationships,
)
from insiphy.preparation.extract import (
    extract_gene,
)
from insiphy.preparation.transcripts import (
    sequence_slice,
    translate_cds,
    transcript_features,
    feature_role,
    _merge_exonic_intervals,
    _annotate_exon,
    _format_intervals,
    _feature_parents,
    _format_attrs,
    _partial_boundary,
    _coding_role,
    _position_role,
    _path_role_record,
    _interval_union_length,
    child_features_for_transcript,
    transcript_sort_key,
    transcript_cds_length,
    select_transcripts,
    introns_from_path,
    exon_features_for_gene,
    load_existing_segments,
    _feature_key,
    _occurrence_id,
)
import math
import csv
import json
import re
