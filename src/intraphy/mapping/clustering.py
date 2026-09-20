"""mapping / clustering: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from intraphy.coding_correspondence import CodingProjectionIndex
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0
from intraphy.coordinates import parse_legacy_blocks
from intraphy.mapping.alternative import alternative_overlap_evidence
from intraphy.mapping.candidate_coordinates import _coordinate_block0
from intraphy.mapping.candidate_coordinates import _format_alignment_blocks
from intraphy.mapping.candidate_coordinates import _genomic_matched_blocks
from intraphy.mapping.fields import EXON_LIKE_ROLES
from intraphy.mapping.match_context import cheap_match_evidence
from intraphy.mapping.match_context import copy_order_context
from intraphy.mapping.match_context import load_distance_table
from intraphy.mapping.match_context import pair_threshold
from intraphy.mapping.match_context import should_align_pair
from intraphy.mapping.match_records import _match_row
from intraphy.mapping.match_records import _projection_compatibility
from intraphy.mapping.match_records import _projection_record
from intraphy.mapping.ordered_chains import _apply_ordered_candidate_chains
from intraphy.mapping.pairwise_matches import match_evidence
from intraphy.mapping.policies import _candidate_sequence_accepted
from intraphy.mapping.policies import graph_components
from intraphy.mapping.policies import inferred_source_label
from intraphy.mapping.policies import known_source_labels
from intraphy.mapping.policies import occurrence_copy_key
from intraphy.mapping.policies import roles_compatible
from intraphy.mapping.policies import sequence_supported_mapping
from intraphy.mapping.short_alignment import _rerun_anchor_bounded_short_candidates
from intraphy.preparation.annotation_index import parse_attributes
from intraphy.storage.values import to_float
from itertools import islice


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
