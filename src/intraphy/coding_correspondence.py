"""coding_correspondence: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import OrderedDict
from collections import defaultdict
from intraphy.preparation.features import FeatureHierarchy
from intraphy import alignment
from intraphy.coding.projection import _aligned_occurrence_pairs
from intraphy.coding.projection import _anchor_metrics
from intraphy.coding.projection import _candidate_json
from intraphy.coding.projection import _coordinate_blocks0
from intraphy.coding.transcripts import _copy_key
from intraphy.coding.transcripts import build_coding_transcript
from intraphy.coding.types import CodingProjectionCandidate
from intraphy.coding.types import CodingProjectionCandidateSet
from intraphy.coding.types import CodingResidueProjection
from intraphy.coding.types import FamilyCodingProjection
import json


class CodingProjectionIndex:
    """Family MSA index with transcript-specific CDS and genome projections."""

    MAX_CACHED_PAIRS = 32
    MAX_CACHED_BASE_PAIRS = 1_000_000

    def __init__(
        self, occurrences, sequences, transcript_paths, raw_features=(), threads=1,
        msa_mode="linsi",
    ):
        self.threads = threads
        self.msa_mode = msa_mode
        self.cache = OrderedDict()
        self.cached_bases = 0
        self.transcripts = {}
        self.by_occurrence = defaultdict(list)
        self.family_alignments = {}
        self.occurrences = {row["occurrence_id"]: row for row in occurrences}
        self.occurrences_by_copy = defaultdict(list)
        for row in occurrences:
            self.occurrences_by_copy[_copy_key(row)].append(row["occurrence_id"])
        by_copy = defaultdict(list)
        for row in raw_features:
            if row.get("ownership") == "target_gene_descendant":
                by_copy[_copy_key(row)].append(row)
        by_transcript = defaultdict(list)
        for row in transcript_paths:
            by_transcript[(*_copy_key(row), row["transcript_id"])].append(row)
        by_id = self.occurrences
        feature_indexes = {key: FeatureHierarchy(rows) for key, rows in by_copy.items()}
        for key, paths in sorted(by_transcript.items()):
            transcript = build_coding_transcript(key, paths, by_id, sequences, by_copy[key[:-1]], feature_index=feature_indexes.get(key[:-1]))
            self.transcripts[key] = transcript
            for occurrence_id in {row["occurrence_id"] for row in paths}:
                self.by_occurrence[occurrence_id].append(key)

    def _family_alignment(self, family_id):
        if family_id in self.family_alignments:
            return self.family_alignments[family_id]
        family_transcripts = [
            transcript for key, transcript in sorted(self.transcripts.items())
            if key[0] == family_id and not transcript.unavailable_reason and transcript.protein
        ]
        aliases_by_protein = defaultdict(list)
        for transcript in family_transcripts:
            aliases_by_protein[transcript.protein].append(transcript.key)
        records = []
        aliases_by_record = {}
        record_by_transcript = {}
        for index, protein in enumerate(sorted(aliases_by_protein), 1):
            record_id = f"coding_{index:06d}"
            aliases = tuple(sorted(aliases_by_protein[protein]))
            records.append((record_id, protein))
            aliases_by_record[record_id] = aliases
            for transcript_key in aliases:
                record_by_transcript[transcript_key] = record_id
        aligned_records = dict(alignment.protein_multiple_alignment(
            tuple(records), mode=self.msa_mode, threads=self.threads,
        )) if records else {}
        expected_ids = {record_id for record_id, _protein in records}
        if set(aligned_records) != expected_ids:
            raise alignment.AlignmentBackendError("family protein MSA returned a different record set")
        lengths = {len(sequence) for sequence in aligned_records.values()}
        if len(lengths) > 1:
            raise alignment.AlignmentBackendError("family protein MSA returned unequal alignment lengths")
        residue_columns = {}
        for record_id, protein in records:
            aligned = aligned_records[record_id]
            if aligned.replace("-", "") != protein:
                raise alignment.AlignmentBackendError("family protein MSA did not preserve input residues")
            residue_columns[record_id] = tuple(
                column0 for column0, amino_acid in enumerate(aligned) if amino_acid != "-"
            )
        projection = FamilyCodingProjection(
            family_id=family_id, mode=self.msa_mode,
            aligned_records=aligned_records,
            aliases_by_record=aliases_by_record,
            record_by_transcript=record_by_transcript,
            residue_columns=residue_columns,
        )
        self.family_alignments[family_id] = projection
        return projection

    def _pair(self, query_key, target_key):
        key = tuple(sorted((query_key, target_key)))
        if key in self.cache:
            self.cache.move_to_end(key)
            pairs, known_columns, aligned_left, aligned_right, _size = self.cache[key]
        else:
            query, target = (self.transcripts[item] for item in key)
            if query.key[0] != target.key[0]:
                return {}, (), "", "", False
            family = self._family_alignment(query.key[0])
            aligned_left = family.aligned_for(query.key)
            aligned_right = family.aligned_for(target.key)
            pairs, known_columns = _aligned_occurrence_pairs(
                query, target, aligned_left, aligned_right,
            )
            size = sum(len(record["positions0"]) for record in pairs.values())
            if size <= self.MAX_CACHED_BASE_PAIRS:
                while self.cache and (len(self.cache) >= self.MAX_CACHED_PAIRS or self.cached_bases + size > self.MAX_CACHED_BASE_PAIRS):
                    _old_key, (_old_pairs, _old_columns, _left, _right, old_size) = self.cache.popitem(last=False)
                    self.cached_bases -= old_size
                self.cache[key] = (pairs, known_columns, aligned_left, aligned_right, size)
                self.cached_bases += size
        inverted = query_key != key[0]
        if not inverted:
            return pairs, known_columns, aligned_left, aligned_right, False
        reversed_pairs = {}
        for (left_occurrence, right_occurrence), record in pairs.items():
            reversed_pairs[(right_occurrence, left_occurrence)] = {
                **record,
                "positions0": [(right0, left0) for left0, right0 in record["positions0"]],
            }
        reversed_columns = tuple(
            (column0, right_aa, left_aa, right_codon, left_codon)
            for column0, left_aa, right_aa, left_codon, right_codon in known_columns
        )
        return reversed_pairs, reversed_columns, aligned_right, aligned_left, True

    def project(self, occurrence_id):
        """Return residue-to-MSA-to-CDS/genome projections for one occurrence."""
        projections = []
        for transcript_key in sorted(self.by_occurrence.get(occurrence_id, ())):
            transcript = self.transcripts[transcript_key]
            if transcript.unavailable_reason:
                continue
            family = self._family_alignment(transcript_key[0])
            columns = family.columns_for(transcript_key)
            for residue_index0, codon in enumerate(transcript.codon_sources):
                if not any(base.occurrence_id == occurrence_id for base in codon):
                    continue
                projections.append(CodingResidueProjection(
                    transcript_key=transcript_key,
                    residue_index0=residue_index0,
                    msa_column0=columns[residue_index0],
                    amino_acid=transcript.protein[residue_index0],
                    cds_bases=codon,
                ))
        return tuple(projections)

    def family_projection(self, family_id):
        """Expose the run-local family MSA and its complete alias mapping."""
        return self._family_alignment(family_id)

    def candidates(self, query_occurrence, target_copy_id):
        """Return local coding evidence for every occurrence in a target copy."""
        query = self.occurrences.get(query_occurrence)
        if query is None:
            return tuple()
        if isinstance(target_copy_id, tuple):
            copy_key = target_copy_id
        else:
            matches = [
                key for key in self.occurrences_by_copy
                if key[0] == query.get("family_id") and key[2] == target_copy_id
            ]
            if len(matches) != 1:
                return tuple()
            copy_key = matches[0]
        return tuple(
            {"target_occurrence_id": target, **self.evidence(query_occurrence, target)}
            for target in sorted(self.occurrences_by_copy.get(copy_key, ()))
        )

    def _competing_occurrences(
        self, query_key, target_key, query_occurrence, target_occurrence,
        query_positions0, target_positions0,
    ):
        """Find alternative occurrence placements supported by compatible paths."""
        competitors = set()
        query_copy = query_key[:3]
        target_copy = target_key[:3]
        target_transcripts = [
            key for key, transcript in self.transcripts.items()
            if key[:3] == target_copy and not transcript.unavailable_reason
        ]
        query_transcripts = [
            key for key, transcript in self.transcripts.items()
            if key[:3] == query_copy and not transcript.unavailable_reason
        ]
        for alternative_target_key in target_transcripts:
            pairs, _columns, _aligned_query, _aligned_target, _inverted = self._pair(
                query_key, alternative_target_key,
            )
            for (other_query, other_target), record in pairs.items():
                if other_query != query_occurrence or other_target == target_occurrence:
                    continue
                other_query_positions0 = {
                    query0 for query0, _target0 in record["positions0"]
                }
                if query_positions0 & other_query_positions0:
                    competitors.add(f"target:{other_target}")
        for alternative_query_key in query_transcripts:
            pairs, _columns, _aligned_query, _aligned_target, _inverted = self._pair(
                alternative_query_key, target_key,
            )
            for (other_query, other_target), record in pairs.items():
                if other_target != target_occurrence or other_query == query_occurrence:
                    continue
                other_target_positions0 = {
                    target0 for _query0, target0 in record["positions0"]
                }
                if target_positions0 & other_target_positions0:
                    competitors.add(f"query:{other_query}")
        return tuple(sorted(competitors))

    def evidence(self, query_occurrence, target_occurrence):
        query_keys = self.by_occurrence.get(query_occurrence, [])
        target_keys = self.by_occurrence.get(target_occurrence, [])
        result = {
            "protein_status": "unavailable",
            "protein_mapping_status": "uncovered",
            "protein_unavailable_reason": "no_CDS_transcript_path",
            "protein_membership_eligible": False,
            "protein_position_eligible": False,
            "protein_hard_observation_eligible": False,
            "protein_candidate_evidence_available": False,
        }
        unavailable, candidates = set(), []
        for query_key in query_keys:
            query = self.transcripts[query_key]
            for target_key in target_keys:
                target = self.transcripts[target_key]
                for side, transcript in (("query", query), ("target", target)):
                    if transcript.unavailable_reason:
                        unavailable.add(f"{side}:{transcript.key[-1]}:{transcript.unavailable_reason}")
                if query.unavailable_reason or target.unavailable_reason:
                    continue
                if not query.coding_lengths.get(query_occurrence) or not target.coding_lengths.get(target_occurrence):
                    unavailable.add("no_CDS_in_requested_occurrence")
                    continue
                result["protein_status"] = "no_aligned_CDS"
                pairs, known_columns, aligned_query, aligned_target, _inverted = self._pair(
                    query_key, target_key,
                )
                record = pairs.get((query_occurrence, target_occurrence))
                if not record:
                    continue
                positions = set(record["positions0"])
                identity = record["aa_matches"] / record["aa_pairs"]
                query_coverage = len(positions) / query.coding_lengths[query_occurrence]
                target_coverage = len(positions) / target.coding_lengths[target_occurrence]
                anchors = _anchor_metrics(
                    record, known_columns, aligned_query, aligned_target,
                    query_occurrence, target_occurrence,
                )
                ordered = sorted(positions)
                position_monotonic = all(
                    right[0] > left[0] and right[1] > left[1]
                    for left, right in zip(ordered, ordered[1:])
                )
                query_positions0 = {query0 for query0, _target0 in positions}
                target_positions0 = {target0 for _query0, target0 in positions}
                candidates.append(CodingProjectionCandidate(
                    query_transcript_key=query_key,
                    target_transcript_key=target_key,
                    coordinate_blocks=_coordinate_blocks0(positions),
                    aa_identity=identity,
                    query_cds_coverage=query_coverage,
                    target_cds_coverage=target_coverage,
                    known_aa_pairs=record["aa_pairs"],
                    blosum62_score=record["blosum62_score"],
                    gap_fraction=anchors["gap_fraction"],
                    left_anchor_pairs=anchors["left_pairs"],
                    right_anchor_pairs=anchors["right_pairs"],
                    left_anchor_score=anchors["left_score"],
                    right_anchor_score=anchors["right_score"],
                    left_anchor_supported=anchors["left_supported"],
                    right_anchor_supported=anchors["right_supported"],
                    terminal_side=anchors["terminal_side"],
                    msa_column_interval=anchors["column_interval"],
                    anchor_resolved=anchors["resolved"],
                    position_monotonic=position_monotonic,
                    competing_occurrences=self._competing_occurrences(
                        query_key, target_key, query_occurrence, target_occurrence,
                        query_positions0, target_positions0,
                    ),
                ))
        result["protein_unavailable_reason"] = ";".join(sorted(unavailable)) or ("NA" if candidates or result["protein_status"] != "unavailable" else "no_CDS_transcript_path")
        if not candidates:
            return result
        candidates = tuple(sorted(
            candidates,
            key=lambda candidate: (
                candidate.query_transcript_key, candidate.target_transcript_key,
                candidate.coordinate_blocks,
            ),
        ))
        candidate_set = CodingProjectionCandidateSet(candidates)
        best = max(candidates, key=lambda candidate: (
            candidate.position_eligible, candidate.anchor_resolved,
            candidate.known_aa_pairs, candidate.blosum62_score,
            candidate.aa_identity,
            max(candidate.query_cds_coverage, candidate.target_cds_coverage),
        ))
        ambiguous = (
            not candidate_set.coordinate_consensus
            or bool(candidate_set.competing_occurrences)
            or any(not candidate.position_monotonic for candidate in candidates)
        )
        mapping_status = (
            "ambiguous_mapping" if ambiguous
            else "resolved_local" if candidate_set.position_eligible
            else "supported_unanchored"
        )
        query_sources = self.transcripts[best.query_transcript_key].source_ids
        target_sources = self.transcripts[best.target_transcript_key].source_ids
        interval = best.msa_column_interval
        # ``protein_status`` remains available for the existing caller. New
        # state construction must use mapping_status and hard eligibility.
        result.update(
            protein_status="ambiguous_transcript_projection" if ambiguous else "supported",
            protein_mapping_status=mapping_status,
            protein_membership_eligible=candidate_set.membership_eligible,
            protein_position_eligible=candidate_set.position_eligible,
            protein_hard_observation_eligible=candidate_set.position_eligible,
            protein_candidate_evidence_available=candidate_set.candidate_evidence_available,
            protein_aa_identity=best.aa_identity,
            protein_query_cds_coverage=best.query_cds_coverage,
            protein_target_cds_coverage=best.target_cds_coverage,
            protein_known_aa_pairs=best.known_aa_pairs,
            protein_blosum62_score=best.blosum62_score,
            protein_gap_fraction=best.gap_fraction,
            protein_left_anchor_pairs=best.left_anchor_pairs,
            protein_right_anchor_pairs=best.right_anchor_pairs,
            protein_left_anchor_score=best.left_anchor_score,
            protein_right_anchor_score=best.right_anchor_score,
            protein_left_anchor_supported=best.left_anchor_supported,
            protein_right_anchor_supported=best.right_anchor_supported,
            protein_terminal_side=best.terminal_side,
            protein_msa_column_start0=interval.start0,
            protein_msa_column_end0=interval.end0,
            protein_msa_column_interval=f"{interval.start0}:{interval.end0}",
            protein_msa_mode=self.msa_mode,
            protein_candidate_mapping_count=len(candidates),
            protein_candidate_coordinate_consensus=candidate_set.coordinate_consensus,
            protein_competing_occurrences=";".join(candidate_set.competing_occurrences) or "NA",
            protein_candidate_details=json.dumps(
                [_candidate_json(candidate) for candidate in candidates],
                sort_keys=True, separators=(",", ":"),
            ),
            protein_candidate_set=candidate_set,
            protein_query_source_features=";".join(query_sources) or "NA",
            protein_target_source_features=";".join(target_sources) or "NA",
            protein_metrics_scope="best_transcript_pair",
            protein_best_query_transcript=best.query_transcript_id,
            protein_best_target_transcript=best.target_transcript_id,
            protein_supporting_transcripts=";".join(
                f"{candidate.query_transcript_id}>{candidate.target_transcript_id}"
                for candidate in candidates
            ),
        )
        if ambiguous:
            return result
        blocks = candidates[0].coordinate_blocks
        result["protein_projected_coordinate_blocks"] = blocks
        result["protein_projected_blocks"] = blocks
        return result


# Public entry points; implementations have a single owner.
from intraphy.coding.projection import (
    _blosum62_score,
    _aligned_occurrence_pairs,
    _coordinate_blocks0,
    _anchor_metrics,
    _candidate_json,
)
from intraphy.coding.transcripts import (
    _copy_key,
    _transcript_features,
    _translation_reason,
    _synthetic_transcript_features,
    build_coding_transcript,
)
from intraphy.coding.types import (
    _STANDARD_CODE,
    _BLOSUM62,
    _EMPTY,
    _ANCHOR_WINDOW,
    _MIN_ANCHOR_PAIRS,
    CodingBaseProjection,
    CodingResidueProjection,
    CodingTranscript,
    FamilyCodingProjection,
    CodingProjectionCandidate,
    CodingProjectionCandidateSet,
)
import json
from intraphy import alignment
