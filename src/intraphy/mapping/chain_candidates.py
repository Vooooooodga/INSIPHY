"""mapping / chain candidates: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.mapping.chain_types import ChainCandidate
from intraphy.mapping.path_membership import _candidate_copy_interval
from intraphy.mapping.path_membership import _candidate_path_memberships
from intraphy.mapping.path_membership import _copy_transcription_bounds
from intraphy.mapping.policies import occurrence_copy_key


def _collect_chain_candidates(occurrences, rows, occurrence_by_id, transcript_paths):
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
    return original_intervals, record_by_candidate, owner_by_candidate, groups, candidate_by_id
