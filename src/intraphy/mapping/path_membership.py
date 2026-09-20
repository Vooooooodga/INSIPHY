"""mapping / path membership: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0
from intraphy.coordinates import local_interval_to_genome
from intraphy.mapping.chain_types import ChainPathMembership
from intraphy.mapping.policies import occurrence_copy_key


def _copy_transcription_bounds(occurrences):
    bounds = {}
    grouped = defaultdict(list)
    for occurrence in occurrences:
        grouped[occurrence_copy_key(occurrence)].append(occurrence)
    for key, rows in grouped.items():
        bounds[key] = (
            min(int(row["start"]) for row in rows),
            max(int(row["end"]) for row in rows),
        )
    return bounds


def _candidate_copy_interval(occurrence, record, side, copy_bounds):
    genomic_blocks = tuple(record.get(f"{side}_genomic_blocks0", ()) or ())
    if genomic_blocks:
        genomic = Interval0(
            min(interval.start0 for interval in genomic_blocks),
            max(interval.end0 for interval in genomic_blocks),
        )
    else:
        start0 = record.get(f"{side}_start0")
        end0 = record.get(f"{side}_end0")
        if start0 in {None, "", "NA"} or end0 in {None, "", "NA"}:
            return None
        try:
            locus = ClosedInterval1(
                int(occurrence["start"]), int(occurrence["end"]),
            ).to_interval0()
            genomic = local_interval_to_genome(
                Interval0(int(start0), int(end0)), locus, occurrence.get("strand"),
            )
        except (KeyError, TypeError, ValueError):
            return None
    if genomic.length == 0:
        return None
    copy_start, copy_end = copy_bounds[occurrence_copy_key(occurrence)]
    if occurrence.get("strand") == "-":
        return Interval0(copy_end - genomic.end0, copy_end - genomic.start0)
    return Interval0(genomic.start0 - (copy_start - 1), genomic.end0 - (copy_start - 1))


def _chain_precedes(left, right):
    return left.query.end0 <= right.query.start0 and left.target.end0 <= right.target.start0


def _candidate_path_memberships(
    query,
    subject,
    transcript_paths,
    reverse=False,
    query_transcript_ids=None,
    subject_transcript_ids=None,
):
    def transcript_tokens(value):
        return {
            token
            for token in str(value or "").replace(",", ";").split(";")
            if token and token != "NA"
        }

    allowed_query_transcripts = transcript_tokens(query_transcript_ids)
    allowed_subject_transcripts = transcript_tokens(subject_transcript_ids)
    by_occurrence = defaultdict(list)
    for path in transcript_paths or ():
        by_occurrence[path.get("occurrence_id")].append(path)
    memberships = []
    for query_path in by_occurrence.get(query.get("occurrence_id"), ()):
        if (
            allowed_query_transcripts
            and query_path.get("transcript_id") not in allowed_query_transcripts
        ):
            continue
        for subject_path in by_occurrence.get(subject.get("occurrence_id"), ()):
            if (
                allowed_subject_transcripts
                and subject_path.get("transcript_id") not in allowed_subject_transcripts
            ):
                continue
            if query_path.get("transcript_id") in {None, "", "NA"}:
                continue
            if subject_path.get("transcript_id") in {None, "", "NA"}:
                continue
            try:
                query_order = int(query_path.get("path_rank", query_path.get("transcript_order")))
                subject_order = int(subject_path.get("path_rank", subject_path.get("transcript_order")))
            except (TypeError, ValueError):
                continue
            query_membership = {
                "path_id": "|".join(
                    str(query_path.get(field, query.get(field, "NA")))
                    for field in ("family_id", "species", "gene_copy_id", "transcript_id")
                ),
                "order": query_order,
                "parent_id": query_path.get("occurrence_id", query.get("occurrence_id", "")),
                "contig": query_path.get("contig", query.get("contig", "NA")),
                "strand": query_path.get("strand", query.get("strand", "NA")),
            }
            subject_membership = {
                "path_id": "|".join(
                    str(subject_path.get(field, subject.get(field, "NA")))
                    for field in ("family_id", "species", "gene_copy_id", "transcript_id")
                ),
                "order": subject_order,
                "parent_id": subject_path.get("occurrence_id", subject.get("occurrence_id", "")),
                "contig": subject_path.get("contig", subject.get("contig", "NA")),
                "strand": subject_path.get("strand", subject.get("strand", "NA")),
            }
            if reverse:
                query_membership, subject_membership = subject_membership, query_membership
            if (
                query_membership["strand"] not in {"+", "-"}
                or subject_membership["strand"] not in {"+", "-"}
                or query_membership["contig"] in {None, "", "NA"}
                or subject_membership["contig"] in {None, "", "NA"}
            ):
                continue
            memberships.append(
                ChainPathMembership(
                    query_path_id=query_membership["path_id"],
                    target_path_id=subject_membership["path_id"],
                    query_order=query_membership["order"],
                    target_order=subject_membership["order"],
                    query_contig=query_membership["contig"],
                    target_contig=subject_membership["contig"],
                    query_strand=query_membership["strand"],
                    target_strand=subject_membership["strand"],
                    query_parent_id=query_membership["parent_id"],
                    target_parent_id=subject_membership["parent_id"],
                )
            )
    return tuple(sorted(set(memberships), key=lambda item: (item.context, item.query_order, item.target_order)))
