"""mapping / reference coverage: explicit implementation ownership."""
from __future__ import annotations

from collections import Counter
from collections import defaultdict
from intraphy.mapping.member_intervals import _block_union_length
from intraphy.mapping.member_intervals import _blocks_overlap
from intraphy.mapping.member_intervals import _hard_position_eligible
from intraphy.mapping.member_intervals import _parse_genomic_blocks


def _assign_reference_coverage(rows, scored, occurrence_by_id):
    """Distinguish complementary fragments, aliases and repeated reference coverage."""
    occurrence_ids = {row["occurrence_id"] for row in rows}
    copy_counts = Counter((row.get("species"), row.get("gene_copy_id")) for row in rows)

    def reference_key(occurrence_id):
        occurrence = occurrence_by_id.get(occurrence_id, {})
        length = int(occurrence.get("end", 0)) - int(occurrence.get("start", 1)) + 1
        copies = copy_counts[(occurrence.get("species"), occurrence.get("gene_copy_id"))]
        direct_matches = sum(occurrence_id in {m.get("query_occurrence_id"), m.get("subject_occurrence_id")}
                             and _hard_position_eligible(m) for m in scored)
        return (-length, -direct_matches, copies, occurrence_id)

    reference_id = min(occurrence_ids, key=reference_key)
    reference = occurrence_by_id.get(reference_id, {})
    reference_length = int(reference.get("end", 0)) - int(reference.get("start", 1)) + 1
    reference_length = max(reference_length, max(
        (_block_union_length(row.get("_interval_records", ()))
         for row in rows if row["occurrence_id"] == reference_id), default=0))
    intervals_by_occurrence = defaultdict(list)
    for match in scored:
        query, subject = match.get("query_occurrence_id"), match.get("subject_occurrence_id")
        if reference_id not in {query, subject} or not _hard_position_eligible(match):
            continue
        member_id = subject if query == reference_id else query
        if member_id not in occurrence_ids:
            continue
        field = "query_genomic_matched_blocks" if query == reference_id else "subject_genomic_matched_blocks"
        genomic = _parse_genomic_blocks(match.get(field))
        if genomic:
            # All comparisons in this pass use coordinates of the same reference occurrence.
            intervals_by_occurrence[member_id].extend(genomic)
    reference_blocks = [block for row in rows if row["occurrence_id"] == reference_id
                        for block in row.get("_interval_records", ())]
    intervals_by_occurrence[reference_id].extend(reference_blocks)
    by_copy = defaultdict(list)
    for row in rows:
        by_copy[(row.get("species"), row.get("gene_copy_id"))].append(row)
    for copy_rows in by_copy.values():
        # Aliases sharing genomic coordinates are the same physical instance.
        instances = []
        for row in sorted(copy_rows, key=lambda x: x["occurrence_id"]):
            groups = [group for group in instances if any(
                _blocks_overlap(row.get("_interval_records", ()), other.get("_interval_records", ()))
                for other in group)]
            if not groups:
                instances.append([row])
            else:
                groups[0].append(row)
                for group in groups[1:]:
                    groups[0].extend(group)
                    instances.remove(group)
        projections = [[block for row in group
                        for block in intervals_by_occurrence[row["occurrence_id"]]]
                       for group in instances]
        covered = _block_union_length([block for group in projections for block in group])
        overlap = sum(_block_union_length(group) for group in projections) - covered
        if not covered:
            continue
        relation = ("repeated_overlap" if overlap > 0 else
                    "complementary_complete" if len(instances) > 1 and covered >= reference_length else
                    "complementary_partial" if len(instances) > 1 else
                    "reference_axis" if reference_id in {row["occurrence_id"] for row in copy_rows} else
                    "single_partial")
        for index, group in enumerate(instances, 1):
            for row in group:
                row.update(reference_occurrence_id=reference_id,
                           reference_length=reference_length,
                           reference_length_source="aligned_reference_coordinates",
                           reference_coverage_relation=relation,
                           reference_covered_bases=covered,
                           reference_overlap_bases=overlap,
                           reference_uncovered_bases=max(0, reference_length-covered))
                if relation == "repeated_overlap":
                    row["repeat_instance_id"] = f"{row['element_id']}.repeat_{index:03d}"
