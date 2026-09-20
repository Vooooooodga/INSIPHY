"""reporting / layout: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.reporting.inputs import biological_role
from insiphy.reporting.inputs import visual_status
from insiphy.reporting.svg import feature_shape
from insiphy.reporting.svg import tip_label_for_group
from insiphy.reporting.svg import to_int
from insiphy.storage.tabular import read_tsv
from pathlib import Path


def transcript_lanes_for_group(rows, path_rows, excluded_occurrences=()):
    excluded_occurrences = set(excluded_occurrences)
    rows_by_occ = {row.get("occurrence_id"): row for row in rows if row.get("occurrence_id")}
    lanes = defaultdict(list)
    seen = set()
    for path in sorted(path_rows, key=lambda item: (item.get("transcript_id", ""), to_int(item.get("path_rank")), item.get("occurrence_id", ""))):
        occ_id = path.get("occurrence_id")
        if not occ_id or occ_id in excluded_occurrences:
            continue
        merged = dict(rows_by_occ.get(occ_id, {}))
        for key, value in path.items():
            if value not in {"", "NA"}:
                merged[key] = value
        lane_id = merged.get("transcript_id") or "observed_segments"
        key = (lane_id, occ_id)
        if key in seen:
            continue
        seen.add(key)
        lanes[lane_id].append(merged)
    used_occ = {occ for _lane, occ in seen}
    leftovers = [row for row in rows if row.get("occurrence_id") not in used_occ]
    if leftovers:
        lanes["unassigned_features" if lanes else "observed_segments"].extend(leftovers)
    if not lanes:
        lanes["observed_segments"] = list(rows)
    return [
        (lane_id, sorted(lane_rows, key=lambda row: (to_int(row.get("transcript_order", row.get("path_rank"))), to_int(row.get("start")), to_int(row.get("end")))))
        for lane_id, lane_rows in sorted(lanes.items(), key=lambda item: item[0])
    ]


def track_groups(input_dir, occurrences, candidate_occurrences=()):
    candidate_occurrences = set(candidate_occurrences)
    grouped = defaultdict(list)
    for row in occurrences:
        grouped[(row.get("species", "NA"), row.get("gene_copy_id", "NA"))].append(row)
    paths_by_group = defaultdict(list)
    for row in read_tsv(Path(input_dir) / "transcript_paths.tsv", optional=True):
        paths_by_group[(row.get("species", "NA"), row.get("gene_copy_id", "NA"))].append(row)
    groups = []
    for key, rows in sorted(grouped.items()):
        starts = [to_int(row.get("start")) for row in rows]
        ends = [to_int(row.get("end")) for row in rows]
        start = min(starts or [0])
        end = max(ends or [start + 1])
        candidate_rows = [row for row in rows if row.get("occurrence_id") in candidate_occurrences]
        observed_rows = [row for row in rows if row.get("occurrence_id") not in candidate_occurrences]
        lanes = transcript_lanes_for_group(observed_rows, paths_by_group.get(key, []), candidate_occurrences)
        if candidate_rows:
            lanes.append(("candidate evidence", sorted(candidate_rows, key=lambda row: (
                to_int(row.get("transcript_order")), to_int(row.get("start")), row.get("occurrence_id", "")
            ))))
        groups.append(
            {
                "key": key,
                "rows": rows,
                "lanes": lanes,
                "start": start,
                "end": end,
                "negative": rows and rows[0].get("strand") == "-",
            }
        )
    return groups


def row_geometry(group, row, left, right_edge):
    start = group["start"]
    end = group["end"]
    span = max(1, end - start + 1)
    segment_start = to_int(row.get("start"))
    segment_end = to_int(row.get("end"), segment_start)
    relative_start = end - segment_end if group["negative"] else segment_start - start
    width = right_edge - left
    x = left + relative_start / span * width
    w = (segment_end - segment_start + 1) / span * width
    return x, w


def connector_pairs(boxes):
    by_group = defaultdict(list)
    for box in sorted(boxes, key=lambda item: (item["row"], item["lane"], item["x"])):
        by_group[box["group"]].append(box)
    groups = list(by_group.values())
    for upper, lower in zip(groups, groups[1:]):
        if upper[0]["species"] == lower[0]["species"]:
            continue
        # Shared isoforms use the facing lanes; distinct split members remain visible.
        outgoing = {box["occurrence_id"]: box for box in upper}
        incoming = {}
        for box in lower:
            incoming.setdefault(box["occurrence_id"], box)
        for left_box in outgoing.values():
            for right_box in incoming.values():
                yield left_box, right_box


def unresolved_without_membership(unresolved, membership_pairs):
    result = {}
    for occurrence_id, records in unresolved.items():
        remaining = [record for record in records
                     if tuple(sorted((occurrence_id, record[2]))) not in membership_pairs]
        if remaining:
            result[occurrence_id] = remaining
    return result


def _candidate_occurrences(occurrences, class_by_occ, status_by_occ):
    candidate_ids = set()
    for row in occurrences:
        occurrence_id = row.get("occurrence_id")
        unit_class = class_by_occ.get(occurrence_id, "context")
        status = visual_status(row, status_by_occ.get(occurrence_id), unit_class)
        if unit_class == "candidate_source" or status == "predicted":
            candidate_ids.add(occurrence_id)
    return candidate_ids


def _species_track_layout(groups, tree, top, lane_h, group_gap, species_gap):
    groups_by_tip = defaultdict(list)
    for group in groups:
        species, copy = group["key"]
        tip = tip_label_for_group(tree, species, copy) if tree else species
        groups_by_tip[tip].append(group)
    if tree:
        leaf_nodes = [node for node in tree.preorder() if node in tree.leaves]
        ordered_tips = [tree.label[node] for node in leaf_nodes]
    else:
        ordered_tips = sorted(groups_by_tip)
    ordered_tips.extend(sorted(set(groups_by_tip) - set(ordered_tips)))
    group_top = {}
    species_y = {}
    cursor = top
    for tip in ordered_tips:
        tip_groups = sorted(groups_by_tip.get(tip, []), key=lambda group: group["key"])
        block_top = cursor
        if not tip_groups:
            cursor += lane_h
        else:
            for group in tip_groups:
                group_top[group["key"]] = cursor
                cursor += 22 + max(1, len(group["lanes"])) * lane_h + group_gap
            cursor -= group_gap
        species_y[tip] = (block_top + cursor) / 2
        cursor += species_gap
    return [group for tip in ordered_tips for group in sorted(groups_by_tip.get(tip, []), key=lambda item: item["key"])], group_top, species_y, cursor


def _tree_y_coordinates(tree, species_y):
    y = {}
    for leaf in tree.leaves:
        y[leaf] = species_y.get(tree.label[leaf], 80.0)

    def assign(node):
        if node in y:
            return y[node]
        child_values = [assign(child) for child in tree.children.get(node, [])]
        y[node] = sum(child_values) / max(1, len(child_values))
        return y[node]

    assign(tree.root)
    return y


def _zoom_windows(group, lanes, left, right, occ_to_element, class_by_occ, status_by_occ):
    span = max(1, group["end"] - group["start"] + 1)
    targets = []
    for lane_id, rows in lanes:
        for row in rows:
            occurrence_id = row.get("occurrence_id", "NA")
            element_id = occ_to_element.get(occurrence_id, "NA")
            unit_class = class_by_occ.get(occurrence_id, "context")
            role = biological_role(row, status_by_occ.get(occurrence_id))
            shape = feature_shape(role)
            if shape == "intron" and lane_id != "candidate evidence":
                continue
            x, width = row_geometry(group, row, left, right)
            if width >= 8:
                continue
            start, end = _oriented_interval(group, row)
            if end <= start:
                continue
            window_span = min(span, max(1000, (end - start) * 12))
            if window_span >= span:
                continue
            center = (start + end) / 2
            window_start = max(0, min(span - window_span, center - window_span / 2))
            targets.append((window_start, window_start + window_span, {occurrence_id}))
    merged = []
    for start, end, ids in sorted(targets, key=lambda item: (item[0], item[1])):
        if merged and start <= merged[-1][1]:
            previous = merged[-1]
            previous[1] = max(previous[1], end)
            previous[2].update(ids)
        else:
            merged.append([start, end, set(ids)])
    return [(start, end, ids) for start, end, ids in merged]


def _oriented_interval(group, row):
    segment_start = to_int(row.get("start"), group["start"])
    segment_end = to_int(row.get("end"), segment_start)
    if group["negative"]:
        return group["end"] - segment_end, group["end"] - segment_start + 1
    return segment_start - group["start"], segment_end - group["start"] + 1


def _row_box_style(row, occurrence_id, element_id, styles, class_by_occ, status_by_occ, match_state="not_assessed"):
    role = biological_role(row, status_by_occ.get(occurrence_id))
    unit_class = class_by_occ.get(occurrence_id, "context")
    style = dict(styles.get(element_id, {"fill": "#D0D0D0", "stroke": "#777777", "stroke_dasharray": "none", "encoding": "color"}))
    style["display_role"] = role
    style["shape"] = "candidate" if unit_class == "candidate_source" else feature_shape(role)
    style["match_state"] = match_state
    return unit_class, role, style
