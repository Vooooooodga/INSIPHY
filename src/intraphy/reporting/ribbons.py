"""reporting / ribbons: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from intraphy.reporting.inputs import visual_status
from intraphy.reporting.layout import _oriented_interval
from intraphy.reporting.layout import _row_box_style
from intraphy.reporting.layout import connector_pairs
from intraphy.reporting.svg import attr
from intraphy.reporting.svg import clipped_ribbon_box
from intraphy.reporting.svg import compact_label
from intraphy.reporting.svg import draw_context_span
from intraphy.reporting.svg import draw_homology_connector
from intraphy.reporting.svg import draw_segment_box
from intraphy.reporting.svg import feature_shape
from intraphy.reporting.svg import svg_text


def draw_correspondence_ribbons(body, element_boxes, styles, matches, membership_ranges=None, membership_evidence_pairs=None):
    membership_ranges = membership_ranges or {}
    # Membership coordinates describe how an established pairwise record maps
    # within an element.  They cannot establish a pair on their own: legacy
    # membership tables may contain shared element IDs for unrelated samples.
    membership_evidence_pairs = membership_evidence_pairs or set(matches)
    membership_connected = set()
    for element_id, boxes in sorted(element_boxes.items()):
        for left_box, right_box in connector_pairs(boxes):
            left_id, right_id = left_box["occurrence_id"], right_box["occurrence_id"]
            key = tuple(sorted((left_id, right_id)))
            match = matches.get(key)
            if match is not None:
                for qs, qe, ts, te in match["blocks"]:
                    if left_id != key[0]:
                        qs, qe, ts, te = ts, te, qs, qe
                    draw_homology_connector(
                        body,
                        clipped_ribbon_box(left_box, qs, qe),
                        clipped_ribbon_box(right_box, ts, te),
                        styles.get(element_id, {}),
                        element_id,
                        match,
                    )
                continue
            if key not in membership_evidence_pairs:
                continue
            left_ranges = membership_ranges.get((element_id, left_id), ())
            right_ranges = membership_ranges.get((element_id, right_id), ())
            left_by_block = defaultdict(list)
            right_by_block = defaultdict(list)
            for block_id, start, end in left_ranges:
                left_by_block[block_id].append((start, end))
            for block_id, start, end in right_ranges:
                right_by_block[block_id].append((start, end))
            shared_blocks = set(left_by_block) & set(right_by_block)
            for block_id in sorted(shared_blocks):
                left_spans = left_by_block[block_id]
                right_spans = right_by_block[block_id]
                if len(left_spans) != len(right_spans):
                    continue
                for (left_start, left_end), (right_start, right_end) in zip(left_spans, right_spans):
                    membership_match = {
                        "match_id": f"membership:{element_id}:{block_id}",
                        "alignment_backend": "element_membership",
                        "correspondence_basis": "element_membership",
                        "projection_field": "matched_blocks",
                        "block_source": "membership_blocks",
                        "connector_kind": "membership_correspondence",
                    }
                    draw_homology_connector(
                        body,
                        clipped_ribbon_box(left_box, left_start, left_end),
                        clipped_ribbon_box(right_box, right_start, right_end),
                        styles.get(element_id, {}),
                        element_id,
                        membership_match,
                    )
                    membership_connected.add(key)
    return membership_connected


def _draw_zoom_panels(body, panels, track_left, track_right, occ_to_element, styles, class_by_occ, status_by_occ, unresolved, matched_occurrences, memberships_by_occurrence, colored_ranges, start_y, lane_h):
    cursor = start_y
    panel_width = track_right - track_left
    for group, window_start, window_end, target_ids in panels:
        visible_lanes = []
        for lane_id, rows in group["lanes"]:
            overlaps = []
            for row in rows:
                start, end = _oriented_interval(group, row)
                if start < window_end and end > window_start:
                    overlaps.append(row)
            if overlaps:
                visible_lanes.append((lane_id, overlaps))
        panel_height = 34 + max(1, len(visible_lanes)) * 22
        span = max(1.0, window_end - window_start)
        body.append(f'<g data-zoom-panel="{attr(group["key"][0] + ":" + group["key"][1])}" data-zoom-start="{window_start:.0f}" data-zoom-end="{window_end:.0f}">')
        body.append(f'<rect x="{track_left:.2f}" y="{cursor:.2f}" width="{panel_width:.2f}" height="{panel_height:.2f}" fill="#FFFFFF" stroke="#777" stroke-width="0.8"/>')
        label = f'Local detail, magnified: {group["key"][0]} {group["key"][1]}, gene-oriented interval {window_start:.0f}-{window_end:.0f} bp (5\u2032 to 3\u2032)'
        body.append(svg_text(track_left + 6, cursor + 14, compact_label(label, 125), 9, weight="bold", fill="#444"))
        for lane_idx, (lane_id, rows) in enumerate(visible_lanes):
            y = cursor + 22 + lane_idx * 22
            body.append(svg_text(track_left - 6, y + 9, compact_label(lane_id, 20), 7, anchor="end", fill="#555"))
            body.append(f'<line x1="{track_left}" y1="{y + 6}" x2="{track_right}" y2="{y + 6}" stroke="#D0D0D0" stroke-width="0.7"/>')
            for row in rows:
                start, end = _oriented_interval(group, row)
                clipped_start, clipped_end = max(start, window_start), min(end, window_end)
                if clipped_end <= clipped_start:
                    continue
                x = track_left + (clipped_start - window_start) / span * panel_width
                width = (clipped_end - clipped_start) / span * panel_width
                occurrence_id = row.get("occurrence_id", "NA")
                element_id = occ_to_element.get(occurrence_id, "NA")
                unit_class, role, style = _row_box_style(
                    row, occurrence_id, element_id, styles, class_by_occ, status_by_occ,
                    "matched_blocks" if occurrence_id in matched_occurrences else (
                        "actual_blocks_missing" if occurrence_id in unresolved else "not_assessed"
                    ),
                )
                state = visual_status(row, status_by_occ.get(occurrence_id, {}), unit_class)
                shape = "candidate" if lane_id == "candidate evidence" else style["shape"]
                style["shape"] = shape
                if shape == "intron" and lane_id != "candidate evidence":
                    draw_context_span(body, x, y, width, role, occurrence_id)
                    continue
                hbox = 11 if shape in {"utr", "noncoding", "candidate"} else 14
                base_style = {
                    "fill": "#D8D8D8", "stroke": "#777777", "stroke_dasharray": "none",
                    "encoding": "neutral", "display_role": role, "shape": shape,
                    "match_state": style.get("match_state", "not_assessed"),
                }
                if lane_id == "candidate evidence" or unit_class == "candidate_source":
                    base_style = style
                draw_segment_box(body, x, y + 2, width, hbox, "NA", base_style, unit_class, state, occurrence_id, lane_id)
                feature_start, _feature_end = _oriented_interval(group, row)
                for membership in memberships_by_occurrence.get(occurrence_id, []):
                    member_element = membership.get("element_id")
                    for block_start, block_end, source in colored_ranges.get((member_element, occurrence_id), []):
                        block_gene_start = feature_start + block_start - 1
                        block_gene_end = feature_start + block_end
                        local_start = max(block_gene_start, window_start)
                        local_end = min(block_gene_end, window_end)
                        if local_end <= local_start:
                            continue
                        overlay_x = track_left + (local_start - window_start) / span * panel_width
                        overlay_width = (local_end - local_start) / span * panel_width
                        overlay_style = dict(styles[member_element])
                        overlay_style.update({"display_role": role, "shape": feature_shape(role), "match_state": source})
                        draw_segment_box(
                            body, overlay_x, y + 2, overlay_width, hbox, member_element, overlay_style,
                            "exon_like", state, occurrence_id, lane_id,
                            block_start, block_end,
                        )
                if occurrence_id in target_ids:
                    body.append(svg_text(x + width + 2, y + 9, compact_label(occurrence_id, 18), 7, fill="#333"))
        body.append("</g>")
        cursor += panel_height + 8
    return cursor
