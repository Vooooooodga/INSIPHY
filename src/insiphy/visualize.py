"""Lightweight SVG visualizations for INSIPHY results."""

from __future__ import annotations

import math
import re
import textwrap
from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape

from .elements import EXON_LIKE_ROLES, collect_element_profiles, element_class_for_occurrence
from .io import read_tsv, write_tsv
from .tree import SpeciesTree


PALETTE = [
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#56B4E9",
    "#F0E442",
    "#999999",
]

PATTERNS = ["diagonal", "dots", "cross", "horizontal", "vertical", "grid", "sparse", "solid"]


def attr(value):
    return escape(str(value), {'"': "&quot;"})


def to_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def compact_label(value, max_len=34):
    text = str(value or "NA")
    if len(text) <= max_len:
        return text
    keep_left = (max_len - 3) // 2
    keep_right = max_len - 3 - keep_left
    return f"{text[:keep_left]}...{text[-keep_right:]}"


def svg_text(x, y, text, size=11, anchor="start", weight="normal", fill="#222"):
    return f'<text x="{x}" y="{y}" font-family="Arial, sans-serif" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}">{escape(str(text))}</text>'


def structural_tree_rows(input_dir, result_dir):
    input_dir = Path(input_dir)
    result_dir = Path(result_dir)
    tree_file = "species_tree.tsv"
    scope_rows = read_tsv(result_dir / "phylogeny_scope.tsv", optional=True)
    for row in scope_rows:
        if row.get("scope") in {"structural_characters", "single_copy_structural_sites"} and row.get("tree_file"):
            tree_file = row["tree_file"]
            break
    rows = read_tsv(input_dir / tree_file, optional=True)
    if rows:
        return rows, tree_file
    for fallback in ["copy_tree.tsv", "gene_tree.tsv", "species_tree.tsv"]:
        rows = read_tsv(input_dir / fallback, optional=True)
        if rows:
            return rows, fallback
    return [], tree_file


def finite_probability(row):
    try:
        value = float(row.get("endpoint_change_probability", "NA"))
    except (TypeError, ValueError):
        return False
    return math.isfinite(value) and 0.0 <= value <= 1.0


def valid_probability_change(row):
    if row.get("structural_change_type") == "posterior_not_reported":
        return False
    conditioning = {part.strip() for part in row.get("conditioning", "").split(";")}
    if "fit_status=success" not in conditioning:
        return False
    return finite_probability(row)


def phylogenetic_change_rows(result_dir):
    result_dir = Path(result_dir)
    if (result_dir / "branch_structural_events.tsv").exists():
        return [
            row for row in read_tsv(result_dir / "branch_structural_events.tsv")
            if row.get("placement_status") in {"required", "possible"}
        ]
    if (result_dir / "structural_changes.tsv").exists():
        rows = read_tsv(result_dir / "structural_changes.tsv")
        rows = [row for row in rows if valid_probability_change(row)]
        return sorted(
            rows,
            key=change_probability,
            reverse=True,
        )
    return read_tsv(result_dir / "event_support_summary.tsv", optional=True)


def change_probability(row):
    for field in (
        "endpoint_change_probability",
        "ctmc_change_probability",
        "stochastic_pr_any_change",
    ):
        try:
            value = float(row.get(field, "NA"))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value) and 0.0 <= value <= 1.0:
            return value
    return 0.0


def change_label(row):
    return row.get("event_type") or row.get("structural_change_type") or row.get("event_class") or "structural_change"


def change_priority(row):
    return ({"required": 2, "possible": 1}.get(row.get("placement_status"), 0), change_probability(row))


def change_symbol(row, x, y, radius=7.0):
    status = row.get("placement_status")
    color = "#D55E00" if any(word in change_label(row) for word in ("loss", "fusion")) else "#0072B2"
    if status:
        fill, opacity = (color if status == "required" else "white"), 1.0
    else:
        opacity = change_probability(row)
        radius *= math.sqrt(opacity)
        fill = color
    title = f"{row.get('site_id', '')}: {change_label(row)}; {status or 'conditional probability'}"
    return (
        f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="{fill}" '
        f'fill-opacity="{opacity:.3f}" stroke="{color}" stroke-width="1.3">'
        f'<title>{escape(title)}</title></circle>'
    )


def tip_label_for_group(tree, species, copy):
    candidates = [f"{species}:{copy}", species, copy]
    for label in candidates:
        if label in tree.leaf_by_label:
            return label
    return candidates[0]


def pattern_defs(styles):
    body = ["<defs>"]
    for _element_id, style in styles.items():
        if style.get("encoding") != "pattern":
            continue
        pid = style["pattern_id"]
        color = style["stroke"]
        kind = style["pattern"]
        body.append(f'<pattern id="{pid}" patternUnits="userSpaceOnUse" width="8" height="8">')
        body.append('<rect width="8" height="8" fill="white"/>')
        if kind == "diagonal":
            body.append(f'<path d="M-2,8 L8,-2 M2,10 L10,2" stroke="{color}" stroke-width="1.4"/>')
        elif kind == "dots":
            body.append(f'<circle cx="2" cy="2" r="1.1" fill="{color}"/><circle cx="6" cy="6" r="1.1" fill="{color}"/>')
        elif kind == "cross":
            body.append(f'<path d="M0,0 L8,8 M8,0 L0,8" stroke="{color}" stroke-width="1.1"/>')
        elif kind == "horizontal":
            body.append(f'<path d="M0,2 L8,2 M0,6 L8,6" stroke="{color}" stroke-width="1.3"/>')
        elif kind == "vertical":
            body.append(f'<path d="M2,0 L2,8 M6,0 L6,8" stroke="{color}" stroke-width="1.3"/>')
        elif kind == "grid":
            body.append(f'<path d="M0,4 L8,4 M4,0 L4,8" stroke="{color}" stroke-width="1.1"/>')
        elif kind == "sparse":
            body.append(f'<circle cx="4" cy="4" r="1.6" fill="none" stroke="{color}" stroke-width="1.1"/>')
        else:
            body.append(f'<rect width="8" height="8" fill="{color}" fill-opacity="0.28"/>')
        body.append("</pattern>")
    body.append("</defs>")
    return "\n".join(body)


def fallback_element_correspondence(input_dir, occurrences):
    homology = read_tsv(Path(input_dir) / "segment_homology.tsv", optional=True)
    evidence = read_tsv(Path(input_dir) / "sequence_synteny_evidence.tsv", optional=True)
    profiles, element_by_homology = collect_element_profiles(homology, occurrences, evidence)
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    rows = []
    for row in homology:
        element_id = element_by_homology.get(row.get("homology_id", ""))
        if not element_id:
            continue
        occ = occ_by_id.get(row.get("occurrence_id"), {})
        rows.append(
            {
                "element_id": element_id,
                "homology_id": row.get("homology_id", "NA"),
                "occurrence_id": row.get("occurrence_id", "NA"),
                "element_class": element_class_for_occurrence(occ, row.get("homology_id", ""), profiles, element_by_homology),
            }
        )
    return rows


def segment_styles(input_dir, result_dir, occurrences=None, correspondence_encoding="color"):
    occurrences = occurrences or []
    element_rows = read_tsv(Path(result_dir) / "element_correspondence.tsv", optional=True)
    if not element_rows:
        element_rows = fallback_element_correspondence(input_dir, occurrences)
    by_occ = {}
    class_by_occ = {}
    status_by_occ = {}
    styled_elements = set()
    for row in element_rows:
        element_class = row.get("element_class", "context")
        if element_class in {"context", "absent"}:
            continue
        by_occ[row["occurrence_id"]] = row["element_id"]
        class_by_occ[row["occurrence_id"]] = element_class
        status_by_occ[row["occurrence_id"]] = row
        if element_class == "exon_like":
            styled_elements.add(row["element_id"])
    styles = {}
    for idx, element_id in enumerate(sorted(styled_elements)):
        dash = "none" if idx % 3 == 0 else "4 2" if idx % 3 == 1 else "1.5 2"
        if correspondence_encoding == "color":
            styles[element_id] = {
                "encoding": "color",
                "fill": PALETTE[idx % len(PALETTE)],
                "stroke": "#111111",
                "pattern": "solid",
                "pattern_id": f"eg_{idx + 1}",
                "stroke_dasharray": "none",
            }
        else:
            styles[element_id] = {
                "encoding": "pattern",
                "fill": f"url(#eg_{idx + 1})",
                "stroke": "#111111",
                "pattern": PATTERNS[idx % len(PATTERNS)],
                "pattern_id": f"eg_{idx + 1}",
                "stroke_dasharray": dash,
            }
    return by_occ, styles, class_by_occ, status_by_occ


def visual_status(row, status_row, unit_class):
    joined = " ".join(
        str(value)
        for value in [
            row.get("role"),
            row.get("presence_status"),
            row.get("boundary_state"),
            row.get("evidence"),
            row.get("path_status"),
            status_row.get("support_type") if status_row else "",
            status_row.get("membership_call") if status_row else "",
            status_row.get("inferred_role") if status_row else "",
        ]
        if value
    ).lower()
    if unit_class == "candidate_source":
        return "sequence_candidate"
    if (
        "predicted_exon_candidate" in joined
        or "predicted_cds" in joined
        or "supports_hidden_segment" in joined
        or "hidden" in joined
        or "predicted" in joined
        or "inferred" in joined
    ):
        return "predicted"
    if "unknown" in joined or "ambiguous" in joined:
        return "unknown"
    if unit_class == "exon_like":
        return "annotated" if row.get("role") in EXON_LIKE_ROLES else "unknown"
    return "context"


def transcript_lanes_for_group(rows, path_rows):
    rows_by_occ = {row.get("occurrence_id"): row for row in rows if row.get("occurrence_id")}
    lanes = defaultdict(list)
    seen = set()
    for path in sorted(path_rows, key=lambda item: (item.get("transcript_id", ""), to_int(item.get("path_rank")), item.get("occurrence_id", ""))):
        occ_id = path.get("occurrence_id")
        if not occ_id:
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


def track_groups(input_dir, occurrences):
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
        groups.append(
            {
                "key": key,
                "rows": rows,
                "lanes": transcript_lanes_for_group(rows, paths_by_group.get(key, [])),
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
    w = max(4, (segment_end - segment_start + 1) / span * width)
    return x, w


def draw_context_span(body, x, y, w, role):
    if role == "intron":
        body.append(f'<line x1="{x:.2f}" y1="{y + 11}" x2="{x + w:.2f}" y2="{y + 11}" stroke="#B8B8B8" stroke-width="1.1" stroke-dasharray="2 2"/>')
        body.append(f'<line x1="{x:.2f}" y1="{y + 7}" x2="{x:.2f}" y2="{y + 15}" stroke="#B8B8B8" stroke-width="0.8"/>')
        body.append(f'<line x1="{x + w:.2f}" y1="{y + 7}" x2="{x + w:.2f}" y2="{y + 15}" stroke="#B8B8B8" stroke-width="0.8"/>')
    else:
        body.append(f'<rect x="{x:.2f}" y="{y + 9}" width="{w:.2f}" height="4" rx="1" fill="#EEEEEE" stroke="#B8B8B8" stroke-width="0.5"/>')


def accepted_ribbon_matches(input_dir, result_dir, occurrences):
    path = Path(result_dir) / "segment_matches.tsv"
    if not path.exists():
        path = Path(input_dir) / "segment_matches.tsv"
    lengths = {
        row["occurrence_id"]: int(row["end"]) - int(row["start"]) + 1
        for row in occurrences
    }
    matches = {}
    conflicts = set()
    for row in read_tsv(path, optional=True):
        if row.get("match_status") != "mapped":
            continue
        query = row.get("query_occurrence_id")
        target = row.get("subject_occurrence_id")
        if query not in lengths or target not in lengths or query == target:
            continue
        key = tuple(sorted((query, target)))
        basis = row.get("correspondence_basis") or "DNA"
        coding_projection = basis == "annotated_CDS_protein"
        projection_field = "protein_projected_blocks" if coding_projection else "projected_reference_blocks"
        # Coding blocks follow transcript orientation; the DNA alignment remains separate.
        strand = "+" if coding_projection else row.get("alignment_strand")
        if strand != "+":
            conflicts.add(key)
            continue
        reference = row.get("projected_reference_occurrence_id", target)
        if reference not in {"", "NA", target}:
            conflicts.add(key)
            continue
        blocks = []
        for token in row.get(projection_field, "").split(";"):
            parsed = re.fullmatch(r"(\d+)-(\d+):(\d+)-(\d+)", token)
            if parsed is None:
                break
            qs, qe, ts, te = map(int, parsed.groups())
            if not (1 <= qs <= qe <= lengths[query] and 1 <= ts <= te <= lengths[target]):
                break
            if qe - qs != te - ts:
                break
            if blocks and (qs <= blocks[-1][1] or ts <= blocks[-1][3]):
                break
            blocks.append((qs, qe, ts, te))
        else:
            if query != key[0]:
                blocks = [(ts, te, qs, qe) for qs, qe, ts, te in blocks]
            blocks = tuple(blocks)
            if key in matches and matches[key]["blocks"] != blocks:
                conflicts.add(key)
            elif key not in matches:
                matches[key] = {
                    "blocks": blocks,
                    "match_id": row.get("match_id", "NA"),
                    "alignment_backend": row.get("alignment_backend", "NA"),
                    "correspondence_basis": basis,
                    "projection_field": projection_field,
                }
            continue
        conflicts.add(key)
    return {key: value for key, value in matches.items() if key not in conflicts}


def clipped_ribbon_box(box, start, end):
    # Boxes are already oriented 5-prime to 3-prime on both genomic strands.
    scale = box["w"] / box["sequence_length"]
    return {
        **box,
        "x": box["x"] + (start - 1) * scale,
        "w": (end - start + 1) * scale,
        "local_start": start,
        "local_end": end,
    }


def draw_homology_connector(body, left_box, right_box, style, element_id, match):
    x1 = left_box["x"] + left_box["w"] / 2
    x2 = right_box["x"] + right_box["w"] / 2
    y1 = left_box["y"]
    y2 = right_box["y"]
    mid = (y1 + y2) / 2
    dash = style.get("stroke_dasharray", "none")
    stroke = "#555555"
    provenance = (
        f'data-match-id="{attr(match["match_id"])}" '
        f'data-alignment-backend="{attr(match["alignment_backend"])}" '
        f'data-correspondence-basis="{attr(match["correspondence_basis"])}" '
        f'data-projection-field="{attr(match["projection_field"])}" '
        f'data-source-occurrence="{attr(left_box["occurrence_id"])}" '
        f'data-target-occurrence="{attr(right_box["occurrence_id"])}" '
        f'data-source-start="{left_box["local_start"]}" data-source-end="{left_box["local_end"]}" '
        f'data-target-start="{right_box["local_start"]}" data-target-end="{right_box["local_end"]}"'
    )
    if left_box["class"] == "exon_like" and right_box["class"] == "exon_like":
        x1a = left_box["x"]
        x1b = left_box["x"] + left_box["w"]
        x2a = right_box["x"]
        x2b = right_box["x"] + right_box["w"]
        body.append(
            f'<path d="M{x1a:.2f},{y1:.2f} C{x1a:.2f},{mid:.2f} {x2a:.2f},{mid:.2f} {x2a:.2f},{y2:.2f} '
            f'L{x2b:.2f},{y2:.2f} C{x2b:.2f},{mid:.2f} {x1b:.2f},{mid:.2f} {x1b:.2f},{y1:.2f} Z" '
            f'fill="#777777" fill-opacity="0.10" stroke="none" data-connector="exon_correspondence" data-element-id="{attr(element_id)}" {provenance}/>'
        )
    body.append(f'<path d="M{x1:.2f},{y1:.2f} C{x1:.2f},{mid:.2f} {x2:.2f},{mid:.2f} {x2:.2f},{y2:.2f}" fill="none" stroke="{stroke}" stroke-width="0.85" stroke-opacity="0.55" stroke-dasharray="{dash}" data-connector="exon_correspondence_line" data-element-id="{attr(element_id)}" {provenance}/>')


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


def draw_correspondence_ribbons(body, element_boxes, styles, matches):
    for element_id, boxes in sorted(element_boxes.items()):
        for left_box, right_box in connector_pairs(boxes):
            left_id, right_id = left_box["occurrence_id"], right_box["occurrence_id"]
            key = tuple(sorted((left_id, right_id)))
            match = matches.get(key)
            if match is None:
                continue
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


def draw_segment_box(body, x, ybox, w, hbox, label, style, unit_class="exon_like", visual_state="annotated", occurrence_id="NA", lane_id="NA"):
    fill = style.get("fill", "#E6E6E6")
    dash = style.get("stroke_dasharray", "none")
    stroke = style.get("stroke", "#111111")
    if visual_state == "predicted":
        fill = "#FFFFFF"
        dash = "3 2"
    elif visual_state == "unknown":
        fill = "#F2F2F2"
        stroke = "#777777"
        dash = "1 2"
    if unit_class == "candidate_source":
        fill = "#FAFAFA"
        stroke = "#666666"
        dash = "2 2"
    body.append(
        f'<rect x="{x:.2f}" y="{ybox}" width="{w:.2f}" height="{hbox}" rx="2" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="0.8" stroke-dasharray="{dash}" '
        f'data-element-id="{attr(label)}" data-occurrence-id="{attr(occurrence_id)}" '
        f'data-unit-class="{attr(unit_class)}" data-visual-status="{attr(visual_state)}" '
        f'data-lane-id="{attr(lane_id)}">'
        f'<title>{escape(str(label))}; {escape(str(unit_class))}; {escape(str(visual_state))}</title></rect>'
    )
    if unit_class == "exon_like" and w > 28:
        label_length = max(4, int(w / 5))
        body.append(svg_text(x + w / 2, ybox + hbox + 11, compact_label(label, label_length), 8, anchor="middle", fill="#333"))


def draw_synteny(input_dir, result_dir, output_dir, correspondence_encoding="color"):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", optional=True)
    occ_to_element, styles, class_by_occ, status_by_occ = segment_styles(input_dir, result_dir, occurrences, correspondence_encoding)
    matches = accepted_ribbon_matches(input_dir, result_dir, occurrences)
    groups = track_groups(input_dir, occurrences)
    width = 1200
    left = 220
    right = 40
    lane_h = 36
    group_gap = 12
    top = 50
    height = top + max(1, sum(22 + len(group["lanes"]) * lane_h + group_gap for group in groups)) + 90
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        pattern_defs(styles),
        svg_text(24, 28, "INSIPHY intragenic synteny", 16, weight="bold"),
    ]
    element_boxes = defaultdict(list)
    pending_boxes = []
    cursor = top
    for ridx, group in enumerate(groups, start=0):
        species, copy = group["key"]
        body.append(svg_text(18, cursor + 15, compact_label(f"{species}  {copy}", 38), 11))
        for lane_idx, (lane_id, rows) in enumerate(group["lanes"]):
            y = cursor + 22 + lane_idx * lane_h
            body.append(svg_text(left - 10, y + 14, compact_label(lane_id, 30), 8, anchor="end", fill="#555"))
            body.append(f'<line x1="{left}" y1="{y + 11}" x2="{width - right}" y2="{y + 11}" stroke="#C8C8C8" stroke-width="1"/>')
            for row in rows:
                x, w = row_geometry(group, row, left, width - right)
                occ_id = row.get("occurrence_id", "NA")
                element_id = occ_to_element.get(occ_id, "NA")
                style = styles.get(element_id, {"pattern_id": "missing", "stroke_dasharray": "none", "stroke": "#999999"})
                role = row.get("role", "segment")
                unit_class = class_by_occ.get(occ_id, "context") if element_id != "NA" else "context"
                state = visual_status(row, status_by_occ.get(occ_id, {}), unit_class)
                if unit_class == "context":
                    draw_context_span(body, x, y, w, role)
                    continue
                ybox = y + 3 if unit_class == "exon_like" else y + 8
                hbox = 16 if unit_class == "exon_like" else 8
                box = {
                    "x": x,
                    "w": w,
                    "y": ybox + hbox / 2,
                    "class": unit_class,
                    "row": ridx,
                    "lane": lane_idx,
                    "group": group["key"],
                    "species": species,
                    "occurrence_id": occ_id,
                    "sequence_length": int(row["end"]) - int(row["start"]) + 1,
                }
                if unit_class == "exon_like" and state == "annotated":
                    element_boxes[element_id].append(box)
                pending_boxes.append((x, ybox, w, hbox, element_id, style, unit_class, state, occ_id, lane_id))
        cursor += 22 + len(group["lanes"]) * lane_h + group_gap
    draw_correspondence_ribbons(body, element_boxes, styles, matches)
    for args in pending_boxes:
        draw_segment_box(body, *args)
    legend_y = height - 58
    if correspondence_encoding == "color":
        legend = "Color marks exon-group membership; ribbons show directly aligned bases; pale dashed boxes mark predicted or uncertain roles."
    else:
        legend = "Texture marks exon-group membership; ribbons show directly aligned bases; gray spans mark introns/context."
    body.append(svg_text(24, legend_y, legend, 10, fill="#555"))
    body.append("</svg>")
    path = output_dir / "intragenic_synteny.svg"
    path.write_text("\n".join(body))
    return path


def tree_layout_depths(tree):
    unit_lengths = any(tree.length[child] is None for _parent, child in tree.edges())
    depth = {tree.root: 0.0}
    for node in tree.preorder():
        for child in tree.children.get(node, []):
            length = 1.0 if unit_lengths else tree.branch_length(child)
            depth[child] = depth[node] + length
    label = "unit branches (missing lengths)" if unit_lengths else "supplied branch lengths"
    return depth, label


def tree_coordinates(tree, depth=None):
    leaves = [node for node in tree.preorder() if node in tree.leaves]
    y = {node: 60 + idx * 42 for idx, node in enumerate(leaves)}

    def assign_y(node):
        if node in y:
            return y[node]
        vals = [assign_y(child) for child in tree.children.get(node, [])]
        y[node] = sum(vals) / max(1, len(vals))
        return y[node]

    assign_y(tree.root)
    if depth is None:
        depth, _layout = tree_layout_depths(tree)
    max_depth = max(depth.values() or [1.0])
    x = {node: 120 + depth[node] / max(max_depth, 1e-9) * 500 for node in depth}
    return x, y


def draw_phylogeny(input_dir, result_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    tree_rows, tree_file = structural_tree_rows(input_dir, result_dir)
    events = phylogenetic_change_rows(result_dir)
    if not tree_rows:
        path = output_dir / "phylogenetic_event_map.svg"
        path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="120"><text x="20" y="40">phylogenetic tree not available</text></svg>\n')
        return path
    tree = SpeciesTree(tree_rows)
    depth, layout = tree_layout_depths(tree)
    x, y = tree_coordinates(tree, depth)
    branch_events = defaultdict(list)
    for row in events:
        scope = row.get("branch_scope", "")
        if "->" in scope and row.get("call_scope", "core_structural_event") == "core_structural_event":
            branch_events[scope].append(row)
    height = max(max(y.values() or [80]) + 90, 100 + min(12, len(events)) * 54)
    width = 1430
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(24, 28, "INSIPHY phylogenetic structural-event map", 16, weight="bold"),
        svg_text(24, 46, f"tree: {tree_file}; layout: {layout}", 10, fill="#555"),
    ]
    for parent, child in tree.edges():
        body.append(f'<line x1="{x[parent]:.2f}" y1="{y[parent]:.2f}" x2="{x[parent]:.2f}" y2="{y[child]:.2f}" stroke="#555" stroke-width="1.2"/>')
        body.append(f'<line x1="{x[parent]:.2f}" y1="{y[child]:.2f}" x2="{x[child]:.2f}" y2="{y[child]:.2f}" stroke="#555" stroke-width="1.2"/>')
        scope = f"{tree.label[parent]}->{tree.label[child]}"
        evs = branch_events.get(scope, [])
        if evs:
            strongest = max(evs, key=change_priority)
            cx = (x[parent] + x[child]) / 2
            cy = y[child]
            body.append(change_symbol(strongest, cx, cy, 8.0))
    for node in tree.preorder():
        if node in tree.leaves:
            body.append(svg_text(630, y[node] + 4, compact_label(tree.label[node], 24), 11))
        else:
            body.append(f'<circle cx="{x[node]:.2f}" cy="{y[node]:.2f}" r="3" fill="#333"/>')
            if node != tree.root:
                body.append(svg_text(x[node] + 6, y[node] - 5, tree.label[node], 9, fill="#555"))
    aside_x = 850
    parsimony = (Path(result_dir) / "branch_structural_events.tsv").exists()
    if events:
        title = "Most-parsimonious branch placements" if parsimony else "Conditional branch probabilities"
    else:
        title = "No valid branch event display"
    body.append(svg_text(aside_x, 58, title, 12, weight="bold"))
    for idx, row in enumerate(sorted(events, key=change_priority, reverse=True)[:12]):
        y0 = 80 + idx * 54
        body.append(change_symbol(row, aside_x + 5, y0 - 5, 5.0))
        support = row.get("placement_status") or f"Pr={change_probability(row):.3f}"
        description = f"{row.get('site_id', '')} | {change_label(row)} | {support}"
        body.append(f'<g><title>{escape(description)}</title>')
        body.append(svg_text(aside_x + 16, y0, compact_label(description, 60), 9))
        body.append("</g>")
        for line_number, line in enumerate(textwrap.wrap(compact_label(row.get("branch_scope", ""), 144), width=72)):
            body.append(svg_text(aside_x + 16, y0 + 13 * (line_number + 1), line, 9, fill="#555"))
    body.append("</svg>")
    path = output_dir / "phylogenetic_event_map.svg"
    path.write_text("\n".join(body))
    return path


def draw_integrated_phylo_synteny(input_dir, result_dir, output_dir, correspondence_encoding="color"):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    tree_rows, tree_file = structural_tree_rows(input_dir, result_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", optional=True)
    events = phylogenetic_change_rows(result_dir)
    path = output_dir / "integrated_phylo_synteny.svg"
    if not tree_rows or not occurrences:
        path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="120"><text x="20" y="40">phylogenetic tree or segment data not available</text></svg>\n')
        return path

    tree = SpeciesTree(tree_rows)
    occ_to_element, styles, class_by_occ, status_by_occ = segment_styles(input_dir, result_dir, occurrences, correspondence_encoding)
    matches = accepted_ribbon_matches(input_dir, result_dir, occurrences)
    groups = track_groups(input_dir, occurrences)
    tree_leaf_order = [node for node in tree.preorder() if node in tree.leaves]
    leaf_order = {tree.label[node]: idx for idx, node in enumerate(tree_leaf_order)}
    ordered_groups = sorted(groups, key=lambda group: (leaf_order.get(tip_label_for_group(tree, group["key"][0], group["key"][1]), 10**6), group["key"][0], group["key"][1]))

    lane_h = 36
    group_gap = 12
    top = 70
    width = 1590
    left_tree = 46
    right_tree = 250
    left_track = 630
    right = 42
    group_top = {}
    cursor = top
    for group in ordered_groups:
        group_top[group["key"]] = cursor
        cursor += 22 + len(group["lanes"]) * lane_h + group_gap
    height = max(cursor + 84, top + 120)
    tip_rows = defaultdict(list)
    for group in ordered_groups:
        species, copy = group["key"]
        first_lane_y = group_top[group["key"]] + 33
        last_lane_y = group_top[group["key"]] + 22 + (len(group["lanes"]) - 1) * lane_h + 11
        tip_rows[tip_label_for_group(tree, species, copy)].append((first_lane_y + last_lane_y) / 2)
    node_y = {}
    for leaf in tree.leaves:
        label = tree.label[leaf]
        vals = tip_rows.get(label, [top + leaf_order.get(label, 0) * (lane_h + group_gap) + 11])
        node_y[leaf] = sum(vals) / len(vals)

    def assign_y(node):
        if node in node_y:
            return node_y[node]
        vals = [assign_y(child) for child in tree.children.get(node, [])]
        node_y[node] = sum(vals) / max(1, len(vals))
        return node_y[node]

    assign_y(tree.root)
    depth, layout = tree_layout_depths(tree)
    max_depth = max(depth.values() or [1.0])
    node_x = {node: left_tree + depth[node] / max(max_depth, 1e-9) * (right_tree - left_tree) for node in depth}

    branch_events = defaultdict(list)
    for row in events:
        scope = row.get("branch_scope", "")
        if "->" in scope and row.get("call_scope", "core_structural_event") == "core_structural_event":
            branch_events[scope].append(row)

    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        pattern_defs(styles),
        svg_text(24, 30, "INSIPHY integrated phylogenetic intragenic synteny", 16, weight="bold"),
        svg_text(24, 52, f"tree: {tree_file}; layout: {layout}", 10, fill="#555"),
        svg_text(left_track, 52, "exon-like gene-internal synteny", 10, fill="#555"),
    ]
    for parent, child in tree.edges():
        body.append(f'<line x1="{node_x[parent]:.2f}" y1="{node_y[parent]:.2f}" x2="{node_x[parent]:.2f}" y2="{node_y[child]:.2f}" stroke="#555" stroke-width="1.2"/>')
        body.append(f'<line x1="{node_x[parent]:.2f}" y1="{node_y[child]:.2f}" x2="{node_x[child]:.2f}" y2="{node_y[child]:.2f}" stroke="#555" stroke-width="1.2"/>')
        scope = f"{tree.label[parent]}->{tree.label[child]}"
        evs = branch_events.get(scope, [])
        if evs:
            strongest = max(evs, key=change_priority)
            cx = (node_x[parent] + node_x[child]) / 2
            cy = node_y[child]
            body.append(change_symbol(strongest, cx, cy))
    for leaf in tree.leaves:
        body.append(f'<line x1="{node_x[leaf]:.2f}" y1="{node_y[leaf]:.2f}" x2="254" y2="{node_y[leaf]:.2f}" stroke="#BBBBBB" stroke-dasharray="2 2"/>')
        body.append(f'<g><title>{escape(tree.label[leaf])}</title>')
        body.append(svg_text(260, node_y[leaf] + 4, compact_label(tree.label[leaf], 24), 10))
        body.append("</g>")
    for node in tree.preorder():
        if node not in tree.leaves and node != tree.root:
            body.append(svg_text(node_x[node] + 6, node_y[node] - 5, tree.label[node], 9, fill="#555"))

    element_boxes = defaultdict(list)
    pending_boxes = []
    for group_idx, group in enumerate(ordered_groups):
        species, copy = group["key"]
        base_y = group_top[group["key"]]
        body.append(svg_text(480, base_y + 15, compact_label(copy, 15), 9, fill="#333"))
        for lane_idx, (lane_id, rows) in enumerate(group["lanes"]):
            y = base_y + 22 + lane_idx * lane_h
            body.append(svg_text(left_track - 10, y + 14, compact_label(lane_id, 18), 7, anchor="end", fill="#666"))
            body.append(f'<line x1="{left_track}" y1="{y + 11}" x2="{width - right}" y2="{y + 11}" stroke="#D0D0D0" stroke-width="1"/>')
            for row in rows:
                x, w = row_geometry(group, row, left_track, width - right)
                occ_id = row.get("occurrence_id", "NA")
                element_id = occ_to_element.get(occ_id, "NA")
                style = styles.get(element_id, {"fill": "#E6E6E6", "stroke_dasharray": "none", "stroke": "#111111"})
                role = row.get("role", "segment")
                unit_class = class_by_occ.get(occ_id, "context") if element_id != "NA" else "context"
                state = visual_status(row, status_by_occ.get(occ_id, {}), unit_class)
                if unit_class == "context":
                    draw_context_span(body, x, y, w, role)
                    continue
                ybox = y + 3 if unit_class == "exon_like" else y + 8
                hbox = 16 if unit_class == "exon_like" else 8
                box = {
                    "x": x,
                    "w": w,
                    "y": ybox + hbox / 2,
                    "class": unit_class,
                    "row": group_idx,
                    "lane": lane_idx,
                    "group": group["key"],
                    "species": species,
                    "occurrence_id": occ_id,
                    "sequence_length": int(row["end"]) - int(row["start"]) + 1,
                }
                if unit_class == "exon_like" and state == "annotated":
                    element_boxes[element_id].append(box)
                pending_boxes.append((x, ybox, w, hbox, element_id, style, unit_class, state, occ_id, lane_id))
    draw_correspondence_ribbons(body, element_boxes, styles, matches)
    for args in pending_boxes:
        draw_segment_box(body, *args)
    legend = (
        "Filled: required in all minimum-change histories; open: possible in some; ribbons: directly aligned bases."
        if (Path(result_dir) / "branch_structural_events.tsv").exists()
        else "Branch symbols show conditional transition probabilities; ribbons show directly aligned bases; introns are gray context."
    )
    if not events:
        legend = "No valid branch event display; ribbons show directly aligned bases; introns are gray context."
    body.append(svg_text(24, height - 32, legend, 10, fill="#555"))
    body.append("</svg>")
    path.write_text("\n".join(body))
    return path


def visualize_results(input_dir, result_dir, output_dir, correspondence_encoding="color"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    synteny = draw_synteny(input_dir, result_dir, output_dir, correspondence_encoding)
    phylogeny = draw_phylogeny(input_dir, result_dir, output_dir)
    integrated = draw_integrated_phylo_synteny(input_dir, result_dir, output_dir, correspondence_encoding)
    rows = [
        {"path": str(synteny), "type": "svg", "description": "Exon-like gene-internal synteny by species and copy"},
        {"path": str(phylogeny), "type": "svg", "description": "Phylogenetic structural event map"},
        {"path": str(integrated), "type": "svg", "description": "Integrated phylogenetic and exon-like synteny map"},
    ]
    write_tsv(output_dir / "visualization_manifest.tsv", rows, ["path", "type", "description"])
    return rows
