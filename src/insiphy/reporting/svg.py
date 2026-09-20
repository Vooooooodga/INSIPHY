"""reporting / svg: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from xml.sax.saxutils import escape
import math


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


EVENT_SIDEBAR_LIMIT = 12


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


def branch_event_summary(events, branch_scope, x, y, size=7):
    total = len(events)
    has_parsimony_status = any("placement_status" in row for row in events)
    status_attributes = ""
    if has_parsimony_status:
        required = sum(row.get("placement_status") == "required" for row in events)
        possible = sum(row.get("placement_status") == "possible" for row in events)
        label = f"{total} total (R{required}/P{possible})"
        title = f"{branch_scope}: {total} events; {required} required; {possible} possible"
        status_attributes = (
            f' data-event-required="{required}" data-event-possible="{possible}"'
        )
    else:
        label = f"{total} total"
        title = f"{branch_scope}: {total} conditional events"
    return (
        f'<g data-branch-event-summary="true" data-branch-scope="{attr(branch_scope)}" '
        f'data-event-total="{total}"{status_attributes}>'
        f'<title>{escape(title)}</title>'
        f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial, sans-serif" font-size="{size}" '
        f'text-anchor="middle" font-weight="bold" fill="#333" paint-order="stroke" '
        f'stroke="white" stroke-width="2.4">{escape(label)}</text></g>'
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


def draw_membership_legend(body, styles, left, right, start_y, per_row=8):
    elements = sorted(styles)
    row_count = math.ceil(len(elements) / per_row) if elements else 0
    item_width = (right - left) / max(1, per_row)
    for index, element_id in enumerate(elements):
        row = index // per_row
        column = index % per_row
        x = left + column * item_width
        y = start_y + row * 15
        style = styles[element_id]
        body.append(
            f'<rect x="{x:.2f}" y="{y - 9:.2f}" width="10" height="10" '
            f'fill="{style["fill"]}" stroke="#222" stroke-width="0.6" data-legend-element="{attr(element_id)}"/>'
        )
        body.append(svg_text(x + 14, y, compact_label(element_id, 20), 8, fill="#333"))
    return row_count


def draw_context_span(body, x, y, w, role, occurrence_id="NA"):
    if role == "intron":
        body.append(f'<line x1="{x:.2f}" y1="{y + 11}" x2="{x + w:.2f}" y2="{y + 11}" stroke="#B8B8B8" stroke-width="1.1" stroke-dasharray="2 2" data-feature-role="intron" data-occurrence-id="{attr(occurrence_id)}"/>')
        body.append(f'<line x1="{x:.2f}" y1="{y + 7}" x2="{x:.2f}" y2="{y + 15}" stroke="#B8B8B8" stroke-width="0.8" data-occurrence-id="{attr(occurrence_id)}"/>')
        body.append(f'<line x1="{x + w:.2f}" y1="{y + 7}" x2="{x + w:.2f}" y2="{y + 15}" stroke="#B8B8B8" stroke-width="0.8" data-occurrence-id="{attr(occurrence_id)}"/>')
    else:
        body.append(f'<rect x="{x:.2f}" y="{y + 9}" width="{max(1.2, w):.2f}" height="4" rx="1" fill="#EEEEEE" stroke="#B8B8B8" stroke-width="0.5" data-feature-role="{attr(role)}" data-occurrence-id="{attr(occurrence_id)}"/>')


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
    color_mode = style.get("encoding") == "color"
    stroke = style.get("fill", "#555555") if color_mode else "#555555"
    ribbon_fill = stroke if color_mode else "#777777"
    provenance = (
        f'data-match-id="{attr(match["match_id"])}" '
        f'data-block-source="{attr(match.get("block_source", "matched_blocks"))}" '
        f'data-block-field="matched_blocks" '
        f'data-alignment-backend="{attr(match["alignment_backend"])}" '
        f'data-correspondence-basis="{attr(match["correspondence_basis"])}" '
        f'data-projection-field="{attr(match["projection_field"])}" '
        f'data-source-occurrence="{attr(left_box["occurrence_id"])}" '
        f'data-target-occurrence="{attr(right_box["occurrence_id"])}" '
        f'data-source-start="{left_box["local_start"]}" data-source-end="{left_box["local_end"]}" '
        f'data-target-start="{right_box["local_start"]}" data-target-end="{right_box["local_end"]}" '
        f'data-correspondence-kind="{attr(match.get("connector_kind", "aligned_match"))}"'
    )
    connector_kind = match.get("connector_kind", "exon_correspondence")
    if left_box["class"] == "exon_like" and right_box["class"] == "exon_like":
        x1a = left_box["x"]
        x1b = left_box["x"] + left_box["w"]
        x2a = right_box["x"]
        x2b = right_box["x"] + right_box["w"]
        body.append(
            f'<path d="M{x1a:.2f},{y1:.2f} C{x1a:.2f},{mid:.2f} {x2a:.2f},{mid:.2f} {x2a:.2f},{y2:.2f} '
            f'L{x2b:.2f},{y2:.2f} C{x2b:.2f},{mid:.2f} {x1b:.2f},{mid:.2f} {x1b:.2f},{y1:.2f} Z" '
            f'fill="{ribbon_fill}" fill-opacity="0.14" stroke="none" data-connector="{attr(connector_kind)}" data-element-id="{attr(element_id)}" {provenance}/>'
        )
    body.append(f'<path d="M{x1:.2f},{y1:.2f} C{x1:.2f},{mid:.2f} {x2:.2f},{mid:.2f} {x2:.2f},{y2:.2f}" fill="none" stroke="{stroke}" stroke-width="0.95" stroke-opacity="0.70" stroke-dasharray="{dash}" data-connector="{attr(connector_kind + "_line")}" data-element-id="{attr(element_id)}" {provenance}/>')


def draw_segment_box(body, x, ybox, w, hbox, label, style, unit_class="exon_like", visual_state="annotated", occurrence_id="NA", lane_id="NA", block_start="NA", block_end="NA"):
    role = str(style.get("display_role", "exon"))
    shape = style.get("shape", "exon")
    match_state = style.get("match_state", "")
    fill = style.get("fill", "#D0D0D0")
    stroke = style.get("stroke", "#555555")
    dash = style.get("stroke_dasharray", "none")
    visible_width = max(1.2, w)
    draw_x = x - (visible_width - w) / 2

    if shape == "utr":
        hbox = min(hbox, 9)
        if style.get("encoding") in {"color", "neutral"}:
            fill = "#FFFFFF"
            stroke = style.get("fill", "#555555") if style.get("encoding") == "color" else "#777777"
    elif shape == "noncoding":
        hbox = min(hbox, 9)
        if style.get("encoding") == "color":
            fill = style.get("fill", "#888888")
            stroke = style.get("fill", "#555555")
        dash = "1.5 1.5" if dash == "none" else dash
    elif shape == "context_exon":
        fill = "#E6E6E6"
        stroke = "#777777"

    if visual_state == "predicted":
        fill = "#FFFFFF"
        dash = "3 2"
    elif visual_state == "unknown":
        fill = "#F2F2F2"
        stroke = "#777777"
        dash = "1 2"
    if unit_class == "candidate_source" or lane_id == "candidate evidence":
        fill = "#FFFFFF"
        stroke = "#666666"
        dash = "2 2"

    body.append(
        f'<rect x="{draw_x:.2f}" y="{ybox}" width="{visible_width:.2f}" height="{hbox}" rx="1" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="0.9" stroke-dasharray="{dash}" '
        f'data-element-id="{attr(label)}" data-occurrence-id="{attr(occurrence_id)}" '
        f'data-unit-class="{attr(unit_class)}" data-visual-status="{attr(visual_state)}" '
        f'data-feature-shape="{attr(shape)}" data-display-role="{attr(role)}" '
        f'data-lane-id="{attr(lane_id)}" data-match-state="{attr(match_state or "not_assessed")}" '
        f'data-block-start="{attr(block_start)}" data-block-end="{attr(block_end)}">'
        f'<title>{escape(str(label))}; {escape(str(role))}; {escape(str(visual_state))}; '
        f'{escape(str(match_state or "matched-block status unavailable"))}</title></rect>'
    )
    if shape == "utr" and style.get("encoding") == "color" and w > 1:
        mid_y = ybox + hbox / 2
        body.append(f'<line x1="{draw_x:.2f}" y1="{mid_y:.2f}" x2="{draw_x + visible_width:.2f}" y2="{mid_y:.2f}" stroke="{stroke}" stroke-width="0.7"/>')
    if match_state and match_state not in {"matched_blocks", "not_assessed"}:
        body.append(svg_text(draw_x + visible_width + 2, ybox + max(8, hbox - 1), "?", 9, weight="bold", fill="#8A4B08"))
    if lane_id == "candidate evidence":
        body.append(svg_text(draw_x + visible_width + 3, ybox + hbox + 8, compact_label(occurrence_id, 24), 7, fill="#333"))
    elif unit_class == "exon_like" and w > 28:
        label_length = max(4, int(w / 5))
        body.append(svg_text(x + w / 2, ybox + hbox + 11, compact_label(label, label_length), 8, anchor="middle", fill="#333"))


def feature_shape(role):
    normalized = str(role).lower().replace("-", "_").replace(" ", "_")
    if "intron" in normalized:
        return "intron"
    if "utr" in normalized or "untranslated" in normalized:
        return "utr"
    if normalized in {"cds", "coding", "coding_exon", "protein_coding_exon", "start_codon", "stop_codon"}:
        return "cds"
    if any(token in normalized for token in ("noncoding", "non_coding", "ncrna", "lncrna", "lnc_rna", "rrna", "trna", "rna_exon")):
        return "noncoding"
    if "exon" in normalized:
        return "exon"
    return "context_exon"


def _draw_tree(body, tree, node_x, node_y, branch_events=None, show_events=False):
    branch_events = branch_events or {}
    for parent, child in tree.edges():
        body.append(f'<line x1="{node_x[parent]:.2f}" y1="{node_y[parent]:.2f}" x2="{node_x[parent]:.2f}" y2="{node_y[child]:.2f}" stroke="#555" stroke-width="1.2" data-tree-edge="vertical"/>')
        body.append(f'<line x1="{node_x[parent]:.2f}" y1="{node_y[child]:.2f}" x2="{node_x[child]:.2f}" y2="{node_y[child]:.2f}" stroke="#555" stroke-width="1.2" data-tree-edge="horizontal"/>')
        if show_events:
            scope = f"{tree.label[parent]}->{tree.label[child]}"
            events = branch_events.get(scope, [])
            if events:
                event = max(events, key=change_priority)
                event_x = (node_x[parent] + node_x[child]) / 2
                event_y = node_y[child]
                body.append(change_symbol(event, event_x, event_y, 7.0))
                body.append(branch_event_summary(events, scope, event_x, event_y - 11, 7))
    for node in tree.preorder():
        if node in tree.leaves:
            body.append(svg_text(node_x[node] + 6, node_y[node] + 4, compact_label(tree.label[node], 24), 9))
        elif node == tree.root:
            body.append(f'<circle cx="{node_x[node]:.2f}" cy="{node_y[node]:.2f}" r="3.2" fill="#111" data-tree-root="true"/>')
        else:
            body.append(f'<circle cx="{node_x[node]:.2f}" cy="{node_y[node]:.2f}" r="2.4" fill="#333"/>')
            body.append(svg_text(node_x[node] + 5, node_y[node] - 4, compact_label(tree.label[node], 18), 7, fill="#555"))
