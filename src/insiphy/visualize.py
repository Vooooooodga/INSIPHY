"""Lightweight SVG visualizations for INSIPHY results."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape

from .elements import collect_element_profiles, element_class_for_occurrence
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
def svg_text(x, y, text, size=11, anchor="start", weight="normal", fill="#222"):
    return f'<text x="{x}" y="{y}" font-family="Arial, sans-serif" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}">{escape(str(text))}</text>'


def structural_tree_rows(input_dir, result_dir):
    input_dir = Path(input_dir)
    result_dir = Path(result_dir)
    tree_file = "species_tree.tsv"
    scope_rows = read_tsv(result_dir / "phylogeny_scope.tsv", optional=True)
    for row in scope_rows:
        if row.get("scope") == "structural_characters" and row.get("tree_file"):
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
    for row in element_rows:
        element_class = row.get("element_class", "context")
        if element_class in {"context", "absent"}:
            continue
        by_occ[row["occurrence_id"]] = row["element_id"]
        class_by_occ[row["occurrence_id"]] = element_class
    styles = {}
    for idx, element_id in enumerate(sorted(set(by_occ.values()))):
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
    return by_occ, styles, class_by_occ


def draw_context_span(body, x, y, w, role):
    if role == "intron":
        body.append(f'<line x1="{x:.2f}" y1="{y + 11}" x2="{x + w:.2f}" y2="{y + 11}" stroke="#B8B8B8" stroke-width="1.1" stroke-dasharray="2 2"/>')
        body.append(f'<line x1="{x:.2f}" y1="{y + 7}" x2="{x:.2f}" y2="{y + 15}" stroke="#B8B8B8" stroke-width="0.8"/>')
        body.append(f'<line x1="{x + w:.2f}" y1="{y + 7}" x2="{x + w:.2f}" y2="{y + 15}" stroke="#B8B8B8" stroke-width="0.8"/>')
    else:
        body.append(f'<rect x="{x:.2f}" y="{y + 9}" width="{w:.2f}" height="4" rx="1" fill="#EEEEEE" stroke="#B8B8B8" stroke-width="0.5"/>')


def draw_homology_connector(body, left_box, right_box, style):
    x1 = left_box["x"] + left_box["w"] / 2
    x2 = right_box["x"] + right_box["w"] / 2
    y1 = left_box["y"]
    y2 = right_box["y"]
    mid = (y1 + y2) / 2
    dash = style.get("stroke_dasharray", "none")
    stroke = "#555555"
    if left_box["class"] == "exon_like" and right_box["class"] == "exon_like":
        x1a = left_box["x"]
        x1b = left_box["x"] + left_box["w"]
        x2a = right_box["x"]
        x2b = right_box["x"] + right_box["w"]
        body.append(
            f'<path d="M{x1a:.2f},{y1:.2f} C{x1a:.2f},{mid:.2f} {x2a:.2f},{mid:.2f} {x2a:.2f},{y2:.2f} '
            f'L{x2b:.2f},{y2:.2f} C{x2b:.2f},{mid:.2f} {x1b:.2f},{mid:.2f} {x1b:.2f},{y1:.2f} Z" '
            'fill="#777777" fill-opacity="0.10" stroke="none"/>'
        )
    body.append(f'<path d="M{x1:.2f},{y1:.2f} C{x1:.2f},{mid:.2f} {x2:.2f},{mid:.2f} {x2:.2f},{y2:.2f}" fill="none" stroke="{stroke}" stroke-width="0.85" stroke-opacity="0.55" stroke-dasharray="{dash}"/>')


def draw_segment_box(body, x, ybox, w, hbox, label, style, unit_class="exon_like"):
    fill = style.get("fill", "#E6E6E6")
    dash = style.get("stroke_dasharray", "none")
    stroke = style.get("stroke", "#111111")
    if unit_class == "candidate_source":
        fill = "#FAFAFA"
        stroke = "#666666"
        dash = "2 2"
    body.append(f'<rect x="{x:.2f}" y="{ybox}" width="{w:.2f}" height="{hbox}" rx="2" fill="{fill}" stroke="{stroke}" stroke-width="0.8" stroke-dasharray="{dash}"/>')
    if unit_class == "exon_like" and w > 28:
        body.append(svg_text(x + w / 2, ybox + hbox + 11, label, 8, anchor="middle", fill="#333"))


def draw_synteny(input_dir, result_dir, output_dir, correspondence_encoding="color"):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", optional=True)
    occ_to_element, styles, class_by_occ = segment_styles(input_dir, result_dir, occurrences, correspondence_encoding)
    grouped = defaultdict(list)
    for row in occurrences:
        grouped[(row.get("species", "NA"), row.get("gene_copy_id", "NA"))].append(row)
    width = 1200
    left = 220
    right = 40
    row_h = 34
    top = 50
    height = top + max(1, len(grouped)) * row_h + 90
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        pattern_defs(styles),
        svg_text(24, 28, "INSIPHY intragenic synteny", 16, weight="bold"),
    ]
    element_boxes = defaultdict(list)
    pending_boxes = []
    for ridx, (key, rows) in enumerate(sorted(grouped.items()), start=0):
        species, copy = key
        rows = sorted(rows, key=lambda row: (row.get("contig", ""), int(row.get("start", "0")), int(row.get("end", "0"))))
        starts = [int(row.get("start", "0")) for row in rows]
        ends = [int(row.get("end", "0")) for row in rows]
        start = min(starts or [0])
        end = max(ends or [start + 1])
        span = max(1, end - start + 1)
        y = top + ridx * row_h
        body.append(svg_text(18, y + 15, f"{species}  {copy}", 11))
        body.append(f'<line x1="{left}" y1="{y + 11}" x2="{width - right}" y2="{y + 11}" stroke="#C8C8C8" stroke-width="1"/>')
        for row in rows:
            x = left + (int(row.get("start", "0")) - start) / span * (width - left - right)
            w = max(4, (int(row.get("end", "0")) - int(row.get("start", "0")) + 1) / span * (width - left - right))
            element_id = occ_to_element.get(row["occurrence_id"], "NA")
            style = styles.get(element_id, {"pattern_id": "missing", "stroke_dasharray": "none", "stroke": "#999999"})
            role = row.get("role", "segment")
            unit_class = class_by_occ.get(row["occurrence_id"], "context") if element_id in styles else "context"
            if unit_class == "context":
                draw_context_span(body, x, y, w, role)
                continue
            ybox = y + 3 if unit_class == "exon_like" else y + 8
            hbox = 16 if unit_class == "exon_like" else 8
            box = {"x": x, "w": w, "y": ybox + hbox / 2, "class": unit_class, "row": ridx}
            element_boxes[element_id].append(box)
            pending_boxes.append((x, ybox, w, hbox, element_id, style, unit_class))
    for element_id, boxes in sorted(element_boxes.items()):
        if len(boxes) < 2:
            continue
        boxes = sorted(boxes, key=lambda item: (item["row"], item["x"]))
        for left_box, right_box in zip(boxes, boxes[1:]):
            draw_homology_connector(body, left_box, right_box, styles.get(element_id, {}))
    for args in pending_boxes:
        draw_segment_box(body, *args)
    legend_y = height - 58
    if correspondence_encoding == "color":
        legend = "Color and label mark exon-like correspondence groups"
    else:
        legend = "Texture, line style, labels and links mark exon-like correspondence groups; gray spans mark introns/context"
    body.append(svg_text(24, legend_y, legend, 10, fill="#555"))
    body.append("</svg>")
    path = output_dir / "intragenic_synteny.svg"
    path.write_text("\n".join(body))
    return path


def tree_coordinates(tree):
    leaves = sorted(tree.leaves, key=lambda node: tree.label[node])
    y = {node: 60 + idx * 42 for idx, node in enumerate(leaves)}

    def assign_y(node):
        if node in y:
            return y[node]
        vals = [assign_y(child) for child in tree.children.get(node, [])]
        y[node] = sum(vals) / max(1, len(vals))
        return y[node]

    assign_y(tree.root)
    depth = {tree.root: 0.0}
    for node in tree.preorder():
        for child in tree.children.get(node, []):
            depth[child] = depth[node] + tree.branch_length(child)
    max_depth = max(depth.values() or [1.0])
    x = {node: 120 + depth[node] / max(max_depth, 1e-9) * 500 for node in depth}
    return x, y


def draw_phylogeny(input_dir, result_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    tree_rows, tree_file = structural_tree_rows(input_dir, result_dir)
    events = read_tsv(Path(result_dir) / "event_support_summary.tsv", optional=True)
    if not tree_rows:
        path = output_dir / "phylogenetic_event_map.svg"
        path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="120"><text x="20" y="40">phylogenetic tree not available</text></svg>\n')
        return path
    tree = SpeciesTree(tree_rows)
    x, y = tree_coordinates(tree)
    branch_events = defaultdict(list)
    other_events = []
    for row in events:
        scope = row.get("branch_scope", "")
        if "->" in scope and row.get("call_scope") == "core_structural_event":
            branch_events[scope].append(row)
        else:
            other_events.append(row)
    height = max(y.values() or [80]) + 90
    width = 1100
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(24, 28, "INSIPHY phylogenetic structural-event map", 16, weight="bold"),
        svg_text(24, 46, f"tree: {tree_file}", 10, fill="#555"),
    ]
    for parent, child in tree.edges():
        body.append(f'<line x1="{x[parent]:.2f}" y1="{y[parent]:.2f}" x2="{x[parent]:.2f}" y2="{y[child]:.2f}" stroke="#555" stroke-width="1.2"/>')
        body.append(f'<line x1="{x[parent]:.2f}" y1="{y[child]:.2f}" x2="{x[child]:.2f}" y2="{y[child]:.2f}" stroke="#555" stroke-width="1.2"/>')
        scope = f"{tree.label[parent]}->{tree.label[child]}"
        evs = branch_events.get(scope, [])
        if evs:
            high = sum(1 for row in evs if row.get("support_tier") == "high")
            cx = (x[parent] + x[child]) / 2
            cy = y[child]
            fill = "#D55E00" if high else "#E69F00"
            if high:
                body.append(f'<path d="M{cx:.2f},{cy - 9:.2f} L{cx + 9:.2f},{cy:.2f} L{cx:.2f},{cy + 9:.2f} L{cx - 9:.2f},{cy:.2f} Z" fill="{fill}" fill-opacity="0.75" stroke="#111" stroke-width="0.8"/>')
            else:
                body.append(f'<rect x="{cx - 7:.2f}" y="{cy - 7:.2f}" width="14" height="14" fill="{fill}" fill-opacity="0.65" stroke="#111" stroke-width="0.8" stroke-dasharray="3 2"/>')
            label = ",".join(sorted({row.get("event_class", "event") for row in evs})[:2])
            body.append(svg_text(cx + 12, y[child] - 7, label, 9, fill="#333"))
    for node in tree.preorder():
        if node in tree.leaves:
            body.append(svg_text(x[node] + 8, y[node] + 4, tree.label[node], 11))
        else:
            body.append(f'<circle cx="{x[node]:.2f}" cy="{y[node]:.2f}" r="3" fill="#333"/>')
            if node != tree.root:
                body.append(svg_text(x[node] + 6, y[node] - 5, tree.label[node], 9, fill="#555"))
    aside_x = 720
    body.append(svg_text(aside_x, 58, "Event support summary", 12, weight="bold"))
    for idx, row in enumerate(events[:12]):
        y0 = 80 + idx * 18
        tier = row.get("support_tier", "qualitative")
        fill = {"high": "#D55E00", "moderate": "#E69F00", "qualitative": "#777"}.get(tier, "#777")
        if tier == "high":
            body.append(f'<path d="M{aside_x + 5},{y0 - 11} L{aside_x + 10},{y0 - 6} L{aside_x + 5},{y0 - 1} L{aside_x},{y0 - 6} Z" fill="{fill}" stroke="#111" stroke-width="0.7"/>')
        elif tier == "moderate":
            body.append(f'<rect x="{aside_x}" y="{y0 - 10}" width="10" height="10" fill="{fill}" stroke="#111" stroke-width="0.7" stroke-dasharray="3 2"/>')
        else:
            body.append(f'<circle cx="{aside_x + 5}" cy="{y0 - 5}" r="5" fill="{fill}" stroke="#111" stroke-width="0.7"/>')
        body.append(svg_text(aside_x + 16, y0, f"{row.get('event_class')} | {row.get('branch_scope')} | {tier}", 9))
    body.append("</svg>")
    path = output_dir / "phylogenetic_event_map.svg"
    path.write_text("\n".join(body))
    return path


def draw_integrated_phylo_synteny(input_dir, result_dir, output_dir, correspondence_encoding="color"):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    tree_rows, tree_file = structural_tree_rows(input_dir, result_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", optional=True)
    events = read_tsv(Path(result_dir) / "event_support_summary.tsv", optional=True)
    path = output_dir / "integrated_phylo_synteny.svg"
    if not tree_rows or not occurrences:
        path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="120"><text x="20" y="40">phylogenetic tree or segment data not available</text></svg>\n')
        return path

    tree = SpeciesTree(tree_rows)
    occ_to_element, styles, class_by_occ = segment_styles(input_dir, result_dir, occurrences, correspondence_encoding)
    grouped = defaultdict(list)
    for row in occurrences:
        grouped[(row.get("species", "NA"), row.get("gene_copy_id", "NA"))].append(row)
    leaf_order = {tree.label[node]: idx for idx, node in enumerate(sorted(tree.leaves, key=lambda node: tree.label[node]))}
    ordered_groups = sorted(grouped.items(), key=lambda item: (leaf_order.get(tip_label_for_group(tree, item[0][0], item[0][1]), 10**6), item[0][0], item[0][1]))

    row_h = 34
    top = 70
    width = 1320
    left_tree = 46
    right_tree = 250
    left_track = 360
    right = 42
    height = top + max(1, len(ordered_groups)) * row_h + 96
    row_y = {key: top + idx * row_h for idx, (key, _rows) in enumerate(ordered_groups)}
    tip_rows = defaultdict(list)
    for (species, copy), y in row_y.items():
        tip_rows[tip_label_for_group(tree, species, copy)].append(y + 11)
    node_y = {}
    for leaf in tree.leaves:
        label = tree.label[leaf]
        vals = tip_rows.get(label, [top + leaf_order.get(label, 0) * row_h + 11])
        node_y[leaf] = sum(vals) / len(vals)

    def assign_y(node):
        if node in node_y:
            return node_y[node]
        vals = [assign_y(child) for child in tree.children.get(node, [])]
        node_y[node] = sum(vals) / max(1, len(vals))
        return node_y[node]

    assign_y(tree.root)
    depth = {tree.root: 0.0}
    for node in tree.preorder():
        for child in tree.children.get(node, []):
            depth[child] = depth[node] + tree.branch_length(child)
    max_depth = max(depth.values() or [1.0])
    node_x = {node: left_tree + depth[node] / max(max_depth, 1e-9) * (right_tree - left_tree) for node in depth}

    branch_events = defaultdict(list)
    for row in events:
        scope = row.get("branch_scope", "")
        if "->" in scope and row.get("call_scope") == "core_structural_event":
            branch_events[scope].append(row)

    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        pattern_defs(styles),
        svg_text(24, 30, "INSIPHY integrated phylogenetic intragenic synteny", 16, weight="bold"),
        svg_text(24, 52, f"tree: {tree_file}", 10, fill="#555"),
        svg_text(left_track, 52, "exon-like gene-internal synteny", 10, fill="#555"),
    ]
    for parent, child in tree.edges():
        body.append(f'<line x1="{node_x[parent]:.2f}" y1="{node_y[parent]:.2f}" x2="{node_x[parent]:.2f}" y2="{node_y[child]:.2f}" stroke="#555" stroke-width="1.2"/>')
        body.append(f'<line x1="{node_x[parent]:.2f}" y1="{node_y[child]:.2f}" x2="{node_x[child]:.2f}" y2="{node_y[child]:.2f}" stroke="#555" stroke-width="1.2"/>')
        scope = f"{tree.label[parent]}->{tree.label[child]}"
        evs = branch_events.get(scope, [])
        if evs:
            high = any(row.get("support_tier") == "high" for row in evs)
            cx = (node_x[parent] + node_x[child]) / 2
            cy = node_y[child]
            fill = "#D55E00" if high else "#777777"
            body.append(f'<path d="M{cx:.2f},{cy - 7:.2f} L{cx + 7:.2f},{cy:.2f} L{cx:.2f},{cy + 7:.2f} L{cx - 7:.2f},{cy:.2f} Z" fill="{fill}" fill-opacity="0.72" stroke="#111" stroke-width="0.7"/>')
    for leaf in tree.leaves:
        body.append(svg_text(node_x[leaf] + 6, node_y[leaf] + 4, tree.label[leaf], 10))
    for node in tree.preorder():
        if node not in tree.leaves and node != tree.root:
            body.append(svg_text(node_x[node] + 6, node_y[node] - 5, tree.label[node], 9, fill="#555"))

    element_boxes = defaultdict(list)
    pending_boxes = []
    for group_idx, ((species, copy), rows) in enumerate(ordered_groups):
        rows = sorted(rows, key=lambda row: (row.get("contig", ""), int(row.get("start", "0")), int(row.get("end", "0"))))
        starts = [int(row.get("start", "0")) for row in rows]
        ends = [int(row.get("end", "0")) for row in rows]
        start = min(starts or [0])
        end = max(ends or [start + 1])
        span = max(1, end - start + 1)
        y = row_y[(species, copy)]
        body.append(svg_text(260, y + 15, copy, 9, fill="#333"))
        body.append(f'<line x1="{left_track}" y1="{y + 11}" x2="{width - right}" y2="{y + 11}" stroke="#D0D0D0" stroke-width="1"/>')
        for row in rows:
            x = left_track + (int(row.get("start", "0")) - start) / span * (width - left_track - right)
            w = max(4, (int(row.get("end", "0")) - int(row.get("start", "0")) + 1) / span * (width - left_track - right))
            element_id = occ_to_element.get(row["occurrence_id"], "NA")
            style = styles.get(element_id, {"fill": "#E6E6E6", "stroke_dasharray": "none", "stroke": "#111111"})
            role = row.get("role", "segment")
            unit_class = class_by_occ.get(row["occurrence_id"], "context") if element_id in styles else "context"
            if unit_class == "context":
                draw_context_span(body, x, y, w, role)
                continue
            ybox = y + 3 if unit_class == "exon_like" else y + 8
            hbox = 16 if unit_class == "exon_like" else 8
            box = {"x": x, "w": w, "y": ybox + hbox / 2, "class": unit_class, "row": group_idx}
            element_boxes[element_id].append(box)
            pending_boxes.append((x, ybox, w, hbox, element_id, style, unit_class))
    for element_id, boxes in sorted(element_boxes.items()):
        if len(boxes) < 2:
            continue
        boxes = sorted(boxes, key=lambda item: (item["row"], item["x"]))
        for left_box, right_box in zip(boxes, boxes[1:]):
            draw_homology_connector(body, left_box, right_box, styles.get(element_id, {}))
    for args in pending_boxes:
        draw_segment_box(body, *args)
    body.append(svg_text(24, height - 32, "Core structural events are marked on tree branches; links mark exon-like correspondence; introns are gray context spans.", 10, fill="#555"))
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
