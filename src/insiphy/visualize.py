"""Lightweight SVG visualizations for INSIPHY results."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape

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


def pattern_defs(styles):
    body = ["<defs>"]
    for hsg, style in styles.items():
        pid = style["pattern_id"]
        color = style["color"]
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


def segment_styles(homology):
    by_occ = {}
    for row in homology:
        by_occ[row["occurrence_id"]] = row["homology_id"]
    styles = {}
    for idx, hsg in enumerate(sorted(set(by_occ.values()))):
        styles[hsg] = {
            "color": PALETTE[idx % len(PALETTE)],
            "pattern": PATTERNS[idx % len(PATTERNS)],
            "pattern_id": f"hsg_{idx + 1}",
            "stroke_dasharray": "none" if idx % 3 == 0 else "4 2" if idx % 3 == 1 else "1.5 2",
        }
    return by_occ, styles


def draw_synteny(input_dir, result_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", optional=True)
    homology = read_tsv(input_dir / "segment_homology.tsv", optional=True)
    occ_to_hsg, styles = segment_styles(homology)
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
            hsg = occ_to_hsg.get(row["occurrence_id"], "NA")
            style = styles.get(hsg, {"pattern_id": "missing", "stroke_dasharray": "none", "color": "#999999"})
            role = row.get("role", "segment")
            ybox = y + 3 if role in {"CDS", "exon", "UTR", "noncoding_exon"} else y + 8
            hbox = 16 if role in {"CDS", "exon", "UTR", "noncoding_exon"} else 8
            fill = f'url(#{style["pattern_id"]})' if hsg in styles else "#E6E6E6"
            body.append(f'<rect x="{x:.2f}" y="{ybox}" width="{w:.2f}" height="{hbox}" rx="2" fill="{fill}" stroke="#111" stroke-width="0.8" stroke-dasharray="{style["stroke_dasharray"]}"/>')
            if w > 22:
                body.append(svg_text(x + w / 2, ybox + hbox + 11, hsg, 8, anchor="middle", fill="#333"))
    legend_y = height - 58
    body.append(svg_text(24, legend_y, "Texture, line style and label mark homologous segment groups; color is auxiliary and colorblind-friendly", 10, fill="#555"))
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
    tree_rows = read_tsv(input_dir / "species_tree.tsv", optional=True)
    events = read_tsv(Path(result_dir) / "event_support_summary.tsv", optional=True)
    if not tree_rows:
        path = output_dir / "phylogenetic_event_map.svg"
        path.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="120"><text x="20" y="40">species_tree.tsv not available</text></svg>\n')
        return path
    tree = SpeciesTree(tree_rows)
    x, y = tree_coordinates(tree)
    branch_events = defaultdict(list)
    other_events = []
    for row in events:
        scope = row.get("branch_scope", "")
        if "->" in scope:
            branch_events[scope].append(row)
        else:
            other_events.append(row)
    height = max(y.values() or [80]) + 90
    width = 1100
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(24, 28, "INSIPHY phylogenetic structural-event map", 16, weight="bold"),
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


def visualize_results(input_dir, result_dir, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    synteny = draw_synteny(input_dir, result_dir, output_dir)
    phylogeny = draw_phylogeny(input_dir, result_dir, output_dir)
    rows = [
        {"path": str(synteny), "type": "svg", "description": "Gene-internal synteny by species and copy"},
        {"path": str(phylogeny), "type": "svg", "description": "Species-tree structural event map"},
    ]
    write_tsv(output_dir / "visualization_manifest.tsv", rows, ["path", "type", "description"])
    return rows
