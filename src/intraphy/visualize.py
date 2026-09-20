"""visualize: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from intraphy.reporting.inputs import colored_membership_ranges
from intraphy.reporting.inputs import confirmed_memberships
from intraphy.reporting.inputs import membership_match_blocks
from intraphy.reporting.inputs import phylogenetic_change_rows
from intraphy.reporting.inputs import ribbon_match_evidence
from intraphy.reporting.inputs import segment_styles
from intraphy.reporting.inputs import structural_tree_rows
from intraphy.reporting.inputs import visual_status
from intraphy.reporting.layout import _candidate_occurrences
from intraphy.reporting.layout import _oriented_interval
from intraphy.reporting.layout import _row_box_style
from intraphy.reporting.layout import _species_track_layout
from intraphy.reporting.layout import _tree_y_coordinates
from intraphy.reporting.layout import _zoom_windows
from intraphy.reporting.layout import row_geometry
from intraphy.reporting.layout import track_groups
from intraphy.reporting.layout import unresolved_without_membership
from intraphy.reporting.ribbons import _draw_zoom_panels
from intraphy.reporting.ribbons import draw_correspondence_ribbons
from intraphy.reporting.svg import EVENT_SIDEBAR_LIMIT
from intraphy.reporting.svg import _draw_tree
from intraphy.reporting.svg import attr
from intraphy.reporting.svg import branch_event_summary
from intraphy.reporting.svg import change_label
from intraphy.reporting.svg import change_priority
from intraphy.reporting.svg import change_probability
from intraphy.reporting.svg import change_symbol
from intraphy.reporting.svg import clipped_ribbon_box
from intraphy.reporting.svg import compact_label
from intraphy.reporting.svg import draw_context_span
from intraphy.reporting.svg import draw_membership_legend
from intraphy.reporting.svg import draw_segment_box
from intraphy.reporting.svg import feature_shape
from intraphy.reporting.svg import pattern_defs
from intraphy.reporting.svg import svg_text
from intraphy.run_result import result_model
from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from intraphy.topology import SpeciesTree
from pathlib import Path
from xml.sax.saxutils import escape
import math
import textwrap


def _draw_synteny_plot(input_dir, result_dir, output_dir, correspondence_encoding, filename, title, show_events):
    input_dir = Path(input_dir)
    result_dir = Path(result_dir)
    output_dir = Path(output_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", optional=True)
    occ_to_element, styles, class_by_occ, status_by_occ = segment_styles(input_dir, result_dir, occurrences, correspondence_encoding)
    matches, unresolved = ribbon_match_evidence(input_dir, result_dir, occurrences)
    membership_ranges = membership_match_blocks(result_dir, occurrences, class_by_occ, status_by_occ)
    memberships_by_occurrence = confirmed_memberships(input_dir, result_dir, occurrences, class_by_occ, status_by_occ)
    colored_ranges = colored_membership_ranges(matches, memberships_by_occurrence, membership_ranges)
    matched_occurrences = set()
    for pair, match in matches.items():
        first, second = pair
        first_elements = {row.get("element_id") for row in memberships_by_occurrence.get(first, [])}
        second_elements = {row.get("element_id") for row in memberships_by_occurrence.get(second, [])}
        if first_elements & second_elements:
            matched_occurrences.update(pair)
        else:
            for occurrence_id, partner_id in ((first, second), (second, first)):
                unresolved[occurrence_id].append((match["match_id"], "no_confirmed_element_membership", partner_id))
    candidate_ids = _candidate_occurrences(occurrences, class_by_occ, status_by_occ)
    groups = track_groups(input_dir, occurrences, candidate_ids)
    tree_rows, tree_file = structural_tree_rows(input_dir, result_dir)
    tree = SpeciesTree(tree_rows) if tree_rows else None

    width = 1600
    tree_left, tree_right = 32, 175
    label_x, group_label_x, lane_label_x = 188, 330, 440
    track_left, track_right = 455, width - 38
    legend_rows = math.ceil(len(styles) / 8) if styles else 0
    lane_h, group_gap, species_gap, top = 36, 10, 13, 82 + legend_rows * 15
    groups, group_top, species_y, cursor = _species_track_layout(groups, tree, top, lane_h, group_gap, species_gap)
    if not groups:
        cursor = top + lane_h
    zoom_panels = []
    for group in groups:
        for start, end, ids in _zoom_windows(group, group["lanes"], track_left, track_right, occ_to_element, class_by_occ, status_by_occ):
            zoom_panels.append((group, start, end, ids))
    zoom_height = sum(34 + max(1, sum(
        any(_oriented_interval(group, row)[0] < end and _oriented_interval(group, row)[1] > start for row in rows)
        for _lane, rows in group["lanes"]
    )) * 22 + 8 for group, start, end, _ids in zoom_panels)
    events = phylogenetic_change_rows(result_dir) if show_events else []
    extra_event_height = min(12, len(events)) * 20 if show_events else 0
    height = max(150, cursor + zoom_height + extra_event_height + 74)
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" data-encoding="{attr(correspondence_encoding)}">',
        '<rect width="100%" height="100%" fill="white"/>',
        pattern_defs(styles),
        svg_text(24, 28, title, 16, weight="bold"),
        svg_text(24, 47, f'tree: {tree_file}; left tree, right gene structure; encoding: {correspondence_encoding}', 9, fill="#555"),
        svg_text(track_left, 62, "homologous element labels", 9, fill="#555"),
    ]
    draw_membership_legend(body, styles, track_left, track_right, 77)

    if tree:
        depth, layout = tree_layout_depths(tree)
        max_depth = max(depth.values() or [1.0])
        node_x = {node: tree_left + depth[node] / max(max_depth, 1e-9) * (tree_right - tree_left) for node in depth}
        node_y = _tree_y_coordinates(tree, species_y)
        branch_events = defaultdict(list)
        if show_events:
            for row in events:
                scope = row.get("branch_scope", "")
                if "->" in scope and row.get("call_scope", "core_structural_event") == "core_structural_event":
                    branch_events[scope].append(row)
        _draw_tree(body, tree, node_x, node_y, branch_events, show_events)
        body.append(svg_text(24, 60, layout, 8, fill="#666"))
        for tip, y in species_y.items():
            body.append(f'<line x1="292" y1="{y:.2f}" x2="{track_left}" y2="{y:.2f}" stroke="#C6C6C6" stroke-width="0.7" stroke-dasharray="2 3" data-guide="species_gene_tracks" data-species="{attr(tip)}"/>')
    else:
        body.append(svg_text(tree_left, top + 10, "Phylogeny unavailable", 9, fill="#777"))
        for species, y in species_y.items():
            body.append(svg_text(label_x, y + 4, compact_label(species, 26), 9))
            body.append(f'<line x1="292" y1="{y:.2f}" x2="{track_left}" y2="{y:.2f}" stroke="#C6C6C6" stroke-width="0.7" stroke-dasharray="2 3"/>')
    body.append(svg_text(track_left, top - 7, "species and transcript structure, oriented 5' to 3'", 9, fill="#555"))

    element_boxes = defaultdict(list)
    pending_segments = []
    pending_overlays = []
    for group_idx, group in enumerate(groups):
        species, copy = group["key"]
        base_y = group_top[group["key"]]
        body.append(svg_text(group_label_x, base_y + 15, compact_label(copy, 20), 9, fill="#333"))
        for lane_idx, (lane_id, rows) in enumerate(group["lanes"]):
            y = base_y + 22 + lane_idx * lane_h
            body.append(
                f'<g data-lane-id="{attr(lane_id)}" data-species="{attr(species)}" '
                f'data-copy="{attr(copy)}" x1="{track_left}" y1="{y + 11}" '
                f'x2="{track_right}" y2="{y + 11}">'
            )
            body.append(f'<line x1="{track_left}" y1="{y + 11}" x2="{track_right}" y2="{y + 11}" stroke="#C8C8C8" stroke-width="0.8"/>')
            body.append("</g>")
            body.append(svg_text(lane_label_x - 10, y + 14, compact_label(lane_id, 24), 8, anchor="end", fill="#555"))
            for row in rows:
                x, w = row_geometry(group, row, track_left, track_right)
                occurrence_id = row.get("occurrence_id", "NA")
                element_id = occ_to_element.get(occurrence_id, "NA")
                unit_class, role, style = _row_box_style(
                    row, occurrence_id, element_id, styles, class_by_occ, status_by_occ,
                    "matched_blocks" if occurrence_id in matched_occurrences else (
                        "membership_blocks" if any(key[1] == occurrence_id for key in membership_ranges) else (
                            "actual_blocks_missing" if occurrence_id in unresolved else "not_assessed"
                        )
                    ),
                )
                state = visual_status(row, status_by_occ.get(occurrence_id, {}), unit_class)
                shape = "candidate" if lane_id == "candidate evidence" else style["shape"]
                style["shape"] = shape
                if shape == "intron" and lane_id != "candidate evidence":
                    draw_context_span(body, x, y, w, role, occurrence_id)
                    continue
                hbox = 10 if shape in {"utr", "noncoding", "candidate"} else 15 if shape == "cds" else 13
                ybox = y + (5 if hbox <= 10 else 3)
                box = {
                    "x": x, "w": w, "y": ybox + hbox / 2, "class": unit_class,
                    "row": group_idx, "lane": lane_idx, "group": group["key"],
                    "species": species, "occurrence_id": occurrence_id,
                    "sequence_length": int(row["end"]) - int(row["start"]) + 1,
                }
                memberships = memberships_by_occurrence.get(occurrence_id, [])
                for membership in memberships:
                    member_element = membership.get("element_id")
                    if unit_class == "exon_like" and state == "annotated":
                        element_boxes[member_element].append({**box, "element_id": member_element})
                    for block_start, block_end, source in colored_ranges.get((member_element, occurrence_id), []):
                        clipped = clipped_ribbon_box(box, block_start, block_end)
                        member_style = dict(styles[member_element])
                        member_style["display_role"] = role
                        member_style["shape"] = feature_shape(role)
                        member_style["match_state"] = source
                        pending_overlays.append((
                            clipped["x"], ybox, clipped["w"], hbox, member_element,
                            member_style, "exon_like", state, occurrence_id, lane_id,
                            block_start, block_end,
                        ))
                base_style = style
                base_label = element_id
                if lane_id != "candidate evidence" and unit_class != "candidate_source":
                    base_style = {
                        "fill": "#D8D8D8", "stroke": "#777777", "stroke_dasharray": "none",
                        "encoding": "neutral", "display_role": role, "shape": shape,
                        "match_state": style.get("match_state", "not_assessed"),
                    }
                    base_label = element_id
                pending_segments.append((x, ybox, w, hbox, base_label, base_style, unit_class, state, occurrence_id, lane_id))
    evidence_pairs = set(matches)
    legacy_pair_records = defaultdict(set)
    legacy_pair_reasons = defaultdict(set)
    for occurrence_id, records in unresolved.items():
        for match_id, reason, partner_id in records:
            key = tuple(sorted((occurrence_id, partner_id)))
            legacy_pair_records[key].add(match_id)
            legacy_pair_reasons[key].add(reason)
    evidence_pairs.update(
        key for key, reasons in legacy_pair_reasons.items()
        if reasons == {"legacy_projection_without_matched_blocks"}
        and len(legacy_pair_records[key]) == 1
    )
    membership_pairs = draw_correspondence_ribbons(
        body, element_boxes, styles, matches, membership_ranges, evidence_pairs
    )
    membership_connected = {occurrence_id for pair in membership_pairs for occurrence_id in pair}
    for args in pending_segments:
        x, ybox, w, hbox, element_id, style, unit_class, state, occurrence_id, lane_id = args
        if occurrence_id in membership_connected and occurrence_id not in matched_occurrences:
            style["match_state"] = "membership_blocks"
    matched_occurrences.update(membership_connected)
    visible_unresolved = unresolved_without_membership(unresolved, membership_pairs)
    for args in pending_segments:
        draw_segment_box(body, *args)
    for args in pending_overlays:
        draw_segment_box(body, *args)

    zoom_start = cursor + 4
    if zoom_panels:
        body.append(svg_text(track_left, zoom_start + 12, "Local detail panels for short or compressed features", 10, weight="bold", fill="#444"))
        zoom_start = _draw_zoom_panels(
            body, zoom_panels, track_left, track_right, occ_to_element, styles, class_by_occ,
            status_by_occ, visible_unresolved,
            matched_occurrences, memberships_by_occurrence, colored_ranges, zoom_start + 18, lane_h,
        )
    missing_count = sum(len(rows) for rows in visible_unresolved.values())
    if correspondence_encoding == "color":
        legend = "Colorblind-safe color: homologous exon membership; CDS: solid block; UTR: open block; noncoding exon: narrow dotted block; intron: thin gray span."
    else:
        legend = "Pattern: homologous exon membership; CDS: solid block; UTR: open block; noncoding exon: narrow dotted block; intron: thin gray span."
    legend += " Candidate evidence has its own dashed lane. ? means no usable matched_blocks; no ribbon is drawn."
    if zoom_panels:
        legend += " Local detail panels are magnified and labeled."
    body.append(svg_text(24, height - 22, legend, 9, fill="#444"))
    if missing_count:
        body.append(f'<text x="{track_left}" y="{height - 40}" font-family="Arial, sans-serif" font-size="8" fill="#8A4B08" data-unresolved-match-count="{missing_count}">? Correspondence records unresolved for ribbon display: {missing_count}; legacy projection fields were not used as ribbon geometry.</text>')
    if show_events and not events:
        body.append(svg_text(track_left, height - 54, "No valid branch event display", 9, fill="#777"))
    elif show_events and events:
        support_note = (
            " R=required and P=possible."
            if any("placement_status" in row for row in events)
            else ""
        )
        body.append(svg_text(
            track_left,
            height - 54,
            f"Branch history marks: {len(events)}; branch labels give event totals.{support_note}",
            9,
            fill="#555",
        ))
    body.append("</svg>")
    path = output_dir / filename
    path.write_text("\n".join(body))
    return path


def draw_synteny(input_dir, result_dir, output_dir, correspondence_encoding="color"):
    return _draw_synteny_plot(
        input_dir, result_dir, output_dir, correspondence_encoding,
        "intragenic_synteny.svg", "IntraPhy intragenic synteny", False,
    )


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
    shown_event_count = min(EVENT_SIDEBAR_LIMIT, len(events))
    height = max(max(y.values() or [80]) + 90, 138 + shown_event_count * 54)
    width = 1430
    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        svg_text(24, 28, "IntraPhy phylogenetic structural-event map", 16, weight="bold"),
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
            body.append(branch_event_summary(evs, scope, cx, cy - 12, 8))
    for node in tree.preorder():
        if node in tree.leaves:
            body.append(svg_text(630, y[node] + 4, compact_label(tree.label[node], 24), 11))
        else:
            body.append(f'<circle cx="{x[node]:.2f}" cy="{y[node]:.2f}" r="3" fill="#333"/>')
            if node != tree.root:
                body.append(svg_text(x[node] + 6, y[node] - 5, tree.label[node], 9, fill="#555"))
    aside_x = 850
    parsimony = result_model(result_dir) == "parsimony"
    if events:
        title = (
            "Most-parsimonious branch placements (bounded preview)"
            if parsimony
            else "Conditional branch probabilities (bounded preview)"
        )
    else:
        title = "No valid branch event display"
    body.append(svg_text(aside_x, 58, title, 12, weight="bold"))
    source_table = "branch_structural_events.tsv" if parsimony else "structural_changes.tsv"
    additional_event_count = len(events) - shown_event_count
    sidebar_summary = f"Valid events: total={len(events)}; shown={shown_event_count}; additional={additional_event_count}."
    sidebar_source_summary = f"Additional events remaining in {source_table}: {additional_event_count}."
    if parsimony:
        sidebar_source_summary += " R=required; P=possible."
    body.append(
        f'<g data-event-sidebar-total="{len(events)}" data-event-sidebar-shown="{shown_event_count}" '
        f'data-event-sidebar-additional="{additional_event_count}" data-event-source="{attr(source_table)}">'
        f'{svg_text(aside_x, 76, sidebar_summary, 9, fill="#555")}'
        f'{svg_text(aside_x, 89, sidebar_source_summary, 9, fill="#555")}</g>'
    )
    for idx, row in enumerate(sorted(events, key=change_priority, reverse=True)[:EVENT_SIDEBAR_LIMIT]):
        y0 = 112 + idx * 54
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
    return _draw_synteny_plot(
        input_dir, result_dir, output_dir, correspondence_encoding,
        "integrated_phylo_synteny.svg", "IntraPhy integrated phylogenetic intragenic synteny", True,
    )


def visualize_results(input_dir, result_dir, output_dir, correspondence_encoding="color", *,
                      layout="target-groups", targets=(), target_manifest=None):
    """Render one target per group by default; legacy overview is explicit.

    Selection affects overlays only. It does not filter annotation or inference.
    """
    if layout == "target-groups":
        from intraphy.reporting.target_views import visualize_target_groups
        return visualize_target_groups(input_dir, result_dir, output_dir,
                                       selectors=targets, target_manifest=target_manifest)
    if layout != "legacy-overview":
        raise ValueError(f"unknown visualization layout: {layout!r}")
    if targets or target_manifest is not None:
        raise ValueError("target selection requires layout=target-groups")
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


# Public entry points; implementations have a single owner.
from intraphy.reporting.inputs import (
    structural_tree_rows,
    finite_probability,
    valid_probability_change,
    phylogenetic_change_rows,
    fallback_element_correspondence,
    segment_styles,
    visual_status,
    _parse_matched_blocks,
    ribbon_match_evidence,
    accepted_ribbon_matches,
    membership_match_blocks,
    _parse_membership_blocks,
    confirmed_memberships,
    colored_membership_ranges,
    biological_role,
)
from intraphy.reporting.layout import (
    transcript_lanes_for_group,
    track_groups,
    row_geometry,
    connector_pairs,
    unresolved_without_membership,
    _candidate_occurrences,
    _species_track_layout,
    _tree_y_coordinates,
    _zoom_windows,
    _oriented_interval,
    _row_box_style,
)
from intraphy.reporting.ribbons import (
    draw_correspondence_ribbons,
    _draw_zoom_panels,
)
from intraphy.reporting.svg import (
    PALETTE,
    PATTERNS,
    EVENT_SIDEBAR_LIMIT,
    attr,
    to_int,
    compact_label,
    svg_text,
    change_probability,
    change_label,
    change_priority,
    change_symbol,
    branch_event_summary,
    tip_label_for_group,
    pattern_defs,
    draw_membership_legend,
    draw_context_span,
    clipped_ribbon_box,
    draw_homology_connector,
    draw_segment_box,
    feature_shape,
    _draw_tree,
)
import math
import re
import textwrap
