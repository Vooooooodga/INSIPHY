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


def draw_context_span(body, x, y, w, role, occurrence_id="NA"):
    if role == "intron":
        body.append(f'<line x1="{x:.2f}" y1="{y + 11}" x2="{x + w:.2f}" y2="{y + 11}" stroke="#B8B8B8" stroke-width="1.1" stroke-dasharray="2 2" data-feature-role="intron" data-occurrence-id="{attr(occurrence_id)}"/>')
        body.append(f'<line x1="{x:.2f}" y1="{y + 7}" x2="{x:.2f}" y2="{y + 15}" stroke="#B8B8B8" stroke-width="0.8" data-occurrence-id="{attr(occurrence_id)}"/>')
        body.append(f'<line x1="{x + w:.2f}" y1="{y + 7}" x2="{x + w:.2f}" y2="{y + 15}" stroke="#B8B8B8" stroke-width="0.8" data-occurrence-id="{attr(occurrence_id)}"/>')
    else:
        body.append(f'<rect x="{x:.2f}" y="{y + 9}" width="{max(1.2, w):.2f}" height="4" rx="1" fill="#EEEEEE" stroke="#B8B8B8" stroke-width="0.5" data-feature-role="{attr(role)}" data-occurrence-id="{attr(occurrence_id)}"/>')


def _parse_matched_blocks(value, query_length, target_length):
    """Parse inclusive occurrence-relative blocks from the explicit match column."""
    blocks = []
    for token in str(value or "").split(";"):
        token = token.strip()
        if not token:
            return None
        parsed = re.fullmatch(r"(\d+)-(\d+):(\d+)-(\d+)", token)
        if parsed is None:
            return None
        qs, qe, ts, te = map(int, parsed.groups())
        if not (1 <= qs <= qe <= query_length and 1 <= ts <= te <= target_length):
            return None
        if qe - qs != te - ts:
            return None
        if blocks and (qs <= blocks[-1][1] or ts <= blocks[-1][3]):
            return None
        blocks.append((qs, qe, ts, te))
    return tuple(blocks) if blocks else None


def ribbon_match_evidence(input_dir, result_dir, occurrences):
    path = Path(result_dir) / "segment_matches.tsv"
    if not path.exists():
        path = Path(input_dir) / "segment_matches.tsv"
    lengths = {
        row["occurrence_id"]: int(row["end"]) - int(row["start"]) + 1
        for row in occurrences
    }
    matches = {}
    conflicts = set()
    unresolved_pairs = set()
    unresolved = defaultdict(list)
    for row in read_tsv(path, optional=True):
        query = row.get("query_occurrence_id")
        target = row.get("subject_occurrence_id")
        if query not in lengths or target not in lengths or query == target:
            continue
        key = tuple(sorted((query, target)))
        match_id = row.get("match_id", "NA")
        if row.get("match_status") != "mapped":
            unresolved_pairs.add(key)
            reason = "ambiguous_match" if row.get("match_status") == "ambiguous" else "unconfirmed_match"
            for occurrence_id in key:
                unresolved[occurrence_id].append((match_id, reason, key[1] if occurrence_id == key[0] else key[0]))
            continue
        basis = row.get("correspondence_basis") or "DNA"
        coding_projection = basis == "annotated_CDS_protein"
        projection_field = "protein_projected_blocks" if coding_projection else "projected_reference_blocks"
        # Coding blocks follow transcript orientation; the DNA alignment remains separate.
        strand = "+" if coding_projection else row.get("alignment_strand")
        if strand != "+":
            unresolved_pairs.add(key)
            for occurrence_id in key:
                unresolved[occurrence_id].append((match_id, "orientation_unresolved", key[1] if occurrence_id == key[0] else key[0]))
            continue
        reference = row.get("projected_reference_occurrence_id", target)
        if reference not in {"", "NA", target}:
            unresolved_pairs.add(key)
            for occurrence_id in key:
                unresolved[occurrence_id].append((match_id, "reference_unresolved", key[1] if occurrence_id == key[0] else key[0]))
            continue
        blocks = _parse_matched_blocks(row.get("matched_blocks"), lengths[query], lengths[target])
        if blocks is None:
            unresolved_pairs.add(key)
            legacy_blocks = _parse_matched_blocks(row.get(projection_field), lengths[query], lengths[target])
            reason = "legacy_projection_without_matched_blocks" if legacy_blocks else "matched_blocks_missing"
            for occurrence_id in key:
                unresolved[occurrence_id].append((match_id, reason, key[1] if occurrence_id == key[0] else key[0]))
            continue
        if query != key[0]:
            blocks = tuple((ts, te, qs, qe) for qs, qe, ts, te in blocks)
        if key in matches and matches[key]["blocks"] != blocks:
            conflicts.add(key)
        elif key not in matches:
            matches[key] = {
                "blocks": blocks,
                "match_id": match_id,
                "alignment_backend": row.get("alignment_backend", "NA"),
                "correspondence_basis": basis,
                "projection_field": projection_field,
                "block_source": "matched_blocks",
            }
    for key in conflicts:
        matches.pop(key, None)
        unresolved_pairs.add(key)
        for occurrence_id in key:
            unresolved[occurrence_id].append(("NA", "conflicting_matched_blocks", key[1] if occurrence_id == key[0] else key[0]))
    for key in unresolved_pairs:
        matches.pop(key, None)
    return matches, unresolved


def accepted_ribbon_matches(input_dir, result_dir, occurrences):
    return ribbon_match_evidence(input_dir, result_dir, occurrences)[0]


def membership_match_blocks(result_dir, occurrences, class_by_occ, status_by_occ):
    """Read occurrence-relative ranges attached to confirmed element memberships."""
    lengths = {row["occurrence_id"]: int(row["end"]) - int(row["start"]) + 1 for row in occurrences}
    occurrence_by_id = {row["occurrence_id"]: row for row in occurrences}
    ranges = defaultdict(list)
    for row in read_tsv(Path(result_dir) / "element_correspondence.tsv", optional=True):
        occurrence_id = row.get("occurrence_id")
        element_id = row.get("element_id")
        if not occurrence_id or occurrence_id not in lengths or not element_id:
            continue
        if row.get("element_class") != "exon_like" or row.get("membership_call", "core_member") not in {"core_member", "resolved_member", "resolved"}:
            continue
        if visual_status(occurrence_by_id[occurrence_id], row, class_by_occ.get(occurrence_id, "exon_like")) in {"predicted", "unknown", "sequence_candidate"}:
            continue
        parsed = _parse_membership_blocks(row.get("matched_blocks", ""), lengths[occurrence_id])
        ranges[(element_id, occurrence_id)].extend(parsed)
    return {key: sorted(set(blocks)) for key, blocks in ranges.items() if blocks}


def _parse_membership_blocks(value, sequence_length):
    tokens = str(value or "").split(";")
    if not tokens or tokens == [""] or tokens == ["NA"]:
        return []
    blocks = []
    current_side = None
    current_match_id = None
    match_block_index = 0
    local_block_index = 0
    for token in tokens:
        token = token.strip()
        tagged = re.fullmatch(r"([^:;]+):(query|subject):(\d+)-(\d+):(\d+)-(\d+)", token)
        aligned = re.fullmatch(r"(\d+)-(\d+):(\d+)-(\d+)", token)
        local = re.fullmatch(r"(\d+)-(\d+)", token)
        if tagged:
            if local_block_index:
                return []
            current_match_id, current_side = tagged.group(1), tagged.group(2)
            query_start, query_end, subject_start, subject_end = map(int, tagged.groups()[2:])
            match_block_index = 0
            start, end = (query_start, query_end) if current_side == "query" else (subject_start, subject_end)
            key = f"{current_match_id}:{match_block_index}"
            match_block_index += 1
        elif aligned and current_side:
            query_start, query_end, subject_start, subject_end = map(int, aligned.groups())
            start, end = (query_start, query_end) if current_side == "query" else (subject_start, subject_end)
            key = f"{current_match_id}:{match_block_index}"
            match_block_index += 1
        elif local and current_side is None:
            if local_block_index:
                return []
            start, end = map(int, local.groups())
            key = f"local:{local_block_index}"
            local_block_index += 1
        else:
            return []
        if not 1 <= start <= end <= sequence_length:
            return []
        blocks.append((key, start, end))
    return blocks


def confirmed_memberships(input_dir, result_dir, occurrences, class_by_occ, status_by_occ):
    rows = read_tsv(Path(result_dir) / "element_correspondence.tsv", optional=True)
    if not rows:
        rows = fallback_element_correspondence(input_dir, occurrences)
    occurrence_by_id = {row["occurrence_id"]: row for row in occurrences}
    by_occurrence = defaultdict(list)
    seen = set()
    for row in rows:
        occurrence_id = row.get("occurrence_id")
        element_id = row.get("element_id")
        if occurrence_id not in occurrence_by_id or not element_id:
            continue
        if row.get("element_class") != "exon_like":
            continue
        if row.get("membership_call", "core_member") not in {"core_member", "resolved_member", "resolved"}:
            continue
        if visual_status(occurrence_by_id[occurrence_id], row, class_by_occ.get(occurrence_id, "exon_like")) in {
            "predicted", "unknown", "sequence_candidate"
        }:
            continue
        key = (element_id, occurrence_id, row.get("matched_blocks", ""))
        if key not in seen:
            seen.add(key)
            by_occurrence[occurrence_id].append(row)
    for occurrence_id in by_occurrence:
        by_occurrence[occurrence_id].sort(key=lambda row: row.get("element_id", ""))
    return by_occurrence


def colored_membership_ranges(matches, memberships_by_occurrence, explicit_membership_ranges):
    ranges = defaultdict(list)
    for key, blocks in explicit_membership_ranges.items():
        ranges[key].extend((start, end, "membership_blocks") for _block_id, start, end in blocks)
    for pair, match in matches.items():
        first, second = pair
        first_elements = {row.get("element_id") for row in memberships_by_occurrence.get(first, [])}
        second_elements = {row.get("element_id") for row in memberships_by_occurrence.get(second, [])}
        for element_id in sorted(first_elements & second_elements):
            for qs, qe, ts, te in match["blocks"]:
                ranges[(element_id, first)].append((qs, qe, "matched_blocks"))
                ranges[(element_id, second)].append((ts, te, "matched_blocks"))
    return {
        key: sorted(set(blocks), key=lambda item: (item[0], item[1], item[2]))
        for key, blocks in ranges.items()
    }


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


def unresolved_without_membership(unresolved, membership_pairs):
    result = {}
    for occurrence_id, records in unresolved.items():
        remaining = [record for record in records
                     if tuple(sorted((occurrence_id, record[2]))) not in membership_pairs]
        if remaining:
            result[occurrence_id] = remaining
    return result


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


def biological_role(row, status_row=None):
    role = str(row.get("role") or "")
    if not role or role.lower() in {"na", "unknown", "segment"}:
        role = str((status_row or {}).get("display_role") or (status_row or {}).get("inferred_role") or "exon")
    return role


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
        "intragenic_synteny.svg", "INSIPHY intragenic synteny", False,
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
            body.append(branch_event_summary(evs, scope, cx, cy - 12, 8))
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
        "integrated_phylo_synteny.svg", "INSIPHY integrated phylogenetic intragenic synteny", True,
    )


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
