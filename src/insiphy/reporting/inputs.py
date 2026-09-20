"""reporting / inputs: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from insiphy.elements import EXON_LIKE_ROLES
from insiphy.elements import collect_element_profiles
from insiphy.elements import element_class_for_occurrence
from insiphy.reporting.legacy import posterior_row
from insiphy.reporting.svg import PALETTE
from insiphy.reporting.svg import PATTERNS
from insiphy.reporting.svg import change_probability
from insiphy.run_result import result_model
from insiphy.storage.tabular import read_tsv
from pathlib import Path
import math
import re


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
    return str(row.get("posterior_available", "")).lower() in {"true", "1", "yes"} and finite_probability(row)


def phylogenetic_change_rows(result_dir):
    result_dir = Path(result_dir)
    if result_model(result_dir) == "parsimony":
        return [
            row for row in read_tsv(result_dir / "branch_structural_events.tsv")
            if row.get("placement_status") in {"required", "possible"}
        ]
    if result_model(result_dir) in {"er-ard", "foreground"}:
        rows = read_tsv(result_dir / "structural_changes.tsv")
        rows = [posterior_row(row) for row in rows]
        rows = [row for row in rows if valid_probability_change(row)]
        return sorted(
            rows,
            key=change_probability,
            reverse=True,
        )
    return read_tsv(result_dir / "event_support_summary.tsv", optional=True)


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
    status = status_row or {}
    values = {str(value).lower() for value in (
        row.get("role"), row.get("presence_status"), row.get("boundary_state"),
        row.get("path_status"), status.get("support_type"),
        status.get("membership_call"), status.get("inferred_role"),
    ) if value is not None}
    if unit_class == "candidate_source":
        return "sequence_candidate"
    predicted = {"predicted_exon_candidate", "predicted_cds", "predicted_exon",
                 "supports_hidden_segment", "sequence_supported_prediction",
                 "predicted", "inferred_exon", "inferred_cds", "hidden"}
    if values & predicted:
        return "predicted"
    uncertain = {"unknown", "ambiguous", "uncertain", "unresolved",
                 "ambiguous_member", "candidate_member", "ambiguous_mapping"}
    if values & uncertain:
        return "unknown"
    if unit_class == "exon_like":
        return "annotated" if row.get("role") in EXON_LIKE_ROLES else "unknown"
    return "context"


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


def biological_role(row, status_row=None):
    role = str(row.get("role") or "")
    if not role or role.lower() in {"na", "unknown", "segment"}:
        role = str((status_row or {}).get("display_role") or (status_row or {}).get("inferred_role") or "exon")
    return role
