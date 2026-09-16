"""Likelihood analysis of single-copy intragenic structural sites."""

from __future__ import annotations

import json
import math
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.optimize import brentq, minimize
from scipy.stats import chi2

from .alignment import AlignmentBackendError, local_alignment_stats
from .io import norm_state, parse_fasta, read_tsv, write_tsv
from .tree import SpeciesTree


EXONIC_ROLES = {"CDS", "exon", "UTR", "noncoding_exon"}
EXON_COMPLETION_CALLS = {
    "hidden_segment_candidate",
    "shifted_splice_site_candidate",
    "joined_exon_candidate",
    "hidden_segment_with_frame_disruption",
}
RATE_MIN = 1e-8
RATE_MAX = 100.0
MULTIPLIER_MIN = 1e-3
MULTIPLIER_MAX = 1e3
PROFILE_DROP_95 = 0.5 * float(chi2.ppf(0.95, 1))


def _fmt(value):
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value):.8g}"


def _validated_tree_rows(path, branch_length_mode):
    rows = read_tsv(path, ["node_id", "parent_id", "label"])
    out = []
    for row in rows:
        item = dict(row)
        if not row.get("parent_id"):
            item["branch_length"] = 0.0
        elif branch_length_mode == "unit":
            item["branch_length"] = 1.0
        else:
            raw = row.get("branch_length") or row.get("length") or row.get("distance")
            if raw in {None, "", "NA"}:
                raise SystemExit(
                    "species_tree.tsv has a branch without length; use --branch-length-mode unit to analyze branch counts"
                )
            length = float(raw)
            if length <= 0:
                raise SystemExit("all non-root branches must have positive lengths")
            item["branch_length"] = length
        out.append(item)
    return out


def _single_copy_families(occurrences):
    copies = defaultdict(set)
    species_by_family = defaultdict(set)
    for row in occurrences:
        key = (row.get("family_id", "NA"), row.get("species", "NA"))
        copies[key].add(row.get("gene_copy_id", "NA"))
        species_by_family[key[0]].add(key[1])
    excluded = []
    valid = set(species_by_family)
    for (family, species), gene_copies in sorted(copies.items()):
        if len(gene_copies) > 1:
            valid.discard(family)
            excluded.append(
                {
                    "family_id": family,
                    "species": species,
                    "copy_count": len(gene_copies),
                    "gene_copy_ids": ";".join(sorted(gene_copies)),
                    "reason": "multiple_gene_copies_in_single_copy_mode",
                }
            )
    return valid, species_by_family, excluded


def _element_site_rows(
    occurrences,
    element_rows,
    completion_rows,
    valid_families,
    species_by_family,
):
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    mapped = defaultdict(lambda: defaultdict(list))
    family_by_element = {}
    element_by_homology = {}
    for row in element_rows:
        if row.get("element_class") != "exon_like":
            continue
        occ = occ_by_id.get(row.get("occurrence_id", ""))
        if not occ or occ.get("family_id") not in valid_families:
            continue
        family = occ["family_id"]
        element = row["element_id"]
        family_by_element[element] = family
        element_by_homology[row.get("homology_id", "")] = element
        mapped[element][occ["species"]].append((occ, row))

    completion_by_element = defaultdict(lambda: defaultdict(list))
    for row in completion_rows:
        element = element_by_homology.get(row.get("homology_id", ""))
        family = row.get("family_id")
        species = row.get("species")
        if element and family in valid_families and species:
            completion_by_element[element][species].append(row)

    rows = []
    for element, family in sorted(family_by_element.items()):
        for species in sorted(species_by_family[family]):
            observations = mapped[element].get(species, [])
            presence_states = set()
            presence_evidence = set()
            states = set()
            evidence = set()
            for occ, membership in observations:
                if membership.get("membership_call", "core_member") != "core_member":
                    presence_evidence.add("ambiguous_correspondence")
                    evidence.add("ambiguous_correspondence")
                    continue
                presence = norm_state(occ.get("presence_status"))
                if presence == "unknown":
                    presence_evidence.add("uncertain_sequence_or_annotation")
                    evidence.add("uncertain_sequence_or_annotation")
                    continue
                presence_states.add(presence)
                presence_evidence.add(
                    "homologous_sequence_observed" if presence == "present" else "explicit_sequence_absence"
                )
                if presence == "absent":
                    continue
                role = occ.get("role", "unknown")
                completion = [
                    row
                    for row in completion_by_element[element].get(species, [])
                    if row.get("gene_copy_id") == occ.get("gene_copy_id")
                ]
                completed_exon = any(
                    row.get("completion_call") in EXON_COMPLETION_CALLS
                    and row.get("inferred_role") in EXONIC_ROLES
                    for row in completion
                )
                if completed_exon:
                    states.add("exonic")
                    evidence.add("sequence_supported_exon_completion")
                else:
                    states.add("exonic" if role in EXONIC_ROLES else "not_exonic")
                    evidence.add("annotated_exon" if role in EXONIC_ROLES else "homologous_non_exonic_sequence")
            for completion in completion_by_element[element].get(species, []):
                call = completion.get("completion_call")
                if call == "supports_true_absence":
                    presence_states.add("absent")
                    presence_evidence.add("sequence_supported_true_absence")
                elif call in EXON_COMPLETION_CALLS and completion.get("inferred_role") in EXONIC_ROLES:
                    presence_states.add("present")
                    presence_evidence.add("sequence_supported_hidden_exon")
                    states.add("exonic")
                    evidence.add("sequence_supported_exon_completion")
            presence_state = next(iter(presence_states)) if len(presence_states) == 1 else "unknown"
            if len(presence_states) > 1:
                presence_evidence.add("conflicting_presence_states")
            rows.append(
                {
                    "family_id": family,
                    "layer": "exon_presence",
                    "site_id": element,
                    "species": species,
                    "state": presence_state,
                    "state_0": "absent",
                    "state_1": "present",
                    "evidence": ";".join(sorted(presence_evidence)) or "no_observation",
                }
            )
            state = next(iter(states)) if len(states) == 1 else "unknown"
            if len(states) > 1:
                evidence.add("conflicting_transcript_states")
            rows.append(
                {
                    "family_id": family,
                    "layer": "exon_role",
                    "site_id": element,
                    "species": species,
                    "state": state,
                    "state_0": "not_exonic",
                    "state_1": "exonic",
                    "evidence": ";".join(sorted(evidence)) or "no_observation",
                }
            )
    return rows


def _junction_site_rows(input_dir, output_dir, occurrences, element_rows, valid_families, species_by_family):
    paths = read_tsv(Path(input_dir) / "transcript_paths.tsv", optional=True)
    if not paths:
        return []
    occ_by_id = {row["occurrence_id"]: row for row in occurrences}
    sequences = parse_fasta(Path(input_dir) / "segment_sequences.fasta")
    elements_by_occ = defaultdict(list)
    for row in element_rows:
        if row.get("element_class") == "exon_like" and row.get("membership_call", "core_member") == "core_member":
            elements_by_occ[row["occurrence_id"]].append(row["element_id"])

    path_groups = defaultdict(list)
    for row in paths:
        occ = occ_by_id.get(row.get("occurrence_id", ""))
        if not occ or occ.get("family_id") not in valid_families:
            continue
        key = (occ["family_id"], occ["species"], occ["gene_copy_id"], row.get("transcript_id", "NA"))
        path_groups[key].append(row)

    observations = defaultdict(lambda: defaultdict(set))
    evidence = defaultdict(lambda: defaultdict(set))
    family_by_site = {}
    element_counts = defaultdict(lambda: defaultdict(list))
    boundary_rows = []
    reference_by_element = {}
    for row in element_rows:
        occurrence_id = row.get("occurrence_id", "")
        sequence = sequences.get(occurrence_id, "")
        element = row.get("element_id", "")
        if sequence and (element not in reference_by_element or len(sequence) > len(reference_by_element[element][1])):
            reference_by_element[element] = (occurrence_id, sequence)
    canonical_positions = defaultdict(list)

    def boundary_site(family, element, left_occurrence, right_occurrence):
        reference_id, reference = reference_by_element.get(element, ("NA", ""))
        left_sequence = sequences.get(left_occurrence, "")
        right_sequence = sequences.get(right_occurrence, "")
        method = "relative_exon_coordinate"
        if reference and left_sequence and right_sequence:
            try:
                left_alignment = local_alignment_stats(left_sequence, reference)
                right_alignment = local_alignment_stats(right_sequence, reference)
                coordinate = int(round((left_alignment.target_end + right_alignment.target_start) / 2))
                method = "projected_sequence_alignment_coordinate"
            except AlignmentBackendError:
                coordinate = int(round(1000 * len(left_sequence) / max(1, len(left_sequence) + len(right_sequence))))
        else:
            left_length = max(1, int(occ_by_id.get(left_occurrence, {}).get("end", 0)) - int(occ_by_id.get(left_occurrence, {}).get("start", 0)) + 1)
            right_length = max(1, int(occ_by_id.get(right_occurrence, {}).get("end", 0)) - int(occ_by_id.get(right_occurrence, {}).get("start", 0)) + 1)
            coordinate = int(round(1000 * left_length / (left_length + right_length)))
        tolerance = max(3, int(round(0.05 * max(1, len(reference))))) if reference else 25
        positions = canonical_positions[(family, element)]
        canonical = next((value for value in positions if abs(value - coordinate) <= tolerance), None)
        if canonical is None:
            canonical = coordinate
            positions.append(canonical)
        return f"JG_{element}_ALN_{canonical}", coordinate, reference_id, method
    for (family, species, _copy, _transcript), path_rows in sorted(path_groups.items()):
        path_rows.sort(key=lambda row: int(row.get("path_rank", "0")))
        previous = None
        previous_occurrence = None
        intron_between = False
        path_elements = []
        for row in path_rows:
            occ_id = row.get("occurrence_id", "")
            role = occ_by_id.get(occ_id, {}).get("role", row.get("role", "unknown"))
            if role == "intron":
                if previous is not None:
                    intron_between = True
                continue
            current_elements = sorted(set(elements_by_occ.get(occ_id, [])))
            if not current_elements:
                previous = None
                previous_occurrence = None
                intron_between = False
                continue
            current = current_elements[0]
            path_elements.append(current)
            if previous is not None:
                if previous == current:
                    site_id, coordinate, reference_id, method = boundary_site(
                        family, current, previous_occurrence, occ_id
                    )
                    boundary_rows.append(
                        {
                            "family_id": family,
                            "species": species,
                            "gene_copy_id": _copy,
                            "transcript_id": _transcript,
                            "site_id": site_id,
                            "element_id": current,
                            "projected_alignment_coordinate": coordinate,
                            "reference_occurrence_id": reference_id,
                            "projection_method": method,
                            "left_occurrence_id": previous_occurrence,
                            "right_occurrence_id": occ_id,
                        }
                    )
                else:
                    site_id = f"JG_{previous}__{current}"
                family_by_site[site_id] = family
                observations[site_id][species].add("present" if intron_between else "absent")
                evidence[site_id][species].add(
                    "annotated_intron_between_homologous_exons" if intron_between else "direct_exonic_adjacency"
                )
            previous = current
            previous_occurrence = occ_id
            intron_between = False
        for element in set(path_elements):
            element_counts[(family, element)][species].append(path_elements.count(element))

    for site_id, family in list(family_by_site.items()):
        if "_ALN_" not in site_id:
            continue
        element = site_id.split("_ALN_", 1)[0].removeprefix("JG_")
        for species, counts in element_counts[(family, element)].items():
            if counts and set(counts) == {1}:
                observations[site_id][species].add("absent")
                evidence[site_id][species].add("single_unsplit_exon_occurrence")

    rows = []
    for site_id, family in sorted(family_by_site.items()):
        for species in sorted(species_by_family[family]):
            values = observations[site_id].get(species, set())
            state = next(iter(values)) if len(values) == 1 else "unknown"
            ev = set(evidence[site_id].get(species, set()))
            if len(values) > 1:
                ev.add("conflicting_transcript_states")
            rows.append(
                {
                    "family_id": family,
                    "layer": "splice_junction",
                    "site_id": site_id,
                    "species": species,
                    "state": state,
                    "state_0": "absent",
                    "state_1": "present",
                    "evidence": ";".join(sorted(ev)) or "no_observation",
                }
            )
    write_tsv(
        Path(output_dir) / "splice_boundary_correspondence.tsv",
        boundary_rows,
        [
            "family_id", "species", "gene_copy_id", "transcript_id", "site_id", "element_id",
            "projected_alignment_coordinate", "reference_occurrence_id", "projection_method",
            "left_occurrence_id", "right_occurrence_id",
        ],
    )
    return rows


def build_structural_site_matrix(input_dir, output_dir):
    occurrences = read_tsv(
        Path(input_dir) / "segment_occurrences.tsv",
        ["occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status"],
    )
    element_rows = read_tsv(Path(output_dir) / "element_correspondence.tsv", optional=True)
    if not element_rows:
        raise SystemExit("element_correspondence.tsv is required before single-copy phylogenetic analysis")
    completion_rows = read_tsv(Path(output_dir) / "annotation_completion_candidates.tsv", optional=True)
    valid, species_by_family, excluded = _single_copy_families(occurrences)
    rows = _element_site_rows(
        occurrences,
        element_rows,
        completion_rows,
        valid,
        species_by_family,
    )
    rows.extend(_junction_site_rows(input_dir, output_dir, occurrences, element_rows, valid, species_by_family))
    sites_with_observed_state_1 = {
        (row["family_id"], row["layer"], row["site_id"])
        for row in rows
        if row["state"] == row["state_1"]
    }
    rows = [
        row for row in rows
        if (row["family_id"], row["layer"], row["site_id"]) in sites_with_observed_state_1
    ]
    return rows, excluded


@lru_cache(maxsize=32768)
def _transition_matrix(gain, loss, branch_length, multiplier=1.0):
    q = np.array(
        [[-gain * multiplier, gain * multiplier], [loss * multiplier, -loss * multiplier]],
        dtype=float,
    )
    return expm(q * branch_length)


def _root_prior(gain, loss, root_presence=None):
    if root_presence is not None:
        return np.array([1.0 - root_presence, root_presence], dtype=float)
    total = gain + loss
    if total <= 0:
        return np.array([0.5, 0.5], dtype=float)
    return np.array([loss / total, gain / total], dtype=float)


def _decode_parameters(theta, model, root_frequency="estimated", root_presence=0.5):
    values = np.exp(np.asarray(theta, dtype=float))
    if model == "ER":
        gain, loss, multiplier, offset = float(values[0]), float(values[0]), 1.0, 1
    elif model == "ARD":
        gain, loss, multiplier, offset = float(values[0]), float(values[1]), 1.0, 2
    elif model == "ARD_FOREGROUND":
        gain, loss, multiplier, offset = float(values[0]), float(values[1]), float(values[2]), 3
    else:
        raise ValueError(f"unknown model: {model}")
    if root_frequency == "stationary":
        rho = gain / max(gain + loss, 1e-300)
    elif root_frequency == "fixed":
        rho = float(root_presence)
    else:
        logit = float(theta[offset])
        rho = 1.0 / (1.0 + math.exp(-logit))
    return gain, loss, multiplier, rho


def _inside_messages(tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence=None):
    inside = {}
    cumulative_scale = {}
    for node in tree.postorder():
        if not tree.children.get(node):
            observed = observations.get(tree.label[node], "unknown")
            vector = np.array([1.0, 1.0], dtype=float)
            if observed == 0:
                vector = np.array([1.0, 0.0], dtype=float)
            elif observed == 1:
                vector = np.array([0.0, 1.0], dtype=float)
            scale = float(vector.sum())
            inside[node] = vector / scale
            cumulative_scale[node] = math.log(scale)
            continue
        vector = np.ones(2, dtype=float)
        child_scale = 0.0
        for child in tree.children[node]:
            multiplier = foreground_multiplier if child in foreground_children else 1.0
            matrix = _transition_matrix(gain, loss, tree.branch_length(child), multiplier)
            vector *= matrix @ inside[child]
            child_scale += cumulative_scale[child]
        scale = max(float(vector.sum()), 1e-300)
        inside[node] = vector / scale
        cumulative_scale[node] = child_scale + math.log(scale)
    prior = _root_prior(gain, loss, root_presence)
    root_term = max(float(prior @ inside[tree.root]), 1e-300)
    return inside, cumulative_scale[tree.root] + math.log(root_term)


def _pattern_log_likelihood(tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence=None):
    _inside, log_likelihood = _inside_messages(
        tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence
    )
    return log_likelihood


def _ascertainment_log_probability(tree, observations, gain, loss, foreground_multiplier, foreground_children, root_presence, mode):
    observed_labels = {label for label, value in observations.items() if value in {0, 1}}
    zero = {label: (0 if label in observed_labels else "unknown") for label in tree.leaf_by_label}
    one = {label: (1 if label in observed_labels else "unknown") for label in tree.leaf_by_label}
    p_zero = math.exp(_pattern_log_likelihood(tree, zero, gain, loss, foreground_multiplier, foreground_children, root_presence))
    excluded = p_zero
    if mode == "variable-only":
        p_one = math.exp(_pattern_log_likelihood(tree, one, gain, loss, foreground_multiplier, foreground_children, root_presence))
        excluded += p_one
    return math.log(max(1e-300, 1.0 - excluded))


def _dataset_log_likelihood(tree, patterns, model, theta, foreground_children, ascertainment, root_frequency, root_presence):
    gain, loss, multiplier, rho = _decode_parameters(theta, model, root_frequency, root_presence)
    total = 0.0
    for observations in patterns:
        value = _pattern_log_likelihood(tree, observations, gain, loss, multiplier, foreground_children, rho)
        if ascertainment in {"observed-at-least-one", "variable-only"}:
            value -= _ascertainment_log_probability(
                tree, observations, gain, loss, multiplier, foreground_children, rho, ascertainment
            )
        total += value
    return total


def _model_bounds(model, root_frequency="estimated"):
    rate_bounds = [(math.log(RATE_MIN), math.log(RATE_MAX))]
    if model in {"ARD", "ARD_FOREGROUND"}:
        rate_bounds.append((math.log(RATE_MIN), math.log(RATE_MAX)))
    if model == "ARD_FOREGROUND":
        rate_bounds.append((math.log(MULTIPLIER_MIN), math.log(MULTIPLIER_MAX)))
    if root_frequency == "estimated":
        rate_bounds.append((-13.8155095579, 13.8155095579))
    return rate_bounds


def _model_starts(model, root_frequency="estimated", root_presence=0.5):
    if model == "ER":
        starts = [[math.log(value)] for value in (0.01, 0.1, 1.0)]
    elif model == "ARD":
        starts = [[math.log(a), math.log(b)] for a, b in ((0.01, 0.1), (0.1, 0.1), (0.1, 1.0), (1.0, 0.1))]
    else:
        starts = [
        [math.log(a), math.log(b), math.log(m)]
        for a, b, m in ((0.01, 0.1, 0.5), (0.1, 0.1, 1.0), (0.1, 1.0, 2.0), (1.0, 0.1, 5.0))
        ]
    if root_frequency == "estimated":
        rho = min(1.0 - 1e-6, max(1e-6, float(root_presence)))
        root_logit = math.log(rho / (1.0 - rho))
        starts = [start + [root_logit] for start in starts]
    return starts


def _profile_interval(objective, optimum, index, bounds, max_log_likelihood, transform=math.exp):
    target = max_log_likelihood - PROFILE_DROP_95
    optimum = np.asarray(optimum, dtype=float)

    def profile_at(fixed):
        free_indices = [idx for idx in range(len(optimum)) if idx != index]
        if not free_indices:
            trial = optimum.copy()
            trial[index] = fixed
            return -float(objective(trial))

        def free_objective(free_values):
            trial = optimum.copy()
            trial[index] = fixed
            trial[free_indices] = free_values
            return objective(trial)

        start = optimum[free_indices]
        free_bounds = [bounds[idx] for idx in free_indices]
        result = minimize(free_objective, start, method="L-BFGS-B", bounds=free_bounds)
        return -float(result.fun)

    def crossing(direction):
        edge = bounds[index][0] if direction < 0 else bounds[index][1]
        points = np.linspace(optimum[index], edge, 18)[1:]
        previous_x = optimum[index]
        previous_value = max_log_likelihood - target
        for point in points:
            value = profile_at(float(point)) - target
            if value <= 0 <= previous_value:
                return brentq(lambda x: profile_at(x) - target, float(point), float(previous_x)), "closed"
            previous_x = float(point)
            previous_value = value
        return None, "open"

    lower, lower_status = crossing(-1)
    upper, upper_status = crossing(1)
    status = "two_sided"
    if lower_status == "open" and upper_status == "open":
        status = "unbounded"
    elif lower_status == "open":
        status = "lower_open"
    elif upper_status == "open":
        status = "upper_open"
    return {
        "low": transform(lower) if lower is not None else None,
        "high": transform(upper) if upper is not None else None,
        "status": status,
        "raw_low": lower,
        "raw_high": upper,
    }


def _profile_support_thetas(objective, optimum, bounds, intervals):
    optimum = np.asarray(optimum, dtype=float)
    support = [optimum.copy()]
    for index, interval in enumerate(intervals):
        for fixed in (interval.get("raw_low"), interval.get("raw_high")):
            if fixed is None:
                continue
            free_indices = [idx for idx in range(len(optimum)) if idx != index]
            trial = optimum.copy()
            trial[index] = fixed
            if free_indices:
                result = minimize(
                    lambda values: objective(np.array([
                        fixed if idx == index else values[free_indices.index(idx)]
                        for idx in range(len(optimum))
                    ])),
                    optimum[free_indices],
                    method="L-BFGS-B",
                    bounds=[bounds[idx] for idx in free_indices],
                )
                if result.success:
                    trial[free_indices] = result.x
            support.append(trial)
    return support


def _observed_information(objective, optimum):
    optimum = np.asarray(optimum, dtype=float)
    n = len(optimum)
    step = 1e-4
    hessian = np.zeros((n, n), dtype=float)
    center = float(objective(optimum))
    for i in range(n):
        unit_i = np.zeros(n)
        unit_i[i] = step
        hessian[i, i] = (objective(optimum + unit_i) - 2 * center + objective(optimum - unit_i)) / (step * step)
        for j in range(i + 1, n):
            unit_j = np.zeros(n)
            unit_j[j] = step
            value = (
                objective(optimum + unit_i + unit_j)
                - objective(optimum + unit_i - unit_j)
                - objective(optimum - unit_i + unit_j)
                + objective(optimum - unit_i - unit_j)
            ) / (4 * step * step)
            hessian[i, j] = hessian[j, i] = value
    eigenvalues = np.linalg.eigvalsh(hessian)
    identifiable = bool(np.all(np.isfinite(eigenvalues)) and np.min(eigenvalues) > 1e-7)
    condition = float(np.max(eigenvalues) / np.min(eigenvalues)) if identifiable else math.inf
    return identifiable and condition < 1e10, eigenvalues, condition


def fit_model(
    tree,
    patterns,
    model,
    foreground_children=frozenset(),
    ascertainment="observed-at-least-one",
    threads=1,
    root_frequency="estimated",
    root_presence=0.5,
):
    if root_frequency not in {"estimated", "stationary", "fixed"}:
        raise SystemExit("root_frequency must be estimated, stationary or fixed")
    if not 0 < float(root_presence) < 1:
        raise SystemExit("root_presence must lie strictly between 0 and 1")
    if ascertainment not in {"observed-at-least-one", "complete-universe", "variable-only"}:
        raise SystemExit("unsupported ascertainment mode")
    if ascertainment == "observed-at-least-one":
        for pattern in patterns:
            if not any(value == 1 for value in pattern.values()):
                raise SystemExit("observed-at-least-one ascertainment requires each structural site to be present in at least one observed tip")
    if ascertainment == "variable-only":
        for pattern in patterns:
            observed = {value for value in pattern.values() if value in {0, 1}}
            if len(observed) < 2:
                raise SystemExit("variable-only ascertainment requires every included site to vary among observed tips")
    bounds = _model_bounds(model, root_frequency)

    def objective(theta):
        return -_dataset_log_likelihood(
            tree, patterns, model, theta, frozenset(foreground_children), ascertainment,
            root_frequency, root_presence,
        )

    starts = _model_starts(model, root_frequency, root_presence)
    if threads > 1 and len(starts) > 1:
        with ThreadPoolExecutor(max_workers=min(threads, len(starts))) as executor:
            results = list(
                executor.map(
                    lambda start: minimize(objective, start, method="L-BFGS-B", bounds=bounds),
                    starts,
                )
            )
    else:
        results = [minimize(objective, start, method="L-BFGS-B", bounds=bounds) for start in starts]
    best = min(results, key=lambda result: float(result.fun))
    theta = np.asarray(best.x, dtype=float)
    log_likelihood = -float(best.fun)
    gain, loss, multiplier, rho = _decode_parameters(theta, model, root_frequency, root_presence)
    boundary = any(
        abs(theta[idx] - lower) < 1e-5 or abs(theta[idx] - upper) < 1e-5
        for idx, (lower, upper) in enumerate(bounds)
    )
    intervals = []
    for idx in range(len(theta)):
        transform = math.exp
        if root_frequency == "estimated" and idx == len(theta) - 1:
            transform = lambda value: 1.0 / (1.0 + math.exp(-value))
        intervals.append(_profile_interval(objective, theta, idx, bounds, log_likelihood, transform))
    identifiable, information_eigenvalues, information_condition = _observed_information(objective, theta)
    support_thetas = _profile_support_thetas(objective, theta, bounds, intervals)
    return {
        "model": model,
        "theta": theta,
        "gain_rate": gain,
        "loss_rate": loss,
        "foreground_multiplier": multiplier,
        "root_presence": rho,
        "log_likelihood": log_likelihood,
        "parameter_count": len(theta),
        "aic": 2 * len(theta) - 2 * log_likelihood,
        "converged": bool(best.success),
        "optimizer_message": str(best.message),
        "boundary": boundary,
        "intervals": intervals,
        "start_count": len(results),
        "identifiable": identifiable,
        "information_eigenvalues": information_eigenvalues,
        "information_condition": information_condition,
        "root_frequency": root_frequency,
        "profile_support_thetas": support_thetas,
    }


def _posterior_messages(tree, observations, gain, loss, multiplier, foreground_children, root_presence=None):
    inside, _log_likelihood = _inside_messages(
        tree, observations, gain, loss, multiplier, foreground_children, root_presence
    )
    prior = _root_prior(gain, loss, root_presence)
    outside = {tree.root: prior / prior.sum()}
    node = {}
    edge = {}
    for current in tree.preorder():
        weights = outside[current] * inside[current]
        node[current] = weights / max(float(weights.sum()), 1e-300)
        for child in tree.children.get(current, []):
            sibling_product = np.ones(2, dtype=float)
            for sibling in tree.children[current]:
                if sibling == child:
                    continue
                sibling_multiplier = multiplier if sibling in foreground_children else 1.0
                sibling_matrix = _transition_matrix(
                    gain, loss, tree.branch_length(sibling), sibling_multiplier
                )
                sibling_product *= sibling_matrix @ inside[sibling]
            child_multiplier = multiplier if child in foreground_children else 1.0
            matrix = _transition_matrix(gain, loss, tree.branch_length(child), child_multiplier)
            parent_context = outside[current] * sibling_product
            joint = parent_context[:, None] * matrix * inside[child][None, :]
            joint /= max(float(joint.sum()), 1e-300)
            edge[(current, child)] = joint
            child_outside = parent_context @ matrix
            outside[child] = child_outside / max(float(child_outside.sum()), 1e-300)
    return node, edge


@lru_cache(maxsize=32768)
def _conditional_transition_count(gain, loss, branch_length, multiplier, start, end, src, dst):
    matrix = _transition_matrix(gain, loss, branch_length, multiplier)
    denominator = max(float(matrix[start, end]), 1e-300)
    q_value = gain * multiplier if (src, dst) == (0, 1) else loss * multiplier

    def integrand(time):
        left = _transition_matrix(gain, loss, time, multiplier)
        right = _transition_matrix(gain, loss, branch_length - time, multiplier)
        return float(left[start, src]) * q_value * float(right[dst, end])

    integral = quad(integrand, 0.0, branch_length, epsabs=1e-9, epsrel=1e-7)[0]
    return max(0.0, integral / denominator)


def _expected_transition_count(joint, gain, loss, branch_length, multiplier, src, dst):
    total = 0.0
    for start in (0, 1):
        for end in (0, 1):
            total += float(joint[start, end]) * _conditional_transition_count(
                gain, loss, branch_length, multiplier, start, end, src, dst
            )
    return total


def _bh_adjust(rows):
    by_test = defaultdict(list)
    for index, row in enumerate(rows):
        try:
            by_test[row["test_id"]].append((index, float(row["p_value"])))
        except (TypeError, ValueError):
            row["q_value"] = "NA"
            row["q_value_method"] = "not_available"
    for entries in by_test.values():
        if len(entries) == 1:
            rows[entries[0][0]]["q_value"] = _fmt(entries[0][1])
            rows[entries[0][0]]["q_value_method"] = "Benjamini-Hochberg_m_equals_1"
            continue
        ordered = sorted(entries, key=lambda item: item[1])
        adjusted = [0.0] * len(ordered)
        running = 1.0
        for rank in range(len(ordered), 0, -1):
            value = min(running, ordered[rank - 1][1] * len(ordered) / rank)
            adjusted[rank - 1] = value
            running = value
        for (index, _p), value in zip(ordered, adjusted):
            rows[index]["q_value"] = _fmt(value)
            rows[index]["q_value_method"] = "Benjamini-Hochberg_within_test_type"


def _read_foreground_children(path, tree):
    if not path:
        return frozenset()
    rows = read_tsv(path)
    label_to_node = {label: node for node, label in tree.label.items()}
    children = set()
    invalid = []
    edge_set = set(tree.edges())
    for row in rows:
        parent = row.get("parent_node") or row.get("parent_id") or ""
        child = row.get("child_node") or row.get("child_id") or ""
        if (parent, child) in edge_set:
            children.add(child)
            continue
        scope = row.get("branch_scope") or row.get("branch") or ""
        if "->" in scope:
            parent_label, child_label = (part.strip() for part in scope.split("->", 1))
            parent_node = label_to_node.get(parent_label)
            child_node = label_to_node.get(child_label)
            if (parent_node, child_node) in edge_set:
                children.add(child_node)
                continue
        invalid.append(scope or f"{parent}->{child}")
    if invalid:
        raise SystemExit("foreground branch file contains unmatched branches: " + ", ".join(invalid))
    if not children:
        raise SystemExit("foreground branch file did not match any branch in species_tree.tsv")
    if len(children) == len(edge_set):
        raise SystemExit("foreground model is unidentifiable when every branch is foreground")
    return frozenset(children)


def infer_single_copy_phylogeny(
    input_dir,
    output_dir,
    model="er-ard",
    foreground_branches=None,
    branch_length_mode="supplied",
    ascertainment="observed-at-least-one",
    threads=1,
    root_frequency="estimated",
    root_presence=0.5,
):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    site_rows, excluded_rows = build_structural_site_matrix(input_dir, output_dir)
    write_tsv(
        output_dir / "excluded_families.tsv",
        excluded_rows,
        ["family_id", "species", "copy_count", "gene_copy_ids", "reason"],
    )
    write_tsv(
        output_dir / "structural_site_matrix.tsv",
        site_rows,
        ["family_id", "layer", "site_id", "species", "state", "state_0", "state_1", "evidence"],
    )
    tree_rows = _validated_tree_rows(input_dir / "species_tree.tsv", branch_length_mode)
    tree = SpeciesTree(tree_rows)
    matrix_species = {row["species"] for row in site_rows}
    missing_tips = sorted(matrix_species - set(tree.leaf_by_label))
    if missing_tips:
        raise SystemExit(
            "species_tree.tsv lacks structural-matrix species: " + ", ".join(missing_tips)
        )
    foreground_children = _read_foreground_children(foreground_branches, tree)
    if model == "foreground" and not foreground_children:
        raise SystemExit("--model foreground requires --foreground-branches")

    grouped = defaultdict(lambda: defaultdict(dict))
    state_labels = {}
    for row in site_rows:
        grouped[(row["family_id"], row["layer"])][row["site_id"]][row["species"]] = row["state"]
        state_labels[(row["family_id"], row["layer"])] = (row["state_0"], row["state_1"])

    fit_rows = []
    test_rows = []
    node_rows = []
    branch_rows = []
    change_rows = []
    for (family, layer), sites in sorted(grouped.items()):
        state_0, state_1 = state_labels[(family, layer)]
        patterns = []
        encoded_sites = []
        for site_id, observations in sorted(sites.items()):
            encoded = {
                species: 0 if state == state_0 else 1 if state == state_1 else "unknown"
                for species, state in observations.items()
            }
            observed_count = sum(value in {0, 1} for value in encoded.values())
            if observed_count >= 2:
                patterns.append(encoded)
                encoded_sites.append((site_id, encoded))
        if not patterns:
            continue
        informative = sum(
            1 for pattern in patterns if len({value for value in pattern.values() if value in {0, 1}}) > 1
        )
        observed_taxa = len(
            {species for pattern in patterns for species, value in pattern.items() if value in {0, 1}}
        )
        if model == "foreground":
            null_fit = fit_model(tree, patterns, "ARD", ascertainment=ascertainment, threads=threads, root_frequency=root_frequency, root_presence=root_presence)
            alternative_fit = fit_model(
                tree,
                patterns,
                "ARD_FOREGROUND",
                foreground_children=foreground_children,
                ascertainment=ascertainment,
                threads=threads,
                root_frequency=root_frequency,
                root_presence=root_presence,
            )
            test_id = "homogeneous_vs_foreground"
        else:
            null_fit = fit_model(tree, patterns, "ER", ascertainment=ascertainment, threads=threads, root_frequency=root_frequency, root_presence=root_presence)
            alternative_fit = fit_model(tree, patterns, "ARD", ascertainment=ascertainment, threads=threads, root_frequency=root_frequency, root_presence=root_presence)
            test_id = "equal_rates_vs_gain_loss"

        for fit in (null_fit, alternative_fit):
            gain_ci = fit["intervals"][0]
            loss_ci = fit["intervals"][0] if fit["model"] == "ER" else fit["intervals"][1]
            multiplier_ci = fit["intervals"][2] if fit["model"] == "ARD_FOREGROUND" else {"low": 1.0, "high": 1.0, "status": "fixed"}
            root_ci = fit["intervals"][-1] if root_frequency == "estimated" else {"low": fit["root_presence"], "high": fit["root_presence"], "status": root_frequency}
            fit_rows.append(
                {
                    "family_id": family,
                    "layer": layer,
                    "model": fit["model"],
                    "n_taxa": observed_taxa,
                    "n_structural_sites": len(patterns),
                    "n_informative_patterns": informative,
                    "gain_rate": _fmt(fit["gain_rate"]),
                    "gain_rate_ci_low": _fmt(gain_ci["low"]),
                    "gain_rate_ci_high": _fmt(gain_ci["high"]),
                    "gain_rate_ci_status": gain_ci["status"],
                    "loss_rate": _fmt(fit["loss_rate"]),
                    "loss_rate_ci_low": _fmt(loss_ci["low"]),
                    "loss_rate_ci_high": _fmt(loss_ci["high"]),
                    "loss_rate_ci_status": loss_ci["status"],
                    "foreground_multiplier": _fmt(fit["foreground_multiplier"]),
                    "foreground_multiplier_ci_low": _fmt(multiplier_ci["low"]),
                    "foreground_multiplier_ci_high": _fmt(multiplier_ci["high"]),
                    "foreground_multiplier_ci_status": multiplier_ci["status"],
                    "root_presence": _fmt(fit["root_presence"]),
                    "root_presence_ci_low": _fmt(root_ci["low"]),
                    "root_presence_ci_high": _fmt(root_ci["high"]),
                    "root_presence_ci_status": root_ci["status"],
                    "root_frequency_mode": root_frequency,
                    "log_likelihood": _fmt(fit["log_likelihood"]),
                    "parameter_count": fit["parameter_count"],
                    "aic": _fmt(fit["aic"]),
                    "converged": str(fit["converged"]).lower(),
                    "parameter_at_boundary": str(fit["boundary"]).lower(),
                    "identifiable": str(fit["identifiable"]).lower(),
                    "information_condition": _fmt(fit["information_condition"]),
                    "optimizer_starts": fit["start_count"],
                    "optimizer_message": fit["optimizer_message"],
                    "ascertainment": ascertainment,
                }
            )

        likelihood_difference = alternative_fit["log_likelihood"] - null_fit["log_likelihood"]
        lrt = 2.0 * likelihood_difference if likelihood_difference >= -1e-7 else None
        estimable = (
            informative >= 1
            and null_fit["converged"]
            and alternative_fit["converged"]
            and not null_fit["boundary"]
            and not alternative_fit["boundary"]
            and null_fit["identifiable"]
            and alternative_fit["identifiable"]
            and lrt is not None
        )
        p_value = float(chi2.sf(max(0.0, lrt), 1)) if estimable else None
        if likelihood_difference < -1e-7:
            test_status = "optimization_failure_alternative_below_null"
        elif estimable:
            test_status = "tested"
        else:
            test_status = "parameters_not_estimable"
        test_row = {
            "family_id": family,
            "layer": layer,
            "test_id": test_id,
            "null_model": null_fit["model"],
            "alternative_model": alternative_fit["model"],
            "null_log_likelihood": _fmt(null_fit["log_likelihood"]),
            "alternative_log_likelihood": _fmt(alternative_fit["log_likelihood"]),
            "lrt_statistic": _fmt(lrt),
            "df": 1,
            "p_value": _fmt(p_value),
            "q_value": "NA",
            "q_value_method": "not_available",
            "test_status": test_status,
            "reference_distribution": "chi_square_df1_asymptotic" if estimable else "not_available",
            "n_taxa": observed_taxa,
            "n_structural_sites": len(patterns),
            "n_informative_patterns": informative,
        }
        test_rows.append(test_row)
        eligible_fits = [fit for fit in (null_fit, alternative_fit) if fit["converged"] and fit["identifiable"]]
        if informative == 0:
            selected_fit = null_fit
        else:
            selected_fit = min(eligible_fits, key=lambda fit: fit["aic"]) if eligible_fits else null_fit

        for site_id, observations in encoded_sites:
            node_posterior, edge_posterior = _posterior_messages(
                tree,
                observations,
                selected_fit["gain_rate"],
                selected_fit["loss_rate"],
                selected_fit["foreground_multiplier"],
                foreground_children,
                selected_fit["root_presence"],
            )
            support_posteriors = []
            for support_theta in selected_fit["profile_support_thetas"]:
                support_gain, support_loss, support_multiplier, support_root = _decode_parameters(
                    support_theta, selected_fit["model"], root_frequency, root_presence
                )
                support_posteriors.append(
                    _posterior_messages(
                        tree, observations, support_gain, support_loss, support_multiplier,
                        foreground_children, support_root,
                    )
                )
            for node_id, probabilities in node_posterior.items():
                for index, probability in enumerate(probabilities):
                    sensitivity = [float(nodes[node_id][index]) for nodes, _edges in support_posteriors]
                    node_rows.append(
                        {
                            "family_id": family,
                            "layer": layer,
                            "site_id": site_id,
                            "node_id": node_id,
                            "node_label": tree.label[node_id],
                            "state": state_0 if index == 0 else state_1,
                            "posterior_probability": _fmt(probability),
                            "profile_probability_low": _fmt(min(sensitivity)),
                            "profile_probability_high": _fmt(max(sensitivity)),
                            "model": selected_fit["model"],
                            "conditioning": "empirical_Bayes_conditional_on_MLE",
                        }
                    )
            for (parent, child), joint in edge_posterior.items():
                branch_multiplier = (
                    selected_fit["foreground_multiplier"] if child in foreground_children else 1.0
                )
                expected_gain = _expected_transition_count(
                    joint,
                    selected_fit["gain_rate"],
                    selected_fit["loss_rate"],
                    tree.branch_length(child),
                    branch_multiplier,
                    0,
                    1,
                )
                expected_loss = _expected_transition_count(
                    joint,
                    selected_fit["gain_rate"],
                    selected_fit["loss_rate"],
                    tree.branch_length(child),
                    branch_multiplier,
                    1,
                    0,
                )
                gain_probability = float(joint[0, 1])
                loss_probability = float(joint[1, 0])
                change_probability = gain_probability + loss_probability
                support_joints = [edges[(parent, child)] for _nodes, edges in support_posteriors]
                gain_sensitivity = [float(value[0, 1]) for value in support_joints]
                loss_sensitivity = [float(value[1, 0]) for value in support_joints]
                change_sensitivity = [gain + loss for gain, loss in zip(gain_sensitivity, loss_sensitivity)]
                for src, dst, probability in (
                    (state_0, state_1, gain_probability),
                    (state_1, state_0, loss_probability),
                ):
                    branch_rows.append(
                        {
                            "family_id": family,
                            "layer": layer,
                            "site_id": site_id,
                            "parent_node": parent,
                            "child_node": child,
                            "parent_label": tree.label[parent],
                            "child_label": tree.label[child],
                            "branch_length": _fmt(tree.branch_length(child)),
                            "from_state": src,
                            "to_state": dst,
                            "endpoint_transition_probability": _fmt(probability),
                            "profile_transition_probability_low": _fmt(min(gain_sensitivity if src == state_0 else loss_sensitivity)),
                            "profile_transition_probability_high": _fmt(max(gain_sensitivity if src == state_0 else loss_sensitivity)),
                            "total_endpoint_change_probability": _fmt(change_probability),
                            "profile_total_change_probability_low": _fmt(min(change_sensitivity)),
                            "profile_total_change_probability_high": _fmt(max(change_sensitivity)),
                            "expected_gain_count": _fmt(expected_gain),
                            "expected_loss_count": _fmt(expected_loss),
                            "model": selected_fit["model"],
                            "conditioning": "empirical_Bayes_conditional_on_MLE",
                        }
                    )
                change_rows.append(
                    {
                        "family_id": family,
                        "layer": layer,
                        "site_id": site_id,
                        "parent_node": parent,
                        "child_node": child,
                        "branch_scope": f"{tree.label[parent]}->{tree.label[child]}",
                        "structural_change_type": "bidirectional_transition_probabilities",
                        "structural_pattern": f"{state_0}<->{state_1}",
                        "endpoint_change_probability": _fmt(change_probability),
                        "gain_endpoint_probability": _fmt(gain_probability),
                        "loss_endpoint_probability": _fmt(loss_probability),
                        "direction_probability": "NA",
                        "expected_gain_count": _fmt(expected_gain),
                        "expected_loss_count": _fmt(expected_loss),
                        "model": selected_fit["model"],
                        "rate_test_status": test_status,
                    }
                )

    _bh_adjust(test_rows)
    write_tsv(
        output_dir / "model_fits.tsv",
        fit_rows,
        [
            "family_id", "layer", "model", "n_taxa", "n_structural_sites", "n_informative_patterns",
            "gain_rate", "gain_rate_ci_low", "gain_rate_ci_high", "gain_rate_ci_status",
            "loss_rate", "loss_rate_ci_low", "loss_rate_ci_high", "loss_rate_ci_status",
            "foreground_multiplier", "foreground_multiplier_ci_low", "foreground_multiplier_ci_high",
            "foreground_multiplier_ci_status", "root_presence", "root_presence_ci_low",
            "root_presence_ci_high", "root_presence_ci_status", "root_frequency_mode", "log_likelihood",
            "parameter_count", "aic", "converged", "parameter_at_boundary", "identifiable",
            "information_condition", "optimizer_starts", "optimizer_message", "ascertainment",
        ],
    )
    write_tsv(
        output_dir / "model_tests.tsv",
        test_rows,
        [
            "family_id", "layer", "test_id", "null_model", "alternative_model", "null_log_likelihood",
            "alternative_log_likelihood", "lrt_statistic", "df", "p_value", "q_value", "q_value_method",
            "test_status", "reference_distribution", "n_taxa", "n_structural_sites", "n_informative_patterns",
        ],
    )
    write_tsv(
        output_dir / "node_state_posteriors.tsv",
        node_rows,
        ["family_id", "layer", "site_id", "node_id", "node_label", "state", "posterior_probability",
         "profile_probability_low", "profile_probability_high", "model", "conditioning"],
    )
    write_tsv(
        output_dir / "branch_transition_posteriors.tsv",
        branch_rows,
        [
            "family_id", "layer", "site_id", "parent_node", "child_node", "parent_label", "child_label",
            "branch_length", "from_state", "to_state", "endpoint_transition_probability",
            "profile_transition_probability_low", "profile_transition_probability_high",
            "total_endpoint_change_probability", "profile_total_change_probability_low",
            "profile_total_change_probability_high", "expected_gain_count", "expected_loss_count", "model", "conditioning",
        ],
    )
    write_tsv(
        output_dir / "structural_changes.tsv",
        change_rows,
        [
            "family_id", "layer", "site_id", "parent_node", "child_node", "branch_scope",
            "structural_change_type", "structural_pattern", "endpoint_change_probability",
            "gain_endpoint_probability", "loss_endpoint_probability", "direction_probability",
            "expected_gain_count", "expected_loss_count", "model", "rate_test_status",
        ],
    )
    write_tsv(
        output_dir / "phylogeny_scope.tsv",
        [{
            "scope": "single_copy_structural_sites",
            "tree_file": "species_tree.tsv",
            "tree_scope": "species_tree",
            "layers": "exon_presence;exon_role;splice_junction",
            "note": "All homologous structural sites within each layer share fitted gain/loss parameters.",
        }],
        ["scope", "tree_file", "tree_scope", "layers", "note"],
    )
    parameters = {
        "analysis_scope": "single-copy",
        "model_test": model,
        "branch_length_mode": branch_length_mode,
        "ascertainment": ascertainment,
        "root_frequency": root_frequency,
        "root_presence_when_fixed": root_presence if root_frequency == "fixed" else None,
        "tree_file": str(input_dir / "species_tree.tsv"),
        "foreground_branches": str(foreground_branches) if foreground_branches else None,
        "threads": max(1, int(threads)),
        "fixed_inputs": ["species_tree", "ortholog_set", "structural_site_states"],
        "estimated_parameters": ["gain_rate", "loss_rate"] + (["root_presence"] if root_frequency == "estimated" else []) + (["foreground_multiplier"] if model == "foreground" else []),
        "excluded_family_count": len({row["family_id"] for row in excluded_rows}),
    }
    (output_dir / "run_parameters.json").write_text(json.dumps(parameters, indent=2, sort_keys=True) + "\n")
    return site_rows, fit_rows, test_rows, change_rows
