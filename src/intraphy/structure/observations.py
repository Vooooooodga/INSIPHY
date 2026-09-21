"""Read-only compatibility between structural observations and candidate states."""
from __future__ import annotations

from itertools import product
import numpy as np
from .types import Catalogue, ExonConfiguration, ExonSpan, ObservationEvidence
from .material import normalize_exons
from .space import StateSpace


def _outside(exons: tuple[ExonSpan, ...], windows: tuple[ExonSpan, ...]) -> tuple[tuple[int, int], ...]:
    """Compare exon structure outside unknown windows without joining across them."""
    parts = []
    for exon in exons:
        remaining = [(exon.start, exon.end)]
        for w in windows:
            cut = []
            for a, b in remaining:
                if b <= w.start or a >= w.end:
                    cut.append((a, b))
                else:
                    if a < w.start:
                        cut.append((a, w.start))
                    if b > w.end:
                        cut.append((w.end, b))
            remaining = cut
        parts.extend(remaining)
    return tuple(parts)


def compatibility(space: StateSpace, observation: ObservationEvidence,
                  view: str = "evidence", selected: ExonConfiguration | None = None) -> np.ndarray:
    if view not in {"evidence", "annotation"}:
        raise ValueError("Observation view must be evidence or annotation")
    if observation.kind == "excluded":
        return np.ones(len(space.states))
    if observation.kind == "coexisting" and selected is None:
        raise ValueError("Coexisting annotated structures require explicit conditional scenarios")
    natives = (selected,) if selected is not None else observation.configurations
    candidates = list(natives)
    if view == "evidence":
        keys = {v.key for v in natives}
        candidates.extend(a.configuration for a in observation.alternatives if a.replaces_key in keys)
    allowed = np.zeros(len(space.states))
    for i, state in enumerate(space.states):
        if any(p is not None and (state.material[k] == 1) != (p == 1)
               for k, p in enumerate(observation.material_presence)):
            continue
        if observation.kind == "unknown":
            allowed[i] = 1
            continue
        windows = observation.unknown_intervals
        if view == "evidence":
            windows = tuple(sorted(set((*windows, *observation.alternative_exons))))
        for native in candidates:
            # Observed absence never selects unintroduced versus subsequently lost.
            native_exons = normalize_exons(space.catalogue, state.material, native.exons)
            if _outside(state.exons, windows) == _outside(native_exons, windows):
                allowed[i] = 1
                break
    return allowed


def observation_scenarios(space: StateSpace, taxa: tuple[str, ...], view: str,
                          max_scenarios: int = 64):
    by_species = {o.species: o for o in space.catalogue.observations}
    coexist = [(s, by_species[s].configurations) for s in taxa
               if s in by_species and by_species[s].kind == "coexisting"]
    count = 1
    for _, choices in coexist:
        count *= len(choices)
    if count > max_scenarios:
        raise ValueError("observation_scenarios_incomplete: too many coexisting structures")
    scenarios = []
    for choices in product(*(cs for _, cs in coexist)):
        selected = dict(zip((s for s, _ in coexist), choices))
        weights = {}
        for s in taxa:
            if s not in by_species:
                weights[s] = np.ones(len(space.states))
            else:
                weights[s] = compatibility(space, by_species[s], view, selected.get(s))
            if not weights[s].any():
                raise ValueError(f"Observation incompatible with candidate catalogue: {s}")
        label = ";".join(f"{s}={config.key}" for s, config in sorted(selected.items())) or "single_structure"
        scenarios.append((label, weights))
    return scenarios
