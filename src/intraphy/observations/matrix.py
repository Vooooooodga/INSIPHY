"""The sole preparation boundary for formal structural observation matrices.

There is no dependency on either inference engine. Files remain optional public
checkpoints, not a mutable message bus between the numeric engines.
"""
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Mapping, Tuple

from intraphy.observations.schema import read_structural_site_matrix, write_structural_site_matrix
from intraphy.storage.tabular import read_tsv


@dataclass(frozen=True)
class ObservationMatrix:
    rows: Tuple[Mapping[str, str], ...]
    excluded: Tuple[Mapping[str, str], ...]
    source: Path
    mode: str
    full_rows: Tuple[Mapping[str, str], ...] = ()
    selection: Mapping = None

    @classmethod
    def snapshot(cls, rows, excluded, source, mode, full_rows=None, selection=None):
        return cls(tuple(MappingProxyType(dict(row)) for row in rows),
                   tuple(MappingProxyType(dict(row)) for row in excluded),
                   Path(source), mode,
                   tuple(MappingProxyType(dict(row)) for row in (rows if full_rows is None else full_rows)),
                   MappingProxyType(dict(selection or {})))

    def as_legacy_tuple(self):
        # These views are read-only; row ownership is no longer shared with callers.
        return self.rows, self.excluded, self.source, self.mode


def prepare_observation_matrix(input_dir, output_dir, source=None, annotation_view="repertoire",
                               *, analysis_range="all", min_callable_fraction=0.70):
    from intraphy.structural_sites import build_structural_site_matrix, structural_matrix_annotation_view
    from intraphy.observations.schema import _validate_structural_site_rows, _normalize_structural_site_contract, normalize_structural_site_row

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    target = output_dir / "structural_site_matrix.tsv"
    if source is not None:
        source = Path(source)
        rows = read_structural_site_matrix(source)
        excluded = read_tsv(source.parent / "excluded_families.tsv", optional=True)
        mode = "frozen_matrix"
    else:
        rows, excluded = build_structural_site_matrix(input_dir, output_dir, annotation_view=annotation_view)
        source, mode = target, "prepared_in_analysis"
    rows = [_normalize_structural_site_contract(normalize_structural_site_row(
        _normalize_structural_site_contract(dict(row))), upgrade_schema=False) for row in rows]
    from .characters import (coordinate_evidence, add_interval_dependencies,
                             validate_sequence_applicability, write_character_catalogue)
    evidence_dir = source.parent if mode == "frozen_matrix" else output_dir
    coordinates = coordinate_evidence(rows, evidence_dir, frozen=mode == "frozen_matrix")
    rows = add_interval_dependencies(rows, coordinates)
    _validate_structural_site_rows(rows)
    validate_sequence_applicability(rows)
    view = structural_matrix_annotation_view(rows)
    if view not in {"view_independent", annotation_view}:
        raise SystemExit(f"requested annotation view does not match the frozen structural-site matrix: requested={annotation_view}, matrix={view}")
    from intraphy.observations.eligibility import ScopePolicy, assess_scope, write_scope_reports
    from intraphy.observations.schema import validate_structural_site_tip_rows
    from intraphy.topology import SpeciesTree
    policy = ScopePolicy(analysis_range, min_callable_fraction)
    tree_path = Path(input_dir) / "species_tree.tsv"
    tree = SpeciesTree(read_tsv(tree_path)) if tree_path.exists() else None
    if tree is not None:
        validate_structural_site_tip_rows(rows, tree.leaf_by_label)
    full, selected, audit, summaries = assess_scope(
        rows, tree.leaf_by_label if tree else None, policy, tree)
    selection = write_scope_reports(output_dir, audit, summaries, policy)
    write_character_catalogue(output_dir, full, coordinates)
    matrix = ObservationMatrix.snapshot(selected, excluded, source, mode, full, selection)
    # The complete background always survives, even if a high-coverage fit is requested.
    write_structural_site_matrix(target, full)
    write_structural_site_matrix(output_dir / "analysis_structural_site_matrix.tsv", selected)
    return matrix


def load_site_matrix(input_dir, output_dir, structural_site_matrix_path, annotation_view,
                     *, analysis_range="all", min_callable_fraction=0.70):
    """Compatibility entry point; new workflows pass ObservationMatrix directly."""
    return prepare_observation_matrix(input_dir, output_dir, structural_site_matrix_path,
                                      annotation_view, analysis_range=analysis_range,
                                      min_callable_fraction=min_callable_fraction).as_legacy_tuple()
