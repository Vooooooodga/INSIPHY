"""The sole preparation boundary for formal structural observation matrices.

There is no dependency on either inference engine. Files remain optional public
checkpoints, not a mutable message bus between the numeric engines.
"""
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Mapping, Tuple

from insiphy.observations.schema import read_structural_site_matrix, write_structural_site_matrix
from insiphy.storage.tabular import read_tsv


@dataclass(frozen=True)
class ObservationMatrix:
    rows: Tuple[Mapping[str, str], ...]
    excluded: Tuple[Mapping[str, str], ...]
    source: Path
    mode: str

    @classmethod
    def snapshot(cls, rows, excluded, source, mode):
        return cls(tuple(MappingProxyType(dict(row)) for row in rows),
                   tuple(MappingProxyType(dict(row)) for row in excluded),
                   Path(source), mode)

    def as_legacy_tuple(self):
        # These views are read-only; row ownership is no longer shared with callers.
        return self.rows, self.excluded, self.source, self.mode


def prepare_observation_matrix(input_dir, output_dir, source=None, annotation_view="repertoire"):
    from insiphy.structural_sites import build_structural_site_matrix, structural_matrix_annotation_view
    from insiphy.observations.schema import _validate_structural_site_rows, _normalize_structural_site_contract, normalize_structural_site_row

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
    _validate_structural_site_rows(rows)
    view = structural_matrix_annotation_view(rows)
    if view not in {"view_independent", annotation_view}:
        raise SystemExit(f"requested annotation view does not match the frozen structural-site matrix: requested={annotation_view}, matrix={view}")
    matrix = ObservationMatrix.snapshot(rows, excluded, source, mode)
    write_structural_site_matrix(target, matrix.rows)
    return matrix


def load_site_matrix(input_dir, output_dir, structural_site_matrix_path, annotation_view):
    """Compatibility entry point; new workflows pass ObservationMatrix directly."""
    return prepare_observation_matrix(input_dir, output_dir, structural_site_matrix_path,
                                      annotation_view).as_legacy_tuple()
