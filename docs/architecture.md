# IntraPhy 0.19.1 architecture

![Implementation ownership](figures/architecture.svg)

Ordinary functions and immutable domain objects, not a generic workflow framework.
One exon or a dependent local exon configuration is the structural model object.

| Responsibility | Owner |
|---|---|
| Native FASTA/GFF coordinates and target selection | `inputs/`, `preparation/`, `structure/native.py` |
| Full supplied-locus MSA and corroboration | `structure/alignment.py`, `corroboration.py` |
| Whole-configuration annotation alternatives | `structure/alternatives.py` |
| Candidate assembly and observation compatibility | `structure/build.py`, `observations.py` |
| Coordinate-identity validation at every entry | `structure/validation.py`, called by `serialization.py` |
| Immutable exons, source payloads and configurations | `structure/types.py`, `material.py` |
| Elementary edits, consequences and opportunity weights | `structure/edits.py` |
| Finite state catalogue and exact geometry preflight | `structure/space.py` |
| Canonical branches and declared source opportunities | `structure/tree_context.py`, `origins.py` |
| Sparse directed edit distances | `structure/paths.py` |
| All-optimal parsimony | `inference/configuration_dp.py`, `configuration_history.py` |
| Finite CTMC and bounded run-local kernels | `inference/configuration_ctmc.py`, `configuration_model.py`, `kernel_cache.py` |
| Pooled scales and whole-gene resampling | `inference/exon_rates.py`, `exon_resampling.py` |
| Read-only result diagrams | `reporting/exon_results.py`, `exon_drawing.py` |
| Computed teaching examples and four current plates | `reporting/exon_guide*.py` |

`analyze` prepares files and calls configuration inference. Explicit schema-2
catalogues bypass alignment, not validation. Pure numerical functions never read
GFF, repair annotation or generate figures. The renderer reads saved results;
only the explicitly synthetic teaching-example generator runs inference.

Legacy V18 P/R/J modules remain explicit baselines; they are not silently used by
configuration inference. Current and legacy catalogue/model identifiers differ.
One entry point cannot assign extra weight to the same coordinates by renaming a
unit. A failure cannot authorize reuse of old successful probability output.

See [audit resolutions](v0191_audit_resolution.md), [model](exon_structure_model.md)
and [validation](v0191_validation.md). Source modules remain under 500 lines.
