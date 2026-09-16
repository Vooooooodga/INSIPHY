# Publication Readiness Review

Date: 2026-09-16

## Decision

The single-copy computational framework is implemented and testable. The
present evidence supports a software prototype and methods development report.
The project still needs multiple curated positive single-copy genes and a
larger conserved-control panel before a general biological performance claim.

## Scientific checks completed

- Gene-level orthology is fixed upstream.
- One gene per species is enforced in the formal path.
- Exon sequence presence, exonic role, and splice-junction presence are
  separated.
- Missing annotation and unknown sequence evidence are marginalized.
- Model parameters are shared across homologous sites within each layer.
- Null and alternative likelihood models are explicit.
- LRT P values are withheld when directions are not estimable or fits reach a
  boundary.
- Ancestral and branch histories are reported as posterior probabilities.
- Molecular mechanism assignment is left to downstream biological analysis.

## Software checks completed

- CLI and Python library use the same implementation.
- Internal alignment uses `PairwiseAligner`; minimap2 and miniprot remain
  optional.
- Independent optimizer starts can use multiple threads.
- Main outputs record model parameters, intervals, tests, posteriors, and run
  settings.
- Figures combine a phylogeny with linked exon-like synteny tracks.
- Color and pattern modes are available.

## Evidence still needed

See `publication_gap.md` for the empirical case set, comparator analyses,
sensitivity checks, and runtime profiling required for a first manuscript.
