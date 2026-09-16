# Changelog

## 0.8.0 - 2026-09-16

- Reworked SVG synteny figures around exon-like biological units instead of
  exposing raw HSG graph IDs as primary visual objects.
- Added link/ribbon-style correspondence between exon-like blocks across
  species/copy tracks.
- Rendered introns and other non-exonic intervals as gray context spans, with
  dashed candidate-source boxes only when connected to exon-like evidence.
- Clarified that HSGs are internal correspondence-evidence clusters, while
  biological modeling and figures should be read through exon-like elements,
  splice boundaries, adjacencies and event tables.

## 0.7.0 - 2026-09-16

- Split event reporting into observable `structural_pattern`,
  `mechanism_hypothesis` and `call_scope` fields.
- Reclassified high-identity paralogous segment matches as ambiguous evidence
  instead of core gene-conversion event calls.
- Added `hsg_phylogenetic_coverage.tsv` to summarize tree-spanning, partial
  and tip-specific HSG support.
- Updated benchmark summaries so core structural precision/recall exclude
  copy-context and ambiguous evidence calls.
- Added pattern-only HSG visualization by default and an integrated
  species-tree plus gene-internal synteny SVG.

## 0.6.0 - 2026-09-16

- Clarified the public method boundary: INSIPHY starts from supplied homologous
  gene/copy sets and analyzes gene-internal structure on a species tree.
- Added tree-distance-aware `progressive_correspondence.tsv` for interpreting
  segment correspondence inside the supplied gene set.
- Added q values, `event_support_summary.tsv` and gene-conversion candidate
  calls to phylogenetic event reporting.
- Added colorblind-friendly SVG visualizations with texture, shape, line style
  and labels as primary encodings.
- Updated accession-level Drosophila real-case manifests for `jingwei`, `Sdic`
  and the RpL32 conserved control.
- Removed the experimental genome-alignment import layer from the package API.

## 0.5.0 - 2026-09-15

- Added optional alignment backend selection for segment correspondence and
  hidden-segment scans: `internal`, `minimap2` and `miniprot`.
- Added candidate prefiltering, role-shift-aware intron/exon compatibility,
  and threaded pair scoring to reduce unnecessary pairwise alignment work.
- Added CTMC parametric bootstrap calibration in `hypothesis_bootstrap.tsv`.
- Added stochastic character mapping summaries in
  `branch_history_posteriors.tsv`.
- Added optional foreground/background structural-rate tests in
  `foreground_tests.tsv`.
- Added benchmark calibration summaries and additional named simulation
  scenarios for exonization, source joining, tandem duplication and
  split/fusion checks.

## 0.4.0 - 2026-09-15

- Added manifest-level `role_hint`, `source_label` and `copy_role` support for
  source/background/derived copy sets supplied by upstream homology workflows.
- Propagated source labels into extracted `segment_occurrences.tsv` and HSG
  `segment_homology.tsv` rows.
- Added source-label inference for derived-copy segments from strongest
  source-copy HSG matches, enabling chimeric source mixture calls for cases
  such as jingwei and Sdic.
- Documented that INSIPHY starts from a supplied homologous gene/copy set and
  does not perform whole-genome orthogroup inference.

## 0.3.0 - 2026-09-15

- Added `hypothesis_tests.tsv` with invariant-model likelihood-ratio tests,
  p values, fitted CTMC/Mk rates and null/alternative AIC/BIC values.
- Added CTMC branch-change probability columns to
  `branch_event_probabilities.tsv` and `candidate_structural_events.tsv`.
- Tightened HSG clustering with role-compatible sequence support and a
  one-segment-per-gene-copy component constraint to reduce transitive
  over-merging on real locus inputs.
- Added `insufficient_observed_tips` reporting for phylogenetic tests with
  fewer than two informative terminal observations.
- Cached CTMC transition matrices for faster real-case runs.
- Documented the statistical model mapping from intragenic structural
  characters to parsimony, CTMC/Mk likelihood and LRT outputs.

## 0.2.0 - 2026-09-15

- Added transcript-aware GFF/GTF extraction with `transcript_paths.tsv` and
  `intron_sites.tsv`.
- Added gapped local alignment helpers, splice motif scoring and frame-status
  reporting for hidden-segment scans.
- Added graph-based HSG correspondence, reciprocal-best calls, membership
  scores and `hsg_graph_edges.tsv`.
- Added copy relationship calls for tandem, same-contig and dispersed/retrocopy
  candidates.
- Added branch-length-aware CTMC/Mk model fitting with rate, AIC and BIC
  outputs in `model_fit.tsv`.
- Expanded simulation and benchmark support for annotation-dropout negative
  controls and detailed event metrics.

## 0.1.0 - 2026-09-15

- Added the INSIPHY package scaffold and CLI.
- Added genome FASTA/GFF extraction for single-gene copies.
- Added first-pass HSG, correspondence, adjacency and copy-context builders.
- Added fixed-tree structural inference for intragenic segment evolution.
- Added jingwei and Sdic curated internal fixtures.
- Added literature review, publication-gap notes, likelihood-like character
  scores, baseline comparisons, simulation and benchmark commands.
- Added real-case preparation commands for annotation inspection, case manifest
  builds and hidden-segment scans.
