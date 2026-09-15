# Changelog

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
- Added jingwei and Sdic curated micro-demos.
- Added literature review, publication-gap notes, likelihood-like character
  scores, baseline comparisons, simulation and benchmark commands.
- Added real-case preparation commands for annotation inspection, case manifest
  builds and hidden-segment scans.
