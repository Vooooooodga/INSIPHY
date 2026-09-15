# Changelog

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
