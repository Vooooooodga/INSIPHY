# Changelog

## 0.12.3 - 2026-09-16

- Apply the same length-aware alignment policy when calculating within-group
  mean exon identity, so medium exon pairs are not forced through the bounded
  internal dynamic-programming backend.

## 0.12.2 - 2026-09-16

- Read the first optimal Biopython alignment by iteration, avoiding integer
  overflow when low-complexity exons admit more optimal paths than can be
  counted in a 64-bit integer.

## 0.12.1 - 2026-09-16

- Preserve literal pipe characters in GFF/GTF feature identifiers and parent
  links, including GenBank WGS transcript identifiers.
- Add an `auto` correspondence aligner that uses exact global alignment for
  short exons, MAFFT for medium sequence pairs, and minimap2 for long pairs.
- Make `auto` the default for case construction and record the backend selected
  for every pairwise comparison.

## 0.12.0 - 2026-09-16

- Made complete exon intervals the primary structural observations, with CDS,
  UTR, GFF phase and reading-frame information retained as attributes.
- Connected sequence-supported annotation completion to the standard run and
  distinguished supported absence from unresolved missing evidence.
- Projected within-exon splice boundaries onto homologous sequence coordinates
  so distinct split/fusion sites remain separate.
- Replaced unconstrained transitive homology merging with tree-ordered,
  coordinate-constrained profile merging and bounded pair scoring.
- Corrected aligned-pair coverage, strand normalization, CDS phase continuity,
  protein queries for miniprot and long-alignment failure behavior.
- Added observed-at-least-one ascertainment, independently estimated root
  frequency, numerical identifiability checks and strict foreground validation.
- Added open-ended profile interval states and profile-likelihood sensitivity
  ranges for empirical-Bayes node and branch probabilities.
- Removed automatic directional event labels and updated tree-aligned,
  colorblind-aware synteny figures.

## 0.11.0 - 2026-09-16

- Made single-copy ortholog analysis the formal default and isolated retained
  multi-copy routines behind `experimental-multicopy`.
- Added an OrthoFinder single-copy orthogroup importer without gene-homology
  re-inference.
- Defined three biological structural-site layers: exon sequence presence,
  exonic role and splice-junction presence.
- Added shared-parameter binary CTMC likelihoods with ER/ARD and homogeneous
  versus foreground model comparisons.
- Added continuous optimization, profile-likelihood intervals, regularity
  checks for LRT P values and BH correction across valid tests.
- Added marginal node-state posteriors, joint branch endpoint posteriors and
  expected directional transition counts.
- Added `structural_site_matrix.tsv`, `model_fits.tsv`, `model_tests.tsv`,
  `node_state_posteriors.tsv`, `branch_transition_posteriors.tsv`,
  `structural_changes.tsv`, `excluded_families.tsv` and `run_parameters.json`.
- Updated phylogenetic figures to use continuous posterior probabilities and
  retained colorblind-aware color and pattern encodings.
- Replaced deprecated Biopython `pairwise2` calls with `PairwiseAligner`.
- Reframed documentation and the real RpL32 demo around the formal
  single-copy model.

## 0.10.0 - 2026-09-16

- Moved possible biological readings out of core event/support tables into
  `interpretation_hints.tsv`; core statistics now report observable structural
  changes and support only.
- Added `progressive_element_correspondence.tsv`,
  `ancestral_element_graph.tsv` and `ancestral_intragenic_paths.tsv`.
- Renamed newly generated internal correspondence components to `HC_*` and
  public internal tables to `internal_homology_*`.
- Added `build-case --copy-tree/--gene-tree` so upstream gene/copy trees can be
  carried into real-case analyses.
- Added real Drosophila demo records for RpL32, jingwei and Sdic, including
  copy-tree positive cases and a conserved-control negative case.
- Tightened internal aligner candidate filtering so intron-intron pairs remain
  context evidence and large internal alignments use a fast approximation.
- Added `insiphy calibrate` for simulation-based operating characteristics:
  false positive rate, power, precision/recall, branch placement accuracy and
  bootstrap behavior.
- Updated documentation to frame simulation calibration as method evaluation,
  separate from real-data event calls.

## 0.9.1 - 2026-09-16

- Changed the default synteny correspondence encoding to color, with
  `--correspondence-encoding pattern` retained for color-independent figures.
- Added `copy_tree.tsv` / `gene_tree.tsv` support for EG presence, EG role,
  EG adjacency and source-mixture phylogenetic inference in multi-copy gene
  sets.
- Kept `copy_multiplicity` on the species tree and added
  `phylogeny_scope.tsv` so each result records the tree scope used by each
  statistical layer.
- Updated Sankoff branch placement to use a global root-to-tip backtrace
  instead of independent node-local state sets, improving copy-lineage event
  placement.
- Updated simulation fixtures to emit copy trees and branch truth for
  copy-lineage structural events.

## 0.9.0 - 2026-09-16

- Promoted exon-like groups (EGs) from display labels to the main public
  correspondence and phylogenetic-statistical objects.
- Added `element_correspondence.tsv` and
  `element_phylogenetic_coverage.tsv`; retained internal homology component
  tables as evidence-graph records.
- Switched primary structural layers to `element_presence`,
  `element_role_state` and `element_adjacency_state`.
- Updated event object ids, synteny graph edges and visual links to use EG ids.
- Expanded simulation scenarios to cover exonization, TE-associated
  exonization, splice-boundary shift, split/fusion, source joining, tandem
  duplication, dispersed processed-copy context, annotation dropout, negative
  control and ambiguous paralogous similarity.
- Added documented real+simulation demo commands while keeping the package
  CLI-first.

## 0.8.0 - 2026-09-16

- Reworked SVG synteny figures around exon-like biological units instead of
  exposing raw internal graph IDs as primary visual objects.
- Added link/ribbon-style correspondence between exon-like blocks across
  species/copy tracks.
- Rendered introns and other non-exonic intervals as gray context spans, with
  dashed candidate-source boxes only when connected to exon-like evidence.
- Clarified that internal homology components are correspondence-evidence
  clusters, while biological modeling and figures should be read through
  exon-like elements, splice boundaries, adjacencies and event tables.

## 0.7.0 - 2026-09-16

- Split event reporting into observable structural-pattern and call-scope
  fields.
- Reclassified high-identity paralogous segment matches as ambiguous evidence
  instead of core gene-conversion event calls.
- Added internal homology phylogenetic coverage to summarize tree-spanning,
  partial and tip-specific component support.
- Updated benchmark summaries so core structural precision/recall exclude
  copy-context and ambiguous evidence calls.
- Added pattern-only correspondence visualization by default and an integrated
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
- Propagated source labels into extracted `segment_occurrences.tsv` and
  `segment_homology.tsv` rows.
- Added source-label inference for derived-copy segments from strongest
  source-copy internal homology matches, enabling chimeric source mixture calls
  for cases such as jingwei and Sdic.
- Documented that INSIPHY starts from a supplied homologous gene/copy set and
  does not perform whole-genome orthogroup inference.

## 0.3.0 - 2026-09-15

- Added `hypothesis_tests.tsv` with invariant-model likelihood-ratio tests,
  p values, fitted CTMC/Mk rates and null/alternative AIC/BIC values.
- Added CTMC branch-change probability columns to
  `branch_event_probabilities.tsv` and `candidate_structural_events.tsv`.
- Tightened internal homology clustering with role-compatible sequence support and a
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
- Added graph-based internal homology correspondence, reciprocal-best calls,
  membership scores and internal graph-edge output.
- Added copy relationship calls for tandem, same-contig and dispersed/retrocopy
  candidates.
- Added branch-length-aware CTMC/Mk model fitting with rate, AIC and BIC
  outputs in `model_fit.tsv`.
- Expanded simulation and benchmark support for annotation-dropout negative
  controls and detailed event metrics.

## 0.1.0 - 2026-09-15

- Added the INSIPHY package scaffold and CLI.
- Added genome FASTA/GFF extraction for single-gene copies.
- Added first-pass internal homology, correspondence, adjacency and
  copy-context builders.
- Added fixed-tree structural inference for intragenic segment evolution.
- Added jingwei and Sdic curated internal fixtures.
- Added literature review, publication-gap notes, likelihood-like character
  scores, baseline comparisons, simulation and benchmark commands.
- Added real-case preparation commands for annotation inspection, case manifest
  builds and hidden-segment scans.
