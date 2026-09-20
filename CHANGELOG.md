# 0.18.0 — structural characters and elementary changes (2026-09-20)

- Remove the previous import namespace and executable alias; require Python 3.10+.
- Add character catalogues, coordinate evidence, dependence metadata and applicability checks.
- Remove formal compound-event summaries; report per-layer minimum character changes and a compatible optimum.
- Block all independent-character likelihood inference when known linked included sites are present.
- Permit explicit unknown tips in independently specified complete catalogues.
- Correct complementary-fragment repeat classification using actual reference coverage.
- Restore conservative CI assertions and add raw FASTA/GFF3 integration regressions.
- Add CLI preflight, explicit output ownership, preserved backups, failure logs and version provenance.
- Split large implementation modules and provide English documentation and a rendered vector methods figure.

No biological-accuracy benchmark or finite-sample statistical calibration is claimed.
Earlier release notes and unchanged historical validation remain in repository history.


### V18 file-input and explanatory-figure completion

- Accept genomic FASTA/GFF files or directories and a supplied rooted tree without
  a hand-written table; exact ortholog-FASTA identifiers select target loci.
- Support combined genomic locus FASTA with distinct species sequence IDs.
- Export paired locus FASTA/GFF with automatic flanks, strand/phase preservation
  and reversible source-coordinate mapping.
- Add optional, explicit AGAT normalization with command/failure provenance.
- Add original tree/structure/ribbon teaching plates and `intraphy explain`;
  keep synthetic guide and data-derived target figures clearly separate.
- Test file contracts, coordinate round-trips, conservative rejection and the
  illustrated parsimony examples. No new biological/statistical calibration.

### Biology-facing method figures and explicit model explanation

- Replace the previous five teaching plates with three focused figures: overview,
  homology inference, and gene structures to a phylogenetic probability model.
- Show the species tree and gene structures together, including sequence loss,
  intron loss, missing annotation and uncertainty rather than only exon splitting.
- Explain compatible-match chain selection and distinguish complementary coverage
  from competing copies. Use named biological characters instead of internal keys.
- Calculate illustrative ancestral-state probabilities and test them against the
  production pruning and posterior algorithms; fixed teaching parameters are explicit.
- Document likelihood construction, missing-state treatment, rate-sharing scope,
  character dependence and conditional inference in `docs/model_bridge.md`.
- Record which primary methods, figure captions and image assets were examined.
- Local full-suite validation: 446 tests passed with real MAFFT and minimap2; this
  is implementation validation, not a biological benchmark or statistical calibration.
