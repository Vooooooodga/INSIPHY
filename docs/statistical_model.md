# Statistical Model Notes

INSIPHY currently separates biological evidence from phylogenetic interpretation.

## Current First-Pass Model

- Segment presence and adjacency use asymmetric transition costs: gains are
  penalized more than losses.
- Role state changes distinguish absent, CDS, exon/UTR and intron/noncoding
  states.
- Source mixture detects whether one copy combines segments assigned to
  different source labels.
- Copy multiplicity summarizes whether a species has one or multiple target
  copies.
- Node probabilities are converted from Sankoff scores by local exponential
  weighting.
- Each structural character also receives a branch-length-aware CTMC/Mk
  likelihood fit. The current implementation performs a small grid search over
  one rate parameter and reports log likelihood, AIC and BIC.
- Annotation dropout is represented as evidence uncertainty and hidden-segment
  support. It does not by itself define segment absence.

## Model Comparison

Two lightweight comparisons are reported:

- strict annotation vs annotation-error model;
- independent character changes vs one compound chimeric/copy event.
- annotation-only, sequence-only and synteny-aware phylogenetic baselines.

These scores are intended as first-pass evidence summaries. A publication-grade
version should add simulation-based calibration, likelihood-style transition
models and robustness checks against annotation incompleteness.

## Current Outputs

- `character_model_scores.tsv`: parsimony score, likelihood score and fitted
  CTMC rate per character.
- `model_fit.tsv`: CTMC/Mk likelihood, AIC, BIC and observed-tip count.
- `branch_event_probabilities.tsv`: branch-level change candidates, including
  branch lengths when supplied in `species_tree.tsv`.
- `benchmark_summary.tsv` and `benchmark_detailed.tsv`: event-level and
  branch-aware benchmark summaries for simulated truth sets.
