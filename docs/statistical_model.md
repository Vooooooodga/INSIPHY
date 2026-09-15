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
  weighting, giving an interpretable approximation rather than a calibrated
  posterior.
- Each structural character also receives a likelihood-like pruning score under
  a simple discrete transition model. This score is useful for ranking and
  debugging, and it should be calibrated by simulation before manuscript-level
  accuracy claims.

## Model Comparison

Two lightweight comparisons are reported:

- strict annotation vs annotation-error model;
- independent character changes vs one compound chimeric/copy event.
- annotation-only, sequence-only and synteny-aware phylogenetic baselines.

These scores are intended as first-pass evidence summaries. A publication-grade
version should add simulation-based calibration, likelihood-style transition
models and robustness checks against annotation incompleteness.
