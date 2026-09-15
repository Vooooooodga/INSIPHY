# Statistical Model Notes

INSIPHY currently separates biological evidence from phylogenetic interpretation.

## Current Statistical Layers

1. **Evidence scoring for HSGs**: sequence identity, coverage, splice motif,
   intron phase, strand, boundary class and local order are combined into
   homologous segment group scores. These are evidence scores for
   correspondence, not phylogenetic p values. Real-case clustering also
   requires role-compatible sequence support and constrains each HSG component
   to contain at most one segment from the same gene copy.
2. **Weighted Sankoff reconstruction**: segment presence, segment role,
   intragenic adjacency, source mixture and copy multiplicity are reconstructed
   on a fixed species tree. Gains, losses and role shifts can have different
   costs, so the resulting branch calls are interpretable event candidates.
3. **Branch-length-aware CTMC/Mk likelihood**: each structural character is
   treated as a discrete trait evolving on the fixed tree. INSIPHY fits a
   one-rate continuous-time Markov model by grid search and reports the
   maximum log likelihood, fitted rate, AIC and BIC.
4. **Invariant-model LRT**: for each character, INSIPHY compares the fitted
   CTMC/Mk model against a no-change model. The test reports
   `lrt_statistic`, `df`, `p_value`, `p_value_method`, null/alternative
   likelihoods and AIC/BIC values.
5. **Branch support**: branch event tables include the parsimony change status
   plus a CTMC posterior approximation for change on that branch.

Annotation dropout is represented as evidence uncertainty and hidden-segment
support. A missing annotation alone is not treated as biological segment loss.

## Model Comparison

Two lightweight comparisons are reported:

- strict annotation vs annotation-error model;
- independent character changes vs one compound chimeric/copy event.
- annotation-only, sequence-only and synteny-aware phylogenetic baselines.

These scores are first-pass evidence summaries. The formal phylogenetic
statistics are in `model_fit.tsv`, `character_model_scores.tsv` and
`hypothesis_tests.tsv`.

## Hypotheses

The default hypothesis test is:

- **H0**: the intragenic character is invariant on the supplied species tree,
  allowing a tiny tip observation error for annotation uncertainty.
- **H1**: the character evolves under a one-rate CTMC/Mk model with transition
  probabilities scaled by branch length.

The LRT uses a boundary-rate chi-square mixture approximation,
`0.5 * chi-square(df=1)`, because the invariant model fixes the transition rate
at zero. With only a few closely related species, the p value is best read as a
calibrated-looking evidence statistic that still needs simulation-based
calibration before publication-level claims.

If fewer than two terminal taxa have observed non-unknown states for a
character, INSIPHY reports `insufficient_observed_tips` and sets the p value to
`NA`.

Future foreground/background tests can compare a one-rate model against a
two-rate model in which a user-specified branch set has its own structural
change rate. This is the direct structural analogue of branch or branch-site
tests used in molecular evolution.

## Current Outputs

- `character_model_scores.tsv`: parsimony score, likelihood score and fitted
  CTMC rate per character.
- `model_fit.tsv`: CTMC/Mk likelihood, AIC, BIC and observed-tip count.
- `hypothesis_tests.tsv`: invariant-model LRT, p value, null/alternative
  likelihoods and model-selection values.
- `branch_event_probabilities.tsv`: branch-level change candidates, including
  branch lengths and CTMC change probabilities when supplied in
  `species_tree.tsv`.
- `benchmark_summary.tsv` and `benchmark_detailed.tsv`: event-level and
  branch-aware benchmark summaries for simulated truth sets.
