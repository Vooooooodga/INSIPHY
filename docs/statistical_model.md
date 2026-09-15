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
   plus a CTMC endpoint posterior approximation for change on that branch.

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

Future foreground/background tests should compare a one-rate model against a
two-rate model in which a user-specified branch set has its own structural
change rate. This is the direct structural analogue of branch or branch-site
tests used in molecular evolution. The user-facing question would be: does a
specified evolutionary branch or clade show an elevated rate of intragenic
structural change compared with the background branches?

## Branch Event Posteriors

The current `branch_event_probabilities.tsv` values are endpoint posterior
summaries. For each branch, INSIPHY estimates the posterior probability of the
parent and child endpoint states under the fitted CTMC and reports
`ctmc_change_probability` when those endpoint states differ.

This is useful for prioritizing branches, but it is not yet stochastic character
mapping. A stochastic-map implementation should sample complete histories along
each branch conditional on observed tips and fitted rates. The expected output
should include:

- `map_sample_count`;
- `posterior_pr_any_change`;
- `posterior_expected_change_count`;
- `posterior_most_frequent_transition`;
- `posterior_transition_probability`;
- `posterior_change_count_low` and `posterior_change_count_high`;
- optional per-sample history tables for debugging.

The branch posterior should be reported per character and then summarized per
biological event class. For example, a `source_mixture` shift to `multi_source`
and an intragenic adjacency joining two source labels on the same branch should
increase confidence in one compound chimeric event.

## Publication Calibration

The current P values should be calibrated by parametric bootstrap before being
used as manuscript-level significance claims. The calibration workflow should:

1. fit the null invariant or one-rate background model;
2. simulate many structural datasets on the same species tree;
3. rerun the same INSIPHY inference;
4. compute the empirical tail probability of the observed LRT statistic;
5. report the empirical P value, Monte Carlo standard error and simulation
   settings.

This is especially important for small trees, sparse observations, missing
annotation and boundary tests where asymptotic chi-square approximations can be
optimistic.

## Current Outputs

- `character_model_scores.tsv`: parsimony score, likelihood score and fitted
  CTMC rate per character.
- `model_fit.tsv`: CTMC/Mk likelihood, AIC, BIC and observed-tip count.
- `hypothesis_tests.tsv`: invariant-model LRT, p value, null/alternative
  likelihoods and model-selection values.
- `branch_event_probabilities.tsv`: branch-level change candidates, including
  branch lengths and CTMC endpoint change probabilities when supplied in
  `species_tree.tsv`.
- `benchmark_summary.tsv` and `benchmark_detailed.tsv`: event-level and
  branch-aware benchmark summaries for simulated truth sets.
