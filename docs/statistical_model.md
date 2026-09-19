# Statistical Model

## Analysis unit and observations

The analysis conditions on an upstream-defined single-copy ortholog set, the supplied rooted species tree, the inferred local sequence correspondences, and a frozen `structural_site_matrix.tsv`. The gene family is the reporting unit; each homologous structural site is an observation. Nucleotide count, transcript count, sequence length, and number of computational runs do not increase the number of structural sites.

The matrix separates three binary layers. `state_0` and `state_1` hold the literal labels for a layer; `state` is one of those labels or `unknown`.

| Layer | State 0 | State 1 |
|---|---|---|
| `exon_presence` | `absent` | `present` |
| `exon_role` | `not_exonic` in the selected supplied-annotation view | `exonic` in the selected supplied-annotation view |
| `splice_junction` | `absent` in the selected transcript view | `present` in the selected transcript view |

Unknown observations permit both states in parsimony and are marginalized in likelihood. Applicability is recorded separately. If the DNA unit is absent, exon role is inapplicable or unknown; it is not coded as `not_exonic`. Exon-role state 0 is conditional on adequate coverage of either the supplied transcript repertoire or the explicitly selected canonical path, as recorded in `annotation_view`. Predicted-only features remain candidate evidence and do not establish an observed state.

The method fits each layer separately. This keeps DNA presence, exon role, and junction presence as separate observations. Junctions from one split/fusion pattern can be linked with `linked_group_id`; they may share a structural cause, so the optional likelihood's conditional independence assumption should be interpreted cautiously. No joint three-layer ancestral transcript is inferred.

## Default: equal-cost maximum parsimony

For a binary site, each tip contributes the singleton set `{0}`, `{1}`, or `{0,1}` for an unknown observation. The cost is zero for the same endpoint state and one for a change. For node `v` and state `s`, the inside recursion is

```text
D_v(s) = sum over children u [ min_t (D_u(t) + c(s,t)) ]
D_tip(s) = 0 when s is allowed at the tip, otherwise infinity
minimum_changes = min_s D_root(s)
```

An outside pass determines all node states and parent-child endpoint pairs compatible with the global minimum. `required_gain` or `required_loss` means the directed change occurs on that branch in every minimum-change history. `possible_gain` or `possible_loss` means it occurs there in at least one minimum-change history. All-optimum sets are reported; the software does not select a single tied history or assign probabilities from the number of tied histories. All-unknown sites are retained and marked as having no observed states. Sites with no observed contrast are distinguished from fully observed conservation.

The branch is one parent-to-child edge of the supplied tree. A branch placement describes a state change under this cost model. A set of equally optimal placements indicates that the observed tips do not uniquely locate the change. It does not provide a P value, a time estimate, or a molecular mechanism. This is a fixed-tree Sankoff-style reconstruction ([Sankoff 1975](https://epubs.siam.org/doi/10.1137/0128004)).

## Optional continuous-time Markov models

For a binary site, the CTMC rate matrix is

```text
Q = [ -g   g ]
    [  l  -l ]
P(t) = exp(Q t)
```

Here `g` is the 0→1 rate, `l` is the 1→0 rate, and `t` is a supplied branch length. Likelihoods are calculated by pruning on the fixed tree. Unknown tips have likelihood vector `(1, 1)`, so both states are integrated over. A separate root-state frequency is estimated by default; stationary and user-fixed root frequencies are explicit alternatives.

Sites within a family-layer share rate parameters. The likelihood is the sum of site log likelihoods, with site-pattern compression preserving site identities and observation masks. A site enters CTMC fitting only when at least two tree tips have explicit binary states and the observed pattern satisfies the selected ascertainment rule. Sites failing either condition remain in the frozen matrix and are counted in the analysis scope, with an exclusion reason. This is a model for binary structural characters. It does not use DNA alignment length as a substitute for the number of structural observations.

### ER versus ARD

The nested comparison is:

| Model | Null or alternative | Rate constraints |
|---|---|---|
| ER | H0 | `g = l = q` |
| ARD | H1 | `g` and `l` estimated separately |

With the same root-frequency treatment in both models, the likelihood-ratio statistic is `2 * (logL_ARD - logL_ER)`. Under regular conditions its reference distribution is asymptotic chi-square with degrees of freedom equal to the difference in free parameters (one for this comparison). It tests whether the data support separate gain and loss rate parameters under these models. It does not test whether a particular gain or loss occurred on a named branch.

### Foreground rate model

The foreground comparison uses branches specified before examining the result:

| Model | Null or alternative | Branch rates |
|---|---|---|
| H0 | one ARD process | `g` and `l` shared across all branches |
| H1 | foreground multiplier | both rates on selected branches are multiplied by `m` |

The null is `m = 1`. Under regular conditions the comparison uses the parameter-count difference, ordinarily one degree of freedom. It tests an overall rate multiplier on the selected branches under this parameterization. Foreground branches must be declared before inspecting the result, must match the supplied tree, and must leave at least one background branch. A significant comparison does not identify a specific event or its cause.

### Fit and test availability

Rates are optimized on the log scale with multiple starting points. The fit and likelihood-ratio test are reported separately. A test is available only when the required models converge, the alternative likelihood is not below the null beyond numerical tolerance, numerical identifiability and boundary diagnostics permit the asymptotic comparison, and the tested contrast has a finite two-sided profile interval containing its estimate. If a parameter is not identifiable, an estimate is boundary-limited, a profile endpoint is open or non-finite, or optimization fails, the affected test is reported unavailable with its reason. The degrees of freedom are based on the actual free-parameter difference.

The current likelihood multiplies site likelihoods conditionally independently. When one `linked_group_id` contains more than one site, those sites can be components of one structural change. The family-layer LRT is then unavailable with `lrt_unavailable_reason=correlated_linked_sites_not_modelled`; parsimony and site-level descriptive outputs remain available.

`total_structural_sites` counts all sites in that family-layer in the frozen matrix. `n_structural_sites` and `included_structural_sites` count only sites admitted to the CTMC after the explicit-tip and ascertainment filters. `n_informative_patterns` is a further subset with an observed binary contrast. These quantities must not be interchanged when reporting sample size.

The chi-square reference is asymptotic. Its accuracy for the small number of structural sites commonly present in a single gene has not been assessed by the method specification; the output records `small_sample_accuracy=unassessed`. A P value is evidence about the stated rate-model comparison under its assumptions; it is not evidence that a particular exon split, fusion, gain, or loss happened. Multiple-testing adjustment does not correct a poor model fit or an inaccurate small-sample reference distribution.

## Ancestral and branch probabilities

For a valid fitted model, inside-outside likelihood messages give node-state probabilities and joint parent/child endpoint probabilities conditional on the tree, observed states, model, and fitted parameter estimates. Expected numbers of transitions along a branch integrate over possible CTMC paths and can exceed the probability that the two endpoints differ. These quantities are distinct.

Parameter-conditional probabilities do not integrate over uncertainty in the alignment, annotation, correspondence, tree, or fitted rates. A probability range is reported only when the parameter profile supports a finite, interpretable range. Non-identifiable or boundary-limited fits and open profile intervals do not receive a seemingly precise posterior sensitivity range. An unavailable range is labelled unavailable; an optimizer's search limit is not treated as a confidence endpoint. Likelihood-test availability and conditional posterior availability are reported separately.

## Ascertainment and the character universe

The likelihood's inclusion rule must describe the same character universe that was supplied to the model. The default `observed-at-least-one` mode is for discovered structural sites selected because at least one observed tip has state 1. For each site, the likelihood is conditioned on that same event using the same observed-tip mask:

```text
log P(data | at least one 1) = log P(data) - log[1 - P(all observed tips are 0)]
```

`variable-only` is for a matrix intentionally restricted to sites with both state 0 and state 1 among observed tips. Its Lewis-style conditional likelihood is

```text
log P(data | variable) = log P(data) - log[1 - P(all 0) - P(all 1)]
```

The denominator is computed over the observed tips for that site. `complete-universe` is for a separately defined candidate catalogue that includes sites with all-zero observations; membership in the catalogue does not turn unobserved species into state 0. Excluded sites and reasons are retained. A site set filtered under one rule cannot be analyzed as if it had been selected under another rule.

These corrections condition on a fixed observation mask. They do not model annotation or assembly failure that depends on the true state. Linked or neighboring sites can also be dependent; the current binary-site likelihood does not jointly model a deletion or one structural event affecting several characters. Separate layers or linked sites must not be combined into a single independent-event count.

## Multiple testing and branch lengths

Benjamini-Hochberg adjustment is applied within each declared test family across valid tests in the run. Missing P values remain missing, and the number of tests included is reported.

In `supplied` mode, every non-root branch length must be finite and nonnegative. A zero-length branch has the identity transition matrix. In `unit` mode, every branch has length one; rates then describe transitions per unit branch and carry no calendar-time interpretation. Rates inherit the unit of the supplied tree lengths.

## Statistical scope and references

The analysis conditions on the upstream single-copy orthology, fixed species tree, local correspondence, annotation-derived states, and ascertainment rule. It does not estimate a gene tree or integrate uncertainty across these inputs. A small number of sites may leave rate parameters or ancestral states weakly identified. Numerical convergence alone does not establish reliable finite-sample P values.

The likelihood menu follows standard phylogenetic practice: specify nested hypotheses and shared parameters, estimate them by maximum likelihood on a fixed tree, and use a likelihood-ratio reference only under its regularity assumptions. See the [PAML manual](https://github.com/abacus-gene/paml/wiki/BASEML) and [HyPhy methods](https://www.hyphy.org/methods/) for examples of problem-specific model menus and their assumptions. INSIPHY's characters are homologous gene-internal structural sites; interpretation remains with the user.
