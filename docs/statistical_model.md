# Statistical Model

## 1. Analysis unit

The formal v0.13 analysis accepts a single-copy ortholog set and a fixed rooted
species tree. Gene-level homology is treated as known input. For each gene
family and structural layer, the optional likelihood analysis treats homologous
sites as conditionally independent observations sharing evolutionary parameters.
The default analysis localizes events using equal-cost parsimony.

The observed tip state (x_{is}) for site (i) and species (s) is one of
`0`, `1`, or `unknown`:

| layer | state 0 | state 1 | biological observation |
|---|---|---|---|
| `exon_presence` | absent | present | homologous exon sequence |
| `exon_role` | not_exonic | exonic | role of an observed homologous sequence |
| `splice_junction` | absent | present | splice junction between homologous units |

An absent state requires explicit sequence evidence. Missing coverage,
ambiguous correspondence, and annotation uncertainty are coded as unknown.

## Default qualitative reconstruction

For each structural site, minimize the sum of state-change costs over the
supplied rooted tree. Staying in the same state costs zero; either direction
of change costs one. Unknown tips allow both states. Inside and outside
minimum-cost messages identify all node states and branch endpoint pairs
compatible with the global minimum. No arbitrary tie-breaking history is
chosen. Pattern compression retains site identifiers and their missing masks.

`required` means that every minimum-cost history has the same directed change
on that branch; `possible` means that at least one minimum-cost history has
that change. These labels express consequences of the parsimony criterion,
not probabilities. Missing observations remain missing in the input table.
Sites with no observed contrast are identified separately from fully observed
conservation. Time and evolutionary rates are not estimated in this mode.

One within-exon splice boundary encodes one split/fusion comparison. The
presence, role and junction layers are not combined into a likelihood or an
independent count of biological events. Remaining within-layer dependence is
an assumption of the optional probability model.

## 2. Optional continuous-time Markov model

Each binary structural site evolves along the supplied tree under

```text
Q = [ -q01   q01 ]
    [  q10  -q10 ]
```

where (q_{01}) is the gain rate and (q_{10}) is the loss rate per unit
branch length. Transition probabilities on a branch of length (t) are
(P(t)=exp(Qt)). The default model estimates the root presence probability
independently of the gain and loss rates. `stationary` and user-supplied
`fixed` root frequencies remain explicit alternatives.

Likelihoods are computed with scaled Felsenstein pruning. At an unknown tip,
the conditional likelihood vector is ((1,1)), which integrates over both
states. Scaling prevents numerical underflow on larger trees.

## 3. Shared-parameter likelihood

For a gene family (g) and layer (l), one parameter set is fitted to all
usable homologous sites:

```text
ln L(theta | X_g,l, T) = sum_i ln P(x_i | theta, T)
```

A site enters fitting when at least two terminal species have observed states.
Sites with fewer than two observations stay in
`structural_site_matrix.tsv` but do not contribute to parameter fitting.

This pooling follows the same statistical principle used by sequence
likelihood methods: model parameters are estimated from a collection of sites,
then site- and branch-specific histories are conditioned on those shared
parameters.

## 4. ER versus ARD

The optional `--model er-ard` nested comparison is:

- H0, `ER`: (q_{01}=q_{10}=q), one free rate parameter;
- H1, `ARD`: (q_{01}) and (q_{10}), two free rate parameters.

With the default estimated root frequency, each model has one additional
parameter: two parameters in ER and three in ARD. The nested comparison still
differs by one parameter. AIC and the reported parameter counts include the
estimated root frequency.

Rates are optimized on the log scale with bounded L-BFGS-B and several starting
points. `model_fits.tsv` reports maximum log likelihood, AIC, convergence,
boundary status, and 95% profile-likelihood intervals.

The likelihood-ratio statistic is

```text
LR = 2 * (ln L_ARD - ln L_ER).
```

The asymptotic `chi-square(df=1)` P value is reported only when both models
converge, the alternative likelihood is at least the null likelihood, the
observed-information matrix has full numerical rank, and no estimate lies on
the optimization boundary.

With no observed tip variation, gain and loss directions cannot be identified.
Boundary estimates also violate the regular chi-square approximation. These
cases receive `test_status=parameters_not_estimable` and `p_value=NA`.

## 5. Foreground model

The foreground analysis compares:

- H0: one ARD process on all branches;
- H1: the same gain/loss rates multiplied by (m) on user-specified branches.

The null value is (m=1), and the LRT uses one degree of freedom when regular
conditions hold. A fitted multiplier above one means the selected branches
have a higher model-based transition rate. Biological causes remain outside
the statistical test.

## 6. Ascertainment

Automatically generated candidate sites exist because state 1 was observed in
at least one species. The default correction therefore uses:

```text
ln P(x_i | at least one state 1) =
ln P(x_i) - ln(1 - P(all 0)).
```

`--ascertainment complete-universe` requires an externally defined
`structural_site_universe.tsv` containing meaningful candidate sites. Unobserved
species are not assigned zero by membership in that catalogue.

`--ascertainment variable-only` applies a Lewis-style Mkv correction:

```text
ln P(x_i | variable) =
ln P(x_i) - ln(1 - P(all 0) - P(all 1)).
```

Every included site must vary among its observed terminal states. For every
mode, the inclusion rule and its likelihood normalizer use exactly the same
set of observed tips. This treats missingness as a fixed observation mask;
it does not model state-dependent annotation or assembly errors.

## 7. Node and branch posteriors

After model selection by AIC, inside-outside messages give marginal
empirical-Bayes probabilities for every node and joint endpoint probabilities
for every branch, conditional on maximum-likelihood parameter estimates:

```text
P(X_parent=a, X_child=b | tip states, fitted model).
```

`node_state_posteriors.tsv` contains node marginals.
`branch_transition_posteriors.tsv` contains:

- posterior probability of each directed endpoint change;
- total endpoint change probability;
- conditional expected numbers of gains and losses along the branch.

The expected count integrates over all CTMC paths conditional on branch
endpoints. It can exceed the endpoint-change probability because an even
number of hidden transitions may return to the starting state.

`structural_changes.tsv` gives a compact branch/site view and retains both
directions. These probabilities are conditional on the fitted parameter set.
Joint parameter sensitivity is currently not estimated; range columns are
`NA` and labelled `sensitivity_not_estimated`. Separate one-dimensional
profile endpoints are not presented as a joint confidence region. Profiles
that reach the numerical search bounds are labelled range-limited, and
optimization failures remain distinguishable from wide statistical intervals.

## 8. Multiple testing

Benjamini-Hochberg adjustment is applied separately within each `test_id`
across all valid family-layer tests in one run. With one valid test, its BH
q value equals its P value.

## 9. Branch lengths

`--branch-length-mode supplied` requires a positive length for every
non-root branch. Estimated rates then use the same unit as the tree, such as
substitutions per site or time.

`--branch-length-mode unit` sets every branch length to one. Rates then mean
expected structural transitions per branch. Unit branches preserve topology
but discard elapsed-time information.

## 10. Interpretation limits

The model conditions on:

- the upstream single-copy ortholog set;
- the supplied species-tree topology and branch lengths;
- the inferred exon correspondence;
- the structural state coding.

Uncertainty in gene orthology and tree topology is not integrated in v0.13.
Site independence is an approximation because neighboring exon and junction
states can be biologically coupled. A small number of exons gives wide
likelihood intervals and limited LRT power. These limitations should be
reported directly in gene-level analyses.

## 11. Relation to established phylogenetic statistics

The default minimum-change reconstruction uses a Sankoff-style cost recursion
on a fixed tree, with outside costs retaining all globally optimal branch
endpoint pairs. The classical fixed-tree minimum-cost formulation is described
by [Sankoff (1975)](https://epubs.siam.org/doi/10.1137/0128004). INSIPHY applies
this framework to observed gene-structure states and does not search for a tree.

The implementation follows the likelihood logic emphasized in the PAML manual:
define explicit null and alternative models, estimate shared parameters by
maximum likelihood on a fixed tree, compare nested models with an LRT when
regularity conditions hold, and treat ancestral reconstructions as
model-conditional probabilities.

The binary state process is related to Pagel-style discrete-trait CTMC models
and Mk/Mkv models. The biological observations here are homologous
gene-structure sites rather than organismal phenotypes.

## 12. Primary output files

- `structural_site_matrix.tsv`
- `model_fits.tsv`
- `model_tests.tsv`
- `node_state_posteriors.tsv`
- `branch_transition_posteriors.tsv`
- `structural_changes.tsv`
- `run_parameters.json`
- `excluded_families.tsv`
