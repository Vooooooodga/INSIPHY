> **V18 legacy documentation.** For the V19 exon configuration model, see
> [exon_structure_model.md](exon_structure_model.md). The binary layers below are
> retained historical baselines, not the current default model.

# Statistical model and interpretation

## Observation and sampling units

A family-layer contains binary structural characters observed at tree tips.
Nucleotide length, exon drawing fragments, transcript number and number of runs
are not replicate structural characters. Presence, exon identity and splice
junctions are fitted separately. Cross-layer totals are not independent-event
counts. A character may undergo repeated change; a mutation may affect multiple
characters. No mutation-event total is estimated.

Unknown tips allow both states in parsimony and have likelihood vector `(1,1)`
in CTMC pruning. DNA-absent exon roles are inapplicable, not observed non-exonic
states. Inapplicability and unknown are distinguished in output metadata even
when both are marginalized by the separate binary model. This is not a joint
model of DNA presence and conditional RNA structure.

## Default equal-cost parsimony

For node v and state s:

```
D_v(s) = sum over children u [min_t {D_u(t) + 1(s != t)}]
```

Tip costs are zero for allowed states and infinity otherwise. An outside pass
finds all states and endpoint pairs attaining the global minimum. Branch lengths
do not affect this equal-cost reconstruction. `required` and `possible` refer
to the minimum-change criterion, not statistical significance or confidence.

`minimum_change_history.tsv` is one reproducible witness of the optimum.
`gene_change_summary.tsv` adds site minimum changes within a layer; it never adds
all possible branch placements. Unknown-only characters contribute no inferred
history. A zero tree-internal minimum does not imply that a structure never
originated before the sampled root.

## Optional CTMC

```
Q = [[-g, g], [l, -l]]
P(t) = exp(Q t)
```

ER constrains `g=l=q`; ARD estimates gain and loss separately. The foreground
alternative multiplies both rates on predeclared branches by `m`, with null
`m=1`. It tests overall turnover, not a separate gain-only effect. The foreground
must leave at least one background branch.

Rates are shared by characters within a family-layer. Root-state frequency is
estimated by default, with stationary and fixed alternatives. Unknown tips are
marginalized by log-space pruning. Pattern compression changes computation only,
not sample size or character identities. At least two explicit tips and the
selected ascertainment criterion are operational admission conditions, not
reliability guarantees.

### Dependence is checked before fitting

When a declared linked group contains multiple included characters, the current
independent-character likelihood is not fitted. Parameter estimates, AIC,
likelihood intervals, LRT P/q values, model selection and posteriors are all
unavailable with `correlated_linked_sites_not_modelled`. The matrix and parsimony
outputs remain available. Removing a compound-event label or a linked-group ID
is not a valid correction for dependence.

The check covers declared linked sites and overlapping aligned intervals; it
cannot detect every shared mutation or dependence process. Identical observed
patterns alone neither prove dependence nor justify character merging.

### Fit and test availability

Models use log-rate optimization, multiple starts, identifiability and boundary
diagnostics, and profile likelihood. A rate-model LRT requires valid null and
alternative fits, nonnegative likelihood improvement within numerical tolerance,
positive parameter-count difference and a finite two-sided profile for the tested
contrast. The usual parameter-count difference is one. Failed diagnostics yield
an unavailable result with a reason, not P=1.

Where admitted, the reference distribution is asymptotic chi-square. The output
explicitly records `small_sample_accuracy=unassessed`. Passing numeric diagnostics
does not establish finite-sample calibration. Single-gene structural data may
contain too little information to estimate rates and root frequency separately.
No formal parameter-bootstrap calibration is claimed in v18.

### Character discovery and ascertainment

`observed-at-least-one` conditions on at least one observed tip being 1:

```
log L_cond = log P(data) - log(1 - P(all observed tips are 0))
```

`variable-only` conditions on both states being observed, removing the all-zero
and all-one probabilities over that same observed-tip mask. `complete-universe`
requires an independently specified catalogue. Explicit unknown tips with missing
masks and reasons are permitted; catalogue membership does not create absence
observations. Observed all-zero characters can be included in an independent
catalogue.

These corrections condition on a fixed mask. They do not model state-dependent
annotation dropout, assembly errors, correspondence discovery or preferential
retention of slowly evolving sequence. Coverage filtering does not solve these
problems. Different discovery schemes must not be compared as though they were
the same observed dataset.

### Posteriors, model choice and counts

For valid fits, the program chooses an eligible model by AIC and reports node
states, joint branch endpoint probabilities and expected transition counts,
conditional on fitted parameters. An endpoint change probability differs from
the probability of any transition on that branch: gain followed by loss can leave
identical endpoints. Expected counts may exceed endpoint probabilities.

The posterior does not integrate uncertainty in orthology, alignment, annotation,
tree, model selection or estimated rates. Finite profile envelopes are parameter
sensitivity summaries, not full joint credible intervals. An optimizer limit is
not a confidence endpoint.

### Branch lengths and multiple tests

Supplied branch lengths must be finite and nonnegative. Rates inherit their unit.
Substitution-length trees do not provide calendar-time rates. Unit branches
produce rates per unit edge and no absolute dates. Zero-length branches have the
identity transition matrix.

BH adjustment applies to valid tests within the declared test family in a run.
Separately run genes require project-level aggregation before a project-wide FDR
claim. Multiple-testing adjustment cannot repair misspecified likelihoods or
uncalibrated P values. Foreground selection or alternative character definitions
must not be chosen by searching for the smallest P value.

See [counting](event_counting.md), [scope](scope_policy.md), and the primary
[references](references.md), especially Sankoff (1975), Lewis (2001) and Malin.
