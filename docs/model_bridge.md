> **V18 legacy documentation.** For the V19 exon configuration model, see
> [exon_structure_model.md](exon_structure_model.md). The binary layers below are
> retained historical baselines, not the current default model.

# From homologous gene regions to a phylogenetic model

This page specifies what the drawings mean mathematically. It is not a new
joint model of gene sequences, exon homology, annotation error and evolution.
The sequence correspondence and species phylogeny are inferred upstream or
provided as input; the current likelihood is conditional on them.

## 1. What is compared?

A homologous interval and an exon are not interchangeable. An exon can contain
several aligned subintervals; homologous DNA can be exonic in one species and
non-exonic in another. A splice position can be comparable even when the internal
sequence of the intron is not alignable. Matching flanking exons does not establish
homology of every intervening nucleotide.

Protein alignments are projected through coding coordinates. Nucleotide matches
and their flanking genomic context support additional correspondence. Candidate
matches must satisfy the implemented order, strand, coverage and ambiguity rules.
These rules supply evidence, not calibrated probabilities of homology.

For an exon missing from a GFF, search the available genomic interval. A supported
DNA match can establish presence without establishing exonic status. A transcript
that explicitly spans the region can support exonic or non-exonic status;
otherwise status remains unknown. Genomic sequence alone cannot distinguish an
unannotated expressed exon from a non-exonic homologous interval. Candidates are
not automatically promoted to confirmed exons. Supported deletion requires
additional positional and alignment evidence; failure to find a match is not
sufficient.

## 2. Each character asks one biological question

| Biological character | Internal field | 1 | 0 | Unresolved / inapplicable |
|---|---|---|---|---|
| Presence of homologous DNA | `exon_presence` | DNA present | Supported DNA absence | Unresolved correspondence or inadequate sequence |
| Exonic status of that DNA | `exon_role` | Exonic in the chosen annotation view | Non-exonic in informative annotation | Unknown if annotation is uninformative; inapplicable if DNA is absent |
| Intron at a homologous position | `splice_junction` | Corresponding junction present | Informative continuous transcript span | Unresolved position/path; inapplicable if the required DNA is absent |

The same column always refers to the same aligned region or position, not to
"exon 2" by rank. Different columns may concern different positions in one gene.
Drawing one interval as several boxes does not create extra characters.

The default annotation view is the set of annotated isoforms. Presence in at least
one such isoform does not measure usage frequency. Missing isoforms are not
biologically disproved. The alternative canonical-transcript view requires an
explicit choice and a separate character matrix.

## 3. A binary continuous-time Markov model

Let `X_k(v)` be the state of character `k` at node `v` of the rooted species
phylogeny `T`. The observations occur at tips. Internal states are unobserved.

For one character class in one gene family, the instantaneous rate matrix is

```text
             to 0    to 1
from 0        -g      g
from 1         l     -l

Q = [[-g, g], [l, -l]]
P(t) = exp(Q t)
```

`g` is the rate of 0-to-1 changes and `l` the rate of 1-to-0 changes. For DNA
presence these are changes of presence, not evidence for a particular sequence
origin mechanism. For exonic status they are conditional changes in annotated
status. For the intron character they describe intron presence at the comparable
position. A single character may change repeatedly, including reversal.

Branch length `t` determines the available evolutionary distance. Rates inherit
its units. Sequence-substitution branch lengths do not yield changes per year.
The model does not estimate a species tree or integrate tree uncertainty.

ER constrains `g = l`; ARD estimates them separately. The foreground model scales
both rates on designated branches by a common multiplier. It does not compare
an exon-shuffling mechanism against an exon-splitting mechanism. In the current
implementation, parameters are shared within a gene family and character class,
not estimated independently for every site and not fitted to a whole transcript
as one multi-state object.

## 4. From observed states to likelihood

For a tip `i`, define the compatibility weights for latent state 0 and 1:

```text
observed 0: (1, 0)
observed 1: (0, 1)
unknown:   (1, 1)
```

Unknown allows either state. These are likelihood weights, not a normalized
50:50 probability assignment. Inapplicability is retained as a separate reason
in the matrix; the separate binary analysis treats it as unobserved, not as 0.
This does not constitute a coupled evolutionary model of DNA and exon status.

For one character, the unconditioned likelihood is

```text
L_k = sum over node states x [
        root_prior(x_root)
        * product over branches (u,v) P_uv(x_u, x_v)
        * product over tips i compatibility_i(x_i)
      ].
```

The implementation evaluates this sum with pruning, in log space. It does not
first pick an ancestral history and then fit rates to that chosen history.
Posterior state probabilities use inside and outside messages. The diagram's
node pies are computed by exhaustive enumeration for an explicit four-taxon
example and tested against these production algorithms.

For an eligible independent-character set `K`, the likelihood factors over
characters. The supported discovery scheme determines whether each term must be
conditioned on ascertainment (for example, at least one observed presence).
Do not apply an unconditional product to an ascertained site set. Corrections
condition on the observed mask; they do not model all causes of annotation loss
or homology-detection failure. See [statistical model](statistical_model.md).

Known dependent characters included together block the independent-character fit,
AIC, intervals, tests and posterior output. Dependency metadata is not a claim
that one physical mutation caused them. Removing that metadata would not make
the characters independent.

## 5. What can be inferred?

* Ancestral states of individual DNA, exonic-status or intron-position characters,
  conditional on their evidence and model.
* Minimal changes and compatible branch placements by maximum parsimony, or
  conditional state/branch-endpoint probabilities and expected transitions by
  the optional likelihood model when estimable.
* Comparisons of gain/loss rate models and a specified foreground rate multiplier,
  subject to validity checks and currently unassessed finite-sample calibration.

A branch-endpoint difference, the probability of at least one transition, and an
expected transition count are different quantities. A branch can begin and end
in the same state after two or more transitions. Separate node marginals do not
form a sampled joint history. Alternative parsimonious placements are not added.

Several structural characters can change in the same gene, and several can be
affected by one mutation. Therefore neither altered exon counts nor sums across
character classes are physical mutation counts. Intron gain can describe exon
splitting without adding a second event. Multiple cutpoints remain multiple
characters; no compound-event interpretation is produced.

## 6. What the current model does not infer

It does not jointly reconstruct complete ancestral transcripts, exon-shuffling
mechanisms, gene duplications, selection or phenotypic causes. It does not assign
probabilities to homology alternatives or annotation errors and integrate them
into the tree likelihood. An incorrect but definite tip call can produce an
apparently confident reconstruction. This is a real limitation, not a graphical
uncertainty that a pie chart solves.

## Reproducible illustration

`intraphy explain --output-dir guide` writes three figures and
`illustrative_model.json`. The probability example fixes `g=0.4`, `l=0.2`, root
presence `2/3`, internal branch lengths `0.4` and terminal lengths `1`. These are
teaching parameters, not fitted estimates. Intron tip states are A=1, B=0, C=1,
D=1. The other two columns demonstrate different questions and are not pooled
into this probability calculation. Every drawing is schematic in distance.

Implementation: `reporting/methods/example.py`, `inference/ctmc.py`,
`inference/posterior.py`, `inference/layer_fitting.py`,
`observations/roles.py`, and `observations/junction_matrix.py`.
