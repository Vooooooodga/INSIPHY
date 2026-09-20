# Method

## Question and inputs

IntraPhy analyzes local gene structure in upstream-defined single-copy orthologs.
The inputs are genomic FASTA, GFF3/GTF gene annotations and a rooted species phylogeny.
Optional ortholog FASTA records select genes; a user-written manifest is not required. The reporting unit is a gene family. The statistical
observation is a structural character. The tree and orthology are conditions of
the analysis; neither is estimated by this program.

## Preserve annotation before inference

The extraction layer retains native genomic coordinates, strand, CDS and UTR
intervals, introns, feature ownership, original attributes and annotated transcript
isoforms. Identical coordinate structures may share computation while their path
identities remain traceable. Repertoire observations indicate use by at least one
annotated isoform. Canonical observations concern the explicitly selected path and
must be analyzed in a separate frozen matrix.

The supplied annotation is not a complete census of all tissues or conditions.
Role state 0 requires informative covering paths in the selected annotation view.
It does not establish that an interval is never transcribed. Annotation dropout
can change an annotation-conditional character even when biological RNA use did
not change. Genomic sequence alone cannot resolve this distinction.

## Local sequence correspondence

For coding sequence, a family protein alignment is projected through codons to
genomic intervals. MAFFT L-INS-i is the default, with E-INS-i as an explicit
alternative. A protein match supports only its actual projected CDS blocks; it
does not extend through unmatched UTR or intronic sequence.

Nucleotide candidates retain scores, aligned blocks, orientation, alternatives
and search completeness. Short feature-bounded alignments are candidates. A hard
position observation requires a unique ordered flank configuration and alignment
within the actual bounded genomic interval. Predicted coding projections do not
become confirmed annotated exons.

Ordered candidate chains retain near-optimal alternatives. A local membership
identifies a homologous unit, a native occurrence and actual matched subintervals.
Partial matches do not expand to the entire parent exon. Complementary ordered
reference coverage is distinguished from repeated coverage of the same reference
positions. Genomically separated fragments are not labeled repeats solely because
they occupy distinct genomic intervals. Transcript aliases at the same physical
interval are not separate repeat instances.

Sequence correspondence, conserved intron position and common ancestry of an
intron insertion are distinct questions. Two introns may occupy corresponding
coding positions without demonstrable homology of their complete internal DNA.
Parallel gains at one position remain possible. The reference coordinate system
is not an ancestral-state assertion.

## Define structural characters

The three separately analyzed layers are:

| Layer | State 0 | State 1 |
|---|---|---|
| `exon_presence` | Corresponding DNA unit demonstrably absent | Corresponding DNA unit present |
| `exon_role` | Non-exonic in adequately covering annotated isoforms | Exonic in the selected supplied annotation |
| `splice_junction` | An informative path continuously spans the corresponding position | A annotated isoform contains the corresponding junction |

Unknown and inapplicable states are retained separately in the observation
metadata. DNA absence makes the corresponding exon role inapplicable. Junction
absence requires positional/path evidence; a missing GFF feature or missing
alignment hit is insufficient.

The character catalogue records identities, labels, discovery rules, counting
units and dependence. The coordinate sidecar records actual evidence. Multiple
plotting blocks do not create additional characters. Different characters with
identical tip patterns are not merged merely because the patterns are equal.
Overlapping aligned sequence intervals are conservatively linked for dependence;
sharing a transcript alone does not imply dependence. Declared linked junctions
remain linked even though compound-event interpretations have been removed.

## Assemble the structural character matrix

The full structural matrix is retained. Default `--analysis-range all` does not
select characters by coverage. The optional high-coverage analysis uses known 0
and 1 entries divided by the complete supplied tree panel. This is a sensitivity
analysis, not an accuracy threshold. Native coordinates and transcript adjacencies
are never changed by character selection. A missing intermediate exon cannot
create a new splice connection between its neighbors.

Frozen-matrix analyses use their archived coordinate sidecar when available;
current annotation is not reread to manufacture coordinates for an old matrix.
Missing sidecars remain explicit limitations. See [scope](scope_policy.md).

## Reconstruct elementary changes

Equal-cost Sankoff reconstruction reports all node states and branch endpoint
pairs compatible with a global minimum. A change is `required` when every optimum
requires that directed change on the branch, and `possible` when at least one
optimum permits it. Neither label is a confidence probability.

One deterministic compatible optimum is provided for inspection. Its total cost
is checked against the dynamic-programming optimum. It is not a sampled history.
Per-layer minimum-change totals do not add alternative possible placements.
Junction gain and the corresponding exon-split description are one character
change. Multiple cutpoints are not summarized as one compound event.

Separate layer reconstructions need not form a valid complete ancestral
transcript. An applicability audit flags incompatible singleton reconstructions;
it does not repair the history or establish a joint model.

## Likelihood and interpretation

Optional binary CTMC models compare gain/loss rates or a foreground multiplier,
conditional on the same frozen matrix. Known dependent characters jointly used
in one family-layer prevent this independent-character analysis before fitting.
Rates, AIC, intervals, LRT and posterior outputs are all unavailable in that case.
See [statistical model](statistical_model.md) for ascertainment and fit diagnostics.

The method reports local structural evidence and conditional reconstruction.
It does not infer mutation-event counts, rearrangement mechanisms, functional
consequences or selection. The [references](references.md) distinguish existing
sequence/structure comparison, character coding and phylogenetic methods from
this implementation's unvalidated claims of biological performance.

The [illustrated model explanation](model_bridge.md) connects biological states
to transition rates, tree likelihood and ancestral-state probabilities.
