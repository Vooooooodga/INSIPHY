# V19 exon configuration model

## Purpose and object

The object is a native exon, or a dependent group of neighboring exons in an
orthologous gene locus. The state is their **ordered configuration**, not an exon
usage measurement. No transcriptome is required. Transcript IDs identify the
provided annotation structures; their counts do not estimate expression.

A configuration contains ordered, non-overlapping exon intervals in a shared
local coordinate system and necessary genomic-material provenance. One exon can
be represented as `[a b]`, another configuration as `[a]--[b]`. The intervals are
not independent characters: the edit acts on the configuration as a whole.
Native genomic coordinates, strand, transcript membership and CDS subintervals
remain in the exon instances and coordinate tables. A display fragment is not a
new biological unit. Unknown annotation is not a deletion.

## Evidence to a local catalogue

1. Resolve genomic FASTA, GFF and the rooted species tree. Orthology and selected
   loci are upstream conditions. Gene-only annotations can remain unknown.
2. Reuse native extraction and coding-family alignment/projection. Align the
   oriented genomic locus windows with MAFFT. Audit exon-to-locus matches using
   minimap2 and exact short-exon copy/orientation checks.
3. Cross-check genomic column correspondence against eligible coding projections.
   A conflict becomes unresolved; it is not resolved by choosing the history with
   fewer events. UTR/noncoding exons do not require a coding projection.
4. Link overlapping whole-exon projections and genomic gap tracts crossing those
   exons into local dependent analyses. A continuous deletion across two exons
   therefore enters one configuration. A flank used only for locating a region
   is not an extra state component.
5. Record actual matched intervals, candidate annotations, native boundaries and
   unknown ranges. Boundaries are homologously projected interbase coordinates;
   an upstream indel shifting raw chromosome coordinates is not a splice shift.

The genomic MSA and its column correspondence are **conditional hypotheses**, not
calibrated posterior homology probabilities. The current implementation does not
integrate alternative alignments. It retains observed copy/orientation conflicts,
qualifies missing sequence with paired flanks, and reports uncertain boundaries.
This does not prove all hidden paralogs or alignment errors have been detected.

The default nucleotide identity, paired-flank length/identity and protein-column
agreement are explicit engineering evidence thresholds. They are recorded in
`exon_evidence_policy.json` and are not universal biological constants. The full
sequence window is not silently trimmed when an alignment resource limit is hit.
An over-budget family remains in the main scope summary as unresolved.

A completely gapped tract is a supported absence only when the spanning positional
alignment and both flanks are qualified. N/ambiguous bases, partial gap coverage
or missing flanks leave the affected material and exon structure unknown. Overlapping
but nonidentical candidate indel tracts are not forced into independent sources.

## Observation states and sensitivity

Each tip has an explicit observation record, distinct from the latent configuration:

* **observed**: a supplied local structure with eligible correspondence;
* **partial**: material/boundary information is incomplete or has plausible alternatives;
* **unknown**: no discriminating local structure is available;
* **coexisting**: multiple distinct local structures actually occur in the supplied annotation;
* **excluded**: unsupported copy/order/correspondence or other scope failure.

For the `annotation` view, exact native local boundaries are conditions. For the
`evidence` view, supported alternative exon/boundary intervals permit additional
configurations. These views are calculated separately, and neither overwrites GFF.
A lower edit count after permitting a candidate does not confirm the candidate.

A missing annotated exon with corresponding DNA is a prediction/uncertainty, not a
confirmed exonization or exon death. True coexistence is not treated as random
misannotation: each actual local configuration is analyzed as an explicit conditional
scenario, with conclusions common to all scenarios reported. V19 does not reconstruct
a joint ancestral transcript repertoire or assign usage probabilities to those scenarios.
CTMC posteriors are withheld for coexisting observations in the current engine.

Unknown tip likelihood weights are one for every allowed state. They are not a
normalized uniform probability prior. Partially unknown ranges permit variation
only there. Masking does not join separated exons or create new splice connections.

## Material identity, introduction and deletion

A maximal qualified insertion/deletion tract receives a source identity. Its lifecycle
is unintroduced, present, or deleted. These are components necessary to distinguish
legally different structural edits, not a parallel DNA presence analysis.

* A source that is present may be deleted.
* A deleted source cannot return through an inverse edge.
* An unintroduced source can enter only under an explicit introduction opportunity.
* Exonic structure may appear/disappear repeatedly on material that remains present.

For each candidate source, V19 explicitly enumerates one introduction opportunity:
either inherited at the root, or available on one specified tree branch. Root-inherited
sources are present at the root; branch-introduced sources are initially unintroduced.
Introduction is a CTMC edit on the designated branch and is not allowed elsewhere.
The opportunity prior is uniform over the root and all tree branches, before observing
tip states. Incompatible opportunities may be removed for computation, but their prior
mass is not redistributed. Multiple sources have a product opportunity prior.

This is a **declared finite source-opportunity model**. It is not an exact stochastic
Dollo immigration model with births integrated over arbitrary continuous positions,
not a calibrated prior on real sequence origin, and not evidence that a mutation
occurred exactly once. Opportunity availability does not force the edit to occur.
Its sensitivity and empirical adequacy require investigation.

Parsimony enumerates those same legal origin scenarios but does not treat the prior
as an edit cost. The root is not assigned the reference species' structure. For CTMC,
the root is uniform over configurations compatible with each origin scenario. This
root prior is explicit, not estimated from a small single-gene data set. Branch lengths
inherit their input units; substitution distances are not years. Unit lengths must be
requested explicitly when real branch lengths are unavailable.

## Elementary edit registry

| Edit | Definition | Counting |
|---|---|---|
| split | Add one declared separator within an exon on retained material | One edit, not a loss plus two births |
| fusion | Remove one separator between adjacent exons with retained intervening material | One edit, not an extra intron-loss count |
| donor/acceptor shift | Change one supported homologous endpoint while retaining the other | Each endpoint change is a basic edit |
| exonization | Introduce an exon at a declared exon candidate on retained material | Not a claim of transcription or function inferred from DNA alone |
| exon_inactivation | Remove an exon from the configuration without removing its DNA | Annotation-conditional unless independently supported |
| dna_insertion | Introduce one defined source tract, possibly splitting an exon or carrying an exon | Structural consequences do not add edits |
| dna_deletion | Delete one continuous source tract and normalize all affected exons | Two affected exons can still be one interval edit |

Precise deletion of intervening intronic DNA can fuse exons in the same deletion
edit. It is not counted again as fusion. Two disjoint spacers disappearing remain
two elementary changes even on one branch; V19 does not infer a compound cDNA event.
A split by internal intronization need not conserve the combined mature exon length.

Retained-source fusion/exonization and material deletion are different paths; they
are not blindly pooled as a single generic reversible absence state. Actual CDS
sequences are joined according to annotation. Internal split-codon bases are retained;
phase 1/2 precise deletion is not rejected. Translation exceptions and unsupported
translation tables are recorded rather than silently translated with code 1. UTRs
are not rejected for lacking ORFs. Coding diagnostics do not change transition rates
or force frameshift paths to zero. No NMD or fitness conclusion is inferred.

Opportunity weights split one type/opportunity's rate among distinct legal target
configurations. Repeating a candidate row does not increase total hazard. Separate
eligible exons, endpoints or source tracts are distinct opportunities. Thus the rate
unit is per eligible declared edit opportunity per unit branch length, and depends
on the declared finite catalogue, not an unobserved whole-genome opportunity universe.
There are no separately estimated NHEJ/MMEJ/NAHR enzyme rates.

## Finite state space

States include all valid non-overlapping configurations on the declared candidate
boundary/exon catalogue, combined with compatible source lifecycle states. Endpoints
include observed/candidate exon boundaries, separator endpoints and material boundaries.
Necessary intermediate ancestral configurations are included; only observed terminal
states or only shortest paths would be insufficient for CTMC.

Enumeration is complete **relative to this explicit finite catalogue only**. It does
not discover every possible ancestral exon, every position in an unalignable intron,
unknown copy origins or every possible boundary. New exons extinct in all sampled tips
may be absent from the catalogue entirely. Candidate alphabet choice is a model condition.

The maximum state count, candidate span budget, source tract budget and origin/scenario
limits are recorded. Exceeding a limit returns `state_space_incomplete` (or an explicit
scenario-limit reason), not a truncated normalized probability. There is no claim of an
infinite-state process or finite-state-projection error bound.

## Generalized parsimony

Shortest directed paths in the elementary-edit graph define branch cost matrices.
Positive unit edit costs are the CLI baseline; the API can evaluate explicit alternative
positive costs. A cost is not a probability. Generalized Sankoff and outside messages
retain all globally optimal node states and branch endpoint pairs. Edit identities on
all shortest branch paths obtain minimum/maximum occurrence bounds.

`required` means every relevant optimal history contains that edit on that branch.
`possible` means at least one does. Neither is a frequency/probability estimate. Marginal
possible placements are not added. One lexically deterministic, jointly compatible
witness is saved and explicitly labeled as representative, not the unique truth.
Separate annotation/coexistence scenarios are not silently replaced by their simplest one.

Per-gene output aggregates minimum edits only within resolved local units. It separately
reports unresolved scope and never labels the sum a number of physical mutations or a
fully reconstructed ancestral gene. Independent local analyses cannot jointly reconstruct
unmodeled dependencies, such as a single mechanism deleting several disjoint introns.

## Finite-state probability model

For a legal edit e mapping state S to S':

```
q(S,S') = sum_e lambda[type(e)] * weight(e | S)
q(S,S) = -sum_(S' != S) q(S,S')
P(t) = exp(Q*t)
```

The observation likelihood sums over all internal configurations using log-space tree
pruning. It also sums over the declared origin scenarios with their prior masses.
It does not fit rates to a previously selected minimum-edit history. Outside messages
compute node marginals and branch endpoint probabilities. Optional marked block matrix
exponentials compute expected edits by type. The no-jump expression gives the probability
of at least one edit. Those quantities differ from different branch endpoints.

No posterior combines uncertainty in the input species tree, homology or arbitrary
annotation error. Probabilities use explicit fixed rate files; no automatic default
biological rates exist. Both the input rate provenance and all assumptions are saved.

## Parameter estimation and tests

`fit-exon-rates` takes a prespecified collection of genes and fixed relative operation
rates. It fits one global scale; a foreground comparison adds one shared branch multiplier.
Every local unit of a gene stays together under the gene bootstrap. Overlapping/duplicate
catalogues are rejected. Conditional independence of non-overlapping local units is an
assumption; a gene bootstrap does not magically prove it or repair all model misspecification.

Multiple starts/bounded optimization, boundaries and local numerical curvature are checked.
At least two genes are required as an implementation floor, not a guarantee of biological
identifiability. Failed bootstrap fits are retained, and intervals are withheld rather
than silently computed on survivors. AIC values are meaningful only for the same fixed
objects, observations, discovery scope and root/origin assumptions.

Automatically built catalogues are labeled `annotation_discovered`. Their fitted scales
are exploratory conditional summaries with no complete discovery/annotation-error correction.
The code never emits asymptotic P values. Foreground Monte Carlo P is allowed only for
an independently declared finite catalogue, compatible missingness masks, successful fits,
and successful refits for every null simulation. The simulated tip emission hides
unintroduced versus deleted source states just as the actual observation does. Partial
observation or coexistence rules not replayed exactly block that test. The P-value's
Monte Carlo resolution is reported. A single declared contrast is not a whole-genome FDR.

This conditional simulation is not raw FASTA/GFF discovery-pipeline calibration. Relabeling
an ascertained catalogue as independent to obtain P values is invalid. Independent biological
benchmarks and discovery-aware operating characteristics are still required for publication.

## Scope not claimed

Local duplications and inversions are detected as correspondence anomalies but their copy
phylogeny/rearrangement histories are not modeled. Exon shuffling, gene fusion, circular or
trans-spliced structures, complete ancestral isoform repertoires, RNA usage, selection and
repair-pathway mechanisms are not inferred. A single reference genome does not establish
species fixation. Low-information or unresolved data must not be called conserved.

## Sources and implementation ownership

The biological specification follows the preceding project plan. Relevant primary methods:
ReSplicer (Bocco & Csűrös 2016, doi:10.1093/gbe/evw157), CESAR2 (Sharma et al. 2017,
doi:10.1093/bioinformatics/btx527), ExOrthist (Márquez et al. 2021,
doi:10.1186/s13059-021-02441-9), and ThorAxe (Zea et al. 2021,
doi:10.1101/gr.274696.120). They motivate coordinate, boundary and exon correspondence
choices, not validation of this specific new configuration model. No third-party algorithm
source was copied. Matrix exponentials/optimization use SciPy. The proposed finite model,
opportunity rules and observation qualification must be evaluated on their own merits.
