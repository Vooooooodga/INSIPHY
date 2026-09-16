# Publication Gap

## Current status

INSIPHY v0.13 is a single-copy method prototype. It accepts an
upstream ortholog set, reconstructs exon-like correspondence from genome
sequence and annotation, builds explicit structural sites, and locates changes
under equal-cost parsimony on a supplied rooted tree. ER/ARD and foreground
CTMC analyses remain optional.

The code can support a methods manuscript description of the model and the
RpL32 conserved-control demonstration. General claims of biological accuracy
still require a broader real-data evaluation.

## Implemented

- OrthoFinder single-copy orthogroup import;
- genome/GFF extraction and annotation-completion candidates;
- exon-like correspondence using sequence and local structural evidence;
- explicit `exon_presence`, `exon_role`, and `splice_junction` matrices;
- all optimal parsimony node states and required/possible branch placements;
- Felsenstein pruning with missing-state marginalization;
- ER/ARD and foreground nested model comparisons;
- profile-likelihood intervals, LRT P values, and BH q values;
- node-state and branch-transition posteriors;
- expected gain/loss counts;
- qualitative branch-placement and optional conditional-probability figures;
- colorblind-aware color and pattern encodings;
- isolated experimental multi-copy code path.

## Required empirical work

1. **Curated positive single-copy cases**
   Add genes with literature-supported exon gain/loss, exonization, and
   split/fusion histories. Confirm every EG correspondence by inspecting the
   underlying genomic alignment and splice boundaries.

2. **Larger conserved controls**
   Analyze a panel of stable single-copy genes across the same species. Report
   the frequency of unsupported changes and ambiguous correspondences.

3. **Cross-clade evaluation**
   Include at least one vertebrate, plant, or other metazoan clade to show how
   sequence divergence affects internal correspondence.

4. **Comparator analyses**
   Compare annotation-only, sequence-only, and full sequence-plus-structure
   correspondence. Report which missing annotations are recovered and which
   exon relationships change.

5. **Sensitivity analyses**
   Evaluate canonical transcript choice, alignment backend, correspondence
   threshold, tree branch lengths, and exclusion of ambiguous EG members.

6. **Independent biological review**
   Have domain experts review structural matrices and branch posteriors without
   relying on software-generated mechanism labels.

7. **Runtime profiling**
   Report elapsed time and peak memory as species count, exon count, and
   candidate pair count increase.

## Statistical limits to report

- Parsimony conditions on equal change costs and the supplied root. It retains
  equally optimal alternatives and cannot reveal unobserved reversals.
- Joint parameter sensitivity in optional probability reconstructions remains
  unimplemented; corresponding ranges are NA. Conditional fitted probabilities
  do not incorporate this uncertainty.
- A single gene usually supplies few exon and junction sites, so rate intervals
  can be wide and LRT power low.
- An invariant gene can support a low overall transition rate but cannot
  identify gain/loss asymmetry.
- Neighboring exon and junction sites are biologically coupled, while the
  current likelihood treats sites as conditionally independent.
- Tree topology, branch lengths, and correspondence are conditioned upon.
- A high branch posterior describes model-conditional event placement and does
  not establish a molecular cause.

## Status of the nine reviewed issues

| Issue | Implemented correction | Remaining scope |
|---|---|---|
| Gene bounds truncate search | Preserve parent/linked/search bounds; extend terminal searches within the declared limit | Sequence beyond the configured interval remains unobserved |
| DNA presence treated as exon identity | Separate sequence presence, annotated role, predicted CDS and boundary conflicts | Transcription and isoform usage are not established by genomic projection |
| Unsupported absence | Require a source/target alignment spanning homologous flanks and a gap covering the expected sequence | Missing anchors or unresolved alignment remain unknown |
| Whole-exon grouping misses split relations | Permit ordered, disjoint projections onto a common reference exon | Event-detection sensitivity still needs positive real cases |
| Heuristic splice coordinates | Project actual boundary bases into shared coordinates; require continuous exon support for boundary absence | Boundary shifts and unaligned termini can remain unmatched |
| Lost orientation | Retain relative strand and compose genomic orientation; exclude opposite-strand splice projections | Inversion-specific phylogenetic models are not implemented |
| Inconsistent alignment measures | Count known paired bases consistently; report gaps, blocks, mode and backend | Algorithms can return different alignments despite common metric definitions |
| Ascertainment mismatch | Retain raw states; condition likelihood on observed masks and selection; require an explicit full catalogue | Dependence of missingness on annotation quality remains outside the model |
| Understated uncertainty | Distinguish numerical failures/search limits; label fitted conditional probabilities; remove surrogate sensitivity ranges | Joint parameter-uncertainty propagation remains unimplemented |

The last correction removes an unsupported precision claim. A complete
probabilistic treatment of parameter uncertainty is still outstanding.

## Multi-copy work deferred

Duplicated genes require a reconciled gene tree, orthology/paralogy-aware
internal correspondence, and explicit treatment of duplication and loss. The
existing experimental code is retained for development. Formal multi-copy
claims should wait until these components have dedicated real-data tests.

## Minimum first-paper claim

A defensible first paper can present:

- a new formalization of gene-internal synteny as homologous exon and junction
  sites;
- sequence-assisted annotation completion inside known orthologs;
- qualitative phylogenetic event placement, with optional shared-parameter
  likelihood models and explicitly conditional branch probabilities;
- several deeply curated positive single-copy examples;
- a larger conserved-control panel;
- transparent limitations and sensitivity analyses.
