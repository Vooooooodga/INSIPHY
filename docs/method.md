# INSIPHY Method Overview

INSIPHY treats each supplied gene copy as an ordered set of biological
structural elements. The primary elements are annotated exons, CDS intervals,
UTRs and sequence-supported candidate exonic regions. Introns are represented
as intervals, splice boundaries, phase and motif context. A non-exonic sequence
interval becomes a primary event object only when sequence and boundary
evidence support a candidate exonization, hidden exon or role-shift history.

The compared gene/copy set is supplied by an upstream homology workflow or a
curated case definition. INSIPHY starts from that set and analyzes structure
inside the homologous genes.

The key biological unit in interpretation is the exon-like structural element.
HSG labels in intermediate files are evidence-graph components used to assemble
correspondence; they are not the named biological unit in figures or manuscript
claims. This distinction matters for introns: an intronic interval can support
splice-boundary conservation, phase compatibility or exonization evidence, but
it is not promoted to a displayed homologous block without role-shift support.

## Evidence Layers

1. **Annotation extraction**: genome annotation supplies exons, CDS intervals,
   UTRs, transcript paths and introns for each target copy.
2. **Sequence-supported completion**: local sequence evidence marks missing
   annotated segments, shifted splice boundaries or joined-segment candidates.
3. **Correspondence evidence graph**: pairwise similarity, boundary
   compatibility, intron phase, splice motif, strand, flanking context and
   local order are combined into graph-based evidence clusters. These clusters
   are promoted to EGs when they correspond to exons, CDS/UTR intervals or
   sequence-supported candidate exonized source intervals. HSG identifiers are
   internal graph labels. Biological interpretation is made through EGs, splice
   boundaries, adjacencies and source-copy context.
4. **Tree-guided progressive interpretation**: segment support is summarized by
   species-tree distance. Close-species support, within-clade support and
   deep-tree support are kept visible in `progressive_correspondence.tsv`.
5. **Intragenic synteny graph**: ordered exon-like elements and their
   adjacencies describe the internal synteny of each gene copy. Introns
   contribute boundary and phase evidence without becoming default homologous
   blocks in the user-facing graph.
6. **Phylogenetic reconstruction**: EG presence, EG role state, adjacency,
   source mixture and copy multiplicity are optimized on the species tree.

The sequence layer can use the internal fallback aligner, minimap2 for
nucleotide segment matching, or miniprot for protein-to-genome style completion
tests when suitable input is provided.

## Evolutionary Questions

INSIPHY asks whether a supplied gene/copy set shows gene-internal structural
change on the species tree:

- gain or loss of a homologous internal segment;
- exonization of intronic or noncoding sequence;
- source joining in a chimeric or duplicated gene;
- split or fusion of adjacent internal segments;
- splice-boundary shift;
- duplicated-copy divergence and copy expansion;
- high-identity paralogous segment matches as ambiguous evidence requiring
  outside support for mechanism-level interpretation;
- annotation gap supported by local sequence evidence.

Candidate events are reported with both low-level state changes and biological
`structural_pattern` labels, such as `exonization_candidate`,
`segment_fusion_or_new_adjacency`, `chimeric_source_join_candidate`,
`copy_duplication_or_expansion` or `ambiguous_paralogous_similarity`.
`mechanism_hypothesis` stores the possible biological explanation, and
`call_scope` separates core structural events from copy-context and ambiguous
evidence.

## Three-Layer Translation

INSIPHY makes the biological-to-statistical translation explicit.

1. **Biological meaning**: a gene-internal event is a change in which internal
   exon-like elements exist, what role they play, where they sit relative to
   neighboring elements, and which source copy they resemble. For example, a
   chimeric gene is represented by adjacent exon-like elements with different
   source labels; copy expansion is represented by multiple related copies in
   the same species; exonization is represented by a non-exonic source interval
   acquiring exon/CDS status.
2. **Mathematical object**: each biological question becomes one or more
   discrete structural characters on a fixed tree. EG presence is
   `present/absent`; EG role is
   `CDS/exon_or_UTR/non_exonic_source/absent`; EG adjacency is
   `present/absent`; source mixture is
   `single_source/multi_source`; copy multiplicity is
   `single_copy/tandem_multi_copy/dispersed_multi_copy/...`.
3. **Statistical computation**: each character is analyzed with a discrete
   phylogenetic model. Sankoff reconstruction places low-cost changes on the
   tree; CTMC/Mk likelihood estimates transition rate; the invariant-model LRT
   asks whether the character supports structural change; bootstrap and
   stochastic mapping calibrate and place the event; foreground/background
   tests ask whether selected branches have elevated structural-change rate.

The real-case outputs keep these layers linked: `candidate_structural_events.tsv`
names the observed structural pattern and mechanism hypothesis, `object_id`
points back to the structural character, and `event_support_summary.tsv` joins
that event to p values, q values, bootstrap evidence and branch-history support
where available.

`element_correspondence.tsv` and `element_phylogenetic_coverage.tsv` are the
main correspondence and tree-coverage tables. `hsg_phylogenetic_coverage.tsv`
remains an internal evidence-coverage table. Figures label displayed exon-like
correspondence groups as `EG_*` and draw introns as gray context spans.

## Statistical Mapping

INSIPHY adapts discrete-character phylogenetic models to intragenic structure:

- **Weighted Sankoff reconstruction** estimates low-cost ancestral histories
  for segment presence, role state, adjacency, source mixture and copy
  multiplicity.
- **Felsenstein pruning / Mk likelihood** fits a branch-length-aware CTMC model
  for each structural character on the fixed species tree.
- **Likelihood-ratio testing** compares an invariant no-change model against a
  one-rate CTMC model and reports p value, fitted rate, AIC and BIC.
- **Benjamini-Hochberg correction** reports q values across structural
  characters.
- **Parametric bootstrap** reports empirical p values for small-tree
  calibration when requested.
- **Stochastic character mapping** samples complete CTMC histories along
  branches and reports posterior summaries for event placement.
- **Foreground/background rate testing** asks whether user-selected branches
  show elevated structural-change rates compared with background branches.

The first publication version focuses on real Drosophila duplicated and
chimeric gene cases plus a conserved negative control. Large-scale calibration
remains a manuscript-readiness task.
