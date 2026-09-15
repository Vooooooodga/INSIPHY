# INSIPHY Method Overview

INSIPHY treats a gene copy as an ordered set of homologous internal segments.
Segments can be annotated exons, CDS intervals, introns, or sequence-supported
hidden candidates. The method reconstructs how these segments and their
adjacencies changed on a fixed species tree.

The method assumes that the compared gene copies have already been selected by
an upstream homology workflow or by a curated case definition. INSIPHY analyzes
the internal structure of that supplied copy set.

## Evidence Layers

1. **Annotation extraction**: genome annotation supplies observed exons, CDS
   intervals, UTRs, transcript paths and introns for each target copy.
2. **Sequence-supported completion**: local sequence evidence can mark a
   missing annotated segment as a hidden-segment candidate or an annotation
   conflict.
3. **Homologous segment grouping**: pairwise segment similarity, boundary
   compatibility, intron phase, splice motif, strand, flanking context and
   local order are combined into graph-based HSG assignments. The sequence
   alignment layer can use the internal fallback aligner, minimap2 for
   nucleotide segment matching, or miniprot for protein-to-genome style
   completion tests when suitable input is provided.
4. **Intragenic synteny graph**: ordered HSG adjacencies describe the internal
   synteny of each gene copy.
5. **Phylogenetic reconstruction**: segment presence, role state, adjacency,
   source mixture and copy multiplicity are optimized on the species tree.

## Evolutionary Questions

INSIPHY asks whether an event can be explained by gene-internal structural
changes:

- gain or loss of a homologous segment;
- exonization of intronic or noncoding sequence;
- fusion of segments from different source loci;
- split or fusion of adjacent internal segments;
- duplicated-copy divergence;
- annotation gap supported by sequence evidence.

Candidate events are reported with both low-level state changes and a biological
`event_class`, such as `exonization_candidate`,
`segment_fusion_or_new_adjacency`, `chimeric_source_join_candidate` or
`copy_duplication_or_expansion`.

The current implementation uses a Sankoff-style discrete character model for
state reconstruction and a branch-length-aware CTMC/Mk likelihood fit for each
structural character. INSIPHY reports fitted event-rate parameters, log
likelihood, AIC/BIC, invariant-model LRT p values, CTMC branch-change
probabilities, bootstrap-calibrated empirical p values when requested,
stochastic-map branch-history summaries, foreground/background rate tests,
candidate branch events and lightweight competing-model scores. The publication
version should expand real accession-level case studies, independent benchmark
sets and large-scale calibration.

## Statistical Method Mapping

INSIPHY adapts four families of phylogenetic methods to intragenic structure:

- **Maximum parsimony / Sankoff reconstruction**: finds low-cost ancestral
  histories for discrete segment states. Here the characters are segment
  presence, segment role, source mixture, adjacency and copy multiplicity.
- **Felsenstein pruning / Mk likelihood**: computes the likelihood of observed
  tip states on a fixed tree under a continuous-time Markov model. Here the
  Mk states are biological structural states such as `present/absent` or
  `CDS/intron_or_noncoding`.
- **Likelihood-ratio testing**: compares a no-change structural model against
  a fitted change model and returns `p_value`, `lrt_statistic` and model
  parameters.
- **Parametric bootstrap**: simulates structural characters under the null on
  the same tree and reports an empirical p value for small-tree calibration.
- **Stochastic character mapping**: samples complete CTMC histories along
  branches, giving posterior summaries for event count and transition type.
- **Ancestral adjacency reconstruction**: whole-genome synteny methods treat
  neighboring genes as phylogenetic characters; INSIPHY applies the same idea
  inside one gene by treating neighboring HSGs as intragenic synteny edges.
