# INSIPHY Method Overview

INSIPHY treats a gene copy as an ordered set of homologous internal segments.
Segments can be annotated exons, CDS intervals, introns, or sequence-supported
hidden candidates. The method reconstructs how these segments and their
adjacencies changed on a fixed species tree.

## Evidence Layers

1. **Annotation extraction**: genome annotation supplies observed exons, CDS
   intervals and introns for each target copy.
2. **Sequence-supported completion**: local sequence evidence can mark a
   missing annotated segment as a hidden-segment candidate or an annotation
   conflict.
3. **Homologous segment grouping**: pairwise segment similarity, boundary
   compatibility and local order are combined into HSG assignments.
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
state reconstruction and a lightweight competing-model score. This gives a
transparent first-pass event history. The planned publication version will add
calibrated likelihood-style model comparison and simulation benchmarks.
