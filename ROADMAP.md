# INSIPHY Roadmap

## Current v0.11 status

- Formal single-copy input from OrthoFinder or a curated ortholog set.
- Sequence-assisted annotation completion inside selected genes.
- Exon-like correspondence with sequence, boundary, phase, and order evidence.
- Exon presence, exonic role, and splice-junction structural layers.
- ER/ARD and foreground binary CTMC models.
- Profile likelihood intervals, LRT, BH adjustment, and branch posteriors.
- Colorblind-aware linked synteny and integrated phylogeny figures.
- Experimental multi-copy implementation retained outside formal claims.

## Next biological work

1. Curate several single-copy genes with published exon gain/loss,
   exonization, or split/fusion histories.
2. Analyze a larger panel of conserved single-copy controls.
3. Compare annotation-only, sequence-only, and full correspondence.
4. Evaluate sensitivity to transcript policy, alignment backend,
   correspondence threshold, and branch lengths.
5. Extend within-element splice-boundary correspondence for partial
   exon split/fusion cases.
6. Add cross-clade real-data examples.

## Next engineering work

1. Profile runtime and peak memory across species and exon counts.
2. Add chunked candidate-pair generation for larger ortholog panels.
3. Add machine-readable JSON schemas for principal output tables.
4. Package accession-level result summaries for releases.
5. Develop reconciled gene-tree statistics before promoting multi-copy
   analysis to formal status.
