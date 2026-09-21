# 0.19.1 audit resolution and acceptance map

This is an implementation correction to the **exon structural evolution** method.
It is not a new exon usage project. It does not assert that all exon evolution can
be identified from genomes or that software tests establish biological accuracy.

| Audit issue | Implementation | Acceptance evidence |
|---|---|---|
| A01 asymmetric annotation alternatives | `structure/alternatives.py`, typed `ConfigurationAlternative`, compatibility of complete candidates | Identical-DNA single-end, double-end, mis-split and mis-fused annotations retain an alternative; changed-cut sequence controls retain a structural change |
| A02 terminal annotation dropout crops available DNA | `alignment.py` searches full supplied/extracted locus and writes `search_windows.tsv` before/after execution | Both missing terminal exons remain in the 1140-base audit input; budget failure records 0 searched bases, not a shorter successful search |
| A03 changed labels double-count a physical region | `structure/validation.py` shared by JSONL reader/writer, preparation/inference and pooled loading | Same region under different IDs, or across files, is rejected; adjacent regions accepted |
| A04 unary tree subdivision changes prior | `tree_context.py`, normalization in parsimony, CTMC and simulation; mapping output | 1.0 versus 0.5+0.5 unobserved branch gives identical likelihood/root posterior; explicit rate boundaries not silently discarded |
| A05 state growth and redundant numerical work | Geometry-count DP, cached index, sparse directed distances, bounded numerical cache, online origin mixture | Four-exon audit case retains all 556 states; small-space dense reference and semigroup checks pass; explicit smaller cap still refuses |
| A06 module insertion split into arbitrary births | `InsertionPayload`, one source opportunity with complete supplied exon payload | Two-source-exon payload introduced in one edit; deleted material cannot return; duplicate payload does not raise hazard |
| A07 legal repeated-gene bootstrap refused | Original collection qualification separated from positive sampled labels | `{g:2,h:0}` reaches weighted numeric fitting; actual failed/boundary fits are still retained |
| Nested alternative below null silently clamped | Null optimum is an explicit alternative-model start; negative improvement outside tolerance is an error status | Numerical nesting check; true optimizer failures are not zero-statistic biological results |
| Explanation figures remained V18 | Four regenerated current plates; old files archived | Current manifest, current model ID, computed-example values and XML checks |
| CTMC ignored by result rendering | Separate parsimony/CTMC SVGs; ancestor glyphs, bars, branch quantities read from saved fields | Altering a saved probability changes its corresponding SVG element; rendering cannot call inference |

## Important interpretation of A01

Whole-structure alternatives originate from another supplied annotated path. To be
allowed they must have compatible known material status, paired exon sequence,
qualified location flanks and matching two-base contexts on either side of every
changed cut. This is a conservative engineering rule, not a universal splice
motif rule or a calibrated homology probability. Noncanonical motifs are not
rejected solely for being noncanonical.

No empty annotation is promoted into evidence that another species has no exon.
Missing-exon candidates are built from positively annotated source exons. Atomic
alternatives belong to a particular native configuration, preserving real
coexisting paths rather than forming arbitrary mixtures of endpoints.

An alternative does not confirm a transcript. Failure to qualify an alternative
does not confirm biological divergence. Both annotation-condition and
candidate-compatible histories remain available; events state their conditions.

## Changes that remain conditional rather than newly solved biology

The following are not implemented or empirically certified by this patch:

- Integration over alternative genomic alignments, calibrated annotation-error
  probabilities, or full species-tree uncertainty.
- Source copy genealogy, arbitrary exon shuffling/inversions, or automatic proof
  that several new exons arrived in one physical insertion.
- Exact continuous-time Dollo immigration. The finite root/branch opportunity
  prior remains explicit, now insensitive to redundant tree nodes.
- All possible extinct ancestral exons or arbitrary boundary positions outside
  the finite candidate catalogue.
- Independent published-case benchmarking or complete raw-discovery operating
  characteristics. The earlier real-case plan is not fabricated as completed.

## Compatibility and reproducibility

`intraphy.exon-configurations/2` and `exon_configuration_v2` prevent old candidate
masks or results being silently relabeled. Rebuild from original genomic inputs,
or construct a new explicit catalogue with a documented biological scope. The
package version is 0.19.1; the Python namespace remains `intraphy` only.

Generalized Sankoff, fixed-rate CTMC and rate-scale fitting share the same
configuration/edit objects. A multiple-exon effect is not a second copy of its
underlying interval edit. Possible placements remain non-additive. Rate units
remain per declared eligible opportunity per branch-length unit, not per year
unless the input branch lengths and model actually have that meaning.

The local delivery contains no operation that pushes a repository. CI files are
ordinary source files for later user-controlled installation/publication.

Mixed foreground/background regimes inside a subdivided unary edge are preserved
as input information but rejected by the current configuration probability engine
and simulator. They are not silently collapsed, and their extra boundary cannot
change the prior between a null and alternative model. Uniform regimes across
redundant segments are merged and remapped normally.
