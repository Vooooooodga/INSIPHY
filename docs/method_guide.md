# Exon structural evolution: current visual guide

`intraphy explain --output-dir guide` generates three scientific plates and one
implementation plate. Add `--png` with CairoSVG installed for preview images.
No external aligner or biological input is needed for these synthetic examples.

## 1. Purpose and range

![Overview](figures/method_overview.svg)

The same tree and exon structures show a split, a fusion, a boundary difference
and DNA deletion. A no-change annotation alternative is not automatically a
confirmed correction. Complex copy/order histories are outside the core model.

## 2. How corresponding structures and alternatives are obtained

![Correspondence](figures/homology_inference.svg)

The diagram follows the implementation: full available genomic locus MSA,
coding-projection corroboration and nucleotide copy/orientation checks. Compatible
one-to-many coverage differs from competing copies. Alternatives are whole
annotated configurations, not freely recombined endpoint masks. This version
still conditions on one qualified alignment; the graphic does not claim joint
alignment inference or automatic use of CESAR2.

## 3. From structures to a probability model

![Model](figures/phylogenetic_model.svg)

States are concrete local exon configurations and arrows are allowed edits.
The tree's ancestral probabilities and branch quantities are computed by the
same installed configuration engine using explicitly illustrative fixed rates.
They are not fitted biological estimates. All example inputs, parameters, history
and normal result figures are provided under `figures/example_*`; numeric data
are also in `figures/illustrative_model.json`.

## Actual data results

`intraphy visualize` emits a parsimony view and, for a CTMC run, a separate
probability view. The latter reads saved ancestral probabilities and distinguishes
endpoint difference, at least one edit and expected edit number. Altering a saved
probability changes the matching SVG element; the renderer does not refit.
Internal-node state marginals do not constitute one jointly reconstructed whole
gene or a full ancestral transcript repertoire.

Solid boxes are annotation, dashed boxes are alternatives, multiple solid lanes
are actually coexisting annotations. Ribbons require actual paired nucleotide
coordinates. Missing DNA is not filled with a fictitious homologous intron.

The old binary-model plates are archived under `figures/legacy_v018/` and remain
available only through the explicitly legacy guide. They do not explain V19.
