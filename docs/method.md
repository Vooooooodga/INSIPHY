# Method Overview

## Scope

INSIPHY studies internal structure within an upstream-defined homologous gene
set. The formal v0.13 scope contains one orthologous gene per species and uses
the species tree. Gene discovery, genome-wide orthology inference, expression
analysis, and molecular-mechanism assignment remain upstream or downstream
tasks.

Multi-copy routines remain available under
`--analysis-scope experimental-multicopy` for continued development. Their
results are outside the formal single-copy method.

## Stage 1: recover observed gene structure

`build-case` extracts each selected gene from genome FASTA and GFF3/GTF.
Canonical-transcript mode chooses a reproducible transcript by CDS or span
length; all-transcript mode retains every annotated transcript.

The extracted representation uses complete exon intervals and contains:

- complete exon intervals with CDS and UTR subinterval attributes;
- intron intervals and splice motifs;
- exon phase and frame status;
- transcript order;
- genomic sequence for every interval.

The search interval includes flanking sequence and parent-linked GFF features
outside the original gene bounds. Terminal searches can extend to a declared
maximum; original annotation bounds remain recorded. DNA correspondence supports
sequence presence. Protein-to-genome splice projection can additionally support
a predicted coding exon, without confirming transcription. A deletion call
requires a source/target genomic alignment spanning both homologous flanks and
a source-only gap covering the expected exon. This comparison uses genomic
mapping with exact local alignment offsets and source anchors from one
transcript path. DNA deletion evidence can be evaluated without a protein
reference. Annotation absence stays unknown.

## Stage 2: infer gene-internal correspondence

The correspondence problem asks which exon-like intervals in different
orthologs derive from the same ancestral sequence unit.

Candidate pairs are evaluated with:

- nucleotide or protein-to-genome alignment;
- aligned coverage and length compatibility;
- left and right local context;
- splice-boundary compatibility;
- exon phase and reading-frame compatibility;
- strand and order consistency;
- reciprocal-best support.

Users may select `internal`, `mafft`, `minimap2`, or `lastz` for correspondence.
MAFFT overlap projection is the default. It trims terminal overhangs from a
pairwise multiple-alignment result; it is not a local alignment algorithm.
Exon/non-exonic comparisons use a separately specified local context mapper
(`--context-aligner minimap2` by default); they are not aligned across entire
long introns with MAFFT. Short candidates without adequate mapper support remain
unknown, a sensitivity limitation recorded with the selected backend.
The internal backend uses Biopython `PairwiseAligner` and refuses sequence
pairs above its declared dynamic-programming limit. miniprot is restricted to
reconstructed protein-to-locus searches.

Pairwise evidence is merged from close to distant species on the supplied
tree. A merge must satisfy direct sequence support across the two profiles;
split fragments must project to ordered, disjoint parts of the same reference
exon. They need not resemble one another.
This prevents a chain of unrelated pairwise matches from forming one group.
`HC_*` and `EG_*` are stable identifiers. Their members remain traceable to
complete exon intervals and aligned coordinates.

Tree-guided progressive summarization first evaluates close relatives, then
within-clade and deeper comparisons. This applies the progressive-alignment
principle inside the supplied gene set. It does not run a whole-genome
alignment and does not search for additional homologous genes.

## Stage 3: construct structural sites

Three observation layers separate distinct biological changes.

### Exon sequence presence

For each EG and species:

- `present`: homologous sequence is observed;
- `absent`: explicit searched absence is available;
- `unknown`: correspondence, coverage, or annotation evidence is inadequate.

Gain/loss in this layer describes appearance or disappearance of the homologous
sequence unit within the gene.

### Exonic role

For each observed homologous sequence:

- `exonic`: CDS, exon, UTR, or noncoding exon;
- `not_exonic`: homologous sequence is present in intronic/non-exonic context;
- `unknown`: sequence role cannot be assigned reliably.

A change in this layer is compatible with exonization or loss of exonic role.
The model reports the role transition and leaves its molecular cause open.

### Splice junction

Transcript paths and homologous flanking units define splice-junction sites:

- `present`: an annotated intron separates the corresponding units;
- `absent`: the units are directly contiguous in an exon path;
- `unknown`: the relevant units cannot both be placed reliably.

Every junction is keyed by actual donor/acceptor bases projected through supported
alignment blocks into common reference coordinates. Unaligned boundaries and
opposite-strand matches do not establish collinear junction correspondence.
Distinct junctions inside one homologous exon block remain distinct.
Unknown correspondence breaks adjacency. Junction-state changes describe the
structural pattern used to study exon splitting and fusion.

## Stage 4: locate structural changes

Default analysis uses equal-cost maximum parsimony on the supplied rooted tree.
Gain and loss each cost one change; unknown observations allow either state.
An inside/outside calculation retains all globally optimal node states and
branch endpoint pairs without enumerating all histories. Repeated observed
patterns reuse the calculation.

`required` means a directed change occurs on that branch in every minimum-change
history. `possible` means it occurs in some such histories. These labels describe
reconstruction under the cost assumptions; they are not probabilities or P values.
Within-exon junction gain/loss describes a split/fusion pattern. Sequence presence
and exonic role remain separate layers and do not constitute independent events.

Root placement and change costs affect reconstruction. The method does not
estimate calendar time or identify unobserved repeated changes. Alternative
equally parsimonious histories remain explicit.

## Optional phylogenetic likelihood

All sites in one family-layer share ER or ARD CTMC parameters. Likelihood is
computed on the fixed species tree with Felsenstein pruning. The software
reports:

- maximum-likelihood gain and loss rates;
- profile-likelihood intervals;
- ER-versus-ARD or homogeneous-versus-foreground LRT;
- node-state empirical-Bayes probabilities conditional on fitted parameters;
- branch endpoint-transition probabilities conditional on fitted parameters;
- expected gain and loss counts.

Ancestral states are probabilistic reconstructions. Several histories can
produce the same terminal pattern, especially with few species or short gene
structures. Joint parameter-uncertainty propagation remains unimplemented:
sensitivity ranges are `NA` with `uncertainty_status=sensitivity_not_estimated`.
Profile searches distinguish optimization failures from numerical search limits.
Ascertainment conditions on the observed-taxon mask. Complete-universe analysis
requires an explicit site catalogue; unobserved states remain unknown.

## Transition reporting

Default tables report branch placement under parsimony, including alternatives.
Optional probability tables retain both directional endpoint transitions.
Mechanistic explanations remain the user's biological analysis.

## Visualization

The synteny panel follows the standard comparative-genomics grammar:
species/gene tracks contain ordered exon blocks, and curves or ribbons connect
homologous EGs between tracks. Introns are thin gray context spans.

The default phylogeny panel marks required changes with solid symbols and
possible changes with hollow symbols. Optional probability figures use continuous
marker size and opacity. The integrated panel aligns tree tips with gene tracks.

Default EG colors use the Okabe-Ito colorblind-aware palette. A
`--correspondence-encoding pattern` mode uses hatching, dots, grids, line
styles, labels, and connecting curves.

## Real-data control

The real demonstrations contain RpL32 from five Drosophila genomes and dsx from
seven bee genomes. Both use genome and annotation data; dsx retains annotated
transcript alternatives. Results are recorded in `real_data_benchmark.md`.
These demonstrations do not establish event-detection precision or recall
without independently curated event truth.
