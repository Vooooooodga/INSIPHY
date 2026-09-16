# Method Overview

## Scope

INSIPHY studies internal structure within an upstream-defined homologous gene
set. The formal v0.11 scope contains one orthologous gene per species and uses
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

The extracted representation contains:

- exon/CDS/UTR intervals;
- intron intervals and splice motifs;
- exon phase and frame status;
- transcript order;
- genomic sequence for every interval.

Sequence-supported intervals can recover a candidate missing exon or a shifted
boundary. Annotation absence alone remains uncertainty. A biological absence
requires sequence-level evidence that the homologous unit is missing from the
searched gene interval.

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

Users may select `internal`, `minimap2`, or `miniprot`. The internal
backend uses Biopython `PairwiseAligner`. Candidate scoring can run in
parallel.

Pairwise evidence is assembled into internal homology components (`HC_*`).
Components supported as exon/CDS/UTR or as a sequence-supported candidate
exonic source are exposed as exon-like groups (`EG_*`). Figures and
phylogenetic models use `EG_*` as the biological unit. Introns remain context
unless they provide explicit evidence for a homologous sequence changing its
exonic role.

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

Junction gain and loss provide the structural pattern used to study exon
splitting and fusion. A specific mechanistic history requires additional
sequence and functional evidence.

## Stage 4: phylogenetic likelihood

All sites in one family-layer share ER or ARD CTMC parameters. Likelihood is
computed on the fixed species tree with Felsenstein pruning. The software
reports:

- maximum-likelihood gain and loss rates;
- profile-likelihood intervals;
- ER-versus-ARD or homogeneous-versus-foreground LRT;
- node-state posterior probabilities;
- branch endpoint-transition posteriors;
- expected gain and loss counts.

Ancestral states are probabilistic reconstructions. Several histories can
produce the same terminal pattern, especially with few species or short gene
structures. Output tables preserve this uncertainty.

## Event vocabulary

The formal output uses a compact structural vocabulary:

- `exon_gain`, `exon_loss`;
- `exonic_role_gain`, `exonic_role_loss`;
- `splice_junction_gain`, `splice_junction_loss`;
- `exon_split_pattern`, `exon_fusion_pattern`.

These names describe modeled changes in observable structure. They do not
assign transposon activity, gene conversion, selection, or other causal
mechanisms.

## Visualization

The synteny panel follows the standard comparative-genomics grammar:
species/gene tracks contain ordered exon blocks, and curves or ribbons connect
homologous EGs between tracks. Introns are thin gray context spans.

The phylogeny panel places posterior structural changes on branches. Marker
size and opacity scale continuously with posterior change probability. The
integrated panel aligns tree tips with gene-structure tracks.

Default EG colors use the Okabe-Ito colorblind-aware palette. A
`--correspondence-encoding pattern` mode uses hatching, dots, grids, line
styles, labels, and connecting curves.

## Real-data control

The current complete single-copy demonstration uses RpL32 orthologs from five
Drosophila genomes. Two well-supported exon groups are conserved in sequence
presence and exonic role. A short D. melanogaster interval has ambiguous
cross-species correspondence and remains unknown. The data contain no
identifiable direction-changing site in the fitted layers, so branch-change
posteriors are near zero and ER/ARD P values are withheld.
