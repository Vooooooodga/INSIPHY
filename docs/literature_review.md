# Literature Review For INSIPHY

INSIPHY is positioned at the intersection of exon orthology, annotation
completion, gene-structure evolution and phylogenetic synteny reconstruction.

## ExOrthist

ExOrthist infers exon orthology groups and uses three evidence classes:
upstream/downstream intron position and phase, exon sequence, and flanking exon
context. INSIPHY adapts this idea inside a single gene or gene family by scoring
segment sequence, coverage, local order, phase compatibility, splice motif,
strand and adjacent segment context. INSIPHY then places the resulting internal
structures on a species tree.

Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC8379844/

## TOGA And CESAR

TOGA integrates structural annotation with orthology inference at genome scale.
CESAR uses coding-exon-aware realignment with reading-frame information. INSIPHY
borrows the principle that annotation should be checked against genome sequence
with explicit uncertainty. The current package separates observed annotation
from sequence-supported hidden-segment, shifted-splice-site and joined-exon
candidates.

References:

- https://www.science.org/doi/10.1126/science.abn3107
- https://github.com/hillerlab/CESAR

## GenePainter

GenePainter maps intron positions and intron phases onto alignments to study
gene-structure conservation and intron gain/loss. INSIPHY uses the same
biological signal at a smaller structural unit: intron/exon/CDS phase, splice
motif and position inform whether internal segments are homologous and whether
role shifts suggest exonization or loss of coding/exonic role.

Reference: https://genepainter.motorprotein.de/

## Phylogenetic Discrete-Character Models

INSIPHY treats intragenic structural states as discrete phylogenetic
characters. The likelihood calculation follows the same logic as Felsenstein's
pruning algorithm: conditional likelihoods are propagated from tips to root on
a fixed tree. The state process is an Mk-style continuous-time Markov chain
similar to discrete morphological-character models, with biological transition
weights for gain, loss and role shift. Pagel-style discrete comparative models
motivate later foreground/background and correlated-character tests; in the
current implementation this is represented by an invariant-model LRT for each
structural character.

References:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC7803665/
- https://academic.oup.com/sysbio/article-abstract/50/6/913/1665006
- https://royalsocietypublishing.org/doi/10.1098/rspb.1994.0006

## DeCoSTAR, edgeHOG And MLGO

Whole-genome structural phylogenetic methods often ask how adjacencies, gene
order and ancestral genome structures changed on a species tree. DeCoSTAR
reconstructs ancestral adjacencies using reconciled phylogenies. edgeHOG scales
ancestral gene order inference through hierarchical orthologous groups. MLGO
uses maximum likelihood for gene-order analysis. INSIPHY translates these
questions to the intragenic level: segments are treated as gene-internal
markers, and neighboring segment pairs form an intragenic synteny graph.

References:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC5441342/
- https://www.nature.com/articles/s41559-025-02818-0
- https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-014-0354-6

## Exon Evolution Concepts

The key event classes INSIPHY reports are derived from exon evolution concepts:
exon shuffling, exonization, intronization or coding-role loss, split/fusion of
segments, and duplicate-copy divergence. The first version reports these as
candidate structural histories with explicit alternative explanations, because
incomplete annotation and paralogy ambiguity can mimic true structural change.

## Real-Case Data Sources

The first real case studies should use FlyBase bulk genome and annotation files
where possible. The current template records FlyBase `FB2026_02` as the default
release placeholder for Drosophila species with FlyBase genomes. D. teissieri
may require Ensembl Metazoa or NCBI annotation if the exact assembly/release is
not represented in the same FlyBase genome set; the template records
`Prin_Dtei_1.1` as the assembly name to check.

Large genome FASTA and GFF/GTF files should stay outside the GitHub repository.
INSIPHY stores only the manifest, stable gene IDs, release labels, local paths
and provenance.
