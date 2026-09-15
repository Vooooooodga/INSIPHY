# Literature Review For INSIPHY

INSIPHY is positioned at the intersection of exon orthology, annotation
completion, gene-structure evolution and phylogenetic synteny reconstruction.

## Review Conclusion

The current method direction is biologically defensible: gene-internal segment
homology should be inferred from genome sequence plus annotation context, and
the resulting structural characters should be interpreted on a species tree.
The main remaining biological risk is under-modeling event classes that can
look similar in a single observed gene copy: exonization, intron gain/loss,
splice-boundary drift, tandem exon duplication, gene conversion, processed
retrocopy insertion and assembly/annotation artifacts can all create partially
overlapping evidence patterns. INSIPHY v0.4.0 covers the core duplicated and
chimeric-gene cases, but the publication version should expose these alternative
explanations explicitly in benchmark tables and event confidence summaries.

## ExOrthist

ExOrthist infers exon orthology groups and uses three evidence classes:
upstream/downstream intron position and phase, exon sequence, and flanking exon
context. INSIPHY adapts this idea inside a single gene or gene family by scoring
segment sequence, coverage, local order, phase compatibility, splice motif,
strand and adjacent segment context. INSIPHY then places the resulting internal
structures on a species tree.

Important design points for INSIPHY:

- ExOrthist analyzes exons inside supplied gene orthogroups and can take
  orthologous pairs from tools such as OrthoFinder or Broccoli.
- Its exon match decision sequentially evaluates intron position/phase,
  exon sequence and flanking exon sequence.
- It supports extra, unannotated exons and bona fide curated exon pairs.
- It clusters pairwise exon homology relationships in a graph and reports
  membership-like support.

INSIPHY should keep the same conceptual separation: upstream tools define the
gene/copy set, and INSIPHY performs gene-internal segment correspondence plus
phylogenetic interpretation. INSIPHY differs in scope because it is a Python
method package and CLI for single-gene or gene-family structural histories,
while ExOrthist is a genome-scale Nextflow pipeline.

References:

- https://github.com/biocorecrg/ExOrthist
- https://doi.org/10.1186/s13059-021-02441-9

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

Stochastic character mapping is the relevant next layer for branch event
posteriors. The current INSIPHY branch table reports endpoint posterior
probabilities for parent-child state changes under the fitted CTMC. A
publication version should also sample complete character histories conditional
on the tips, because a branch can contain zero, one or multiple changes even
when the endpoint states are the same. Those samples would provide
posterior-like summaries such as `Pr(change on branch)`, expected number of
changes, most frequent transition type and credible intervals across stochastic
maps.

References:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC7803665/
- https://academic.oup.com/sysbio/article-abstract/50/6/913/1665006
- https://royalsocietypublishing.org/doi/10.1098/rspb.1994.0006
- https://pmc.ncbi.nlm.nih.gov/articles/PMC1403802/

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

The literature suggests the following event inventory should guide benchmarks:

- gene duplication and post-duplication exon-intron divergence;
- retroposition or processed-copy insertion, often with intron loss;
- chimeric new gene formation by exon shuffling or illegitimate recombination;
- exonization of intronic, transposed-element or intergenic sequence;
- intron gain, intron loss, intron sliding and intronization of exonic sequence;
- tandem exon duplication and partial exon duplication;
- exon splitting, exon fusion and splice-boundary shifts;
- alternative-splicing turnover and constitutive-to-alternative exon shifts;
- pseudogenization, frame disruption and copy collapse;
- gene conversion or concerted evolution among close paralogs;
- annotation dropout, fragmented assemblies and unresolved paralogy.

INSIPHY v0.4.0 directly models segment presence, role state, adjacency state,
source mixture and copy multiplicity. It partially captures split/fusion through
adjacency changes and joined-segment candidates. It does not yet explicitly
model transposable-element origin, intron sliding as its own class, gene
conversion among paralogs, isoform-specific alternative-splicing turnover or
tree/topology uncertainty.

References:

- https://academic.oup.com/mbe/article/17/9/1294/994535
- https://pubmed.ncbi.nlm.nih.gov/10958846/
- https://pubmed.ncbi.nlm.nih.gov/22253233/
- https://link.springer.com/article/10.1186/1745-6150-7-11
- https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0004680
- https://rnajournal.cshlp.org/content/13/10/1603

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

## Benchmark Availability

There does not appear to be a mature public benchmark specifically for
single-gene intragenic synteny histories with known branch-level exonization,
source joining, copy expansion and annotation dropout truth. Existing resources
cover adjacent problems: orthology benchmarking, exon orthology tools such as
ExOrthist, curated gene-family studies such as tetraspanins, and well-known
young/chimeric gene cases in Drosophila. INSIPHY therefore needs a composite
benchmark:

- curated positive cases with literature-supported histories, starting with
  `jingwei` and `Sdic`;
- curated negative controls among conserved single-copy genes;
- duplicated-gene case sets with known exon-intron divergence;
- simulation under known histories for power, false positives and branch
  placement accuracy;
- annotation-dropout perturbations of real or simulated annotations;
- ablation baselines: annotation-only, sequence-only, exon-position-only and
  source-label-free models.
