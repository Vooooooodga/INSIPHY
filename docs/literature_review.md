# Biological and Statistical Rationale

## Conclusions adopted in v0.12

The literature supports three separable tasks:

1. establish homology among exon-like sequence units inside an upstream-defined
   ortholog set;
2. code observable sequence presence, exonic role, and splice-junction states;
3. fit explicit evolutionary models on a fixed phylogeny and retain
   uncertainty in ancestral histories.

This separation prevents annotation absence from being equated with biological
loss and prevents a structural pattern from being promoted automatically to a
molecular mechanism.

## Exon and intron evolution

Gene structures can change through exon gain/loss, exonization and
intronization, intron gain/loss, splice-site movement, exon duplication,
shuffling, splitting, and fusion. Several mechanisms can yield the same
terminal exon-intron pattern. Comparative sequence and local structure can
identify the changed unit, while transposon origin, selection, gene conversion,
or expression consequences usually require additional data.

The v0.12 layers reflect this:

- sequence presence addresses gain/loss of a homologous unit;
- exonic role addresses recruitment or loss of an exon role while sequence
  remains present;
- splice-junction presence addresses segmentation changes compatible with
  split/fusion histories.

Intronic sequence enters exon homology only when it is homologous to an
exon-like unit in another species or receives sequence-supported exon
completion evidence. Routine introns remain boundary and spacing context.

Relevant reviews and models include:

- Keren, Lev-Maor and Ast, 2010, *Alternative splicing and evolution:
  diversification, exon definition and function*.
- Irimia and Roy, 2014, *Origin of spliceosomal introns and alternative
  splicing*.
- Carmel et al., 2007, probabilistic reconstruction of intron gain and loss:
  https://doi.org/10.1186/1471-2148-7-192
- Csűrös, Rogozin and Koonin, 2011, comparison of likelihood, MCMC, and Dollo
  approaches to ancestral intron reconstruction:
  https://doi.org/10.1371/journal.pcbi.1002150

## Internal correspondence

Exon correspondence cannot rely on annotation labels alone. Useful evidence
includes sequence similarity, aligned coverage, exon phase, splice motifs,
flanking exon context, strand, order, and agreement among close relatives.
One-to-many and many-to-one alignments must remain possible because one
ancestral sequence interval can be partitioned into multiple descendant
exons, or multiple ancestral intervals can become one exon.

Progressive Cactus demonstrates the general value of ordering difficult
homology comparisons by a guide phylogeny. INSIPHY applies that principle
locally within supplied orthologous genes. It does not run whole-genome Cactus
and does not infer new gene families.

## Phylogenetic likelihood

Felsenstein's pruning algorithm supplies the likelihood foundation: sum over
unobserved internal states while conditioning on a tree and a transition
model. Pagel's discrete-character framework and Lewis's Mk/Mkv model provide
direct precedents for finite-state CTMC analysis and variable-site
ascertainment.

The PAML manual contributes several operational principles used here:

- parameters are shared across many sites rather than fitted independently to
  every observed pattern;
- nested hypotheses are defined by parameter constraints;
- likelihood-ratio tests compare optimized log likelihoods;
- branch and site models require a priori model specification;
- ancestral reconstructions are conditional on the fitted model and tree;
- boundary estimates and weak information require cautious inference.

PAML itself models nucleotide, codon, or amino-acid substitution. INSIPHY uses
the same likelihood discipline for binary homologous gene-structure sites.

Core references:

- Felsenstein, 1981, *Evolutionary trees from DNA sequences: a maximum
  likelihood approach*.
- Pagel, 1994, *Detecting correlated evolution on phylogenies*.
- Lewis, 2001, *A likelihood approach to estimating phylogeny from discrete
  morphological character data*.
- Yang, 2007, *PAML 4: phylogenetic analysis by maximum likelihood*.
- PAML documentation and example control files:
  https://github.com/abacus-gene/paml

## Genome synteny concepts translated to one gene

Whole-genome synteny methods ask whether homologous units retain order,
orientation, adjacency, and ancestral linkage across a tree. At gene scale,
the homologous units are exon-like sequence intervals and splice junctions.
The analogous questions are:

- which exon-like units are conserved;
- whether their order and orientation remain stable;
- where a unit appears, disappears, or changes exonic role;
- where a splice junction appears or disappears;
- which branches have a higher structural transition rate.

Whole-genome microsynteny has also been used directly as phylogenetic
information, as in *Whole-genome microsynteny-based phylogeny of angiosperms*
(Zhao et al., 2021):
https://doi.org/10.1038/s41467-021-23665-0

The linked-block visualization used for conserved chromosome synteny provides
the visual grammar for INSIPHY exon tracks. The tree and homologous connections
are shown together, while posterior structural changes are placed on branches.

Recent synteny software such as Synolog emphasizes scalable orthology,
multi-genome synteny clusters, and linked visual outputs:
https://doi.org/10.64898/2026.04.07.717040
INSIPHY begins after gene orthology has been supplied and operates at the
within-gene structural scale.

## Statistical boundary

The formal model estimates observable structural transitions. It does not
infer:

- causative transposable elements;
- adaptive selection;
- gene conversion;
- expression or isoform abundance;
- pathway consequences;
- gene-level orthology;
- species-tree uncertainty.

These questions can use INSIPHY outputs as structured evidence in a broader
study.

## Empirical benchmark strategy

A single public benchmark with known branch-level histories for homologous
exons is currently unavailable. Evaluation should therefore use:

- deeply curated single-copy genes with established exon/intron histories;
- conserved single-copy controls;
- ortholog sets from multiple clades and annotation releases;
- manual sequence-level review of every proposed gain/loss or split/fusion;
- sensitivity to transcript choice, aligner, branch lengths, and
  correspondence thresholds;
- comparison with annotation-only and sequence-only correspondence.

The current RpL32 control establishes end-to-end execution and conservative
handling of invariant structure. It does not establish general sensitivity or
specificity.
