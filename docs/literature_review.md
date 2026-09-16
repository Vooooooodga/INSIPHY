# Biological And Statistical Rationale

INSIPHY is positioned at the intersection of exon evolution, annotation
completion, gene-structure comparison and discrete-character phylogenetics.

## Review Conclusion

The method direction is biologically defensible: gene-internal segment homology
should be inferred from genome sequence plus annotation context, and the
resulting structural characters should be interpreted on a species tree. The
main remaining biological risk is event ambiguity. Exonization, intron
gain/loss, splice-boundary drift, tandem exon duplication, gene conversion,
processed-copy insertion and assembly/annotation artifacts can create
overlapping evidence patterns in one observed gene copy. INSIPHY reports these
as candidate structural histories with explicit alternative explanations.

## Gene-Internal Homology Evidence

The following evidence classes are important for homologous internal segment
assignment:

- local nucleotide or protein-to-genome sequence similarity;
- alignment coverage and segment length ratio;
- flanking segment context;
- intron position and intron phase;
- splice donor/acceptor motif compatibility;
- strand and frame status;
- local order inside each gene copy;
- copy role and source label when the case is a duplicated or chimeric gene.

INSIPHY keeps the upstream homology step separate. OrthoFinder or an equivalent
workflow supplies the gene/copy set; INSIPHY then evaluates internal structure
inside that set.

## Annotation Completion

Genome annotation can miss short exons, shifted splice sites, noncanonical
transcripts or fragmented gene models. INSIPHY therefore separates observed
annotation from sequence-supported evidence. Hidden-segment candidates,
shifted-splice-site candidates and joined-segment candidates are retained as
evidence classes and then interpreted on the species tree.

## Phylogenetic Structural Inference

INSIPHY treats intragenic structural states as discrete phylogenetic
characters. The likelihood calculation follows the Felsenstein pruning logic on
a fixed species tree. Structural states evolve under Mk-style CTMC models, and
branch histories are summarized with stochastic character mapping when
requested.

The implemented statistics answer different questions:

- Sankoff reconstruction asks which ancestral structural states minimize event
  cost.
- CTMC/Mk likelihood estimates a structural-change rate on the tree.
- The invariant-model LRT asks whether a character is better explained by
  structural change than by no change.
- Bootstrap empirical p values calibrate the LRT on small trees.
- Stochastic mapping estimates event placement along branches.
- Foreground/background tests ask whether specified branches have elevated
  structural-change rates.

## Event Inventory

The benchmark and event vocabulary should cover:

- gene duplication and post-duplication exon-intron divergence;
- retroposition or processed-copy insertion, often with intron loss;
- chimeric new gene formation by source joining;
- exonization of intronic, transposed-element or intergenic sequence;
- intron gain, intron loss, intron sliding and intronization;
- tandem exon duplication and partial exon duplication;
- exon splitting, exon fusion and splice-boundary shifts;
- alternative-splicing turnover where annotation supports isoform differences;
- pseudogenization, frame disruption and copy collapse;
- high-similarity paralogous segments with mechanism ambiguity;
- annotation dropout, fragmented assemblies and unresolved paralogy.

INSIPHY v0.7.0 directly models segment presence, role state, adjacency state,
source mixture and copy multiplicity. It reports splice-boundary shifts,
segment fusion and TE-associated exonization as candidate structural patterns.
High-identity paralogous segment matches are reported as ambiguous evidence
because intragenic structure alone cannot distinguish gene conversion, recent
duplication and unresolved paralogy.

## Progressive Alignment Idea

Progressive whole-genome alignment uses a species tree to organize alignment
across evolutionary distances. INSIPHY adapts that idea at gene-internal scale:
segment correspondence is interpreted by tree distance within the supplied
homologous gene set. This produces close-species, within-clade and deep-tree
support classes without running whole-genome alignment or searching for new
gene-level homologs.

Cactus remains relevant as a conceptual example of tree-guided progressive
alignment. INSIPHY does not require Cactus output as input.

## Visualization Principles

Gene-structure figures should show exon/intron organization, homologous segment
blocks and event placement on the species tree. INSIPHY therefore generates:

- a gene-internal synteny map by species and copy;
- a species-tree event map with structural-event support summaries.

The figures use texture, labels, line styles and shapes before color, so the
main interpretation remains available to colorblind readers and in black-white
printing.

## Benchmark Availability

A mature public benchmark for branch-level intragenic synteny histories is not
available as a single resource. INSIPHY therefore needs a composite benchmark:

- curated positive cases with literature-supported histories, starting with
  `jingwei` and `Sdic`;
- conserved negative controls such as RpL32;
- duplicated-gene case sets with known exon-intron divergence;
- simulation under known histories for power, false positives and branch
  placement accuracy;
- annotation-dropout perturbations of real or simulated annotations;
- ablation baselines: annotation-only, sequence-only, exon-position-only and
  source-label-free models.
