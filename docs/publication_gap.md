# Publication Gap

INSIPHY v0.4.0 is a runnable method prototype. It now includes transcript-aware
extraction, graph-based segment correspondence, source/copy-role propagation,
branch-length-aware CTMC/Mk fitting and invariant-model LRT p values. A
manuscript-grade method still requires stronger evidence in five areas.

## Current Strengths

- The input model matches the project aim: genome FASTA, GFF/GTF annotation and
  a species tree.
- The package already separates annotation completion, homologous segment
  correspondence and phylogenetic structural inference.
- The outputs report candidate events in biological terms, including
  exonization, new adjacency, source joining and copy expansion, and they now
  include likelihood parameters and p values for structural characters.
- Small jingwei and Sdic demos run as package fixtures, and the accession-level
  NCBI Gene locus demos now recover source-mixture evidence for both cases.
- Annotation-dropout negative controls are represented in the simulator and
  benchmark layer.
- The package is CLI-first and does not require Nextflow. Nextflow/Slurm runs
  are server-side formal execution records for internal real-case demos.

## Required Before Publication

1. **Real case studies**: complete accession-level `jingwei` and `Sdic`
   analyses using documented genome and annotation versions, then add at least
   one conserved negative control and one additional duplicated/chimeric case.
2. **Statistical calibration**: compare Sankoff scores, CTMC likelihoods, LRT
   p values and simulated ground truth; add parametric bootstrap for LRT P
   values.
3. **Baseline comparisons**: quantify gains over annotation-only,
   sequence-only and intron/exon-position-only approaches.
4. **Robustness tests**: simulate missing annotation, fragmented gene models,
   tandem duplicates and ambiguous paralogy.
5. **Evidence reporting**: provide manuscript tables for event class, branch,
   support, model score, alternative explanation and input evidence.
6. **Branch histories**: add stochastic character mapping or a clearly labeled
   approximation, with posterior probability of event placement on each branch.

## Current Implementation Gap

The package now includes simulation and benchmark commands, including a first
annotation-dropout negative control. They should be expanded to cover many
trees, multiple event rates, foreground/background branch scenarios, real
accession-level cases and independent simulation generators before claims about
accuracy are made.

## Benchmark Plan

No single public benchmark appears to cover this exact task. The publication
benchmark should therefore combine:

- curated literature-positive cases: `jingwei`, `Sdic`, and one or more young
  Drosophila chimeric/duplicate genes;
- conserved negative controls: single-copy genes with stable exon-intron
  organization across close species;
- duplicated-gene structure-divergence cases from published exon-intron
  divergence studies;
- simulated truth sets for exonization, source joining, exon split/fusion,
  retrocopy-like intron loss, tandem duplication and annotation dropout;
- ablations that remove annotation completion, source labels, adjacency
  evidence or phylogenetic modeling.

Primary metrics:

- event-level precision, recall and F1;
- branch placement accuracy;
- false positive rate under conserved negative controls;
- robustness to missing exons and fragmented annotation;
- calibration of reported P values and branch posterior probabilities.

## Biological Gaps To Track

The current event vocabulary should be extended or explicitly labeled for:

- splice-boundary shifts and intron sliding;
- tandem exon duplication and partial exon duplication;
- transposable-element-associated exonization;
- gene conversion among close paralogs;
- isoform-specific alternative-splicing turnover;
- processed retrocopy evidence from intron loss and insertion context;
- tree and branch-length uncertainty.
