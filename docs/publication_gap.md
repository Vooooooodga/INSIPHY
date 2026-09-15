# Publication Gap

INSIPHY v0.2.0 is a runnable method prototype. It now includes transcript-aware
extraction, graph-based segment correspondence and branch-length-aware CTMC/Mk
fitting. A manuscript-grade method still requires stronger evidence in four
areas.

## Current Strengths

- The input model matches the project aim: genome FASTA, GFF/GTF annotation and
  a species tree.
- The package already separates annotation completion, homologous segment
  correspondence and phylogenetic structural inference.
- The outputs report candidate events in biological terms, including
  exonization, new adjacency, source joining and copy expansion.
- Small jingwei and Sdic demos run in CI.
- Annotation-dropout negative controls are represented in the simulator and
  benchmark layer.

## Required Before Publication

1. **Real case studies**: complete accession-level jingwei and Sdic analyses
   using documented genome and annotation versions.
2. **Statistical calibration**: compare Sankoff scores, likelihood-like
   character scores and simulated ground truth.
3. **Baseline comparisons**: quantify gains over annotation-only,
   sequence-only and intron/exon-position-only approaches.
4. **Robustness tests**: simulate missing annotation, fragmented gene models,
   tandem duplicates and ambiguous paralogy.
5. **Evidence reporting**: provide manuscript tables for event class, branch,
   support, model score, alternative explanation and input evidence.

## Current Implementation Gap

The package now includes simulation and benchmark commands, including a first
annotation-dropout negative control. They should be expanded to cover many
trees, multiple event rates, real accession-level cases and independent
simulation generators before claims about accuracy are made.
