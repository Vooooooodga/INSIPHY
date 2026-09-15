# Publication Gap

INSIPHY v0.1.0 is a runnable method scaffold. It is useful for demonstrating
the method logic, but a manuscript-grade method requires stronger evidence in
four areas.

## Current Strengths

- The input model matches the project aim: genome FASTA, GFF/GTF annotation and
  a species tree.
- The package already separates annotation completion, homologous segment
  correspondence and phylogenetic structural inference.
- The outputs report candidate events in biological terms, including
  exonization, new adjacency, source joining and copy expansion.
- Small jingwei and Sdic demos run in CI.

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

The package now includes simulation and benchmark commands, but these are
minimal. They should be expanded to cover many trees, multiple event rates and
controlled annotation dropout before claims about accuracy are made.
