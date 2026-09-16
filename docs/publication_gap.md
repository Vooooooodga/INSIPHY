# Publication Gap

INSIPHY v0.10.0 is a runnable method prototype. It includes transcript-aware
extraction, graph-based segment correspondence, source/copy-role propagation,
branch-length-aware CTMC/Mk fitting, invariant-model LRT p values,
BH q values, parametric-bootstrap calibration, stochastic-map branch summaries,
event support tables, internal correspondence coverage summaries,
colorblind-aware SVG output, optional pattern encoding, copy/gene tree support
for multi-copy structural characters, optional foreground/background
structural-rate tests, progressive EG summaries and simulation operating
characteristic summaries. A manuscript-grade method still requires stronger
evidence in five areas.

## Current Strengths

- The input model matches the project aim: genome FASTA, GFF/GTF annotation,
  a species tree and optional copy/gene tree for multi-copy structural histories.
- The package already separates annotation completion, homologous segment
  correspondence and phylogenetic structural inference.
- The outputs report observable structural changes separately from optional
  interpretation hints. They include likelihood parameters, p values,
  empirical bootstrap p values and branch posterior summaries for structural
  characters.
- Accession-level real-case manifests are available for Drosophila `jingwei`,
  `Sdic` and a conserved-control direction.
- Annotation-dropout negative controls are represented in the simulator and
  benchmark layer.
- The package is CLI-first. Workflow orchestration and cluster execution
  records remain outside the method package.

## Required Before Publication

1. **Real case studies**: complete accession-level `jingwei`, `Sdic` and RpL32
   conserved-control analyses using documented genome and annotation versions,
   then add at least one additional duplicated/chimeric case.
2. **Statistical calibration**: compare Sankoff scores, CTMC likelihoods, LRT
   p values, bootstrap p values, stochastic-map posterior summaries and
   simulated ground truth across many trees, event rates and annotation-error
   settings.
3. **Baseline comparisons**: quantify gains over annotation-only,
   sequence-only and intron/exon-position-only approaches.
4. **Robustness tests**: simulate missing annotation, fragmented gene models,
   tandem duplicates and ambiguous paralogy.
5. **Evidence reporting**: provide manuscript tables for event class, branch,
   support, model score, alternative explanation and input evidence.
6. **Branch histories**: evaluate stochastic character mapping calibration and
   branch placement accuracy under known simulated histories.

## Current Implementation Gap

The package includes simulation, calibration and benchmark commands, bootstrap
summaries, event support summaries and named event scenarios. They should be
expanded to cover many trees, multiple event rates, real accession-level cases
and independent simulation generators before claims about broad accuracy are
made.

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
- high-similarity paralogous segments that may reflect gene conversion, recent
  duplication or unresolved paralogy;
- isoform-specific alternative-splicing turnover;
- processed retrocopy evidence from intron loss and insertion context;
- tree and branch-length uncertainty.
