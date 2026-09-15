# Publication Readiness Review

Date: 2026-09-15

## Conclusion

INSIPHY v0.4.0 has the correct core frame for a first method paper: upstream
gene/copy homology is supplied by tools or curation, INSIPHY infers
gene-internal homologous segment groups from genome sequence and annotation, and
then analyzes structural characters on a species tree. The remaining gap is
evidence strength, not the project concept. The method needs larger real-case
evaluation, simulation-based statistical calibration and clearer branch-history
posteriors before manuscript-level claims.

## Biological Coverage Check

Covered directly:

- conserved or missing homologous intragenic segments;
- role shifts among CDS, exon/UTR, intron/noncoding and absence;
- intragenic adjacency gain/loss as segment fusion, split or rearrangement
  candidates;
- source mixture and source joining for chimeric genes;
- copy multiplicity changes for duplicated genes;
- annotation dropout and hidden sequence-supported segments;
- source/background/derived copy roles in curated or upstream-defined cases.

Partially covered:

- exon split/fusion through adjacency changes and joined-segment evidence;
- processed-copy or retrocopy-like events through copy relationship classes and
  intron-loss-like patterns;
- splice-boundary shifts through hidden-segment and phase/splice evidence;
- many-to-many segment relationships through HSG graph structure.

Missing or under-modeled:

- intron sliding as a separate event class;
- tandem exon duplication and partial exon duplication as explicit event types;
- transposable-element origin of exonized sequence;
- gene conversion or concerted evolution among close paralogs;
- isoform-specific alternative-splicing turnover;
- uncertainty in species tree topology and branch lengths;
- foreground/background branch-rate tests for user-defined evolutionary
  hypotheses.

## Statistical Status

Implemented:

- Sankoff reconstruction for structural characters;
- branch-length-aware one-rate CTMC/Mk likelihood;
- invariant-model LRT with boundary-rate p-value approximation;
- AIC/BIC and fitted transition rate;
- endpoint posterior probabilities for parent-child state changes on branches;
- simulated event benchmark summaries.

Needed before publication:

- parametric bootstrap for empirical LRT P values;
- stochastic character mapping or a clearly labeled approximation for branch
  histories;
- foreground/background two-rate CTMC tests;
- posterior/event summaries aggregated across correlated characters;
- calibration plots showing false positive rate and power under simulation;
- sensitivity to annotation dropout, branch lengths, thresholds and paralogy
  ambiguity.

## Benchmark Direction

A ready-made benchmark for branch-level intragenic synteny histories is unlikely
to exist. INSIPHY should use a composite benchmark:

- curated positives: `jingwei`, `Sdic`, and additional young duplicated or
  chimeric Drosophila genes;
- conserved negatives: close-species single-copy genes with stable exon-intron
  structure;
- published duplicate-gene exon-intron divergence cases;
- simulated truth sets for each event class;
- ablation baselines that remove sequence, annotation, adjacency, source labels
  or phylogeny.

Benchmark metrics should include event precision/recall/F1, branch placement
accuracy, calibrated P-value behavior, posterior calibration and robustness to
annotation incompleteness.

## Software Boundary

The user-facing software should remain a modern CLI/Python package:

```bash
insiphy build-case --manifest manifest.tsv --species-tree species_tree.tsv --output-dir case
insiphy run --input-dir case --output-dir results
insiphy benchmark --input-dir case --output-dir results
```

Nextflow and Slurm are appropriate for server-side formal runs and large
project records on R730. They should stay outside the distributed package and
outside the user-facing method definition.

## Immediate Next Work

1. Add stochastic-map or bootstrap-calibrated branch event support.
2. Expand the simulator to generate event-specific truth sets across many
   trees, rates and annotation dropout patterns.
3. Add real-case manifests for at least one more curated duplicated/chimeric
   gene and one conserved negative control.
4. Write manuscript-ready tables that connect every event call to sequence
   evidence, structural evidence, phylogenetic support and alternative
   explanations.
