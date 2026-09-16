# Publication Readiness Review

Date: 2026-09-16

## Conclusion

INSIPHY v0.9.1 has the correct core frame for a first method paper: upstream
gene/copy homology is supplied by tools or curation, INSIPHY infers
gene-internal exon-like structural correspondence from genome sequence and annotation, and
then analyzes structural characters on the relevant phylogenetic tree. Single-copy
families can use the species tree; multi-copy structural characters use a
copy/gene tree when provided, while copy multiplicity remains species-level.
The remaining gap is evidence strength and benchmark scale. The method needs
larger real-case evaluation, simulation-based statistical calibration and
clearer branch-history posteriors before manuscript-level claims.

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
- many-to-many segment relationships through the internal correspondence graph.

Missing or under-modeled:

- intron sliding as a separate event class;
- tandem exon duplication and partial exon duplication as explicit event types;
- transposable-element origin of exonized sequence;
- high-similarity paralogous segments with mechanism ambiguity;
- isoform-specific alternative-splicing turnover;
- uncertainty in species tree, copy/gene tree topology and branch lengths;
- larger conserved-control sets for user-defined evolutionary hypotheses.

## Statistical Status

Implemented:

- Sankoff reconstruction for structural characters;
- branch-length-aware one-rate CTMC/Mk likelihood;
- invariant-model LRT with boundary-rate p-value approximation;
- AIC/BIC and fitted transition rate;
- endpoint posterior probabilities for parent-child state changes on branches;
- simulated event benchmark summaries.

Needed before publication:

- larger calibration sets for empirical LRT P values;
- calibration plots for stochastic character mapping branch histories;
- posterior/event summaries across correlated structural characters;
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

Cluster-scale execution records can be maintained outside the distributed
package. The user-facing method definition remains the CLI/Python API above.

## Immediate Next Work

1. Run accession-level `jingwei`, `Sdic` and RpL32 control analyses end to end.
2. Expand the simulator to generate event-specific truth sets across many
   trees, rates and annotation dropout patterns.
3. Add real-case manifests for at least one more curated duplicated/chimeric
   gene and one conserved negative control.
4. Write manuscript-ready tables that connect every event call to sequence
   evidence, structural evidence, phylogenetic support and alternative
   explanations.
