# Real-Data Demonstration

## RpL32 conserved control

The first formal single-copy demonstration uses RpL32 orthologs from:

- Drosophila melanogaster, GCF_000001215.4;
- Drosophila simulans, GCF_016746395.2;
- Drosophila erecta, GCF_003286155.1;
- Drosophila yakuba, GCF_016746365.2;
- Drosophila teissieri, GCF_016746235.2.

Genome FASTA and GFF3 files are stored outside the repository under
`/data/db/genome`. The repository contains accession-level manifests and the
species tree.

## Complete user path

```bash
insiphy build-case \
  --manifest examples/real_cases/rpl32_control/manifest.tsv \
  --species-tree examples/real_cases/rpl32_control/species_tree.tsv \
  --output-dir work/rpl32_control \
  --aligner auto \
  --threads 4

insiphy run \
  --input-dir work/rpl32_control \
  --output-dir results/rpl32_control \
  --analysis-scope single-copy \
  --model er-ard \
  --branch-length-mode unit \
  --threads 4

insiphy visualize \
  --input-dir work/rpl32_control \
  --result-dir results/rpl32_control \
  --output-dir figures/rpl32_control
```

The command above demonstrates the RpL32 case with unit branch lengths. The
bee `dsx` analysis uses the supplied positive branch lengths; both cases were
also fitted with estimated and stationary root frequencies. The formal
Nextflow/Slurm run is recorded under
`/data/projects/intragenic_structure/results/20260916_145605_insiphy` with
INSIPHY v0.12.4, four CPUs per case and color-encoded correspondence figures.
Nextflow schedules the demonstration on this server; it is not a dependency of
the INSIPHY command-line package.

## RpL32: conserved control

The five-species correspondence has conserved core exon sequences, with
unresolved terminal-exon correspondence and one variable splice-junction
pattern among the three junctions included in model fitting. All three fitted
exon-presence patterns have zero informative tip differences. Seven sequence
hits overlapping existing annotated exons remain ambiguous; they are not
called new exons.

For exon presence and exon role, `model_tests.tsv` reports
`test_status=parameters_not_estimable` and `p_value=NA`. The junction model
comparison is sensitive to the root assumption: the estimated-root fit is
not estimable, whereas the stationary-root fit has P=0.2891. A single
variable junction pattern and uncertainty in transcript annotation do not
establish a branch-specific change.

This control does not measure sensitivity to known exon gains or losses.

## Bee dsx: alternative exon organization

Seven bee orthologs from a supplied single-copy orthogroup retain all
annotated transcripts. Twenty-four exon-presence sites enter fitting, with no
informative tip differences; exon gain/loss-rate asymmetry therefore has no
interpretable P value. Five of 21 splice-junction sites have variable tip
patterns. For junctions, the equal-rates (ER) versus asymmetric-rates (ARD)
likelihood-ratio test gives P=0.4524 with estimated root frequency and
P=0.3675 with stationary root frequency. These values provide no evidence
for unequal junction gain and loss rates in this gene and set of taxa.

Three alignments to sequence outside annotated exons are reported as hidden
segment candidates; 167 other queried combinations remain ambiguous. An
aligned genomic sequence alone cannot establish that it is spliced as an
exon. The candidates require examination of splice boundaries and annotation
coverage before any biological claim. The analysis makes no mechanistic or
unique ancestral-event assignment.

These gene-level likelihoods are conditional on the supplied tree, inferred
homologous elements and annotated isoforms. Adjacent exon and junction sites
are correlated; tests pooling them assume conditional independence. Missing
annotations and unknown states limit inference from the apparent absence of
an exon. The figures show correspondence and conditional branch probabilities,
not verified ancestral histories.

## Next real-data set

The next benchmark should contain single-copy genes with independently
documented structural changes. Candidate cases must satisfy:

- one ortholog per species;
- assembly and annotation versions available;
- enough species to distinguish alternative branch placements;
- sequence-level evidence for the altered exon or junction;
- a published history that can be reviewed independently of INSIPHY.

Jingwei and Sdic remain useful for future multi-copy development. Their
duplication histories place them outside the current formal single-copy
benchmark.
