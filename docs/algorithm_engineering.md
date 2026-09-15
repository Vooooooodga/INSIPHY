# Algorithm And Engineering Notes

INSIPHY keeps the method package focused on single genes or supplied homologous
gene-copy sets. Upstream tools such as OrthoFinder, OMA, Broccoli or curated
case manifests define the gene/copy set. INSIPHY then solves three internal
problems: complete missing structure from genome sequence, infer homologous
segments inside the supplied set, and reconstruct structural histories on a
fixed species tree.

## Alignment Strategy

The package has a small internal dynamic-programming fallback so tests and toy
demos run without external software. Publication-scale runs should prefer
mature aligners for the expensive sequence layer:

- `minimap2`: nucleotide segment and local genomic interval matching. It uses
  minimizer seeding, chaining and base-level alignment, and provides CIGAR-like
  PAF tags useful for structural evidence.
- `miniprot`: protein-to-genome alignment with splicing and frameshift states.
  It is useful when a protein sequence from an annotated source copy is mapped
  back to a target genome interval to find hidden or shifted coding segments.
- GenePainter/CESAR concepts: intron phase, intron position and coding-frame
  constraints are treated as structural evidence, even when their code is not
  embedded in INSIPHY.

The current CLI exposes `--aligner internal|minimap2|miniprot`. External tools
are optional and discovered on `PATH`; their source code is not vendored.

## Candidate Filtering

A naive all-by-all segment comparison has O(n²) candidate pairs per family, and
each exact dynamic-programming alignment can cost O(L1 x L2) time and memory.
INSIPHY v0.5 therefore applies cheap filters before expensive alignment:

- same family;
- different gene copy for HSG graph edges;
- minimum length ratio;
- available sequence;
- role compatibility that still keeps intron/exon role-shift candidates.

This design keeps likely exonization and intronization cases in the candidate
set while reducing alignments among impossible pairs. Pair scoring can run with
`--threads`, and external aligners are called with one thread per pair to avoid
oversubscribing CPU cores.

## Phylogenetic Algorithms

INSIPHY uses discrete structural characters:

- segment presence;
- segment role;
- intragenic adjacency;
- source mixture;
- copy multiplicity.

For each character it runs:

- weighted Sankoff parsimony for low-cost ancestral states;
- Felsenstein pruning under a branch-length-aware CTMC/Mk model;
- likelihood-ratio testing against an invariant no-change model;
- optional parametric bootstrap under the invariant null;
- optional stochastic character mapping by CTMC uniformization;
- optional foreground/background CTMC rate comparison.

The stochastic mapping implementation follows the standard SIMMAP/phytools
logic: fit a CTMC rate matrix, condition on observed tips, sample endpoint
states on each branch, and sample complete CTMC paths along branches. The bridge
sampler uses uniformization, which is appropriate for the small state spaces in
intragenic structural characters.

## Memory And Runtime Defaults

The default package mode stays lightweight:

- `--aligner internal`;
- `--threads 1`;
- `--bootstrap-replicates 0`;
- `--stochastic-maps 0`.

Recommended real-case settings for a small gene family are:

```bash
insiphy run \
  --input-dir case_dir \
  --output-dir result_dir \
  --bootstrap-replicates 200 \
  --stochastic-maps 200 \
  --seed 7
```

For larger families, increase `--threads` during `derive-tables` or
`build-case`, inspect `alignment_backend_report.tsv`, and record runtime,
memory and candidate-filtering ratios in the benchmark notes.

## Reference Code And Methods

Key methods reviewed for v0.5:

- ExOrthist: exon orthology from supplied gene orthogroups using intron
  position/phase, exon sequence and flanking exon context.
- GenePainter: mapping intron positions and phases onto protein alignments.
- CESAR/TOGA: coding-exon-aware realignment and genome-scale annotation
  transfer ideas.
- miniprot/minimap2/minisplice: modern sequence alignment, spliced alignment
  and splice-site scoring backends.
- SIMMAP/phytools: stochastic character mapping for discrete traits on
  phylogenies.

The implementation borrows algorithmic ideas and output semantics, while the
INSIPHY codebase remains an independent Python package.
