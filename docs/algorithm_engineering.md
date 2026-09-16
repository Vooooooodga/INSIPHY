# Algorithm And Engineering Notes

INSIPHY is a single-gene or small gene-family method package. Upstream tools
define the gene/copy set. INSIPHY then solves three internal problems:
sequence-supported annotation completion, homologous segment correspondence and
phylogenetic structural inference.

## Alignment Strategy

The package keeps a small internal dynamic-programming fallback for reproducible
local tests. Real accession-level runs should prefer mature local aligners:

- `minimap2`: nucleotide segment and local genomic interval matching.
- `miniprot`: protein-to-genome alignment with splice and frameshift states.

The CLI exposes `--aligner internal|minimap2|miniprot`. External tools are
discovered on `PATH`; their code is not vendored.

INSIPHY does not run whole-genome alignment. The relevant idea from progressive
genome alignment is tree-guided ordering of evidence: segment support from close
species is evaluated first, and deeper support is interpreted in the context of
the supplied species tree. This behavior is represented in
`progressive_correspondence.tsv`.

## Candidate Filtering

A naive all-by-all segment comparison has O(n^2) candidate pairs per family, and
exact dynamic programming can cost O(L1 x L2) time and memory. INSIPHY filters
before alignment:

- same family;
- different gene copy;
- minimum length ratio;
- available sequence;
- role compatibility that keeps intron/exon role-shift candidates.

Pair scoring can run with `--threads`. External aligners are called with one
thread per pair to avoid oversubscribing CPU cores.

## Phylogenetic Algorithms

INSIPHY uses discrete structural characters:

- segment presence;
- segment role;
- intragenic adjacency;
- source mixture;
- copy multiplicity.

For each character it runs weighted Sankoff reconstruction, CTMC/Mk likelihood,
invariant-model LRT, optional parametric bootstrap, optional stochastic
character mapping and optional foreground/background CTMC rate comparison.

The stochastic mapping implementation fits a CTMC rate matrix, conditions on
observed tips, samples endpoint states and samples complete paths by
uniformization. The state spaces are small, which keeps the implementation
tractable for single-gene analyses.

## Code Quality Policy

The implementation favors clear scientific code over defensive scaffolding:

- user-facing input files get explicit required-field checks;
- internal invariants fail directly when broken, so bugs are visible;
- no broad exception swallowing around biological inference;
- no speculative wrappers, duplicated fallback systems or unused abstraction
  layers;
- functions are named after biological or statistical operations;
- comments explain algorithmic intent, not obvious assignments;
- output tables stay flat and inspectable.

This policy is meant to reduce maintenance cost and avoid code that looks
machine-generated: repeated guard clauses, generic helper layers, excessive
normalization and vague error recovery are removed unless they protect a real
user input boundary.

## Memory And Runtime Defaults

The default package mode stays lightweight:

- `--aligner internal`;
- `--threads 1`;
- `--bootstrap-replicates 0`;
- `--stochastic-maps 0`.

Recommended real-case settings for a small gene family:

```bash
insiphy build-case \
  --manifest manifest.tsv \
  --species-tree species_tree.tsv \
  --output-dir case_dir \
  --aligner minimap2 \
  --threads 4

insiphy run \
  --input-dir case_dir \
  --output-dir result_dir \
  --bootstrap-replicates 200 \
  --stochastic-maps 200 \
  --seed 7
```

For larger families, increase `--threads` during table derivation, inspect
`alignment_backend_report.tsv`, and report runtime, memory and candidate counts
in benchmark notes.
