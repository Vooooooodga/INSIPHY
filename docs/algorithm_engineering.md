# Algorithm Engineering

## Architecture

INSIPHY is a Python package with explicit modules for extraction, alignment,
annotation completion, correspondence, structural-site construction,
phylogenetic likelihood, and visualization. CLI commands call the same library
functions exposed to Python users. The package does not embed a workflow engine
or scheduler.

Formal single-copy statistics live in `structural_phylogeny.py`. The older
multi-copy implementation remains isolated behind
`--analysis-scope experimental-multicopy`.

## Alignment

The internal backend uses Biopython `PairwiseAligner` with affine gaps.
`minimap2` and `miniprot` are optional external backends. All backends return
the same `AlignmentStats` interface: identity, coverage, score, coordinates,
CIGAR, and backend.

Candidate generation applies length and role compatibility before dynamic
programming. Large internal comparisons use a bounded fallback to prevent a
single long interval pair from allocating an unbounded dynamic-programming
matrix.

For (n) interval occurrences, unrestricted pair generation is (O(n^2)).
Species, copy, length, and role filters reduce the number of aligned pairs.
External aligners remain preferable for long genomic intervals.

## Correspondence graph

Pairwise matches form a sparse graph. Connected components provide internal
homology candidates, while reciprocal-best status, local order, splice
boundaries, and one-copy constraints limit transitive over-merging. Graph
storage is (O(V+E)).

Public `EG_*` identifiers are assigned only to components containing an
annotated exon-like member or sequence-supported candidate exonic source.
Context introns remain available for boundary calculations without becoming
displayed homologous blocks.

## Likelihood

For (S) structural sites, (N) tree nodes, and two states, one likelihood
evaluation is (O(SN)). The two-state transition matrices are cached by rate,
branch length, and foreground multiplier.

Scaled pruning stores one length-two vector per node. Site processing is
sequential in memory, so the core likelihood memory cost is (O(N)) plus the
input site matrix.

Optimization uses log-transformed positive parameters, bounded L-BFGS-B, and
multiple deterministic starts. `--threads` parallelizes independent optimizer
starts. Profile-likelihood intervals remain serial because each profile point
depends on a constrained nested optimization.

## Posterior calculation

Inside-outside messages compute node marginals and parent-child joint
posteriors in (O(N)) per site. Expected directional transition counts use
conditional CTMC path integrals. These calculations run after one family-layer
model has been selected, avoiding repeated parameter fitting per site.

## Data and memory discipline

- FASTA/GFF extraction operates one selected gene at a time.
- Large all-genome indexes and databases remain external to the package.
- Result tables are plain TSV/JSON and can be processed incrementally.
- Unknown observations use likelihood marginalization and do not create
  imputed sequence states.
- Output directories contain only declared result tables and figures.

## Code policy

The implementation favors small typed concepts and direct library calls:

- one representation for tree rows and branch lengths;
- one alignment result type across backends;
- one structural-site matrix feeding every formal model;
- one likelihood engine for ER, ARD, and foreground variants;
- explicit errors for violated biological assumptions;
- no automatic biological mechanism labels;
- no silent correction of malformed single-copy input;
- no hidden thresholds for statistical significance or branch support.

Validation belongs at biological interfaces: copy count, species-tree labels,
positive branch lengths, foreground edges, and ascertainment assumptions.
Internal helpers rely on those established invariants instead of repeating
defensive checks at every line.

## Reproducibility

`run_parameters.json` records analysis scope, model, branch-length mode,
ascertainment, foreground file, thread count, and tree path.
`excluded_families.tsv` records copy-count violations.
`model_fits.tsv` records optimizer starts, convergence messages, parameter
boundaries, and confidence intervals.

The test suite covers extraction, alignment, single-copy state coding, ER/ARD,
foreground rates, node and branch posteriors, visualization, OrthoFinder
import, and the retained experimental multi-copy path.
