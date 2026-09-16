# Algorithm Engineering

## Architecture

INSIPHY is a Python package with explicit modules for extraction, alignment,
annotation completion, correspondence, structural-site construction,
phylogenetic likelihood, and visualization. CLI commands call the same library
functions exposed to Python users. The package does not embed a workflow engine
or scheduler.

Shared biological observations live in `structural_sites.py`. Default
single-copy event placement lives in `parsimony.py`; optional probability
models remain in `structural_phylogeny.py`. The older
multi-copy implementation remains isolated behind
`--analysis-scope experimental-multicopy`.

## Alignment

The internal backend uses Biopython `PairwiseAligner` with affine gaps.
`minimap2` and `miniprot` are optional external backends. All backends return
the same `AlignmentStats` interface: identity, coverage, score, coordinates,
CIGAR, strand, supported paired blocks and backend. Protein and nucleotide
identity retain explicit different units.

Candidate generation groups by gene family and excludes intron-intron pairs.
Short split-exon fragments are not removed by a fixed whole-exon length ratio.
Exon/exon correspondence uses explicit overlap projections. Exon/non-exonic
comparisons use a separately selected local mapper, avoiding whole-intron
MAFFT alignments. The internal
backend calls the declared Biopython dependency; no hand-written dynamic
programming fallback or silent replacement of a failed external backend is
used.

For (n) interval occurrences, unrestricted pair generation is (O(n^2)).
Species, copy and role restrictions reduce the number of aligned pairs.
External aligners remain preferable for long genomic intervals.

## Correspondence graph

Accepted pairwise matches form a graph for tree-ordered constrained merging.
Direct matches and disjoint projections onto a common reference establish
merge compatibility. Split fragments need not resemble one another; a chain
of similarity without these coordinates is insufficient. Graph storage is
(O(V+E)); the accepted graph can still be dense in the worst case.

Public `EG_*` identifiers are assigned only to components containing an
annotated exon-like member or sequence-supported candidate exonic source.
Context introns remain available for boundary calculations without becoming
displayed homologous blocks.

## Parsimony and likelihood

For binary states, inside/outside parsimony messages give globally compatible
node states and branch endpoint pairs without enumerating all tied histories.
With U distinct observation/missingness patterns and N nodes, message
calculation takes O(UN); producing individual site outputs additionally costs
O(SN). No local node tie-breaking is used for branch calls.

For (S) structural sites, (N) tree nodes, and two states, one likelihood
evaluation is O(UN) after exact observation-pattern compression. Repeated
patterns retain multiplicity weights. Binary transition probabilities use
closed forms with stable exponential subtraction. Structural-site independence
is still a biological assumption; compression does not establish it.

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
exact matrix-reward calculations. These calculations run after one family-layer
model has been selected, avoiding repeated parameter fitting per site.

## Data and memory discipline

- FASTA extraction uses existing uncompressed FASTA indexes when available;
  otherwise interval reads stream the source without retaining a chromosome.
- Original annotation bounds and expanded search bounds are recorded separately.
- GFF gene hierarchies use bounded in-memory caching; no reference files change.
- Pairwise diagnostics are streamed; accepted relations remain in memory for
  correspondence inference. Worst-case dense accepted graphs remain quadratic.
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

The current repair is exercised through the explicitly requested real-data
demo workflow. Simulation, checksum and additional smoke-test runs are outside
this execution. Existing unit tests have not been rerun for v0.13.
