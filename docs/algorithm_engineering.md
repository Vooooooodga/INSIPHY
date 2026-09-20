# Algorithm Engineering

## Architecture

IntraPhy is a Python package with library modules and CLI commands over the same implementation. The package analyzes prepared gene sets; it does not embed a scheduler or workflow engine.

Core modules:

- extraction and case building;
- annotation evidence and bounded locus search;
- alignment adapters;
- correspondence graph construction;
- structural-site matrix construction;
- parsimony;
- optional CTMC likelihood;
- visualization.

The formal path is single-copy. Multi-copy code remains isolated under `experimental-multicopy`.

## One validation boundary

The code should validate biological assumptions at interface boundaries:

- every formal family has one gene per species;
- IDs map exactly to gene/transcript/protein annotation records;
- the tree is rooted, connected, acyclic, and matches species labels;
- foreground branches exist in the tree;
- branch lengths are finite and compatible with the selected mode;
- ascertainment settings match the available candidate universe.

After those checks, internal functions rely on clear invariants. Repeating broad defensive guards inside every helper makes the model harder to read and can hide invalid biological input. Invalid assumptions should fail at the boundary with a concrete message.

## Alignment adapters

IntraPhy should use mature tools for mature alignment problems:

- MAFFT for exon-level multiple/pair projection where appropriate;
- minimap2 for genomic local context and long nucleotide intervals;
- miniprot for protein-to-genome coding projections;
- Biopython `PairwiseAligner` for bounded short intervals.

Each adapter returns a common alignment-statistics record with identity, coverage, strand, coordinates, paired blocks and backend identity. If a requested backend fails or is unavailable, the result is reported as unavailable for that evidence class. Hidden replacement by another algorithm is avoided.

IntraPhy does not try to become a new general-purpose alignment package. The method contribution is the structural interpretation and phylogenetic use of the evidence.

## Correspondence graph

For `n` intervals, unrestricted pair generation is `O(n^2)`. Candidate filters by family, role, size, species distance and local context reduce practical work. Accepted relations form a graph with `O(V + E)` storage, where `E` can still be dense in difficult genes.

Tree-guided merging proceeds from close to distant comparisons. Merges require sequence or projection compatibility, not mere transitive similarity. Split/fusion evidence uses ordered coordinate projections. Repeated hits to the same reference region remain repeated or ambiguous evidence.

## Structural-site construction

The structural matrix is the single source for formal models. Sequence presence, exonic role and splice junctions are separate layers. Predicted roles remain predicted evidence. Candidate source intervals can support sequence presence or role uncertainty without becoming confirmed exon homology.

This design keeps statistical meaning stable across commands and figures.

## Parsimony complexity

For `U` distinct observed patterns and `N` tree nodes, binary parsimony messages cost `O(U N)`. Producing per-site outputs costs `O(S N)` after pattern compression. The implementation keeps all globally optimal node and branch endpoint states, avoiding local tie-breaking.

## Likelihood complexity

For `U` compressed patterns and `N` nodes, one binary likelihood evaluation costs `O(U N)`. `_inside_log_messages` stores length-two log-likelihood vectors per node. Binary child-state sums use `np.logaddexp`; child contributions are added in log space, and the root prior is combined by the same stable log-sum-exp operation.

The posterior pass also propagates log-space messages. Prefix and suffix sums of child messages provide each child's sibling contribution in constant time after linear preparation at its parent. Thus a pattern's inside/outside passes cost `O(E)` over `E` tree edges, including multifurcations, with `O(N)` working message storage. Parameter optimization uses deterministic starts and bounded L-BFGS-B.

Independent optimizer starts can use multiple threads. Profile-likelihood scans are more expensive because every profile point requires a constrained fit.

## Memory and I/O

- FASTA intervals are read from existing sources without rewriting reference files.
- Search bounds and original gene bounds are recorded separately.
- Pairwise diagnostics can be streamed.
- Accepted correspondence relations stay in memory for graph construction.
- Result files are TSV/JSON/SVG and can be inspected with ordinary tools.
- Output directories should contain declared outputs only.

## Code policy

The project follows these concrete policies:

- one structural-site matrix feeds all formal models;
- one tree representation is used by parsimony and likelihood;
- one alignment-statistics interface wraps external tools;
- mature libraries and external aligners are preferred over new local algorithms;
- no silent fallback when a requested backend fails;
- no arbitrary threshold retuning to make a demo look cleaner;
- no automatic mechanism labels in core statistical outputs;
- no invented probabilities for invalid fitted models;
- no replacement of missing biological evidence by absence.

Thresholds should be named parameters or documented defaults. Changing them for a benchmark or manuscript figure requires a recorded reason and rerun.

## Visualization engineering

Visualization is a consumer of biological state, not a second inference engine. It must preserve:

- confirmed exon homology ribbons only for confirmed exon-like members;
- separate visual style for predicted, candidate and unknown states;
- separate transcript lanes for alternative paths;
- colorblind-friendly defaults and optional pattern encoding;
- valid probability symbols only when the selected fitted posterior is valid.

## Reproducibility status

Final v0.14.0 run `20260918_121100_insiphy` passed 136 selected formal regression tests under Slurm `61625` ([report](/data/projects/intragenic_structure/results/20260918_121100_insiphy/regression_tests.txt)), including four new role-conflict tests. All ten re-inference/visualization tasks and five comparison tasks completed using unchanged cases/evidence; all six core tables agree in each case pair. Audits establish recovered Hdac3 structure, partial rec8 recovery, unresolved spo5, incomplete RpL32 coverage and two final dsx role contrasts. Final dsx fits emit no node/branch posteriors. The 132-test run is historical. SVG semantic and selected endpoint checks are complete; no rendered inspection was performed. Broader scientific limitations remain in the [real-data benchmark](real_data_benchmark.md).
