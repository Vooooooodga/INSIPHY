> **V18 legacy documentation / retained reference.** The default V19 model and
> commands are described in [the README](../README.md) and
> [exon_structure_model.md](exon_structure_model.md). Do not interpret the old
> independent-layer analyses as exon configuration inference.

# Validation and remaining scientific work

## File-input and visual-guide completion

Local evidence is archived in
[`validation/v018/file_input_completion`](../validation/v018/file_input_completion/).
The complete suite ran **442 tests: no failures, errors or skips**, with real
MAFFT v7.505 and minimap2 2.26-r1175. The installed wheel completed 11 CLI commands,
including no-manifest preparation, portable locus export, inference and graphics.
Its result gallery contained 27 data-derived target SVGs and five separately
labelled teaching SVGs. Five additional standalone guide SVGs were generated.

The wheel was imported from a separate installation and run with `-I` outside
the checkout. Runtime dependencies were reused from the existing environment;
this was not a fresh dependency-resolution test. The first wheel attempt lacked
Biopython, was corrected by supplying that dependency path, and is retained in
the evidence directory. Remote CI results must be checked independently.

The AGAT tests mock subprocess execution to verify command/failure contracts.
A real AGAT executable was not run. The guide's small parsimony examples are
checked against the algorithm, but are not biological or statistical validation.
No scalability claim is added by this completion.

## Distinct evidence levels

Unit tests verify stated contracts and numerical identities. Integration tests
exercise real external programs. Synthetic raw-input tests start from FASTA/GFF3
without supplied homology or character answers. Prepared historical examples
instead test downstream behavior conditional on provided answers. These categories
must be reported separately.

The v18 release adds tests for character identity, coverage-based dependence,
unknown/inapplicable observations, conservative likelihood availability,
minimum-change witnesses, safe output handling and raw-input splice differences.
Both positive evidence and evidence-insufficient controls are retained.

Current release results and actual versions are recorded under `validation/v018/`.
Historical v16/v17 validation and demo directories are unchanged provenance, not
new evidence. Missing optional programs must be reported as skipped or untested;
they must never be replaced by fake executables to claim integration coverage.

## What passing tests does not establish

V18 does not claim a benchmarked biological accuracy across taxa or all 28 scope
cases. It does not claim finite-sample calibration of LRTs, interval coverage,
posterior calibration, or independence of every apparently unlinked character.
The raw four-taxon example is an integration regression, not a simulation study.

## Publication validation still required

Use independently designed raw-sequence/annotation simulations and curated
biological cases to measure correspondence errors, false absence, positional
error, erroneous required changes, direction error, recall and abstention rate.
Separate true structural evolution from annotation deletion, assembly gaps,
repeats and alignment uncertainty. Remove truth sidecars before analysis.

For rate-model claims, measure null false-positive rate, power, interval coverage,
parameter bias and unavailable-fit frequency at realistic character counts, trees
and missingness patterns. Include state-dependent annotation failure and linked
characters rather than assuming a fixed mask fully describes discovery.

Compare methods only on matching inferential targets. An exon-orthology tool,
a coding-annotation projector and an intron-history model solve different parts
of the problem. Sensitivity to reference choice, transcript view, sequence
thresholds and tree perturbation must be reported. Unknown results are not
correct negatives and cannot silently inflate accuracy.
