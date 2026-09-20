# Software architecture

The package uses explicit imports and ordinary functions. There is no dynamic
source loader, new workflow framework or runtime refactoring dependency.
Production Python modules are kept at or below 500 physical lines; the release
check fails on an oversized file.

| Responsibility | Owners |
|---|---|
| Annotation indexing and extraction | `preparation/` |
| External tools, formats, pairwise/short alignments | `aligners/` |
| Coding coordinates and protein correspondence | `coding/`, `coding_correspondence.py` |
| Candidate context, chaining, local membership and reference coverage | `mapping/` |
| Nucleotide search, protein projection and candidate classification | `evidence/` |
| Observation semantics, character catalogue, masks and applicability | `observations/` |
| Validated rooted topology | `topology.py` |
| Sankoff, CTMC, fitting, profiles, diagnostics and posterior export | `inference/` |
| Formal run dispatch | `phylogeny.py`, `parsimony.py`, `structural_phylogeny.py` |
| Primary sequence-to-result order | `workflow.py` |
| CLI parsing, validation, output ownership and provenance | `commands/`, `cli.py` |
| Target figures and reporting | `reporting/` |
| Explicitly nonformal multicopy code | `experimental/` |
| Synthetic raw inputs and historical verification utilities | `verification/` |

`workflow.run_all` performs initial correspondence, new sequence evidence,
annotation-candidate classification, explicit final correspondence, and inference.
The observation layer creates the matrix once. Numeric engines do not change
homology, create absent tips or reinterpret free-text notes as scientific evidence.

Important refactored owners include `evidence/sequence_search.py`,
`mapping/ordered_chains.py`, `mapping/path_evaluation.py`,
`mapping/chain_qualification.py`, `mapping/reference_coverage.py`,
`inference/layer_fitting.py`, `inference/model_comparison.py`, and
`inference/posterior_export.py`. Patch or modify a function at its implementation
owner. Re-exporting a function does not make a private monkey-patch propagate to
all modules.

The character catalogue, coordinate sidecar and frozen matrix preserve inference
inputs. Correspondence alternatives are evidence, not posterior history samples.
Ancestral consistency is audited without constructing a complete transcript.
A failed execution is marked failed and cannot silently display a previous run.

See the [methods figure](figure.md) for scientific data flow. The figure separates
biological evidence, character construction and conditional reconstruction rather
than presenting implementation files as biological stages.
