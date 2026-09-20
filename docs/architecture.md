# IntraPhy 0.16 architecture

## Data flow

`workflow.run_all`: prepared inputs → initial correspondence → sequence evidence → annotation completion → **explicit** final correspondence → one `ObservationMatrix` → selected inference engine → `RunResult` → report/plots.

An existing `annotation_completion_candidates.tsv` in the output directory is not implicitly read by correspondence. Pass `annotation_rows=` or `annotation_completion_path=` deliberately. Independent stage commands still consume documented checkpoint files. This release does not claim that all legacy TSV checkpoints have disappeared.

## Ownership

| Package/module | Responsibility |
|---|---|
| `storage/` | Plain/gzip TSV and FASTA I/O; scalar boundary conversions |
| `preparation/` | Parsed annotation/index, transcript selection and extraction, copy context |
| `aligners/` | Candidate types, subprocess backends, formats, short/pairwise alignment, projection and legacy codecs |
| `mapping/` | Candidate qualification, coordinate evidence, ordered chains, membership refinement |
| `coding/`, `coding_correspondence.py` | CDS/codon coordinates and family protein correspondence |
| `evidence/`, `annotation.py` | Local evidence search/projection and explicit completion workflow |
| `observations/`, `structural_sites.py` | Schema and the unique biological observation-building layer |
| `topology.py` | Validated rooted topology independent of models |
| `inference/` | Pure Sankoff/CTMC kernels, fitting, profiles, diagnostics and posterior calculations |
| `parsimony.py`, `structural_phylogeny.py` | Formal run drivers and compatible helper exports |
| `phylogeny.py` | Explicit formal/experimental dispatch |
| `reporting/`, `visualize.py` | Descriptive summaries, typed result availability, view preparation, SVG primitives |
| `experimental/` | Legacy multicopy mathematics, explicitly outside formal single-copy scope |
| `verification/` | Distributed simulation/benchmark commands; oracle fixtures remain labelled as such |
| `commands/parser.py`, `cli.py` | Existing argparse definitions and thin dispatch |
| `workflow.py`, `run_result.py` | Stage order, run parameters and completed-result ownership |

Compatibility modules deliberately re-export existing helpers. Function implementation ownership is singular, while public old import paths remain available. Private monkey-patches must target the canonical lookup site. Public function imports are preserved; assigning a new value to a re-exported private name is not an API promise.

## Contracts

`ObservationMatrix` snapshots scalar TSV rows into read-only mappings. It is prepared once by formal dispatch and shared by the chosen engine. Frozen input is validated and checkpointed, not rebuilt from GFF. Existing schema-v3 semantics and legacy schema readers remain intact.

`RunResult` identifies the selected model, scope, tree and actual artifacts. A `running` manifest prevents rendering an older successful result after a failed new run. Without a manifest, unambiguous legacy tables remain readable; ambiguous formal result tables require explicit selection. Free English evidence/conditioning text does not determine current result validity. The pre-0.15 posterior adapter remains an explicit compatibility boundary.

`AnnotationIndex` parses one resource, retaining duplicate IDs, multi-Parent records and antisense neighbors. ID/Parent/token and closed-coordinate spatial indexes are reused across genes. Spatial queries preserve source order. The bounded process-local cache holds at most two resources and uses path/mtime/size for lifecycle invalidation, not content hashes. Input resources should remain immutable during a run; call `clear_annotation_cache()` after deliberate in-place replacement preserving stat metadata.

`FeatureHierarchy` is shared when projecting all transcripts of a copy. Full codon-to-genome correspondence is retained; an affine coordinate shortcut is not applied across introns.

`AlignmentCandidate`/`AlignmentCandidateSet` and the compatibility `AlignmentStats` definitions have one owner in `aligners/types.py`. Legacy conversions are isolated. Some existing backend interfaces still return Stats for API compatibility; this release does not remove all conversion costs. Candidate search completeness, scoring scales, coordinate blocks and alternative mappings remain separate facts.

## Near-optimal chains

A context graph is built once for exact/near-optimal evaluation. Membership means a candidate belongs to at least one admissible near-optimal path. Mandatory membership is evaluated by excluding the candidate and checking whether another path still reaches the score threshold. The best two paths remain a bounded display/count aid, not proof that a node is mandatory. Independent transcript contexts are also included in ambiguity aggregation. The correctness-oriented exclusion pass can cost O(V(V+E)); future optimization needs equivalence tests.

## Deliberately not introduced

No new core Pydantic, pandas, Typer or IntervalTree dependency; no global framework, no pipeline-manager hierarchy, no working-directory changes for concurrency, no checksum pipeline. The in-memory spatial index is a dependency-free baseline, not a claim of universally optimal interval-tree complexity.

No new biological threshold, substitution model or statistical calibration is claimed. Existing real demo outputs were retained as historical material and were not regenerated here. Formal statistical validity still requires its separate scientific validation programme.
