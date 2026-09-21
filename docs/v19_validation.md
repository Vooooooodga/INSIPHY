> Historical **0.19.0** validation record. Current release: [0.19.1 validation](v0191_validation.md).

# V19 implementation validation

This record distinguishes implementation checks from biological accuracy and
finite-sample statistical calibration. It records the local pre-publication
validation only. Remote publication was not performed in this environment.
The updated CI workflow is configured for fresh Python 3.10/3.13 validation;
its result must be checked after an actual push, not inferred from this record.

## Local checks (2026-09-21)

| Check | Observed result |
|---|---|
| Full regression and new unit/integration suite | 516 tests passed, no failures/errors/skips; 107.738 seconds |
| Independently retained raw-input scenario runner | 20/20 scenarios passed |
| Installed-wheel command runner | 11 commands passed, 5 readable SVG files |
| Installed namespace and entry points | `intraphy` and `python -m intraphy` report 0.19.0; removed namespace unavailable |
| Source layout and Python 3.10 syntax | Passed; production/test/tool modules at most 500 lines |
| Real alignment tools | MAFFT v7.505; minimap2 2.26-r1175 |

Local Python was 3.13.5. The wheel was installed outside the source checkout.
The local isolated import environment explicitly reused installed numeric and
Biopython dependency directories because network access was unavailable; it
was **not** a fresh dependency resolution/installation. The GitHub Actions
matrix installs dependencies and the wheel in a new environment on each runner.

Early local failures were retained in the working validation archive: the first
wheel environment could not import Biopython; after explicitly providing its
numerical dependencies, a smoke-runner assertion used an outdated probability
status label. The final runner checks the actual, longer conditional-status
label rather than dropping the check. A gene-only-locus integration failure
was fixed by preserving that locus as structurally unknown, not inventing exons.
These failures were not rewritten as successful runs.

## What is checked

* Original V18 regression semantics, retained only for explicitly selected legacy
  analyses; V19 defaults are separately tested.
* Generalized Sankoff results versus exhaustive ancestral-state enumeration:
  all optimal node states/branch pairs, non-additive alternative placements,
  globally compatible representative histories and unknown observations.
* Finite-state pruning versus enumeration; valid generators, zero-length branches,
  the semigroup identity, conditional state probabilities, marked transition
  counts and the distinction between endpoint change and at least one edit.
* Split/fusion and insertion/deletion accounting, including one interval deletion
  affecting multiple exons and one insertion with a split as its consequence.
* No reintroduction of a deleted material source through an inverse edge;
  ancestor intermediate configurations and explicit refusal at state-space caps.
* Per-opportunity weights are unchanged by duplicate catalogue entries or drawing
  fragments. Model state indices never stand in for native exon identifiers.
* Pooled scale fitting, whole-gene resampling, retained failed bootstrap draws,
  invalid foregrounds and withholding Monte Carlo significance for an
  annotation-discovered catalogue. These checks are not an empirical power/FDR
  calibration and do not validate the discovery process.
* Raw FASTA/GFF/tree tests distinguish real sequence changes from annotation-only
  damage and missing sequence. Truth files are never analysis inputs.

## Raw scenarios

Conserved structure; inserted separator; intronization; exact fusion at coding
phases 0, 1 and 2; one-exon deletion; continuous multi-exon deletion; annotation
dropout; an assembly gap; donor and acceptor changes; a neutral upstream indel;
local duplication; local inversion; a negative-strand locus; genuinely coexisting
annotated structures; UTR; a microexon; and a noncanonical splice motif.

The outcome includes conservative unresolved cases. In particular, duplicate
or reverse-copy evidence is not forced into ordinary exon loss. Plausible
boundary reannotation may remove a change under the evidence-compatible view;
this is reported separately from the annotation-conditional history.

## Optional tools

CESAR 2.0 and AGAT adapters have interface/error/provenance tests, not real-tool
execution in this local validation. CESAR predictions remain separate candidate
evidence and do not automatically modify the structure catalogue or source GFF.
No third-party source or profile is redistributed. The standard raw-input
workflow actually executes MAFFT and minimap2.

## Not established by these results

No independent empirical exon-structure accuracy benchmark has been completed.
No general operating range for divergence, repeats or annotation quality has
been calibrated. No valid one-to-one map from structural edits to molecular
mutation mechanisms is claimed. No transcript usage, RNA-seq requirement or
selection test has been added. Inference is conditional on the declared local
catalogue, alignment, annotations, tree and operation model.

## Reproduce

```bash
python -m pip install -c requirements-ci.txt . build
python tools/check_source_layout.py
python -I tools/verify_namespace.py
python -m unittest discover -s tests -v
python tools/validate_v19.py --output-dir validation/current/v19-raw
python -m build
python -m venv /tmp/intraphy-v19-wheel
/tmp/intraphy-v19-wheel/bin/python -m pip install -c requirements-ci.txt dist/*.whl
PATH=/tmp/intraphy-v19-wheel/bin:$PATH /tmp/intraphy-v19-wheel/bin/python -I tools/verify_namespace.py
/tmp/intraphy-v19-wheel/bin/python tools/smoke_validate.py --output-dir validation/current/wheel-smoke
```

The raw runner retains input files, separate truth, per-case logs and results,
versions, failures and a machine-readable `validation.json`. The installed-wheel
runner retains every command and its exit status. CI uploads these products even
when an earlier validation stage fails.
