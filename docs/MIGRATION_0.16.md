# Migration to 0.16.0

This distribution was refactored from the user-supplied INSIPHY-main.zip (0.15.0 source). It is a complete source tree and locally built wheel, not a patch loader or a GitHub release pushed on your behalf.

## Installation

From the unzipped source directory:

```sh
python -m pip install .
insiphy --version
PYTHONPATH=src python -m unittest discover -s tests -v
```

Python >=3.9 is retained. This session executed Python 3.13 only; 3.9 syntax was parsed and CI now declares 3.9/3.11/3.13 jobs. MAFFT and other selected external aligners are system executables, not bundled by the wheel. No substitute aligner is silently inserted when MAFFT is absent.

## Existing scripts

CLI flags and commands remain. Old Python import paths re-export canonical implementations. Do not copy only one old top-level .py file into this package: the responsibility packages must be installed together. Dataclass module paths have moved; previously pickled private objects are not a supported interchange format. TSV/JSON and declared CLI inputs remain the compatibility surface.

## Deliberate behaviour corrections

- Correspondence no longer silently consumes previous output completion records. To reuse them, pass the explicit keyword input.
- An output directory with conflicting legacy formal results is not guessed by the renderer. New runs record `run_result.json`; incomplete runs are not rendered as prior successes.
- Free discussion words such as “predicted” and “unknown” in evidence prose no longer determine current visual state. Use the structured fields.
- Empty half-open intervals do not overlap other intervals. Malformed TSV width, duplicate headers, invalid evidence JSON and missing explicitly supplied trees produce input errors.
- Near-optimal chain ambiguity includes paths beyond the first two and alternatives across distinct transcript contexts. A graph with no connecting required-anchor path has an explicit status.
- A jointly resolved complementary protein correspondence is no longer vetoed by incomplete DNA-candidate enumeration. It retains its protein score and protein coordinates, while raw DNA score/CIGAR/blocks remain unchanged. Both strand regressions exercise this correction.
- `source_label` and `copy_role` must come from the manifest. A name such as jingwei or Sdic no longer implies its historical role.
- `baseline_comparison.tsv` is now explicitly “not_a_method_comparison”, without hand-written method reward scores.
- Benchmark identities retain site/object identity. Missing truth is unavailable, not a zero-performance evaluation. Older truth records lacking object identity remain legacy/coarse inputs and cannot validate per-site accuracy.

## Tests and private mocks

Several original tests now patch the canonical implementation owner rather than a re-export. Two candidate-context expectations and the legacy baseline assertion were intentionally updated to the corrected contract. Ten external integration tests explicitly skip when MAFFT is not installed; with MAFFT they execute. No failed biological assertion was replaced by a pass.

## Remaining work is not hidden

The legacy simulation/calibration route remains explicitly experimental; it does not validate the formal single-copy CTMC. Backend Stats compatibility and some TSV checkpoints remain. Full candidate-type elimination, optional optimized interval-tree backend, additional plotting view objects and deeper decomposition of long evidence functions can be considered later, with profiling and equivalence tests. They are not claimed complete merely because code was moved.
