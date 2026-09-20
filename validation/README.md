# Validation history

**0.17 current evidence is in [v017/README.md](v017/README.md).** The files listed below are retained 0.16 records, unless explicitly suffixed v017. No historical biological demo was rerun in this delivery.

# Local validation records

- `regression_tests.log`: 297 discovered tests, 287 executed successfully, 10 explicit MAFFT-dependent integration skips.
- `baseline_tests.log`: original uploaded source baseline, 268 tests; 10 missing-MAFFT errors and 2 strand subtest assertion failures. This is not a v0.16 result.
- `cli_smoke_results.json`: actual CLI commands for parsimony, ER/ARD and foreground; 1/2-thread tables compared byte-for-byte; SVGs parsed as XML. This uses a synthetic four-species fixture, not biological validation.
- `module_ownership.tsv`: where all original 23 modules' responsibilities now live. Public-module line counts exclude extracted implementations and must not be described as total code reduction.
- `architecture_before/` and `architecture_after/`: static AST inventory; no content hashes or executed package code. Revision fields refer to the local working repository state, not a remote GitHub release.
- `install_validation.json`: packaging and isolated project-install checks, when completed. Numeric dependencies are inherited from the session environment and no all-platform claim is made.

The retained `demo_results/` files are upstream historical outputs, not new v0.16 runs. External-aligner integration and real cases require their executables/data and were not repeated in this session. The old calibration fixture is not statistical calibration of the formal model.
