# IntraPhy 0.17 actual validation

- Source: 365 tests, 355 passed, 10 explicit missing-MAFFT skips.
- Installed wheel: same results. It is a separate project venv, with existing third-party numeric dependencies exposed; this is not a fresh online dependency installation. The initial dependency-path failure log is retained.
- Added 68 tests: 56 scope/candidate/visual contracts and 12 native-input extraction cases.
- CLI: parsimony, ER/ARD and foreground, 1 and 2 threads; 14 core result tables identical. 99 grouped SVGs parsed as XML.
- Four-taxon prepared-input gallery: 24 SVGs under all/high-coverage. Four final history images were rendered to PNG and inspected.
- Runtime Python 3.13; 3.9 syntax parse only. Other interpreters not executed.
- No new real biological dataset run or statistical calibration. Missing external tools are not replaced with fake executables.

`validation.json` records scope, dependency versions and changed modules. Full logs accompany it. Older files in the parent validation directory are historical 0.16 evidence.
