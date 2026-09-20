# IntraPhy maintenance contract

- Scientific state construction belongs in `observations/`; codecs parse/validate, inference consumes a frozen matrix, reporting displays explicit results.
- Public modules preserve imports/CLI. See `docs/architecture.md` for canonical owners; edit the owner, not a compatibility export.
- No gene-name-based biological roles. Source/copy roles come from input manifests.
- Unknown, absent, not-applicable, annotation-conditional and predicted evidence remain distinct. Preserve negative-strand, codon/phase and alternative-path tests.
- Do not replace near-optimal candidate enumeration with one longest path. Do not replace Sankoff with a matrix exponential.
- External tools use per-task temporary directories and subprocess `cwd`; do not change process cwd. Do not introduce per-layer executors or content-hash audits.
- Keep formal single-copy, experimental multicopy and verification fixtures separate. Do not describe legacy calibration as calibration of the formal model.
- First run targeted tests; before delivery run `PYTHONPATH=src python -m unittest discover -s tests -v`. Missing MAFFT yields explicit integration-test skips, not biological success.
- Keep task output short: changed contract, actual tests, unresolved limits. Do not repeatedly dump whole files or copy long plans into this file.
