# Installation and external programs

Use Python 3.10 or later in an isolated environment. The package dependencies are
NumPy, SciPy, Biopython and NetworkX. `requirements-ci.txt` records a reproducible
numeric baseline with Python-version markers. It is a validation environment,
not a claim that every permitted upstream release behaves identically.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -c requirements-ci.txt .
python -m intraphy --version
intraphy inspect-aligners
```

For editable development, use `python -m pip install -c requirements-ci.txt -e .`.
The distribution, console command and import are named `intraphy`. Uninstall any
separately installed older distribution from the environment; v18 does not
provide a forwarding package or old console alias.

## External tools

MAFFT performs the coding family alignment and default exon-pair alignments.
minimap2 is used for supported nucleotide searches. Both are required by the
standard raw-input example. On Ubuntu, install `mafft` and `minimap2` with apt.
The baseline CI image pins Ubuntu 24.04 package versions and records both package
and executable versions. External-tool changes require a new integration test.

miniprot is optional for protein-to-genome candidate projection. Install its real
executable on PATH before selecting `--evidence-aligner miniprot`; minimap2 remains
required for the separate nucleotide channel. LASTZ is optional and must be
installed when explicitly selected. Optional backends are not claimed validated
by tests that only run MAFFT and minimap2.

`intraphy inspect-aligners` reports resolved paths and executable versions.
`environment.json` also records Python and numeric-library versions for primary
commands. The program does not inspect secrets or dump the process environment.

## Clean package verification

Build a wheel using `python -m build`. Install the wheel in a new environment and
run from outside the source checkout. Check `intraphy --version`,
`python -m intraphy --version`, and the complete raw-input example. The source tree
must not be providing imports accidentally through PYTHONPATH.

The release's validation report distinguishes local source tests, wheel smoke
tests and fresh dependency installation in CI. Historical logs are not evidence
that an optional executable or a different Python interpreter was tested.


## Optional AGAT

An AGAT installation on PATH enables `intraphy normalize-annotation`. The required
executable is `agat_convert_sp_gxf2gxf.pl`. Follow the official
[AGAT installation instructions](https://agat.readthedocs.io/en/latest/).
The core package does not install AGAT or change annotations automatically.
Real MAFFT/minimap2 integration and AGAT mocked-process interface tests are
reported separately. No real-AGAT validation is implied by a passing unit suite.
