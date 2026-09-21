# IntraPhy 0.19.1 validation record

This is a local source/wheel validation, **not a GitHub release or CI run**.
The delivered 0.19.0 source is the baseline. No remote branch was changed.

| Check | Observed outcome |
|---|---|
| Full regression, unit and integration suite | 563 passed; 0 failures, errors or skips; 118.0 seconds |
| Newly added audit/figure regressions | 47 tests, included in that full-suite total |
| Preserved original raw structural examples | 20/20 passed |
| Additional raw audit cases | 9/9 passed |
| Installed wheel, outside source directory | 11 commands passed; 18 SVG files parsed |
| Source-layout guard | Passed; maximum production module 495 physical lines |
| Python 3.10 syntax | Passed under 3.13 AST parser; NOT a 3.10 runtime execution |
| Actual interpreter | 3.13.5 |
| Actual native tools | MAFFT 7.505; minimap2 2.26-r1175 |

Numeric dependencies were NumPy 2.3.5, SciPy 1.17.0,
Biopython 1.86, NetworkX 3.6.1.
The wheel was installed as a package in a separate venv and run with `python -I`
from outside the source directory. Existing dependency directories were explicitly
linked into that venv. This is not a completely fresh dependency installation.
The first isolated package check exposed an old hard-coded version assertion and
missing reused Biopython dependencies; both were corrected and the complete wheel
check rerun. Diagnostic records are preserved rather than claimed as successes.

## Audit outcomes

Identical-DNA single-end, double-end, erroneous split and erroneous fusion now
have zero minimum **candidate-compatible** edits while their original-annotation
histories retain the expected differences. Missing first/last exon annotations
retain all 1140 bases available in the synthetic input. Signal-changing positive
controls and precise DNA deletion are not blanket-masked into no-change results.
These are implementation tests; matching or mutated motifs are not independent
proof of biological splicing.

A repeated physical region under a new unit ID is rejected across input routes.
An unobserved 1.0 branch and its 0.5+0.5 subdivision give equal log likelihood
(-2.753705983974205 in the source audit example) and identical root probabilities.
A declared two-exon source module has a one-insertion path. The four-exon audit
space contains all 556 states and 4803 elementary edges under the default 1024 cap;
an explicit 256 cap still refuses rather than truncating and renormalizing.

The legal repeated-gene draw `g,g` reaches numerical fitting. In the recorded
example it is genuinely nonidentifiable by the curvature diagnostic; no interval
is manufactured from it. This is a different, appropriate outcome from refusing
the draw merely because only one unique label was sampled.

Current explanation values are computed from saved example inputs and parameters.
Read-only CTMC rendering responds to changed node/branch probabilities, rejects
invalid probability vectors and old-model histories, and does not call inference.
Current and legacy explanation plates are separated.

## Reproduce

```bash
python -m unittest discover -s tests -v
python tools/validate_v19.py --output-dir validation/current/raw
python tools/validate_v0191.py --output-dir validation/current/audit
python tools/check_source_layout.py
python tools/render_method_figures.py --output-dir validation/current/guide --png
python tools/smoke_validate.py --output-dir validation/current/wheel-smoke
```

Run the final command with the installed wheel's interpreter outside the source
root. PNG previews need CairoSVG. MAFFT/minimap2 must be real available binaries
for the raw tests; a run that skips them is not equivalent to this validation.

## Not validated by these results

No independent published-case biological benchmark, general finite-sample
operating-characteristic calibration or actual AGAT/CESAR2 integration run was
completed here. No exon-usage or transcriptome analysis was added. The model still
conditions on qualified alignment, the finite candidate catalogue and explicit
rates/priors. Complex copy genealogy, inversion/shuffling histories and a full
continuous-time origin process remain outside this correction release.

The local delivery includes machine-readable summary, complete logs, raw inputs,
separate generating truth, actual inference outputs and source-patch verification.
The 20+9 raw cases overlap with the unit suite and must not be counted as separate
independent biological validations.
