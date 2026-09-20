# Command-line workflow

## Primary commands

Use `example`, `check`, `build-case`, `run`, `infer-phylogeny`, and `visualize`.
`import-orthofinder` prepares a supplied orthogroup. `inspect-aligners` and
`inspect-annotation` provide diagnostics. Each command has `--help`.

The raw example in the README is the recommended first run. Relative resource
paths are interpreted relative to the manifest, not the current working directory.
The tree may be Newick at preparation; the prepared tree is stored as TSV.
No shell-specific quoting is needed for normal paths. Paths with spaces must be
quoted by the shell in the usual way.

## Outputs and failures

Primary outputs have an execution manifest, a log and environment provenance.
A nonempty directory is refused by default. `--force` is limited to recognized
IntraPhy output directories and preserves the previous directory as a timestamped
backup. Output must not contain its own input, be a symlink, or be the root, home
or current working directory. An existing lock is never silently removed.

After an interrupted run, inspect `execution.json`, `run_result.json` and
`intraphy.log` before removing a stale `.intraphy.lock`. Use a new output directory
for the next run. V18 does not claim automatic checkpoint resume.

Normal errors return a nonzero code with an actionable message. `--debug` prints
a Python traceback. `--quiet` suppresses routine terminal progress but retains
file logging. Unavailable statistical inference is a valid, explicitly recorded
scientific outcome and is distinct from a failed command.

## Frozen matrices and views

Run parsimony first and reuse its `structural_site_matrix.tsv` explicitly for
other models. Retain adjacent catalogue and coordinate files. A canonical
sensitivity analysis must generate its own matrix in a separate directory;
canonical and repertoire matrices cannot be mixed.

`--analysis-range all` is the default. `high-coverage` and
`--min-callable-fraction` define an optional sensitivity subset. Known absence is
a called state. Unknown does not mean absence. The coverage threshold is not a
homology threshold and never alters coordinates or creates transcript links.

## Experimental and developer commands

The old multicopy model, simulator and calibration commands remain explicitly
experimental. They do not calibrate the formal single-copy model. Their historical
outputs may contain legacy event categories and are not formal v18 conclusions.
Nonzero bootstrap and stochastic-map flags are rejected for the formal model.
Low-level extraction and correspondence commands are intended for experienced
users inspecting intermediate evidence, not required steps in the basic workflow.
