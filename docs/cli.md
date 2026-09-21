> **V18 legacy documentation / retained reference.** The default V19 model and
> commands are described in [the README](../README.md) and
> [exon_structure_model.md](exon_structure_model.md). Do not interpret the old
> independent-layer analyses as exon configuration inference.

# Command-line workflow

The standard interface uses genomic FASTA, matching GFF/GTF and a rooted tree.
Use `--orthologs` to select genes from whole-genome annotations. No hand-written
manifest is required. See [input examples and coordinate handling](inputs.md).

```bash
intraphy check --fasta genomes --gff annotations --orthologs families \
  --species-tree tree.nwk
intraphy build-case --fasta genomes --gff annotations --orthologs families \
  --species-tree tree.nwk --output-dir work --threads 8
intraphy run --input-dir work --output-dir results --threads 8
intraphy visualize --input-dir work --result-dir results --output-dir figures
```

`check` validates matching file names, selected locus IDs, tree tips and genomic
coordinate bounds without requiring external aligners. It reports resolved loci;
it does not infer orthology. `extract-loci` is an optional coordinate-safe genomic
FASTA/GFF export. `normalize-annotation` is an explicit optional AGAT adapter.

## Primary commands

Use `example`, `check`, `build-case`, `run`, `infer-phylogeny`, and `visualize`.
`explain` renders the standalone synthetic methods guide without genomic inputs.
`import-orthofinder` prepares a supplied orthogroup. `inspect-aligners` and
`inspect-annotation` provide diagnostics. Each command has `--help`.

The raw example in the README is the recommended first run. Direct file paths
are relative to the current working directory. In the optional table interface,
resource paths are relative to the manifest.
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
