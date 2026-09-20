# Migration to 0.18.0

This is a breaking namespace and interpretation update.

The only command and import package is `intraphy`; the prior namespace and console
alias are removed. Use a clean environment or uninstall a separately installed
old distribution. Python 3.10 is the minimum supported interpreter.

Formal parsimony no longer writes `compound_structural_events.tsv`. Each mapped
junction remains an elementary character; split/fusion descriptors do not add
counts. Use `gene_change_summary.tsv`, `character_catalogue.tsv` and
`minimum_change_history.tsv`. Do not count possible branch rows as events.

Known linked characters now prevent the entire independent-character CTMC
analysis before fitting. Old outputs that report rates, AIC or posteriors despite
blocked linked-site LRTs must not be reused as v18 results.

The default scope remains all characters. The optional coverage fraction has no
calibrated biological cutoff. Full coordinates and observations are preserved.
Independent catalogues may retain explicitly unknown tips instead of fabricating
zero states.

The v17 three-test CI regression is repaired by restoring conservative assertions
consistent with the resolved-flank evidence requirement. No alignment threshold
was weakened. Raw-input tests now run separately from historical prepared fixtures.
A complementary-fragment repeat-classification defect is corrected using actual
reference coverage and distinct genomic instances.

Large implementations are separated into explicit owner modules. Existing
`intraphy` entry modules re-export functions with normal imports; unit-test patches
must target the canonical runtime lookup. There is no old-name forwarding package,
metaprogramming-based loader or runtime source extraction.

Historical demo outputs and validation logs are retained as historical records.
They have not been silently relabeled or regenerated as new biological evidence.
The English documentation describes v18 only; obsolete planning notes are available
in repository history.


## Direct genomic inputs and teaching figures

The standard interface now accepts `--fasta` and `--gff` files/directories plus
`--species-tree`. `--orthologs` identifies upstream family members by exact GFF
IDs. Combined genomic FASTA with species-specific sequence IDs is supported.
The generated `input_targets.tsv` is provenance output. `--manifest` remains
optional for pipeline use. No old project-name executable or import is restored.
`extract-loci` writes coordinate-consistent genomic FASTA/GFF with automatic
flanks. AGAT normalization is opt-in. See `docs/inputs.md` and `docs/figure.md`.
The standalone `explain` command renders five original methods figures; each
result gallery links to that separate teaching guide. No scaling benchmark or
statistical calibration is implied by the interface completion.
