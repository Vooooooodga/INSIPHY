# Input Format

## Upstream assumption

Formal analysis starts from one single-copy orthologous gene per species. The
ortholog set can come from OrthoFinder or curation. INSIPHY does not infer
gene-level orthology.

## OrthoFinder importer

`import-orthofinder` reads `Orthogroups/Orthogroups.tsv` and selects one
orthogroup. The genome resource manifest is tab-delimited and requires:

```text
species	genome_fasta	annotation_file
```

Optional provenance fields are `assembly`, `annotation`, `source_url`,
`release`, and `notes`. Species names must match normalized OrthoFinder
column names.

The importer accepts a Newick or INSIPHY TSV species tree. Newick branch
lengths are retained. Default parsimony uses the rooted topology alone.
Optional likelihood analysis needs lengths or `--branch-length-mode unit`.

## Case manifest

`build-case` reads:

```text
case_id	species	family_id	gene_id	gene_copy_id	genome_fasta	annotation_file
```

Formal single-copy input requires exactly one `gene_copy_id` for every
`family_id, species` pair. Recommended provenance fields are:

```text
assembly	annotation	source_url	release	notes
```

Paths may point to plain or gzip-compressed FASTA and GFF3/GTF files.

## Species tree

`species_tree.tsv` requires:

```text
node_id	parent_id	label	branch_length
```

- Exactly one row has an empty `parent_id` and defines the root.
- Every non-root branch length must be positive in `supplied` mode.
- Every species in the structural matrix must match one leaf `label`.
- Internal labels must be unique when they are referenced by foreground files.

Branch lengths may represent time or substitutions per site. The fitted
structural rates inherit that unit.

## Foreground branches

The foreground file may use either:

```text
parent_id	child_id
```

or:

```text
branch_scope
parent_label->child_label
```

Every listed edge must occur in `species_tree.tsv`.

## Prepared case tables

`build-case` and `derive-tables` generate:

- `segment_occurrences.tsv`: observed or sequence-supported intervals;
- `segment_sequences.fasta`: interval sequences;
- `gene_loci.fasta`: oriented search-window sequences;
- `gene_loci.tsv`: original annotation bounds, search bounds and genome source;
- `segment_matches.tsv`: pairwise alignment and context scores;
- `segment_homology.tsv`: internal homology-component membership;
- `transcript_paths.tsv`: ordered transcript paths;
- `intron_sites.tsv`: intron boundary, motif, and phase evidence;
- `physical_adjacencies.tsv`: neighboring intervals within each gene;
- `sequence_synteny_evidence.tsv`: annotation and sequence evidence;
- `copy_context.tsv`: upstream copy metadata retained for provenance.

`segment-correspondence` generates `element_correspondence.tsv`. Important
fields are:

- `element_id`: exon-like homologous unit (`EG_*`);
- `occurrence_id`: observed member;
- `element_class`: `exon_like`, `candidate_source`, or `absent`;
- `membership_score`: correspondence support;
- `membership_call`: `core_member` or `ambiguous_member`.

Only `core_member` observations enter formal tip-state coding. Ambiguous
members remain in the evidence output and are coded as unknown.

## Transcript policy

`--transcript-policy canonical` chooses one transcript using
`--canonical-rule longest_cds` by default. `all` retains every annotated
transcript. Conflicting states among retained transcripts are coded as unknown
for the affected structural site.

## Branch-length and ascertainment flags

- `--branch-length-mode supplied`: use positive lengths from the input tree.
- `--branch-length-mode unit`: replace every non-root length with one.
- `--ascertainment observed-at-least-one`: default for sites discovered because
  at least one species contains the structure.
- `--ascertainment complete-universe`: use only when the input defines
  meaningful all-zero candidates.
- `--ascertainment variable-only`: use only when both constant patterns were
  deliberately excluded.

Ascertainment and branch-length modes apply only to the optional probability
analysis. `--model parsimony` is the default. Gene extraction accepts `--flank`
and `--max-extension` in bp; these control search extent and do not establish
biological confidence.

`--ascertainment complete-universe` requires `structural_site_universe.tsv`
inside the input directory. Its required columns are `family_id`, `layer`,
and `site_id`. A site-only row declares the candidate set; it contributes no
absence observation. Optional `species` and `state` columns provide explicit
curated observations, with `state_0`, `state_1`, and `evidence` recommended.
Explicit observations replace the generated observation at the same
family/layer/site/species key. Missing species remain unknown. The declared
universe must be justified independently of which sites happen to be present.

`--aligner` selects exon correspondence. `--context-aligner` selects local
exon/non-exonic mapping. `run --evidence-aligner` selects supplementary genomic
or protein evidence. MAFFT, minimap2 and miniprot serve different alignment
questions; these choices are recorded with their output evidence.

## Experimental multi-copy input

`copy_tree.tsv` and `gene_tree.tsv` remain supported only under
`--analysis-scope experimental-multicopy`. The formal v0.13 single-copy
statistics ignore these files.
