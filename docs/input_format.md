# Input formats

## Raw gene manifest

Required TSV columns are `species`, `family_id`, `gene_id`, `genome_fasta` and
`annotation_file`. Resource paths may be absolute or relative to the manifest.
Optional `gene_copy_id` defaults to `gene_id`; `case_id` defaults to `family_id`.
Source and copy-role fields are optional metadata and are never inferred from
gene names. Formal input contains one distinct gene locus per species and family.

Genome FASTA names must match GFF3/GTF sequence IDs. Annotation must preserve
feature parentage and strand. All selected supplied transcripts are retained;
`--transcript-policy canonical` is an explicit preparation restriction.
Search flanks and extension bounds do not redefine the annotated gene boundary.

## Rooted tree

A Newick tree can be supplied to `build-case` or `import-orthofinder`. Species
labels must match the manifest. The prepared `species_tree.tsv` uses
`node_id`, `parent_id`, `label` and optional `branch_length`. The root has no parent.
Topology is validated independently of the inference method. CTMC requires finite
nonnegative non-root branch lengths, unless `--branch-length-mode unit` is explicit.
Parsimony does not require branch lengths for its cost.

## OrthoFinder

The resource manifest uses `species`, `genome_fasta` and `annotation_file`.
Supply an upstream orthogroup and the species tree:

```bash
intraphy import-orthofinder \
  --orthofinder-dir OrthoFinder/Results_run \
  --orthogroup OG0001 --genome-manifest genomes.tsv \
  --species-tree species.nwk --output-dir prepared/OG0001
```

Exact annotation ID resolution determines distinct gene loci. Multiple isoform
IDs at one locus do not become extra gene copies. Missing or ambiguous locus
resolution is recorded as an exclusion. This import does not infer orthology.

## Frozen structural matrix

Use a generated `structural_site_matrix.tsv` rather than fabricating negative
states. Schema-v3 rows identify `family_id`, `layer`, `site_id`, `species`,
`state`, `state_0`, `state_1`, annotation view, discovery rule, observation mask,
applicability and evidence. Each character requires exactly one row per tree tip.
`linked_group_id` records character dependence, not a compound-event conclusion.

An independent complete-universe catalogue can contain observed all-zero sites
and explicit unknown tips with missing masks and reasons. Unknown tips must not
be replaced with zero to satisfy a validator. A discovered positive-only catalogue
cannot be relabeled independent after filtering.

Internal tables for occurrences, transcript paths, matches, element membership
and position projection retain source coordinates. Their definitions and readers
are in `observations/schema.py`, `mapping/fields.py` and `storage/`. Stable
character identity is conditional on frozen input and algorithm settings; IDs
are not permanent biological accession numbers.
