# Input Format

## Upstream assumption

Formal v0.14 analysis starts from one single-copy orthologous gene per species. The ortholog set can come from OrthoFinder or manual curation. INSIPHY analyzes internal structure after gene-level homology has been supplied.

## OrthoFinder importer

`import-orthofinder` reads an OrthoFinder result directory and one orthogroup. The resource manifest is tab-delimited:

```text
species	genome_fasta	annotation_file
```

Recommended additional fields:

```text
assembly	annotation	source_url	release	notes
```

Species names must match normalized OrthoFinder column names. Member IDs are mapped exactly to annotation IDs and common attributes such as gene, transcript, protein, Parent, Alias and Dbxref-derived tokens. Multiple transcript or protein members mapping to the same gene locus count once. Single-copy eligibility requires exactly one distinct mapped gene locus per requested species; unresolved IDs, zero loci or multiple loci are recorded as exclusions.

## Case manifest

`build-case` reads:

```text
case_id	species	family_id	gene_id	gene_copy_id	genome_fasta	annotation_file
```

Formal single-copy input requires one `gene_copy_id` for every `family_id, species` pair. Paths may point to plain or gzip-compressed FASTA and GFF3/GTF files.

## Species tree

`species_tree.tsv` requires:

```text
node_id	parent_id	label	branch_length
```

- Exactly one root row has an empty `parent_id`.
- Every non-root node points to an existing parent.
- Leaf labels match species in the structural matrix.
- Node IDs and leaf labels are unique.
- In supplied mode, every non-root branch length must be finite and nonnegative. Zero is retained as a true zero-length branch with identity transition matrix; missing lengths are an input error. Unit mode explicitly assigns length one.

Newick input is accepted by the importer and converted to this table.

## Transcript policy

The current method description treats all annotated transcript paths as the default biological repertoire. In that mode:

- all transcript paths are preserved;
- identical genomic structures can be deduplicated;
- alternative splice boundaries remain separate observations;
- a structure used by any annotated transcript contributes repertoire evidence;
- conflicts among paths are retained as alternative or unknown states rather than forced into one consensus.

Canonical transcript selection remains useful for focused analyses and compatibility testing. When used, it should be stated explicitly.

## Prepared case tables

`build-case` and `derive-tables` generate:

- `segment_occurrences.tsv`: observed, predicted, candidate, or absence intervals;
- `segment_sequences.fasta`: interval sequences;
- `gene_loci.tsv`: original gene bounds, search bounds, strand, and source metadata;
- `gene_loci.fasta`: oriented search-window sequences;
- `transcript_paths.tsv`: ordered transcript paths;
- `intron_sites.tsv`: intron boundary, motif, and phase evidence;
- `physical_adjacencies.tsv`: neighboring intervals within each gene;
- `segment_matches.tsv`: pairwise alignment and context scores;
- `segment_homology.tsv`: internal evidence components;
- `sequence_synteny_evidence.tsv`: annotation-completion and sequence evidence;
- `copy_context.tsv`: retained upstream copy metadata;
- `raw_gene_features.tsv`: original annotation features overlapping the declared search window, including feature type, coordinates, parent IDs, attributes and `ownership` (`target_gene_descendant` or `overlapping_context`). Targeted retention tests passed and the first real runs produced this table. Retaining a feature does not establish its homology or add a statistical layer.

## Alignment blocks used by figures

Ribbons read direct `match_status=mapped` pairs from `segment_matches.tsv`. A result-directory table takes precedence over the prepared input table. Both endpoints must be confirmed exon-like observations. Complete annotation boxes retain their full extent; ribbon endpoints use only accepted aligned intervals.

- `correspondence_basis=annotated_CDS_protein` selects `protein_projected_blocks`; other matches use `projected_reference_blocks`.
- Blocks use `qstart-qend:tstart-tend`, separated by semicolons. Coordinates are 1-based inclusive nucleotide positions local to each occurrence, oriented 5-prime to 3-prime on either genomic strand.
- Protein blocks are already transcript-oriented. DNA blocks require positive relative alignment orientation. Missing, conflicting or reverse blocks do not produce a ribbon; protein matches never fall back to DNA blocks.
- Each SVG ribbon records the match ID, correspondence basis, projection field and endpoint coordinates. Transitive group membership without a direct match retains membership color but has no base-correspondence ribbon.

## Element correspondence

`element_correspondence.tsv` is the primary correspondence table. Important fields:

- `element_id`: stable homologous unit label such as `EG_*`;
- `homology_id`: internal evidence component;
- `occurrence_id`: member occurrence;
- `element_class`: `exon_like`, `candidate_source`, `absent`, or `context`;
- `display_role`: observed or inferred display role;
- `support_type`: evidence source such as annotation, sequence candidate, prediction, or absence support;
- `membership_call`: `core_member`, `ambiguous_member`, or equivalent status;
- `inferred_role`: confirmed or predicted role used by the structure layer;
- `predicted_role`: predicted annotation class such as `CDS`.

`predicted_exon_candidate`, `inferred_role=predicted_CDS`, and `predicted_role=CDS` indicate prediction evidence. Formal state coding and visualization keep them distinct from confirmed exonic role.

## Structural matrix

`structural_site_matrix.tsv` stores one row per family/layer/site/species observation.

Layers:

- `exon_presence`;
- `exon_role`;
- `splice_junction`.

States are `state_0`, `state_1`, or `unknown` under that layer's definition. Unknown means the current input cannot distinguish states.

## Output interface changes in v0.14

The public extant-summary files are:

- `observed_element_tree_coverage.tsv`;
- `observed_intragenic_paths.tsv`.

The old `ancestral_element_graph.tsv` and `ancestral_intragenic_paths.tsv` names are removed from the v0.14 public interface. Complete ancestral transcript graphs are outside the current formal output.

## Branch-length and ascertainment flags

- `--branch-length-mode supplied`: use finite nonnegative lengths from the input tree; zero gives identity transitions and missing lengths are errors.
- `--branch-length-mode unit`: use one for every non-root branch.
- `--ascertainment observed-at-least-one`: default for discovered sites.
- `--ascertainment complete-universe`: use an explicit candidate catalogue.
- `--ascertainment variable-only`: use when constant patterns were intentionally excluded.

Ascertainment affects optional likelihood analysis. Default parsimony uses the topology and observed/unknown states.

## Experimental multi-copy input

`copy_tree.tsv` and `gene_tree.tsv` remain accepted under `--analysis-scope experimental-multicopy`. Formal v0.14 single-copy statistics do not use them.
