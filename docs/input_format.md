# Input and Output Format

## Required biological inputs

The formal first-version analysis requires:

1. A supplied single-copy ortholog set: one gene copy for each requested species in each family.
2. Genome FASTA and existing GFF3/GTF annotation for each species.
3. A user-supplied rooted species tree whose leaves match the species labels.

Gene-level homology is an upstream input. `import-orthofinder` can read an OrthoFinder result and identify one mapped gene locus per requested species. A curated manifest is also accepted. IntraPhy performs gene-internal sequence correspondence after this relationship is supplied; it does not perform genome-wide ortholog discovery or infer a gene tree in the formal first-version analysis. Transcriptome or expression data are not required or consumed.

## OrthoFinder resource manifest

The resource manifest is tab-delimited:

```text
species  genome_fasta  annotation_file
```

Recommended metadata columns are `assembly`, `annotation`, `source_url`, `release`, and `notes`. Species names must match the normalized OrthoFinder columns. Gene, transcript, and protein identifiers are mapped to annotation IDs and common parent/alias attributes. Multiple transcript or protein IDs mapping to one gene locus count as one locus. A family is eligible for formal single-copy analysis when exactly one distinct gene locus maps to every requested species. Unresolved IDs, zero mapped loci, or multiple loci are reported as exclusions.

## Case manifest

`build-case` accepts a tab-delimited manifest with these required columns:

```text
case_id  species  family_id  gene_id  gene_copy_id  genome_fasta  annotation_file
```

Each `family_id, species` pair has one `gene_copy_id`. FASTA and annotation paths may be plain or gzip-compressed. The gene identifiers and copy assignments are supplied by the user or upstream analysis.

## Rooted species tree

`species_tree.tsv` uses:

```text
node_id  parent_id  label  branch_length
```

There is one root row with an empty `parent_id`; every other node names an existing parent. Node identifiers and leaf labels are unique. Leaf labels match the input species. In supplied-length CTMC mode, every non-root length is finite and nonnegative. A zero length is retained as zero; unit-length mode explicitly assigns one to each non-root branch. Newick input can be imported into this table.

## Transcript and annotation policy

The default `transcript_policy=all` retains every annotated transcript path. `annotation_view=repertoire` reports a feature as exonic when at least one supplied path uses it. This is conditional on the supplied annotation and its coverage. It does not establish transcript use in every tissue, condition, or sex. A path that does not cover a queried interval cannot support state 0 there. `canonical` is an explicit alternate analysis and records the selected transcript and selection rule.

Original feature type, coordinates, phase, parent IDs, attributes, and transcript membership are retained for supplied CDS, UTR, coding and noncoding exons, introns, and other annotated feature classes. Feature overlap is allowed when different paths or annotation classes describe the same sequence interval. Predicted coding/exon features are kept as annotation-completion candidates. They do not overwrite the input annotation or independently establish an observed exon-role state.

The original gene-row interval, the union of explicitly associated child features, and the declared genome search interval are distinct. Candidate extension outside an annotation boundary is limited to the searched interval and recorded with its evidence and search completeness. Failure to find sequence past an unsearched or incomplete boundary is unknown. A sequence-absence observation requires a resolved homologous location, ordered flanking anchors, and searched intervening sequence.

## Correspondence evidence and matched blocks

Prepared case and analysis directories contain the applicable subset of:

- `gene_loci.tsv`, `gene_loci.fasta`: original locus bounds, search bounds, strand, and source metadata;
- `transcript_paths.tsv`, `segment_occurrences.tsv`, `segment_sequences.fasta`: provided transcript paths and original or candidate intervals;
- `segment_matches.tsv`, `segment_homology.tsv`: candidate sequence alignments, scoring scheme, candidate alternatives, and correspondence status;
- `element_correspondence.tsv`: element-to-occurrence membership, parent features, transcript IDs, actual matched blocks, coverage relation, repeat instance, and resolution status;
- `sequence_synteny_evidence.tsv`, `annotation_completion_candidates.tsv`: sequence presence, predicted role, local annotation conflicts, and search status;
- `raw_gene_features.tsv`: original features in the declared search window, including feature class, parentage, attributes, and ownership;
- `structural_site_matrix.tsv`: the frozen observation matrix shared by parsimony and optional likelihood analysis.

An aligned block describes only the sequence portion actually matched. It does not imply coverage of the full parent feature. Block records identify their coordinate convention, sequence orientation, and the query/target occurrence. Protein-derived blocks project MSA residues through codons/CDS to genomic sequence; DNA blocks derive from the nucleotide alignment. Missing block coordinates remain unavailable. A legacy match with only a parent interval cannot be expanded into fabricated aligned blocks.

Candidate mappings, alternative hits, score scheme, completeness status, and ambiguity are retained. A complementary 1:n relation has ordered projections to distinct portions of a shared reference interval; repeated overlap maps multiple instances to the same portion. Partial and unresolved relations remain explicit. A resolved hard membership is keyed by element, source occurrence, and actual matched subinterval. One parent occurrence can contribute multiple local memberships. `resolved_local` protein correspondence can provide a hard local observation only when its separate eligibility field is true; candidate, unanchored, ambiguous, and unresolved records create no hard observation.

Short-DNA evidence has two stages. A feature-bounded alignment is a candidate only. When two flanks map uniquely, in order, and on the same transcript path, the intervening genomic interval is extracted and realigned. Only this anchor-bounded result can establish a hard nucleotide membership or position observation. Missing, reversed, cross-path, non-unique, or incompletely enumerated flanks retain candidate status.

## Structural site matrix, schema version 3

The matrix has exactly one row for every `family_id`, `layer`, `site_id`, and tree-tip combination. Its schema is versioned in `schema_version`. The required state labels are stored literally in `state_0` and `state_1`; `state` contains one of those labels or `unknown`. An unobserved tip remains present as `state=unknown`, `observation_mask=missing`, with a non-empty `observation_reason`.

| Field | Meaning |
|---|---|
| `family_id`, `site_id`, `species` | Family, homologous structural character, and tree-tip observation |
| `layer` | `exon_presence`, `exon_role`, or `splice_junction` |
| `state`, `state_0`, `state_1` | Observed state, its layer-specific zero label, and its layer-specific one label |
| `evidence` | Human-readable evidence summary retained for compatibility |
| `schema_version` | Matrix schema version; current value is `3` |
| `annotation_view` | `view_independent` for sequence presence, or the selected `repertoire`/`canonical` transcript view for role and junction observations |
| `applicability` | `applicable`, `inapplicable`, or `undetermined` for the queried character |
| `observation_reason` | Why the state is observed, unknown, or inapplicable |
| `transcript_scope` | Whether the row summarizes a gene-locus repertoire or a transcript-level scope |
| `parent_feature_ids` | Original annotation features underlying the observation |
| `member_interval_ids` | Sequence intervals/members assigned to the structural unit |
| `evidence_ids` | Traceable alignment, annotation, or candidate records supporting the row |
| `discovery_rule` | Rule by which the site entered the candidate character set |
| `discovery_species` | Species or source record that caused discovery, when applicable |
| `observation_mask` | Whether the tip is observed or missing for this site |
| `linked_group_id` | Related sites that may be structurally dependent, such as cuts in one 1:n pattern |
| `site_kind` | Character subtype, such as element presence or within-exon boundary |
| `observation_source` | Origin of the observation, normally genome sequence and supplied annotation |
| `annotation_completeness` | Recorded completeness assessment for the local annotation evidence |
| `confidence_flag` | Qualitative evidence status; not a posterior probability |
| `conclusion_flag` | Qualitative status/reason for the observation; not a statistical test result |

The three layer definitions are:

| Layer | `state_0` | `state_1` |
|---|---|---|
| `exon_presence` | `absent` | `present` |
| `exon_role` | `not_exonic` in the selected supplied-annotation view | `exonic` in the selected supplied-annotation view |
| `splice_junction` | `absent` in the selected transcript view | `present` in the selected transcript view |

Presence, role, and junction states remain separate. For example, aligned homologous DNA can be present while exon role is unknown or annotation-conditionally non-exonic. A prediction remains candidate evidence. Unknown state and inapplicability have distinct metadata; absent DNA does not become an exon-role zero.

Within one site, `state_0`, `state_1`, `linked_group_id`, `site_kind`, `discovery_rule`, and `annotation_view` are invariant across tips. Schema-v3 import rejects a site with a missing tree tip, an extra species, a duplicated tip row, inconsistent site labels, or inconsistent invariant metadata. The tree-tip set is taken from the supplied rooted tree.

The model input uses this single frozen matrix. `observation_mask` is `observed`, `missing`, `generated_all_zero`, or `explicit_all_zero`; recognized legacy missing-mask labels are consumed as missing observations. Other values are rejected. A missing observation always has `state=unknown`; generated or explicit all-zero rows represent declared character-universe observations and do not stand in for absent data. The observation mask and ascertainment rule determine which tips contribute to each conditional likelihood. Matrix sites remain present when all states are unknown; the unavailable reason is reported separately.

The default repertoire run writes the matrix once. Parsimony, ER/ARD, and foreground analyses consume that exact file through `--structural-site-matrix`. A canonical sensitivity analysis generates a separate matrix with `annotation_view=canonical`; mixing repertoire and canonical rows or reusing one view under the other label is rejected.

## Foreground branch input

`--foreground-branches` accepts a tab-delimited table with either node IDs:

```text
parent_node  child_node
```

or labelled edges:

```text
branch_scope
parent_label -> child_label
```

Every row must match one edge in `species_tree.tsv`. An empty set, an unmatched edge, or a set containing every branch is rejected because it does not define an estimable foreground contrast.

## Coordinates and legacy data

GFF3/GTF coordinates and user-facing genomic `start/end` fields are 1-based closed intervals with `start <= end`. Alignment-internal intervals use 0-based half-open coordinates where their field names specify `start0/end0`. The gene's genomic strand and the relative strand of an alignment are separate quantities. A splice cut is a position between bases and is stored as a cut coordinate, not as an arbitrary exon-end base.

Legacy matrices are read through the compatibility adapter. Presence rows can recover their view from the layer definition. Legacy role and junction rows require an unambiguous `transcript_scope`; rows without one are rejected because canonical and repertoire observations have different biological meanings. When an older correspondence table lacks actual matched-block coordinates, the block-level homology is unavailable; the full occurrence interval is not substituted. Coordinate formats for each serialized block field are documented alongside that field.

## Analysis outputs and boundaries

Default qualitative outputs include `structural_site_matrix.tsv`, `node_structural_states.tsv`, `branch_structural_events.tsv`, `structural_site_summary.tsv`, and `compound_structural_events.tsv`. Correspondence outputs include `element_correspondence.tsv`, `splice_boundary_correspondence.tsv`, `observed_element_tree_coverage.tsv`, and `observed_intragenic_paths.tsv`. `alignment_backend_report.tsv` records the selected alignment mode and backend; the exact executable version comes from container or execution provenance. `phylogeny_scope.tsv` records matrix scope, total and included site counts, annotation view, and exclusion or unavailability reasons. Optional CTMC outputs include `model_fits.tsv`, `model_tests.tsv`, `node_state_posteriors.tsv`, and `branch_transition_posteriors.tsv`.

For a `1↔2` junction contrast, `branch_structural_events.tsv` reports intron gain/loss in `event_type` and split/fusion in `structural_relation`. For `1↔n` with `n>=3`, `compound_structural_events.tsv` reports `required` only when all cutpoints co-occur on the same species, gene copy, and transcript path; otherwise the summary is `possible_non_joint`.

The first-version analysis is intended for supplied single-copy orthologs with at least some alignable gene-internal sequence and interpretable annotation paths. Large clusters of highly similar microexons may have unresolved instance mappings; affected intervals retain candidate evidence and an explicit unresolved status. Whole-gene multi-copy homology, transcript usage without transcriptome evidence, and regulatory or molecular-mechanism inference are outside this input model. Such features can remain in descriptive records when supplied, while formal inference is limited to sites with defensible correspondence and state definitions.

See [Method Overview](method.md) for correspondence and role semantics and [Statistical Model](statistical_model.md) for parsimony, optional CTMC models, and ascertainment.
