# INSIPHY Input Format

INSIPHY can start from genome FASTA plus GFF/GTF annotation, or from prepared
TSV tables.

## FASTA/GFF Extraction

`extract-gene` extracts one gene copy from genome sequence and annotation:

```bash
PYTHONPATH=src python3 -m insiphy.cli extract-gene \
  --genome genome.fa \
  --annotation annotation.gff3 \
  --gene-id GeneA \
  --family-id family_a \
  --species SpeciesA \
  --gene-copy-id SpeciesA_GeneA \
  --output-dir work/family_a \
  --transcript-policy canonical
```

Run it once per species/copy. Use `--append` after the first copy to add more
species or paralogous copies to the same directory. `--transcript-policy all`
retains all annotated transcript paths; the default `canonical` uses the
longest-CDS transcript as the copy summary.

`derive-tables` creates first-pass homology, match, adjacency and copy-context
tables from extracted segments:

```bash
PYTHONPATH=src python3 -m insiphy.cli derive-tables \
  --input-dir work/family_a \
  --identity-threshold 0.7 \
  --distance-table species_distance.tsv
```

The optional distance table may contain `species1`, `species2` and `distance`
fields. Distance classes `short`, `medium` and `long` adjust correspondence
stringency in the ExOrthist style.

## Required Tables For Inference

- `species_tree.tsv`: rooted species tree with `node_id`, `parent_id`, `label`.
- `segment_occurrences.tsv`: observed or inferred gene-internal segments.
- `segment_homology.tsv`: homologous segment group membership.
- `physical_adjacencies.tsv`: copy-specific neighboring segment pairs.
- `segment_matches.tsv`: pairwise segment correspondence scores.
- `copy_context.tsv`: copy status per species and copy.
- `transcript_paths.tsv`: transcript-specific ordered paths through segment
  occurrences.
- `intron_sites.tsv`: intron interval, phase and splice motif evidence.
- `copy_relationships.tsv`: copy-pair geometry and tandem/dispersed calls.
- `sequence_synteny_evidence.tsv`: annotation support, hidden candidates or
  annotation conflicts.

The demo directories provide complete examples for all tables.

## Real Case Manifest

`build-case` reads a tab-delimited manifest with at least these fields:

```text
case_id	species	family_id	gene_id	gene_copy_id	genome_fasta	annotation_file
```

Recommended additional fields are:

```text
assembly	annotation	source_url	release	notes
```

The command writes:

- `case_provenance.tsv`
- `case_build_report.tsv`
- extracted INSIPHY input tables

Example:

```bash
PYTHONPATH=src python3 -m insiphy.cli build-case \
  --manifest examples/real_cases/jingwei/manifest.tsv \
  --species-tree examples/real_cases/jingwei/species_tree.tsv \
  --output-dir work/jingwei_case
```

`inspect-annotation` can be used before manifest finalization to discover gene
IDs from symbols or aliases:

```bash
PYTHONPATH=src python3 -m insiphy.cli inspect-annotation \
  --annotation annotation.gff3 \
  --alias-file examples/real_cases/jingwei/gene_aliases.tsv \
  --species Drosophila_yakuba \
  --case-id jingwei \
  --output-dir work/jingwei_inspect
```

`scan-hidden-segments` performs gapped local alignment of source segments
against a target gene interval FASTA and reports identity, coverage, CIGAR-like
alignment, splice motif score and frame status. It is intended for curated
local intervals or flank windows, not whole-genome searches.
