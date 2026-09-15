# INSIPHY Input Format

INSIPHY can start from genome FASTA plus GFF/GTF annotation, or from prepared
TSV tables.

## FASTA/GFF Extraction

`extract-gene` extracts one gene copy from genome sequence and annotation:

```bash
PYTHONPATH=src python -m insiphy.cli extract-gene \
  --genome genome.fa \
  --annotation annotation.gff3 \
  --gene-id GeneA \
  --family-id family_a \
  --species SpeciesA \
  --gene-copy-id SpeciesA_GeneA \
  --output-dir work/family_a
```

Run it once per species/copy. Use `--append` after the first copy to add more
species or paralogous copies to the same directory.

`derive-tables` creates first-pass homology, match, adjacency and copy-context
tables from extracted segments:

```bash
PYTHONPATH=src python -m insiphy.cli derive-tables \
  --input-dir work/family_a \
  --identity-threshold 0.7
```

## Required Tables For Inference

- `species_tree.tsv`: rooted species tree with `node_id`, `parent_id`, `label`.
- `segment_occurrences.tsv`: observed or inferred gene-internal segments.
- `segment_homology.tsv`: homologous segment group membership.
- `physical_adjacencies.tsv`: copy-specific neighboring segment pairs.
- `segment_matches.tsv`: pairwise segment correspondence scores.
- `copy_context.tsv`: copy status per species and copy.
- `sequence_synteny_evidence.tsv`: annotation support, hidden candidates or
  annotation conflicts.

The demo directories provide complete examples for all tables.
