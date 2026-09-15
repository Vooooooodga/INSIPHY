# Sdic Real Case Template

This template records the expected INSIPHY inputs for an accession-level Sdic
case study. Fill the `TBD` fields after genome and annotation files are curated
on the server.

Target biology:

- Sdic family in D. melanogaster;
- source relationship to `Annexin B10` and `sw`;
- copy-number ambiguity and intron-derived segment recruitment.

Suggested steps:

```bash
PYTHONPATH=src python -m insiphy.cli inspect-annotation \
  --annotation /data/db/genome/TBD/annotation.gff3.gz \
  --alias-file examples/real_cases/sdic/gene_aliases.tsv \
  --species Drosophila_melanogaster \
  --case-id sdic \
  --output-dir work/sdic_inspect

PYTHONPATH=src python -m insiphy.cli build-case \
  --manifest examples/real_cases/sdic/manifest.tsv \
  --species-tree examples/real_cases/sdic/species_tree.tsv \
  --output-dir work/sdic_case
```
