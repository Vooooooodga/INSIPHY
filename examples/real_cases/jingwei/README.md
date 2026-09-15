# jingwei Real Case Template

This template records the expected INSIPHY inputs for an accession-level
jingwei case study. Fill the `TBD` fields after genome and annotation files are
curated on the server.

Target biology:

- derived `jgw/jingwei` copy in the yakuba/teissieri lineage;
- source relationship to `Adh` and `yande/yellow emperor`;
- candidate sequence-supported hidden or role-shifted segments.

Suggested steps:

```bash
PYTHONPATH=src python -m insiphy.cli inspect-annotation \
  --annotation /data/db/genome/TBD/annotation.gff3.gz \
  --alias-file examples/real_cases/jingwei/gene_aliases.tsv \
  --species Drosophila_yakuba \
  --case-id jingwei \
  --output-dir work/jingwei_inspect

PYTHONPATH=src python -m insiphy.cli build-case \
  --manifest examples/real_cases/jingwei/manifest.tsv \
  --species-tree examples/real_cases/jingwei/species_tree.tsv \
  --output-dir work/jingwei_case
```
