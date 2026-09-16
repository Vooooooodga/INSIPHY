# RpL32 Conserved Control

This accession-level control uses the local NCBI RefSeq Drosophila genome
FASTA/GFF files under `/data/db/genome`.

Expected signal:

- conserved internal segment structure;
- low source-mixture support;
- low copy-expansion support;
- no high-confidence exonization or source-joining event.

Example:

```bash
PYTHONPATH=src python3 -m insiphy.cli build-case \
  --manifest examples/real_cases/rpl32_control/manifest.tsv \
  --species-tree examples/real_cases/rpl32_control/species_tree.tsv \
  --output-dir work/rpl32_control_case \
  --aligner minimap2 \
  --threads 4

PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir work/rpl32_control_case \
  --output-dir results/rpl32_control \
  --bootstrap-replicates 200 \
  --stochastic-maps 200 \
  --seed 7

PYTHONPATH=src python3 -m insiphy.cli visualize \
  --input-dir work/rpl32_control_case \
  --result-dir results/rpl32_control \
  --output-dir results/rpl32_control_figures
```
