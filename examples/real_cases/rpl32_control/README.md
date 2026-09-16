# RpL32 Conserved Control

This accession-level control uses the local NCBI RefSeq Drosophila genome
FASTA/GFF files under `/data/db/genome`.

Expected signal:

- conserved exon sequence presence and exonic role;
- near-zero branch transition probabilities;
- no ER/ARD P value when all reliable sites are invariant.

Example:

```bash
PYTHONPATH=src python3 -m insiphy.cli build-case \
  --manifest examples/real_cases/rpl32_control/manifest.tsv \
  --species-tree examples/real_cases/rpl32_control/species_tree.tsv \
  --output-dir work/rpl32_control_case \
  --aligner auto \
  --threads 4

PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir work/rpl32_control_case \
  --output-dir results/rpl32_control \
  --analysis-scope single-copy \
  --model er-ard \
  --branch-length-mode supplied \
  --threads 4

PYTHONPATH=src python3 -m insiphy.cli visualize \
  --input-dir work/rpl32_control_case \
  --result-dir results/rpl32_control \
  --output-dir results/rpl32_control_figures
```
