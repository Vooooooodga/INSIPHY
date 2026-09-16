# Sdic Real Case

This accession-level case uses local Drosophila NCBI RefSeq genome FASTA/GFF
files under `/data/db/genome`.

Expected signal:

- derived Sdic copies in D. melanogaster;
- source relationship to `Annexin B10` and `sw`;
- copy expansion, source joining and intragenic adjacency changes.

Example:

```bash
PYTHONPATH=src python3 -m insiphy.cli build-case \
  --manifest examples/real_cases/sdic/manifest.tsv \
  --species-tree examples/real_cases/sdic/species_tree.tsv \
  --copy-tree examples/real_cases/sdic/copy_tree.tsv \
  --output-dir work/sdic_case \
  --aligner minimap2 \
  --threads 4

PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir work/sdic_case \
  --output-dir results/sdic \
  --bootstrap-replicates 200 \
  --stochastic-maps 200 \
  --foreground-branches examples/real_cases/sdic/foreground_branches.tsv \
  --seed 7

PYTHONPATH=src python3 -m insiphy.cli visualize \
  --input-dir work/sdic_case \
  --result-dir results/sdic \
  --output-dir results/sdic_figures
```
