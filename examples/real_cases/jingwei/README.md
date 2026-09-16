# jingwei Real Case

This accession-level case uses local Drosophila NCBI RefSeq genome FASTA/GFF
files under `/data/db/genome`.

Expected signal:

- derived `jgw/jingwei` copy in the yakuba/teissieri lineage;
- source relationship to `Adh` and `yande/ymp`;
- source joining and possible role-shifted segments in the derived copy.

Example:

```bash
PYTHONPATH=src python3 -m insiphy.cli build-case \
  --manifest examples/real_cases/jingwei/manifest.tsv \
  --species-tree examples/real_cases/jingwei/species_tree.tsv \
  --copy-tree examples/real_cases/jingwei/copy_tree.tsv \
  --output-dir work/jingwei_case \
  --aligner minimap2 \
  --threads 4

PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir work/jingwei_case \
  --output-dir results/jingwei \
  --bootstrap-replicates 200 \
  --stochastic-maps 200 \
  --foreground-branches examples/real_cases/jingwei/foreground_branches.tsv \
  --seed 7

PYTHONPATH=src python3 -m insiphy.cli visualize \
  --input-dir work/jingwei_case \
  --result-dir results/jingwei \
  --output-dir results/jingwei_figures
```
