# Real Case Template

This directory documents how a real INSIPHY case study should be recorded
without committing large genome files to the repository.

Use one row per species and gene copy in `manifest.tsv`. Genome FASTA and
annotation files should live in the project or database area on the server, not
inside the GitHub repository.

The expected workflow is:

1. curate species, assembly version, annotation version and target gene IDs;
2. run `insiphy extract-gene` once per species/copy;
3. run `insiphy derive-tables`;
4. add or review `species_tree.tsv`;
5. run `insiphy run`;
6. review `candidate_structural_events.tsv`, `model_comparison.tsv`,
   `baseline_comparison.tsv` and `benchmark_summary.tsv` if simulated truth is
   available.
