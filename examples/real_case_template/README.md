# Real Case Template

This directory documents how a real INSIPHY case study should be recorded
without committing large genome files to the repository.

Use one row per species and gene copy in `manifest.tsv`. Genome FASTA and
annotation files should live in the project or database area on the server, not
inside the GitHub repository.

The expected workflow is:

1. curate species, assembly version, annotation version and target gene IDs;
2. add or review `species_tree.tsv`;
3. run `insiphy build-case --manifest manifest.tsv --species-tree species_tree.tsv`;
4. run `insiphy run` with bootstrap, stochastic maps and foreground branches
   when needed;
5. run `insiphy visualize`;
6. review `candidate_structural_events.tsv`, `event_support_summary.tsv`,
   `case_summary.tsv` and `baseline_comparison.tsv`.
