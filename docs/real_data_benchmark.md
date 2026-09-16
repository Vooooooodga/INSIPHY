# Real-Data Demonstration

## RpL32 conserved control

The first formal single-copy demonstration uses RpL32 orthologs from:

- Drosophila melanogaster, GCF_000001215.4;
- Drosophila simulans, GCF_016746395.2;
- Drosophila erecta, GCF_003286155.1;
- Drosophila yakuba, GCF_016746365.2;
- Drosophila teissieri, GCF_016746235.2.

Genome FASTA and GFF3 files are stored outside the repository under
`/data/db/genome`. The repository contains accession-level manifests and the
species tree.

## Complete user path

```bash
insiphy build-case \
  --manifest examples/real_cases/rpl32_control/manifest.tsv \
  --species-tree examples/real_cases/rpl32_control/species_tree.tsv \
  --output-dir work/rpl32_control \
  --aligner internal \
  --threads 4

insiphy run \
  --input-dir work/rpl32_control \
  --output-dir results/rpl32_control \
  --analysis-scope single-copy \
  --model er-ard \
  --branch-length-mode supplied \
  --threads 4

insiphy visualize \
  --input-dir work/rpl32_control \
  --result-dir results/rpl32_control \
  --output-dir figures/rpl32_control
```

## v0.11 result

The correspondence stage recovered two core exon-like groups spanning all five
species. Their `exon_presence` states are all present, and their
`exon_role` states are all exonic. A short D. melanogaster interval forms a
low-support group and remains unknown in the formal matrix.

No analyzed site varies among observed tips. Gain/loss asymmetry is therefore
not identifiable. `model_tests.tsv` reports
`test_status=parameters_not_estimable` and `p_value=NA` for these layers.
The ER-based branch endpoint-change probabilities are near zero.

This is the expected qualitative result for a conserved control: strong
correspondence, stable structure, and no supported branch-specific change. It
does not measure sensitivity to true exon gains or losses.

## Next real-data set

The next benchmark should contain single-copy genes with independently
documented structural changes. Candidate cases must satisfy:

- one ortholog per species;
- assembly and annotation versions available;
- enough species to distinguish alternative branch placements;
- sequence-level evidence for the altered exon or junction;
- a published history that can be reviewed independently of INSIPHY.

Jingwei and Sdic remain useful for future multi-copy development. Their
duplication histories place them outside the current formal single-copy
benchmark.
