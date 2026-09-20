# RpL32 v0.11 Real-Data Demo

This directory contains the formal single-copy outputs generated from the five
Drosophila assemblies listed in `examples/real_cases/rpl32_control/manifest.tsv`.

## Result

- `EG_0001` and `EG_0003` each contain one core exon-like member from all five
  species.
- `EG_0006` is a short D. melanogaster-only interval with ambiguous membership;
  its formal tip states remain unknown.
- The two core EGs are present and exonic in every species.
- `JG_EG_0001__EG_0003` is present in every species.
- No analyzed site varies among observed tips.
- ER/ARD and foreground model tests report
  `test_status=parameters_not_estimable` and `p_value=NA`.
- The largest endpoint-change posterior is approximately `4.2e-13`.

`foreground_model_fits.tsv` and `foreground_model_tests.tsv` use the terminal
D. melanogaster branch as the prespecified foreground. The invariant structure
contains no information for estimating a foreground rate difference.

The `figures` directory contains the default colorblind-aware exon-synteny,
phylogeny, and integrated SVG outputs.
