# Real-Data Benchmark Plan

The first benchmark uses accession-level Drosophila genome FASTA and GFF files.
The cases are young or duplicated genes with described structural histories,
plus one conserved negative control.

## Priority Cases

1. **jingwei**
   - Expected biology: chimeric new gene involving `Adh`-derived sequence and
     `yande/ymp` source sequence, with coding recruitment of previously
     noncoding or intronic sequence.
   - INSIPHY target signal: source mixture, source-joining adjacency and
     role-state change on the branch where the derived copy appears.

2. **Sdic**
   - Expected biology: young duplicated/chimeric gene family involving
     `Annexin B10` and `sw`, with copy-number ambiguity and derived coding
     structure.
   - INSIPHY target signal: multi-source structure, copy expansion and
     intragenic adjacency changes.

3. **RpL32 conserved control**
   - Expected biology: conserved ribosomal protein gene with stable
     exon-intron organization across close Drosophila species.
   - INSIPHY target signal: low support for source mixing, exonization and
     structural novelty.

## Available Local Data

The current R730 data store contains the following complete genome/annotation
pairs:

- `/data/db/genome/Drosophila_melanogaster/GCF_000001215.4/`
- `/data/db/genome/Drosophila_simulans/GCF_016746395.2/`
- `/data/db/genome/Drosophila_erecta/GCF_003286155.1/`
- `/data/db/genome/Drosophila_yakuba/GCF_016746365.2/`
- `/data/db/genome/Drosophila_teissieri/GCF_016746235.2/`

The repository stores only manifests, species trees and provenance notes. Large
genome FASTA and GFF files stay under `/data/db/genome`.

## Full User Scenario

Each real case should exercise:

- `build-case` from real FASTA/GFF manifest;
- `derive-tables` with `--threads` and at least one external local aligner when
  available;
- `run` with CTMC/LRT, BH q values, bootstrap, stochastic maps and foreground
  branches;
- `visualize` with colorblind-friendly SVG output;
- `benchmark` for curated positive cases and negative control expectations.

## Benchmark Questions

For each case, report:

- supplied copies and source loci;
- hidden or shifted segments supported by sequence;
- exon-like correspondence groups conserved, gained or role-shifted;
- branch with strongest event support;
- core structural, copy-context and ambiguous-evidence call counts;
- LRT p value, q value, empirical bootstrap p value and fitted CTMC rate;
- stochastic-map `Pr(any change)` and expected change count on key branches;
- alternative explanations: annotation dropout, fragmented assembly, paralogy
  ambiguity, high-similarity paralogous segments or weak sequence support.

## v0.10 Demo Run

The current package was exercised on real local Drosophila inputs under:

```text
/scratch/projects/intragenic_structure/insiphy_v010_real_demo
```

Settings:

- `build-case` from accession-level genome FASTA/GFF manifests.
- `--threads 2` with the internal aligner, because minimap2/miniprot were not
  available on `PATH` in this environment.
- `run` with `--bootstrap-replicates 10`, `--stochastic-maps 10` and foreground
  branch files for positive cases.
- `visualize` with default color encoding; RpL32 was also rendered with
  pattern encoding.
- `benchmark` against curated truth tables for jingwei and Sdic.

Observed results:

- **RpL32 control**: no hidden-segment candidates, no source-join candidates,
  no multi-source tips, no branch-event candidates and no copy-context
  candidates. This is the desired conserved-control behavior.
- **jingwei**: structural characters used `copy_tree.tsv`; one core event was
  called, `chimeric_origin_or_source_mixing`, on
  `adh_yak_dup->Drosophila_yakuba:Dyak_jgw`. The LRT p value was
  `8.43996e-06`, BH q value `0.000118159`, CTMC branch-change probability
  `0.973651`, and stochastic-map `Pr(any change)` was `1` with 10 maps.
  Benchmark against the copy-tree truth table gave precision, recall and branch
  accuracy of `1`.
- **Sdic**: structural characters used `copy_tree.tsv`; the expected
  `chimeric_origin_or_source_mixing` event was called on
  `sdic_root->sdic_cluster` with branch accuracy `1`. One additional
  `coding_or_exonic_role_loss` core call was made on an AnxB10 lineage, so core
  benchmark precision was `0.5` and recall was `1`. This extra call is retained
  as an unresolved real-data signal requiring closer annotation and orthology
  review.
- **Calibration smoke run**: two replicates each of `negative_control`,
  `exonization` and `source_join` completed. The negative control had
  false-positive rate `0`; exonization and source-join scenarios had
  mean precision/recall/F1 of `1` in this tiny smoke run. These values show
  command behavior, not final publication calibration.

Important limitations:

- The internal aligner is a fallback; publication runs should use minimap2 or
  another documented local aligner for scalable segment matching.
- Only 10 bootstrap/stochastic-map replicates were used for the real demo.
  Empirical p values therefore have coarse Monte Carlo resolution.
- The Sdic extra role-loss call shows that curated truth tables must separate
  expected focal events from additional lineage-specific structural variation.

## Acceptance For v0.10 Real Run

A v0.10 real run is acceptable when it can:

- run from accession-level genome FASTA/GFF plus a manifest and species tree;
- recover expected qualitative event classes for `jingwei` and `Sdic`;
- include the RpL32 negative control or a clearly documented equivalent
  OrthoFinder-supported conserved control;
- save all INSIPHY output tables needed for biological interpretation;
- generate black-and-white-readable, colorblind-friendly synteny and event-map
  figures with one correspondence encoding mode per figure;
- document unresolved uncertainty.

Publication-level claims still require a larger case set and quantified
false-positive rate.
