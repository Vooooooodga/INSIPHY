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
- HSGs conserved, gained or role-shifted;
- branch with strongest event support;
- LRT p value, q value, empirical bootstrap p value and fitted CTMC rate;
- stochastic-map `Pr(any change)` and expected change count on key branches;
- alternative explanations: annotation dropout, fragmented assembly, paralogy
  ambiguity, gene conversion or weak sequence support.

## Acceptance For v0.6 Real Run

A v0.6 real run is acceptable when it can:

- run from accession-level genome FASTA/GFF plus a manifest and species tree;
- recover expected qualitative event classes for `jingwei` and `Sdic`;
- include the RpL32 negative control or a clearly documented equivalent
  OrthoFinder-supported conserved control;
- save all INSIPHY output tables needed for biological interpretation;
- generate black-and-white-readable, colorblind-friendly synteny and event-map
  figures;
- document unresolved uncertainty.

Publication-level claims still require a larger case set and quantified
false-positive rate.
