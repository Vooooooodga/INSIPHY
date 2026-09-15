# Real-Data Benchmark Plan

The first real-data benchmark should focus on Drosophila young or duplicate
genes because the expected histories are well described in the literature and
the species are close enough for gene-internal synteny to be informative.

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

3. **sphinx or another Drosophila young duplicate/chimeric gene**
   - Expected biology: a published young gene with clear genome-sequence
     origin and close-species distribution.
   - INSIPHY target signal: independent validation outside the two bundled
     examples.

4. **Conserved negative control**
   - Candidate genes: stable, single-copy Drosophila housekeeping genes such as
     `RpL32` or `Act5C`, chosen only after confirming orthology and annotation
     quality in the selected genome set.
   - INSIPHY target signal: low false-positive rate for segment gain/loss,
     exonization and source mixing.

## Data Sources

Recommended primary sources:

- FlyBase genome FASTA and GFF/GTF releases for Drosophila species.
- NCBI or Ensembl Metazoa for species or assemblies missing from FlyBase.
- Published supplementary tables for duplicated-gene exon-intron divergence
  studies when a specific independent benchmark case is selected.

Large genome and annotation files should remain outside the GitHub repository.
Case manifests in `examples/real_cases/` should record local paths, source
release, assembly, annotation version, gene IDs, source labels and copy roles.

## Benchmark Questions

For each real case, report:

- which supplied copies and source loci were analyzed;
- which hidden segments were sequence-supported;
- which HSGs were conserved, gained or role-shifted;
- which branch has the highest posterior support for the event;
- LRT p value, empirical bootstrap p value and fitted CTMC rate;
- stochastic-map `Pr(any change)` and expected change count on key branches;
- alternative explanations: annotation dropout, fragmented assembly, paralogy
  ambiguity, gene conversion or weak sequence support.

## Acceptance For v0.5 Real Demo

A v0.5 real demo is acceptable when it can:

- run from a case manifest using only genome sequence, annotation and a species
  tree;
- reproduce the expected qualitative event class for jingwei and Sdic;
- include at least one negative-control run with no strong source-mixing or
  exonization signal;
- save all INSIPHY output tables needed for biological interpretation;
- document data sources and unresolved uncertainty.

This benchmark remains a first real-world evaluation. Publication-level claims
still require a larger case set and a quantified false-positive rate.
