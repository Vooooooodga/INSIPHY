# INSIPHY Roadmap

## Current v0.4.0 Status

INSIPHY v0.4.0 is a runnable method prototype. It includes:

- genome FASTA + GFF/GTF extraction for target gene copies, transcripts,
  intron sites and transcript paths;
- segment occurrence, sequence, adjacency, copy-context and copy-relationship
  tables;
- graph-based homologous segment grouping, reciprocal correspondence scoring
  and membership calls;
- source/copy-role propagation for source/background/derived-copy cases;
- fixed-tree reconstruction of intragenic structural characters with Sankoff
  and branch-length-aware CTMC/Mk scoring;
- invariant-model likelihood-ratio tests with p values, fitted rates, AIC and
  BIC for structural characters;
- simple baselines, simulated benchmarks and annotation-dropout negative
  controls;
- jingwei and Sdic curated micro-demos;
- unit tests for demos and genome/annotation preprocessing.

## Publication-Grade Milestones

1. **Real-data case studies**
   - Build accession-level jingwei and Sdic analyses from curated genome and
     annotation versions.
   - Record all genome/annotation sources and target gene identifiers.

2. **Statistical calibration**
   - Calibrate CTMC/Mk event support with simulation under known histories.
   - Compare rate models for segment presence, role state, adjacency, source
     mixture and copy multiplicity.
   - Add parametric bootstrap for invariant-model LRT p values.
   - Add stochastic character mapping or an explicitly labeled endpoint
     posterior approximation for branch event placement.
   - Add foreground/background branch-set tests for user-defined evolutionary
     hypotheses.

3. **Benchmarking**
   - Compare against annotation-only, exon-orthology-only and simple
     intron/exon gain-loss baselines.
   - Report event recall, precision, branch placement accuracy and robustness
     to incomplete annotation.
   - Build a composite benchmark because there is no mature public benchmark
     for branch-level intragenic synteny histories.

4. **Method robustness**
   - Stress-test hidden-segment detection, shifted splice sites and joined
     exons on accession-level data.
   - Expand handling of tandem duplicates, many-to-many segment mappings,
     processed pseudogene candidates and ambiguous paralogy.
   - Track splice-boundary shifts, intron sliding, tandem exon duplication,
     transposable-element exonization and gene conversion as explicit
     alternatives where evidence is insufficient.

5. **Release readiness**
   - Add continuous integration after GitHub upload.
   - Add example data provenance, command transcripts and versioned releases.
   - Prepare manuscript figures and result tables from real case studies.
   - Keep the distributed software as a CLI/Python package. Server-side
     Nextflow/Slurm workflows are internal formal run records and should stay
     outside the package.
