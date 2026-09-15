# INSIPHY Roadmap

## Current v0.2.0 Status

INSIPHY v0.2.0 is a runnable method prototype. It includes:

- genome FASTA + GFF/GTF extraction for target gene copies, transcripts,
  intron sites and transcript paths;
- segment occurrence, sequence, adjacency, copy-context and copy-relationship
  tables;
- graph-based homologous segment grouping, reciprocal correspondence scoring
  and membership calls;
- fixed-tree reconstruction of intragenic structural characters with Sankoff
  and branch-length-aware CTMC/Mk scoring;
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

3. **Benchmarking**
   - Compare against annotation-only, exon-orthology-only and simple
     intron/exon gain-loss baselines.
   - Report event recall, precision, branch placement accuracy and robustness
     to incomplete annotation.

4. **Method robustness**
   - Stress-test hidden-segment detection, shifted splice sites and joined
     exons on accession-level data.
   - Expand handling of tandem duplicates, many-to-many segment mappings,
     processed pseudogene candidates and ambiguous paralogy.

5. **Release readiness**
   - Add continuous integration after GitHub upload.
   - Add example data provenance, command transcripts and versioned releases.
   - Prepare manuscript figures and result tables from real case studies.
