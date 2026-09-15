# INSIPHY Roadmap

## Current v0.1.0 Status

INSIPHY v0.1.0 is a runnable method scaffold. It includes:

- genome FASTA + GFF/GTF extraction for target gene copies;
- first-pass segment occurrence, sequence, adjacency and copy-context tables;
- homologous segment grouping and correspondence scoring;
- fixed-tree reconstruction of intragenic structural characters;
- likelihood-like character scores, simple baselines and simulated benchmarks;
- jingwei and Sdic curated micro-demos;
- unit tests for demos and genome/annotation preprocessing.

## Publication-Grade Milestones

1. **Real-data case studies**
   - Build accession-level jingwei and Sdic analyses from curated genome and
     annotation versions.
   - Record all genome/annotation sources and target gene identifiers.

2. **Statistical calibration**
   - Add likelihood-style transition models for segment presence, role state,
     adjacency, source mixture and copy multiplicity.
   - Calibrate event support with simulation under known histories.

3. **Benchmarking**
   - Compare against annotation-only, exon-orthology-only and simple
     intron/exon gain-loss baselines.
   - Report event recall, precision, branch placement accuracy and robustness
     to incomplete annotation.

4. **Method robustness**
   - Improve hidden-segment detection with local realignment and splice motif
     scoring.
   - Add explicit handling of tandem duplicates, many-to-many segment mappings
     and ambiguous paralogy.

5. **Release readiness**
   - Add continuous integration after GitHub upload.
   - Add example data provenance, command transcripts and versioned releases.
   - Prepare manuscript figures and result tables from real case studies.
