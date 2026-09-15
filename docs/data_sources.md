# Data Sources And Demo Scope

This package currently ships two curated micro-demos. They are small method
fixtures built from published biological patterns and hand-written toy
sequences, so they exercise the software logic without claiming to reproduce a
complete accession-level analysis.

## jingwei demo

Biological pattern represented:

- a duplicated gene with segments derived from `Adh` and `yande`;
- coding recruitment of a segment that is intronic or noncoding in the source
  context;
- a hidden segment candidate supported by sequence and local synteny evidence.

Literature background:

- Long and Langley, 1993, Science: origin of jingwei by exon shuffling.
- Long et al., 1999, Genetics: jingwei evolution and recruited sequence.
- Wang et al., 2000, Molecular Biology and Evolution:
  https://doi.org/10.1093/oxfordjournals.molbev.a026413
- Zhang et al., 2004, PNAS: https://doi.org/10.1073/pnas.0407066101

## Sdic demo

Biological pattern represented:

- a duplicated gene family derived from `Annexin B10` and `sw`;
- multi-source gene structure in the derived copy;
- an intron-derived segment whose coding/exonic role is inferred from sequence
  and structure evidence;
- copy-number ambiguity in the derived species.

Literature background:

- Nurminsky et al., 1998, Nature: https://doi.org/10.1038/25126
- Yeh et al., 2012, G3: https://pmc.ncbi.nlm.nih.gov/articles/PMC3277543/
- Zhao et al., 2023, Communications Biology:
  https://doi.org/10.1038/s42003-023-05427-4

## Method Assumptions In The Demo Files

The demo TSVs separate three evidence layers:

1. observed or sequence-completed segment occurrences;
2. homologous segment group assignments and pairwise correspondence support;
3. fixed species-tree structural characters for presence, role, adjacency,
   source mixture, and copy multiplicity.

This keeps annotation completion separate from phylogenetic inference. The
phylogenetic module can then compare whether a hidden segment is better treated
as an annotation gap or as a real structural loss/gain, and whether a gene copy
is better explained by independent segment changes or by one compound
chimeric/copy event.

## Real-Data Benchmark Direction

The first real-data benchmark should remain Drosophila-focused. Jingwei and
Sdic are the positive controls already represented by bundled method fixtures.
The next case should add an independent young duplicate or chimeric gene such
as sphinx, followed by a conserved single-copy negative control. The benchmark
plan and acceptance criteria are recorded in `docs/real_data_benchmark.md`.

Large FASTA/GFF/GTF files are intentionally excluded from the repository.
Manifests should point to local project inputs and record source release,
assembly, annotation, source labels and copy roles.
