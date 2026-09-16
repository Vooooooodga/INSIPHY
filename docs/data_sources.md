# Data Sources

INSIPHY real-case examples are Drosophila-focused because the selected gene
histories are biologically interpretable and the available species are close
enough for gene-internal synteny to be informative.

## jingwei

Expected biological pattern:

- a chimeric/duplicated gene involving `Adh`-derived sequence and `yande/ymp`
  source sequence;
- coding recruitment of sequence that is intronic or noncoding in the source
  context;
- source mixture and source-joining adjacency in the derived copy.

Selected literature background:

- Long and Langley, 1993, Science: origin of jingwei by exon shuffling.
- Long et al., 1999, Genetics: jingwei evolution and recruited sequence.
- Wang et al., 2000, Molecular Biology and Evolution:
  https://doi.org/10.1093/oxfordjournals.molbev.a026413
- Zhang et al., 2004, PNAS: https://doi.org/10.1073/pnas.0407066101

## Sdic

Expected biological pattern:

- a young duplicated/chimeric gene family involving `Annexin B10` and `sw`;
- multi-source gene structure in derived copies;
- copy-number expansion and intragenic adjacency changes.

Selected literature background:

- Nurminsky et al., 1998, Nature: https://doi.org/10.1038/25126
- Yeh et al., 2012, G3: https://pmc.ncbi.nlm.nih.gov/articles/PMC3277543/
- Zhao et al., 2023, Communications Biology:
  https://doi.org/10.1038/s42003-023-05427-4

## RpL32 Conserved Control

Expected biological pattern:

- conserved ribosomal protein gene;
- stable exon-intron structure among close Drosophila species;
- low support for source joining, copy expansion and exonization.

This control is used to estimate false-positive behavior in the first
accession-level benchmark.

## Local Genome Data

The R730 local data store currently contains:

- `/data/db/genome/Drosophila_melanogaster/GCF_000001215.4/`
- `/data/db/genome/Drosophila_simulans/GCF_016746395.2/`
- `/data/db/genome/Drosophila_erecta/GCF_003286155.1/`
- `/data/db/genome/Drosophila_yakuba/GCF_016746365.2/`
- `/data/db/genome/Drosophila_teissieri/GCF_016746235.2/`

Large FASTA/GFF/GTF files are excluded from the repository. Case manifests
record local paths, assembly, annotation source, release, source labels and
copy roles.
