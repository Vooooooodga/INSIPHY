# Data Sources

The existing demonstration resources include Drosophila and yeast. The
mammalian cases below are literature-backed candidates for focused real-data
tests of intron loss/exon joining, splice-boundary displacement and within-gene
repeated-exon homology. They are resource pointers and do not establish that
files have been downloaded, event coordinates have been transferred to current
assemblies, or an upstream analysis has established single-copy membership.
The v0.15 five-case formal run remains pending; case roles and acceptance fields
are defined in [Real-Data Demonstration](real_data_benchmark.md).

## Deferred Experimental-Multicopy Resources

Jingwei and Sdic involve gene duplication, chimeric origin, or copy-number
expansion. They remain resources for `experimental-multicopy` development and
are excluded from the formal v0.15 single-copy benchmark.

### jingwei

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

### Sdic

Expected biological pattern:

- a young duplicated/chimeric gene family involving `Annexin B10` and `sw`;
- multi-source gene structure in derived copies;
- copy-number expansion and intragenic adjacency changes.

Selected literature background:

- Nurminsky et al., 1998, Nature: https://doi.org/10.1038/25126
- Yeh et al., 2012, G3: https://pmc.ncbi.nlm.nih.gov/articles/PMC3277543/
- Zhao et al., 2023, Communications Biology:
  https://doi.org/10.1038/s42003-023-05427-4

## Formal Single-Copy Benchmark Sources

### RpL32 Conserved Control

Expected biological pattern:

- conserved ribosomal protein gene;
- stable exon-intron structure among close Drosophila species;
- low support for source joining, copy expansion and exonization.

This case is a qualitative conservation control. It records whether the
observed, adequately covered structure remains conserved and does not estimate
a false-positive rate.

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

## T6 Mammalian Case Sources

### First-party biological and method papers

| Case | Source and exact location | What it supports | What it does not establish |
|---|---|---|---|
| FDPS | [Sharma et al. 2017, CESAR 2.0, Fig. 1E and Results](https://academic.oup.com/bioinformatics/article/33/24/3985/4095639), DOI `10.1093/bioinformatics/btx527` | Human exons 5-7 correspond to one large rat exon after two introns are absent in the rat comparison; mouse is shown without this joined structure | Independent branch placement on a broader tree, current GFF coordinates, or the proposed molecular mechanism |
| FDPS | [Coulombe-Huntington & Majewski 2007](https://doi.org/10.1101/gr.5703406), [Exalign comparative table](https://academic.oup.com/nar/article/36/8/e47/2410452) | Mammalian intron-loss context; the Exalign table records human `NM_002004`, mouse `NM_134469`, rat `NM_031840` for FDPS | Versioned current transcript IDs or lift-over to selected RefSeq assemblies; paper intron numbering should be reconciled with exon flanks |
| MAMSTR | [Sharma et al. 2017, Fig. 1D and Results](https://academic.oup.com/bioinformatics/article/33/24/3985/4095639), DOI `10.1093/bioinformatics/btx527` | Cow MAMSTR exon 3 has a reported GG acceptor mutation and a 30 bp acceptor shift; this is an example of comparative splice-site recovery | Exact genomic accession/coordinates, independent RNA validation, or a uniquely resolved branch history |
| PKM | [Wang et al. 2012, Figs. 1A-C](https://doi.org/10.1093/jmcb/mjr030) | Human exon 9 and 10 sequence homology; M1/M2 mutually exclusive paths; human minigene and RT-PCR evidence | Origin, timing or branch placement of the exon duplication across mammals |
| PKM | [Echigoya et al. 2008, equine M1/M2 study](https://pubmed.ncbi.nlm.nih.gov/18602015/) | Independent evidence for M1/M2 forms in horse | Cross-species genomic homology map or gene copy assessment |

The exon 9/10 pair in PKM is a within-gene duplicate relationship. It is
compatible with testing a single-copy orthologous gene family, provided each
species contributes one protein-coding PKM locus and the two internal exon
instances remain distinct in the structural representation.

### RefSeq assembly and gene-record candidates

NCBI assembly accessions below were selected from official current RefSeq
assembly/annotation records as viewed on 2026-09-18. The FASTA and GFF files
are available through each linked NCBI Datasets/Assembly record. This table is
not a download manifest. No files were downloaded for T6. Current assembly
accessions are pinned here for reproducibility; the event coordinates still
require extraction from the exact downloaded FASTA/GFF and independent
comparison with the cited figure.

| Use | Species | Assembly | RefSeq annotation | NCBI Gene record | Current genomic accession/interval checked? |
|---|---|---|---|---|---|
| FDPS, MAMSTR, PKM | *Homo sapiens* | GRCh38.p14, [`GCF_000001405.40`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) | `RS_2025_08` | [FDPS 2224](https://www.ncbi.nlm.nih.gov/gene/2224); [MAMSTR 284358](https://www.ncbi.nlm.nih.gov/gene/284358); [PKM 5315](https://www.ncbi.nlm.nih.gov/gene/5315) | Gene records available; FDPS/PKM event intervals not checked. MAMSTR gene bounds recorded in `real_positive_cases.md`; exon-3 cuts not checked. |
| FDPS, MAMSTR, PKM | *Mus musculus* | GRCm39, [`GCF_000001635.27`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/) | `RS_2024_02` | [Fdps 110196](https://www.ncbi.nlm.nih.gov/gene/110196); [Mamstr 74490](https://www.ncbi.nlm.nih.gov/gene/74490); [Pkm 18746](https://www.ncbi.nlm.nih.gov/gene/18746) | Gene records available; event intervals not checked. |
| FDPS, PKM | *Rattus norvegicus* | GRCr8, [`GCF_036323735.1`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036323735.1/) | `RS_2024_02` | [Fdps 83791](https://www.ncbi.nlm.nih.gov/gene/83791); [Pkm 25630](https://www.ncbi.nlm.nih.gov/gene/25630) | Gene records available; the rat FDPS joined exon and PKM exon 9/10 coordinates have not been checked in this assembly. |
| FDPS | *Canis lupus familiaris* | ROS_Cfam_1.0, [`GCF_014441545.1`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_014441545.1/) | Annotation Release 106 | [FDPS 480129](https://www.ncbi.nlm.nih.gov/gene/480129) | Gene record and candidate assembly available; focal intron positions not checked. |
| MAMSTR, PKM | *Bos taurus* | ARS-UCD2.0, [`GCF_002263795.3`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/) | `RS_2024_12` | [MAMSTR 505540](https://www.ncbi.nlm.nih.gov/gene/505540); [PKM 512571](https://www.ncbi.nlm.nih.gov/gene/512571) | Official whole-gene intervals available; MAMSTR exon-3 shift and PKM exon 9/10 coordinates not checked. |
| MAMSTR, PKM | *Equus caballus* | TB-T2T, [`GCF_041296265.1`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_041296265.1/) | `RS_2024_12` | [MAMSTR 100054356](https://www.ncbi.nlm.nih.gov/gene/100054356); [PKM 100063687](https://www.ncbi.nlm.nih.gov/gene/100063687) | Gene records available; focal exon coordinates not checked. Equine M1/M2 transcript publication is an external evidence source. |

These assembly records are candidates for new FASTA and GFF inputs. The old
accessions cited in the papers are sequence or gene records, not necessarily
the selected genome-assembly versions. Never copy their coordinates into a
current manifest without mapping the flanking sequence and recording the
accession.version on which the coordinates were checked.

### Required pre-analysis checks

1. Obtain genome FASTA and RefSeq GFF3 for each selected assembly from its
   official NCBI Assembly/Datasets record. Obtain the corresponding protein
   set for OrthoFinder from that same annotation release.
2. Check each candidate orthogroup contains one protein-coding locus per
   species; review gene-tree topology and inspect additional genome hits,
   pseudogenes and internal exon repeats separately. NCBI/Ensembl ortholog
   labels alone do not satisfy the upstream single-copy input requirement;
   qualify the set with OrthoFinder or a documented equivalent analysis.
3. Record GeneID, transcript/protein accession.version, assembly accession,
   annotation release, chromosome accession.version, strand and exact event
   coordinates in the case manifest after the feature-to-sequence mapping is
   reviewed.
4. Keep published event expectations in the separate truth/scoring record.
   The analysis manifest must not contain expected branches, event labels or
   manually selected states.

### Resources still to stage

No T6 input files are staged. Deduplicating across the three cases yields six
assembly/annotation pairs to obtain if all candidates proceed: human
`GCF_000001405.40`, mouse `GCF_000001635.27`, rat `GCF_036323735.1`, dog
`GCF_014441545.1`, cow `GCF_002263795.3` and horse `GCF_041296265.1`. Required
files per selected assembly are genomic FASTA, RefSeq GFF3 and protein FASTA
for OrthoFinder. The dog and horse assemblies can be dropped if the final
trees do not include those taxa. Existing server-local availability of these
six mammalian datasets has not been inspected as part of T6.
