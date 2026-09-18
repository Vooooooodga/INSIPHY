# Real Positive Case Preparation

This document records genome resources and gene identifiers for INSIPHY
real-data evaluation. Published event expectations are external scoring metadata,
documented here and in the separate `truth_events.tsv` files linked below.
Analysis receives the resource manifest and tree, without those scoring files.

## Ready Metadata: 2026-09-18

The parent workflow reports completed direct downloads of the four yeast
assemblies and D. ananassae. The three other HDAC3 species use existing project
genome/GFF resources. Gene, transcript and CDS records below were read directly
from these exact local GFF files with `zgrep`; annotation providers and assembly
versions come from their headers. Original database files were only read.

| Case | Manifest | Tree | Separate external scoring file | Metadata state |
|---|---|---|---|---|
| spbc29a10_02 | [manifest](../examples/real_cases/spbc29a10_02/manifest.tsv) | [tree](../examples/real_cases/spbc29a10_02/species_tree.tsv) | [scoring](../examples/real_cases/spbc29a10_02/truth_events.tsv) | Four exact GFF GeneIDs; all resource paths supplied |
| spog_00055 | [manifest](../examples/real_cases/spog_00055/manifest.tsv) | [tree](../examples/real_cases/spog_00055/species_tree.tsv) | [scoring](../examples/real_cases/spog_00055/truth_events.tsv) | Four exact GFF GeneIDs; all resource paths supplied |
| hdac3 | [manifest](../examples/real_cases/hdac3/manifest.tsv) | [tree](../examples/real_cases/hdac3/species_tree.tsv) | [scoring](../examples/real_cases/hdac3/truth_events.tsv) | Four exact GFF GeneIDs; all required genome/GFF paths supplied |

These files can now be passed to `build-case`. This metadata task did not run
case construction, inference, tests or additional downloads. Parent analyses
can proceed independently of the deferred Bap170 and BCKDK preparation.

The manifests supply curated gene choices and neutral resource metadata.
No upstream OrthoFinder run or independent single-copy membership analysis was
performed in this task. The raw GFF reads establish the identities and feature
hierarchies in the selected releases; the supplied homology grouping remains
an upstream assumption. Keep each scoring file in the example directory and
outside the generated analysis input directory. No expected branch, event,
ancestral state, `role_hint`, `source_label` or `copy_role` is supplied in these
three manifests.

## Early Download Rows

The following NCBI RefSeq URLs were confirmed by HTTP HEAD on 2026-09-16. They
are ready to be copied into the project download table by the parent workflow
owner. No file was downloaded during this preparation step.

| case_id | species | assembly | assembly_name | genome_gz_size | gff_gz_size | genome_url | gff_url | protein_url |
|---|---|---|---|---:|---:|---|---|---|
| spbc29a10_02 | Schizosaccharomyces_pombe | GCF_000002945.1 | ASM294v2 | 3989633 | 1387873 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_protein.faa.gz |
| spbc29a10_02 | Schizosaccharomyces_japonicus | GCF_000149845.2 | SJ5 | 3492571 | 771207 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/845/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/845/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/845/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_protein.faa.gz |
| spbc29a10_02 | Schizosaccharomyces_octosporus | GCF_000150505.1 | SO6 | 3574252 | 766584 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/505/GCF_000150505.1_SO6/GCF_000150505.1_SO6_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/505/GCF_000150505.1_SO6/GCF_000150505.1_SO6_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/505/GCF_000150505.1_SO6/GCF_000150505.1_SO6_protein.faa.gz |
| spbc29a10_02 | Schizosaccharomyces_cryophilus | GCF_000004155.1 | SCY4 | 3651991 | 798150 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/155/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/155/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/155/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_protein.faa.gz |
| spog_00055 | Schizosaccharomyces_pombe | GCF_000002945.1 | ASM294v2 | 3989633 | 1387873 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_protein.faa.gz |
| spog_00055 | Schizosaccharomyces_japonicus | GCF_000149845.2 | SJ5 | 3492571 | 771207 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/845/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/845/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/845/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_protein.faa.gz |
| spog_00055 | Schizosaccharomyces_octosporus | GCF_000150505.1 | SO6 | 3574252 | 766584 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/505/GCF_000150505.1_SO6/GCF_000150505.1_SO6_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/505/GCF_000150505.1_SO6/GCF_000150505.1_SO6_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/505/GCF_000150505.1_SO6/GCF_000150505.1_SO6_protein.faa.gz |
| spog_00055 | Schizosaccharomyces_cryophilus | GCF_000004155.1 | SCY4 | 3651991 | 798150 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/155/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/155/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/155/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_protein.faa.gz |
| hdac3 | Drosophila_ananassae | GCF_017639315.1 | ASM1763931v2 | 64115269 | 5546242 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/315/GCF_017639315.1_ASM1763931v2/GCF_017639315.1_ASM1763931v2_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/315/GCF_017639315.1_ASM1763931v2/GCF_017639315.1_ASM1763931v2_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/315/GCF_017639315.1_ASM1763931v2/GCF_017639315.1_ASM1763931v2_protein.faa.gz |

## Official Yeast Alternative URLs

The four yeast assemblies also resolve through the official NCBI RefSeq fungi
tree. These URLs point to the same RefSeq assembly versions and file sizes as
the `genomes/all/GCF/...` URLs above. They are suitable alternate source rows if
the all-assembly path is temporarily slow.

| species | assembly | genome_url | gff_url | protein_url |
|---|---|---|---|---|
| Schizosaccharomyces_pombe | GCF_000002945.1 | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_pombe/all_assembly_versions/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_pombe/all_assembly_versions/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_pombe/all_assembly_versions/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_protein.faa.gz |
| Schizosaccharomyces_japonicus | GCF_000149845.2 | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_japonicus/all_assembly_versions/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_japonicus/all_assembly_versions/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_japonicus/all_assembly_versions/GCF_000149845.2_SJ5/GCF_000149845.2_SJ5_protein.faa.gz |
| Schizosaccharomyces_octosporus | GCF_000150505.1 | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_octosporus/all_assembly_versions/GCF_000150505.1_SO6/GCF_000150505.1_SO6_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_octosporus/all_assembly_versions/GCF_000150505.1_SO6/GCF_000150505.1_SO6_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_octosporus/all_assembly_versions/GCF_000150505.1_SO6/GCF_000150505.1_SO6_protein.faa.gz |
| Schizosaccharomyces_cryophilus | GCF_000004155.1 | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_cryophilus/all_assembly_versions/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_cryophilus/all_assembly_versions/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_cryophilus/all_assembly_versions/GCF_000004155.1_SCY4/GCF_000004155.1_SCY4_protein.faa.gz |

## Drosophila and Mammal Candidate URLs

These rows were also confirmed by HTTP HEAD. They are useful after the yeast and
D. ananassae staging steps.

| case_id | species | assembly | assembly_name | genome_gz_size | gff_gz_size | protein_gz_size | genome_url | gff_url | protein_url |
|---|---|---|---|---:|---:|---:|---|---|---|
| bap170 | Drosophila_melanogaster | GCF_000001215.4 | Release_6_plus_ISO1_MT | 44235956 | 8688226 | 9785243 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz |
| bap170 | Drosophila_ananassae | GCF_017639315.1 | ASM1763931v2 | 64115269 | 5546242 | 7948587 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/315/GCF_017639315.1_ASM1763931v2/GCF_017639315.1_ASM1763931v2_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/315/GCF_017639315.1_ASM1763931v2/GCF_017639315.1_ASM1763931v2_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/315/GCF_017639315.1_ASM1763931v2/GCF_017639315.1_ASM1763931v2_protein.faa.gz |
| bap170 | Drosophila_persimilis | GCF_003286085.1 | DperRS2 | 57355118 | 4037154 | 5806890 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/286/085/GCF_003286085.1_DperRS2/GCF_003286085.1_DperRS2_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/286/085/GCF_003286085.1_DperRS2/GCF_003286085.1_DperRS2_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/286/085/GCF_003286085.1_DperRS2/GCF_003286085.1_DperRS2_protein.faa.gz |
| bap170 | Drosophila_pseudoobscura | GCF_009870125.1 | UCI_Dpse_MV25 | 49527736 | 5766771 | 7688435 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/870/125/GCF_009870125.1_UCI_Dpse_MV25/GCF_009870125.1_UCI_Dpse_MV25_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/870/125/GCF_009870125.1_UCI_Dpse_MV25/GCF_009870125.1_UCI_Dpse_MV25_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/870/125/GCF_009870125.1_UCI_Dpse_MV25/GCF_009870125.1_UCI_Dpse_MV25_protein.faa.gz |
| bckdk | Homo_sapiens | GCF_000001405.40 | GRCh38.p14 | 972898531 | 78190483 | 28489875 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz |
| bckdk | Pan_troglodytes | GCF_028858775.2 | NHGRI_mPanTro3-v2.1_pri | 940198818 | 36808605 | 22084551 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/858/775/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri/GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_protein.faa.gz |
| bckdk | Mus_musculus | GCF_000001635.27 | GRCm39 | 834262758 | 43813395 | 23647097 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_protein.faa.gz |
| bckdk | Rattus_norvegicus | GCF_036323735.1 | GRCr8 | 875872536 | 37442693 | 20138563 | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.fna.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.gff.gz | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_protein.faa.gz |

## Exact Annotation Records

Coordinates below are one-based, inclusive, on the named accession.version.
They describe the gene feature in the downloaded assembly, including annotated
UTRs where present. Transcript and protein accessions record the supplied GFF
hierarchy; no RNA experiments or transcriptome sequence files are inputs.

### spbc29a10_02

| Species | Assembly | GeneID | GFF gene ID | Sequence accession:start-end:strand | Transcript | Protein |
|---|---|---|---|---|---|---|
| S. pombe | GCF_000002945.1 | 2540637 | gene-SPOM_SPBC29A10.02 | NC_003423.3:2536150-2538254:+ | NM_001021957.3 | NP_596047.1 |
| S. japonicus | GCF_000149845.2 | 7052533 | gene-SJAG_04271 | NW_011627861.1:2734877-2737107:- | XM_002175352.2 | XP_002175388.1 |
| S. octosporus | GCF_000150505.1 | 25033683 | gene-SOCG_04721 | NW_013185621.1:201288-203392:- | XM_013165313.1 | XP_013020767.1 |
| S. cryophilus | GCF_000004155.1 | 25034398 | gene-SPOG_00066 | NW_013185626.1:150371-152505:- | XM_013167573.1 | XP_013023027.1 |

The S. pombe symbol is `spo5`, with `old_locus_tag=SPBC29A10.02`; all four
products are annotated as meiotic RNA-binding protein 1. Numeric identifiers
match `Dbxref=GeneID:` in the source GFF. These coordinates belong to ASM294v2
and the other assembly versions listed here; they must not be replaced by
coordinates from a newer NCBI Gene assembly without changing the resources.

### spog_00055

| Species | Assembly | GeneID | GFF gene ID | Sequence accession:start-end:strand | Transcript | Protein |
|---|---|---|---|---|---|---|
| S. pombe | GCF_000002945.1 | 2540458 | gene-SPOM_SPBC29A10.14 | NC_003423.3:2568938-2570993:+ | NM_001021970.3 | NP_596059.1 |
| S. japonicus | GCF_000149845.2 | 7050558 | gene-SJAG_03926 | NW_011627861.1:3468273-3470644:- | XM_002175012.2 | XP_002175048.1 |
| S. octosporus | GCF_000150505.1 | 25033695 | gene-SOCG_04733 | NW_013185621.1:170866-172975:- | XM_013165325.1 | XP_013020779.1 |
| S. cryophilus | GCF_000004155.1 | 25034387 | gene-SPOG_00055 | NW_013185626.1:118532-120891:- | XM_013167562.1 | XP_013023016.1 |

The S. pombe symbol is `rec8`. Its gene/mRNA records carry `partial=true`
and `start_range=.,2568938`; this is uncertainty about the annotated 5-prime
end. Retain those flags and evaluate internal correspondence separately from
terminal completeness. `SPOG_00055` belongs to S. cryophilus in this GFF,
consistent with [NCBI Gene 25034387](https://www.ncbi.nlm.nih.gov/gene/25034387).

The S. pombe header attributes its annotation to PomBase. The other three
yeast headers attribute annotation to the INSDC submitter. They supply no
numbered annotation release; the manifests record this explicitly alongside
the assembly version, instead of claiming an unspecified current release.

### hdac3

| Species | Assembly | GeneID | GFF gene ID | Sequence accession:start-end:strand | Transcript | Protein |
|---|---|---|---|---|---|---|
| D. ananassae | GCF_017639315.1 | 6500490 | gene-LOC6500490 | NC_057927.1:3929065-3930814:+ | XM_001953300.4 | XP_001953335.1 |
| D. melanogaster | GCF_000001215.4 | 44446 | gene-Dmel_CG2128 | NT_033777.3:5472810-5474600:- | NM_143721.4 | NP_651978.2 |
| D. simulans | GCF_016746395.2 | 6726814 | gene-LOC6726814 | NC_052523.2:1975669-1977178:- | XM_002102185.4 | XP_002102221.1 |
| D. yakuba | GCF_016746365.2 | 6535413 | gene-LOC6535413 | NC_052530.2:3762902-3764434:- | XM_002096026.4 | XP_002096062.1 |

The GFF headers specify NCBI Annotation Release 102 for D. ananassae,
FlyBase Release 6.54 for D. melanogaster, NCBI Release 103 for D. simulans and
NCBI Release 102 for D. yakuba. D. melanogaster carries
`exception=dicistronic gene`; its HDAC3 mRNA/CDS hierarchy is listed above.
The other three GFF gene names use `LOC` identifiers, while their mRNA/CDS
product is histone deacetylase 3. The manifest uses the stable numeric GeneID.

All paths use `/data/db/genome/<species>/<assembly>/genome.fa.gz` and
`annotation.gff3.gz`, matching their actual download placement. The optional
protein-file fields are empty for the three existing Drosophila comparators;
the method extracts sequence from the genome and supplied annotation.

### Deferred Cases

Bap170 and BCKDK metadata are not finalized in this update. Their candidate
download tables remain available. Published expectations belong to their
separate scoring files. BCKDK locus-interval retrieval research is deferred;
no mammalian download or input preparation is claimed here.

## Figure-to-Sequence Correspondence: 2026-09-18

The three focal published intervals have now been manually located in the
modern sequences. This comparison used the official figure content and existing
`case/segment_sequences.fasta`, `case/segment_occurrences.tsv` and
`case/intron_sites.tsv` under
`/data/projects/intragenic_structure/results/20260916_214900_insiphy/<case>/threads_1/analysis/`.
Assembly and annotation versions are those in Exact Annotation Records above.
All coordinates below are one-based, inclusive; sequence strings follow the
transcript strand. These are external scoring observations, not inference inputs.
No alignment threshold or analysis input was changed to obtain correspondence.

### Hdac3: Fig. S4

The 62 bp third intron of D. ananassae `XM_001953300.4` is
`NC_057927.1:3929894-3929955:+` in `GCF_017639315.1`. Its complete sequence
matches [Farlow et al. 2010, Fig. S4B](https://doi.org/10.1371/journal.pgen.1000819.s004):

```text
GTGGGTATATATGAACCGCCAAGTTGCTTAGGCTACTAATTTGTACTTTCTTCTATTTCCAG
```

The flanking coding sequence also matches the figure. The upstream CDS contains
717 nt, placing the cut after codon 239. The modern D. melanogaster continuous
coding sequence across this cut matches Fig. S4C. The corresponding uninterrupted
cuts in the existing modern comparators are:

| Species | Accession | Consecutive bases flanking the cut, in transcript order | Strand | Evidence |
|---|---|---|---|---|
| D. melanogaster | NT_033777.3 | 5473476 / 5473475 | - | Fig. S4C and modern sequence |
| D. simulans | NC_052523.2 | 1976337 / 1976336 | - | Modern sequence correspondence |
| D. yakuba | NC_052530.2 | 3763580 / 3763579 | - | Modern sequence correspondence |

Fig. S4C includes D. pseudoobscura, not D. simulans or D. yakuba. Its
D. pseudoobscura sequence was not checked against a modern local resource in
this audit. The four-species analysis must keep the root-state ambiguity
described below separate from successful identification of this cut.

### rec8: Fig. 3B

The complete 39 bp underlined S. cryophilus sequence in
[Zhu and Niu 2013, Fig. 3B](https://doi.org/10.1371/journal.pone.0061683.g003)
and its flanks match `SPOG_00055`, transcript `XM_013167562.1`, in
`GCF_000004155.1`. The region lies at
`NW_013185626.1:118741-118779:-`, positions 1278-1316 of the existing
`Schizosaccharomyces_cryophilus_Scry_SPOG_00055_001_exon` sequence:

```text
GTACGATTTATCATACGTATTTCTTACTTACCTACTTAG
```

The corresponding introns in the other three species also match the original
figure sequences, including the complete S. pombe intron. High-resolution
inspection resolved an initial ambiguity in reading the figure thumbnail;
no discrepancy at this site was established.

| Species | Modern interval | Length | Annotated role |
|---|---|---:|---|
| S. cryophilus | NW_013185626.1:118741-118779:- | 39 bp | Coding exon |
| S. pombe | NC_003423.3:2570676-2570715:+ | 40 bp | Fourth intron in transcript order |
| S. octosporus | NW_013185621.1:171120-171162:- | 43 bp | Fourth intron in transcript order |
| S. japonicus | NW_011627861.1:3468615-3468653:- | 39 bp | Fourth intron in transcript order |

Score both homologous sequence retention and the exon/intron role difference.
The S. pombe 5-prime `partial=true` annotation remains unchanged. On negative
strands, use transcript order or genomic coordinates for intron identity;
the initial case table numbered introns by ascending genomic coordinate.

### spo5: Fig. 3A

The complete S. pombe intron and its flanks match
[Zhu and Niu 2013, Fig. 3A](https://doi.org/10.1371/journal.pone.0061683.g003).
It is the 55 bp second intron of `NM_001021957.3`,
`NC_003423.3:2536843-2536897:+`, in `GCF_000002945.1`:

```text
GTATGATACACAGCTTATGAAAAAACTAATACTTACCAAACTTTTCGTCTTTTAG
```

The figure's corresponding coding regions in the three other species were
located using the displayed sequences and flanking anchors:

| Species | Modern interval | Length | Annotated role |
|---|---|---:|---|
| S. pombe | NC_003423.3:2536843-2536897:+ | 55 bp | Second intron |
| S. cryophilus | NW_013185626.1:151778-151824:- | 47 bp | Coding exon |
| S. octosporus | NW_013185621.1:202664-202711:- | 48 bp | Coding exon |
| S. japonicus | NW_011627861.1:2736356-2736406:- | 51 bp | Coding exon |

The separate 48 bp first intron in S. pombe is outside this figure-specific
positive label; no additional historical claim is assigned to it here.

### Recovery Status

In the first iteration, `20260916_214900_insiphy`, all three published events
remained automatically unrecovered despite completed jobs. The manual
figure-to-sequence correspondence above establishes the independent site-level
expectations used to assess subsequent runs. Those coordinates and the separate
`truth_events.tsv` scoring files remain unchanged.

Read-only focal-site review of the repaired run, `20260918_111929_insiphy`,
supports the following observations. Counts describe splice-junction states
at the specific homologous position, not gene-wide exon counts.

| Case | Focal correspondence | Junction present / absent / unknown | Accepted recovery scope |
|---|---|---|---|
| hdac3 | `D470_A471` on the D. yakuba reference exon; exact D. ananassae 62 bp intron `NC_057927.1:3929894-3929955:+` | 1 / 3 / 0 | Junction contrast recovered: ananassae present; melanogaster, simulans and yakuba absent |
| spog_00055 / rec8 | `D1277_A1317` brackets the exact 39 bp cryophilus coding interval `NW_013185626.1:118741-118779:-` | 1 / 1 / 2 | Partial junction contrast recovered: octosporus present; cryophilus absent; pombe and japonicus unknown |
| spbc29a10_02 / spo5 | Published 55 bp intron `NC_003423.3:2536843-2536897:+` remains in the input | No focal junction site constructed | Unrecovered; no focal exon/intron role contrast or junction contrast |

Hdac3 retains two equally parsimonious one-change histories: a `split` on
`root -> Drosophila_ananassae`, or a `fusion` on
`root -> melanogaster_subgroup`. Both have `possible` placement; the recovered
four-species contrast does not uniquely determine a Dana gain.

For rec8, the equally parsimonious alternatives are a junction gain on the
octosporus terminal branch and a junction loss on the cryophilus terminal
branch. Direct correspondence between the 39 bp coding interval and homologous
intronic DNA remains unresolved in the automatic output, and the corresponding
exon/intron role contrast has not been recovered. The accepted result is
therefore partial recovery of the published structural difference.

All three positive cases completed with both 1 and 16 threads. Six core tables
agree within each pair in the completed baseline and the final inference-only
run `20260918_121100_insiphy`. That final run repairs a dsx role conflict using
unchanged sequence/correspondence evidence; these three focal outcomes persist.
The comparison does not establish equality of every output or validate the
biological accuracy of unresolved observations.
See [Real-Data Demonstration](real_data_benchmark.md) for the current automatic
results and execution status. Published truth remains separate from analysis
inputs, and no significant model comparison or unique historical direction is
claimed from these observations.

## Tree Rationale

The two yeast cases use `(japonicus,(pombe,(octosporus,cryophilus)))`, the
four-species restriction of [Zhu and Niu 2013, Figure 1](https://doi.org/10.1371/journal.pone.0061683.g001).
Their paper also used additional fungal outgroups for historical direction;
these four-species examples do not reproduce that complete sampling.

HDAC3 uses `(ananassae,(yakuba,(melanogaster,simulans)))`, the restriction of
the [Drosophila 12 Genomes Consortium 2007 tree](https://doi.org/10.1038/nature06341)
to the four available resources. This is a supplied species-tree topology
under the assumed orthology grouping, not a tree inferred from HDAC3 in this
task. Opposing character states in the two root daughter clades can permit
equally parsimonious ancestral states and transition directions. Root states
and expected directions must remain unconstrained by the scoring labels.

All three trees have zero root length and unit non-root branch lengths for
qualitative parsimony. The units are not calibrated ages or substitutions;
any likelihood analysis on these lengths would be conditional on the arbitrary
unit-length assumption. No rate estimation or significance result is claimed.

## Current Limits

The three active cases have complete resource and gene metadata and confirmed
manual correspondence for the focal published intervals above. All three were
unrecovered in `20260916_214900_insiphy`. In `20260918_111929_insiphy`, Hdac3 has
a complete focal junction contrast but two equally parsimonious directions;
rec8 has a two-species focal junction contrast, two unknown species, and no
automatic homologous-intronic-DNA or exon/intron role contrast; spo5 remains
unrecovered. These outcomes must be scored separately as site correspondence,
observed state coverage, role recovery and directional branch placement.
No overall accuracy estimate, significant P value or uniquely resolved event
direction is established by this audit. The six-table 1-thread versus 16-thread
comparisons passed. Bap170 and BCKDK modern positive-site correspondence
remains deferred. Current automatic results are tracked in
[Real-Data Demonstration](real_data_benchmark.md).

The supplied gene grouping has not undergone an independent genome-wide
single-copy assessment in this task. Retain all equally supported histories
when limited sampling leaves event direction uncertain. Original annotation
flags remain part of the input evidence, and published expectations must not
be used to repair or relabel it. No new data, RNA inputs or analysis jobs were
introduced during the figure-to-sequence audit.
