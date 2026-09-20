# SPOG_00055 Genomic Case

The manifest supplies the rec8 gene group in four Schizosaccharomyces species.
All four numeric GeneIDs and their gene/mRNA/CDS records were read directly from
the assembly-matched local GFF3 files. Download completion was reported by the
parent workflow; no analysis was run during this metadata update.

| Species | Manifest GeneID | GFF gene ID | Annotated transcript |
|---|---|---|---|
| S. pombe | 2540458 | gene-SPOM_SPBC29A10.14 | NM_001021970.3 |
| S. japonicus | 7050558 | gene-SJAG_03926 | XM_002175012.2 |
| S. octosporus | 25033695 | gene-SOCG_04733 | XM_013165325.1 |
| S. cryophilus | 25034387 | gene-SPOG_00055 | XM_013167562.1 |

Use `manifest.tsv` and `species_tree.tsv` as the two inputs to `build-case`.
The `gene_id` values resolve through the GFF `Dbxref=GeneID:` attributes.
The source files remain immutable at their flat database paths. Transcript
accessions identify provided annotations; no transcriptome data are required.
These are curated gene mappings supplied to the method; an OrthoFinder analysis
has not been performed for this case in this preparation.

In the selected ASM294v2 annotation, S. pombe rec8 has `partial=true` and
`start_range=.,2568938` on the gene/mRNA records. Retain this annotation
uncertainty when interpreting its 5-prime end. S. cryophilus uses the historical
`SPOG` locus prefix in this assembly.

The rooted topology is `(japonicus,(pombe,(octosporus,cryophilus)))`, following
[Zhu and Niu 2013, Figure 1](https://doi.org/10.1371/journal.pone.0061683.g001).
Every non-root branch has unit length for qualitative parsimony. These values
carry no divergence-time or substitution-rate estimates.

Published expectations are exclusively in [truth_events.tsv](truth_events.tsv)
for external scoring. Keep that file outside the generated analysis input
directory and pass no event label, expected branch or ancestral state to
`build-case` or `run`. Assembly coordinates and preparation status are in
[the dataset record](../../../docs/real_positive_cases.md).
