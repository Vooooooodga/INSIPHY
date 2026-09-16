# Real-Data Demonstration

## RpL32 conserved control

The first formal single-copy demonstration uses RpL32 orthologs from:

- Drosophila melanogaster, GCF_000001215.4;
- Drosophila simulans, GCF_016746395.2;
- Drosophila erecta, GCF_003286155.1;
- Drosophila yakuba, GCF_016746365.2;
- Drosophila teissieri, GCF_016746235.2.

Genome FASTA and GFF3 files are stored outside the repository under
`/data/db/genome`. The repository contains accession-level manifests and the
species tree.

## Complete user path

```bash
insiphy build-case \
  --manifest examples/real_cases/rpl32_control/manifest.tsv \
  --species-tree examples/real_cases/rpl32_control/species_tree.tsv \
  --output-dir work/rpl32_control \
  --aligner mafft \
  --context-aligner minimap2 \
  --flank 1000 \
  --max-extension 10000 \
  --threads 4

insiphy run \
  --input-dir work/rpl32_control \
  --output-dir results/rpl32_control \
  --analysis-scope single-copy \
  --model parsimony \
  --evidence-aligner miniprot \
  --threads 4

insiphy visualize \
  --input-dir work/rpl32_control \
  --result-dir results/rpl32_control \
  --output-dir figures/rpl32_control
```

These commands show the package's user-facing stages. The server schedules them
with Slurm and an external Nextflow workflow. The package has no Nextflow
dependency. Version 0.13.0 uses the existing dependency container with an explicit
source overlay and four CPUs per case. Default parsimony uses the rooted tree;
optional ER/ARD analysis uses unit branch lengths for RpL32 and supplied lengths
for dsx, with estimated root frequency. No simulations or checksum runs were used.

Completed application run `20260916_202436_insiphy` is stored under
`/data/projects/intragenic_structure/results/20260916_202436_insiphy`.
It reuses the prepared cases archived in run `20260916_164410_insiphy`; genome
extraction and exon comparisons were completed in `20260916_162953_insiphy`.
The current run recomputes annotation evidence, both phylogenetic analyses, and
three SVG figures after repairing deletion-spanning genomic alignment.

## RpL32: conserved control

Three homologous exon groups and two splice boundaries have observations in
all five species. Sequence presence, exon role and both boundaries are conserved
in this matrix. Equal-cost parsimony requires zero changes at these sites.
No annotation-completion query remains for this group.

All three optional ER/ARD comparisons report `parameters_not_estimable` and
`p_value=NA`. With no observed state contrast, these data do not identify
gain/loss asymmetry. The embedded-null optimization start gives the same fitted
likelihood under ER and ARD. Earlier boundary patterns and P values from v0.12
are historical outputs based on different correspondence coding.

This control does not measure sensitivity to known exon gains or losses.

## Bee dsx: alternative exon organization

Seven bee orthologs from OrthoFinder group OG0009654 retain all annotated
transcripts: Apis cerana, Apis mellifera, Bombus ignitus, Bombus pascuorum,
Bombus terrestris, Frieseomelitta varia and Tetragonisca angustula. Genome and
annotation accessions are in the case manifest. No expression or sex-specific
transcript data enter this analysis.

In run `20260916_202436_insiphy`, the raw matrix contains 42 exon-associated
sequence groups, their 42 role observations, and 57 projected junction sites.
These group counts include alternative annotated intervals; they are not the
number of exons in a single transcript. Unknown observations are frequent:

| Layer | Sites | Observed cells | Unknown cells | Sites observed in all seven species |
|---|---:|---:|---:|---:|
| Sequence presence | 42 | 95 | 199 | 1 |
| Exonic role | 42 | 90 | 204 | 0 |
| Splice junction | 57 | 84 | 315 | 0 |

No site has a supported state-0/state-1 contrast in this run. Parsimony requires
zero changes under the observed constraints, with fully missing sites labeled
uninformative. This leaves dsx evolutionary history unresolved; the large number
of unknowns precludes a conclusion of whole-gene structural conservation.

Of 152 supplementary comparisons, 142 remain ambiguous, five suggest a boundary
conflict, four support sequence within already annotated exons, and one predicts
a missing coding interval in B. terrestris (NC_063273.1:2476398-2476441, minus
strand). This 44 bp interval is a protein-projection candidate, not a confirmed
transcribed exon or a dated exon-gain event. Its support score summarizes
alignment identity/coverage and has no posterior-probability interpretation.

Only sites with at least two observed species enter optional likelihood fitting:
24 presence sites, 24 role sites and 18 junction sites. Each layer has zero
observed state contrasts; all three model tests report parameters not estimable
and P=NA. The previously reported P=0.4524 used the older junction coding and
does not describe these repaired observations.

These gene-level likelihoods are conditional on the supplied tree, inferred
homologous elements and annotated isoforms. Adjacent exon and junction sites
are correlated; tests pooling them assume conditional independence. Missing
annotations and unknown states limit inference from the apparent absence of
an exon. Default figures show correspondence and parsimony placement categories.
Optional fitted probabilities condition on estimated parameters; their joint
parameter-sensitivity ranges remain unestimated.

## Deletion evidence and computational limits

The repaired deletion search uses minimap2 on the genomic interval spanning
homologous flanking exons. Both flanks must occur in one supported alignment,
with a query-only gap covering the expected exon and contiguous target bases
beside that gap. Source flanks follow one annotated transcript path. The search
can run without a protein reference, including for noncoding exons.

In the current 152 dsx comparisons, 72 lack a source flanking pair and 72 lack
an observed target homologous flank. Eight comparisons reach genomic alignment:
six have insufficient alignment identity, one lacks joint projection of both
flanks, and one has no gap spanning the expected exon. None supports a deletion.
The earlier internal alignment-size errors are absent from these results.
The structural matrix and annotation-completion counts remain unchanged.

Positive annotation completion in this run uses miniprot. For 122 comparisons,
the chosen reference has no usable protein context; this run does not perform
a supplementary DNA-presence search for those queries. Their unresolved status
therefore reflects both available evidence and this analysis choice. Noncoding
sequence recovery remains an important real-data evaluation requirement.

Each case received four CPUs and 12 GB RAM. Nextflow reports task execution
times of 9.9 s for RpL32 and 64 s for dsx, with peak resident memory of 107.1 MB
and 140.8 MB respectively. These measurements cover annotation completion,
inference and drawing from prepared cases; they exclude the earlier extraction
and pairwise alignments. Both Slurm jobs completed with exit status zero.

Figures in each case's `figures/` directory are `intragenic_synteny.svg`,
`phylogenetic_event_map.svg`, and `integrated_phylo_synteny.svg`. The default
color encoding includes homologous exon connections; the event panel displays
only the changes permitted by the reported reconstruction.

## Next real-data set

The next benchmark should contain single-copy genes with independently
documented structural changes. Candidate cases must satisfy:

- one ortholog per species;
- assembly and annotation versions available;
- enough species to distinguish alternative branch placements;
- sequence-level evidence for the altered exon or junction;
- a published history that can be reviewed independently of INSIPHY.

Jingwei and Sdic remain useful for future multi-copy development. Their
duplication histories place them outside the current formal single-copy
benchmark.
