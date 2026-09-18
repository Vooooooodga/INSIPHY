# Method Overview

## Scope

INSIPHY studies gene-internal structure within an upstream-defined homologous gene set. The formal v0.14 path starts from one orthologous gene per species and a supplied rooted species tree. Gene discovery, genome-wide orthology inference, expression analysis, and molecular-mechanism assignment belong outside the current formal model.

Multi-copy routines remain available under `--analysis-scope experimental-multicopy` for future development. Current manuscript-level claims should use the repaired single-copy path.

## Biological units

The method distinguishes four related concepts:

- **raw feature**: an original annotation feature overlapping the target search window, including exon, CDS, UTR, noncoding exon, intron, nested RNA, and regulatory/context features when supplied in the annotation;
- **occurrence**: one observed or sequence-supported interval in one species and gene copy;
- **element**: a homologous gene-internal unit such as an exon-like block, candidate non-exonic source, or explicit absence record;
- **structural site**: a binary observation derived from elements or splice junctions and used by parsimony or likelihood models.

`EG_*` labels are stable element identifiers. They are not a claim that every member is a confirmed coding exon. Element class and evidence fields carry that biological status.

## Three core observations

1. **Sequence presence**
   A homologous DNA unit is present, explicitly absent, or unknown. Absence requires sequence evidence from the searched interval and flanking context.

2. **Exonic role**
   A present homologous sequence has confirmed exonic role, confirmed non-exonic role, predicted exonic role, or unresolved role. Protein projection and `predicted_CDS` are recorded as prediction evidence.

3. **Splice junction**
   Two homologous units are separated by an intron, directly joined in one exon path, or unknown. Junction homology uses projected donor/acceptor coordinates, not a length-ratio proxy.

## Four operational stages

### Stage 1: recover observed structure

`build-case` extracts target loci from genome FASTA and GFF3/GTF. The v0.14 design treats all transcript paths as the default biological repertoire. Identical genomic structures across isoforms can be deduplicated; distinct splice boundaries and transcript-specific paths remain visible.

The extracted representation records:

- complete exon intervals;
- CDS and UTR subintervals;
- intron intervals and splice motifs;
- phase and frame attributes;
- transcript order and path membership;
- original annotation bounds and expanded search bounds.

`raw_gene_features.tsv` preserves annotation feature types, coordinates, parent relationships and attributes within the search window, distinguishing target-gene descendants from overlapping context. Targeted retention tests passed and the first real runs produced this table. Additional feature types need their own correspondence evidence and observation model before they can enter formal phylogenetic analysis.

### Stage 2: supplement incomplete annotation

Within the target gene and declared flanking/search window, INSIPHY records evidence for hidden or incomplete annotation:

- DNA evidence supports sequence presence;
- miniprot or protein-to-genome projection can support a predicted CDS-like interval;
- boundary conflicts mark disagreement between annotation and projected structure;
- deletion evidence requires ordered flanks and searched intervening sequence.

Predicted records do not overwrite the input GFF. They enter downstream tables with evidence status and role qualifiers.

### Stage 3: infer internal correspondence

The correspondence problem asks which gene-internal intervals across species represent the same structural unit. Pairwise evidence uses:

- sequence identity and aligned coverage;
- exact boundary projection;
- local left/right context;
- strand and order consistency;
- phase and reading-frame compatibility where applicable;
- tree-guided close-to-distant merging inside the supplied ortholog set.

Split/fusion relations are represented by ordered complementary projections. A `1:n` split can involve more than two descendants when each part maps to a distinct, ordered part of the same reference unit. Repeated full-overlap hits are retained as repeats or ambiguous evidence.

INSIPHY uses mature alignment tools and adapter code. It does not introduce a new general-purpose alignment algorithm.

### Stage 4: locate changes on the tree

The default analysis is equal-cost maximum parsimony on the supplied rooted tree. For every structural site, INSIPHY retains all globally optimal node states and branch endpoint pairs.

- `required`: every minimum-change history places that directed change on the branch.
- `possible`: at least one minimum-change history places that directed change on the branch.

These terms describe branch placement under the cost model. They do not provide a P value.

Optional ER/ARD and foreground likelihood models estimate shared rates within a family-layer and compare nested models when regular conditions hold. These models can help describe rate asymmetry or foreground rate shifts, but the current first-version emphasis is qualitative branch placement and transparent uncertainty.

## Classic branch units and gene-structure events

Classical systematics often assigns events to branches of a fixed tree. INSIPHY follows that framework, but the characters are gene-internal structural sites. A branch event therefore means:

```text
one structural site changes state between the parent and child endpoints of one tree branch
```

It does not directly name a molecular mechanism. For example, a junction gain/loss pattern can describe an exon split/fusion structure, while the cause may involve mutation, annotation uncertainty, lineage-specific transcript usage, or additional molecular processes.

## Visualization

The default figure uses horizontal gene tracks, complete annotated exon boxes and ribbons clipped to directly aligned bases between confirmed exon-like observations. Membership colors remain on the complete boxes. Protein-supported correspondence selects `protein_projected_blocks`, limiting the ribbon to coding bases without extending across UTR. Missing direct matches produce no ribbon. Transcript alternatives are drawn as separate lanes. Shared occurrences use the facing lanes of neighboring species for correspondence ribbons; distinct split members retain their connections to complementary intervals. Candidate source, predicted exon, and unknown states use distinct visual styles and do not create confirmed homology ribbons.

The integrated figure aligns the supplied tree with gene tracks. Parsimony symbols show required or possible branch placements. Optional probability symbols are drawn only when the selected fitted posterior is valid.

Both tree figures use supplied non-root branch lengths when all are known, preserving zero lengths. If any non-root length is missing, every branch is drawn with unit length and the figure labels this topological layout. A missing root length has no effect. This display rule leaves the input tree and optional CTMC branch-length requirements unchanged.

## Current empirical status

The final v0.14.0 assessment passed 136 selected formal regression tests under Slurm `61625` and completed ten re-inference/visualization tasks plus five comparisons from unchanged cases/evidence. Six core tables agree in every pair. Hdac3's focal contrast is recovered with ambiguous direction; rec8 is partial, spo5 unrecovered, and RpL32/dsx retain incomplete observations. Final dsx has two role contrasts, nonidentifiable/boundary-limited fits and no node/branch posteriors. SVG semantic and selected endpoint checks are complete; rendered inspection and broader biological performance assessment remain outstanding. The [benchmark record](real_data_benchmark.md) separates this assessment from the historical 97- and 132-test iterations and v0.13 demonstrations.
