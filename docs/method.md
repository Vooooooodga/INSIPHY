# Method Overview

## Scope and input

INSIPHY compares gene-internal structure in a supplied set of single-copy orthologous genes. The user supplies one gene copy per species, genome FASTA, existing GFF3/GTF annotation, and a rooted species tree. Gene-level homology is an upstream result; an OrthoFinder result can be imported, or an equivalent curated manifest can be provided. INSIPHY does not search the genome for orthologous genes or infer a gene tree in this first-version scope.

The analysis asks which internal sequence units correspond across the supplied genes, what their observed annotated roles and splice boundaries are, and which structural-state changes are supported on the fixed tree. It reports evidence and model results. It does not assign molecular mechanisms or claim that a statistical result proves a particular event occurred.

## Biological observations

The method carries original features, inferred correspondences, and phylogenetic observations as separate records. An element identifier tracks a proposed homologous unit; its biological label remains exon, intron, CDS, UTR, noncoding exon, or another feature type supplied by the annotation. An `EG_*` identifier is stable only for a frozen input, parameter set, and correspondence algorithm version.

Three observation layers are analyzed separately:

| Layer | State 0 | State 1 | What the comparison describes |
|---|---|---|---|
| `exon_presence` | homologous DNA unit absent | homologous DNA unit present | presence of the corresponding sequence interval |
| `exon_role` | not exonic in the selected supplied-annotation view | exonic in the selected supplied-annotation view | annotation-conditional exon identity of DNA that corresponds across species |
| `splice_junction` | the homologous cut is absent from the selected transcript view | the homologous cut is present | an intron at a corresponding position between or within mapped units |

DNA presence and exon role answer different questions. A DNA alignment can support sequence presence while its role remains unknown. A coding projection or another predicted feature is retained as a candidate and does not by itself establish an annotated exon role. A role contrast is coded only where the homologous DNA interval is sufficiently covered and the supplied transcript paths inform that interval.

The input annotation may contain CDS, 5′/3′ UTR, noncoding exons, introns, and other sequence-feature classes. Their original type, coordinates, parentage, and transcript membership are retained. A feature class enters formal phylogenetic inference only when a homologous unit and a clear state definition are available. Overlapping annotations may describe different transcript paths or feature classes at the same DNA interval.

## Evidence construction

### 1. Preserve transcript paths and annotation context

The default `transcript_policy=all` retains every supplied transcript path. Identical genomic structures may be deduplicated for computation while their transcript identifiers remain traceable. `annotation_view=repertoire` summarizes whether a structure is used by any supplied path. A `canonical` view restricts role and junction observations to the explicitly selected canonical path and records its selection rule. Sequence-presence observations are view independent. Repertoire and canonical are separate view-specific matrices; canonical is a sensitivity analysis and is never substituted into a frozen repertoire analysis.

Role states are conditional on the annotation supplied by the user. For example, role state 0 means the corresponding interval is non-exonic in the provided, adequately covering paths; it does not establish that the interval is never transcribed in any tissue or condition. If no supplied path covers the interval, its role is unknown. Annotation and search boundaries are kept distinct so a feature extending past a gene row, or a candidate found in a declared flanking search interval, remains traceable.

### 2. Build coding sequence coordinates

For coding evidence, each transcript CDS is assembled in transcript orientation, including negative-strand paths. A family-level protein multiple sequence alignment provides shared residue columns. Residue columns are projected through codons and CDS bases back to genomic intervals, retaining the source feature and transcript for each block. Identical proteins can share alignment computation while their transcript paths remain distinct. Protein correspondence is local: a common MSA column alone does not establish an unambiguous genomic mapping; local sequence and positional anchors, competing placements, and actual projected blocks are retained.

L-INS-i and E-INS-i are available for the family protein alignment. The output records the selected mode and backend. The exact executable version belongs to the container or execution provenance and is not inferred from the mode name. This uses a family-wide coding coordinate system; the supplied species tree remains the tree for evolutionary inference. For mixed CDS/UTR exons, protein evidence covers only the projected CDS blocks. UTR sequence needs nucleotide correspondence evidence or remains uncovered. The use of progressive alignment ideas is limited to sequence correspondence; INSIPHY does not run a whole-genome alignment.

### 3. Align short DNA intervals with local anchors

Short intronic, noncoding, and boundary intervals use two stages. A feature-bounded search first records candidate alignments with their scoring scheme, orientation, score, candidate count, completeness status, and aligned blocks. It cannot create a hard membership. When both flanks map uniquely, in the expected order and on the same transcript path, the intervening genomic interval is extracted and realigned. Only this anchor-bounded result can support a hard nucleotide membership or position observation. The current bounded nucleotide candidate scheme `nt_blastn_v1` uses match `+2`, mismatch `-3`, gap-open `-7`, and gap-extension `-2`. Optimal ties are retained where enumerated. Ambiguous or incomplete flank configurations remain candidate evidence.

### 4. Resolve ordered correspondences and block coverage

Candidate alignments are connected into monotone chains when their blocks are compatible in query and target transcript order. Chaining compares scores only within the same molecule type and scoring scheme. Structural role and splice conservation do not add score, so the correspondence step does not prefer a candidate merely because it preserves an exon pattern. Near-optimal mappings and incomplete candidate enumeration remain explicit.

A resolved hard membership is identified by the element, the source occurrence, and the actual matched subinterval blocks. One parent occurrence can therefore contribute several local memberships. A partial alignment supports only its covered subinterval and never expands to the full parent feature. Candidate and unresolved mappings remain review records and create no hard phylogenetic observation. For one-to-many relations, projections that cover distinct, ordered, complementary portions of the same reference interval are classified as complementary. Projections that repeatedly cover the same reference portion are classified as repeated overlap and remain separate repeat instances or ambiguous mappings. These sequence-coverage relationships define the structural correspondence; parsimony subsequently evaluates mapped junction states on the tree. Multiple junctions within one 1:n relation share a linked-group identifier so their dependence is visible.

An absence call requires a resolved homologous location, ordered flanking anchors, and evidence that the intervening sequence was searchable. Missing annotation, an unsearched boundary, an ambiguous repeated hit, or an incomplete search yields unknown rather than absence.

## Phylogenetic analysis

The default analysis uses equal-cost maximum parsimony on the supplied rooted species tree. For each site, it retains every node state and branch endpoint pair found in all globally minimum-change histories. A directed change is `required` on a branch when every optimal history assigns it there, and `possible` when at least one optimal history does. Ties are reported as sets; the method does not choose one arbitrary history or convert counts of equally parsimonious histories into probabilities. The method follows the fixed-tree minimum-change framework of [Sankoff (1975)](https://epubs.siam.org/doi/10.1137/0128004).

Tree changes are reported at the resolution of the observed layer: sequence-unit gain/loss, annotation-conditional exon-role gain/loss, or intron/junction gain/loss. For a `1↔2` relation, the junction character reports `event_type=intron_gain` or `intron_loss`, while `structural_relation` reports `split` or `fusion`. These two columns describe the same observed structural contrast at different resolutions.

For a `1↔n` relation with `n>=3`, all component cutpoints must co-occur in one species, one gene copy, and one transcript path before the summary can report a `required` compound structural event. Component changes that lack this joint path evidence are retained as individual sites and summarized as `possible_non_joint`. `compound_structural_events.tsv` records this distinction. These outputs describe structural state changes; they do not assign a molecular mechanism or an absolute time.

The repertoire matrix is generated once and frozen. Parsimony, ER/ARD, and foreground CTMC consume that same schema-v3 file; optional likelihood runs receive it through `--structural-site-matrix`. A canonical sensitivity analysis first creates its own frozen matrix and is kept in a separate result directory.

Optional ER/ARD and foreground CTMC analyses are described in [Statistical Model](statistical_model.md). They estimate rates and compare explicit rate models conditional on the supplied tree, structural matrix, and ascertainment rule. Their P values concern those model comparisons; they do not establish that a particular branch event occurred.

## Interpretation limits

The formal first-version scope is strongest for single-copy orthologs with at least some locally alignable gene-internal sequence and usable annotation paths. Extreme clusters of near-identical microexons can have multiple equally plausible mappings; the affected interval is reported as unresolved while independent regions remain analyzable. Whole-gene multi-copy families are outside this analysis scope because gene-copy correspondence requires additional inference.

Genome sequence and annotation alone do not establish tissue-specific transcript use, sex-specific expression, regulatory mechanism, or functional consequences. INSIPHY reports sequence, supplied-annotation, and tree-based structural evidence for users to interpret with other data. Complex rearrangements and feature classes without a defensible homology/state definition may be retained descriptively without entering the formal tree model.

Relevant biological background includes comparative studies of [splice-boundary evolution](https://academic.oup.com/gbe/article/8/8/2340/2198117), a review of [intron biology](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2023.1150212/full), and work on [cross-species transcript paths](https://pmc.ncbi.nlm.nih.gov/articles/PMC8327911/). These citations provide biological context; the state definitions and inference rules used here are specified above.
