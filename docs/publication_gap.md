# Publication Gap

## Current status

INSIPHY v0.12 is a runnable single-copy method prototype. It accepts an
upstream ortholog set, reconstructs exon-like correspondence from genome
sequence and annotation, builds explicit structural sites, and fits ER/ARD or
foreground CTMC models on a fixed species tree.

The code can support a methods manuscript description of the model and the
RpL32 conserved-control demonstration. General claims of biological accuracy
still require a broader real-data evaluation.

## Implemented

- OrthoFinder single-copy orthogroup import;
- genome/GFF extraction and annotation-completion candidates;
- exon-like correspondence using sequence and local structural evidence;
- explicit `exon_presence`, `exon_role`, and `splice_junction` matrices;
- Felsenstein pruning with missing-state marginalization;
- ER/ARD and foreground nested model comparisons;
- profile-likelihood intervals, LRT P values, and BH q values;
- node-state and branch-transition posteriors;
- expected gain/loss counts;
- continuous-probability phylogenetic visualization;
- colorblind-aware color and pattern encodings;
- isolated experimental multi-copy code path.

## Required empirical work

1. **Curated positive single-copy cases**
   Add genes with literature-supported exon gain/loss, exonization, and
   split/fusion histories. Confirm every EG correspondence by inspecting the
   underlying genomic alignment and splice boundaries.

2. **Larger conserved controls**
   Analyze a panel of stable single-copy genes across the same species. Report
   the frequency of unsupported changes and ambiguous correspondences.

3. **Cross-clade evaluation**
   Include at least one vertebrate, plant, or other metazoan clade to show how
   sequence divergence affects internal correspondence.

4. **Comparator analyses**
   Compare annotation-only, sequence-only, and full sequence-plus-structure
   correspondence. Report which missing annotations are recovered and which
   exon relationships change.

5. **Sensitivity analyses**
   Evaluate canonical transcript choice, alignment backend, correspondence
   threshold, tree branch lengths, and exclusion of ambiguous EG members.

6. **Independent biological review**
   Have domain experts review structural matrices and branch posteriors without
   relying on software-generated mechanism labels.

7. **Runtime profiling**
   Report elapsed time and peak memory as species count, exon count, and
   candidate pair count increase.

## Statistical limits to report

- A single gene usually supplies few exon and junction sites, so rate intervals
  can be wide and LRT power low.
- An invariant gene can support a low overall transition rate but cannot
  identify gain/loss asymmetry.
- Neighboring exon and junction sites are biologically coupled, while the
  current likelihood treats sites as conditionally independent.
- Tree topology, branch lengths, and correspondence are conditioned upon.
- A high branch posterior describes model-conditional event placement and does
  not establish a molecular cause.

## Multi-copy work deferred

Duplicated genes require a reconciled gene tree, orthology/paralogy-aware
internal correspondence, and explicit treatment of duplication and loss. The
existing experimental code is retained for development. Formal multi-copy
claims should wait until these components have dedicated real-data tests.

## Minimum first-paper claim

A defensible first paper can present:

- a new formalization of gene-internal synteny as homologous exon and junction
  sites;
- sequence-assisted annotation completion inside known orthologs;
- shared-parameter phylogenetic likelihood models and branch posteriors;
- several deeply curated positive single-copy examples;
- a larger conserved-control panel;
- transparent limitations and sensitivity analyses.
