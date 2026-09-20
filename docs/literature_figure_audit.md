# Literature and figure-design audit

The guide uses original vector drawings. It does not reproduce or relabel another
method's artwork. Literature was used to check biological vocabulary, comparison
units and visual reasoning, not to claim validation of IntraPhy.

## Sources and what was examined

| Primary source | Material examined | Design or modeling lesson | What is not imported |
|---|---|---|---|
| Marquez et al., ExOrthist (2021), DOI 10.1186/s13059-021-02441-9 | Full methods and Figure 1/related captions: protein alignment, intron positions, neighboring exons, reciprocal matches | Separate evidence for exon correspondence from downstream interpretation; genomic context resolves otherwise similar matches | ExOrthist's scores, calibrated thresholds or full exon definitions are not claimed as this program's implementation |
| Zea et al., ThorAxe (2021), DOI 10.1101/gr.274696.120 | Full text, algorithm description and Figures 1-2 captions: aligned subexons, genomic order, evolutionary splicing graph | Use corresponding subintervals rather than forcing complete exons into one-to-one boxes | An evolutionary splicing graph is not itself an ancestral transcript model; ThorAxe's clustering/alignment algorithm is not our chain algorithm |
| Csűrös, Malin (2008), DOI 10.1093/bioinformatics/btn226 | Methods/features and intron gain/loss model description | Aligned intron positions become binary characters on a rooted tree; ancestral states and rates have separate interpretations | Malin supports heterogeneity models beyond the current IntraPhy implementation |
| Wang, Devos and Bennetzen (2014), DOI 10.1371/journal.pgen.1004843 | Full text, Methods, Figure 1 caption and Figure 1 image | Put the tree and corresponding gene structures in the same view; inheritance, gain, loss and recurrent change are distinguishable questions | A minimal history is not certainty about mutation mechanism; taxon sampling can alter inferred direction |
| Wertheim et al., RELAX (2015), DOI 10.1093/molbev/msu400 | Model description and Figures 1-2 captions | Draw the biological meaning of parameters and the hypothesis contrast, rather than a software flowchart | Codon selection, omega and relaxation parameters do not describe the IntraPhy model |
| Kowalczyk et al., RERconverge (2019), DOI 10.1093/bioinformatics/btz468 | Methods summary and Figure 1 workflow description | Keep tree-based biological comparison visible alongside the quantities being estimated | Relative protein rates are not structural gain/loss rates |
| TOGA2 (2026), bioRxiv DOI 10.64898/2026.06.30.735536, version 1 | Full preprint: exon alignment/projection, missing annotation and reference-absent introns | Reference-relative structural differences do not establish evolutionary direction | Preprint is not peer-reviewed; IntraPhy does not implement its SpliceAI integration or inherit its benchmarks |

For ExOrthist, ThorAxe and RELAX, the accessible full text and/or figure captions
were inspected; not all high-resolution figure assets could be retrieved. Figure
1 of Wang et al. was directly viewed. These distinctions matter: a caption read
is not a claim that an inaccessible image was inspected. The previous review's
code findings were checked against the present IntraPhy owner modules listed
below rather than treated as proof that methods are identical.

## Resulting hierarchy

1. **Overview:** one orthologous gene family, three biological character classes,
   then their ancestral states/branch changes. The species tree and gene diagrams
   stay together. No stand-alone tutorial panel on the parsimony principle.
2. **Correspondence algorithm:** protein and DNA matches, genomic context, an
   ordered-chain graph, complementary coverage versus competing matches, then
   missing-annotation outcomes. A split/fused pair illustrates one mapping case,
   not the method's exclusive purpose or a direction of evolution.
3. **Probability model:** the same gene diagrams and species become explicit
   matrix rows and named biological columns; an intron character illustrates
   a reversible two-state process; actual computed ancestral probabilities appear
   on the same tree beside the genes.

Detailed assumptions and notation belong in [the model explanation](model_bridge.md)
and figure captions, not six dense explanatory panels. Original coordinates and
annotation provenance remain in the software's data tables.

## Terminology

Use *rooted species phylogeny*, *gene/transcript annotation*, *annotated isoforms*,
*homologous sequence interval*, *intron at a homologous position*, *structural
character matrix*, and *ancestral-state reconstruction*. Internal implementation
terms such as frozen observations or annotated paths need not be displayed to
biological users. Existing schema keys remain stable and are defined explicitly.
Neither *exon identity* nor a bare letter such as `b` should stand in for an
unexplained biological question.

## Code cross-check

`coding_correspondence.py` and `coding/projection.py`: coding alignment projection.
`aligners/short.py`: bounded affine-gap local/global alignment and incomplete
optimal-alignment enumeration handling.
`mapping/chain_graph.py`: order/strand-compatible DAG and best/near-optimal chains.
`observations/roles.py`: DNA and annotation-conditional status.
`observations/junction_matrix.py`: comparable splice positions.
`inference/ctmc.py` and `inference/posterior.py`: transition matrices, pruning,
ancestral probabilities and expected transitions.
`inference/layer_fitting.py`: dependence restriction before independent-site fit.

The illustrated probability values are checked by enumeration versus pruning and
inside/outside messages. This validates the displayed calculation, not biological
accuracy or statistical calibration of the entire method.
