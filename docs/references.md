# Primary methods and conceptual references

These references motivate comparisons and terminology. They do not validate
IntraPhy's operating characteristics or establish that every implementation rule
is a published standard.

- Sankoff D. (1975). Minimal mutation trees of sequences.
  SIAM Journal on Applied Mathematics. https://doi.org/10.1137/0128004
  Fixed-tree minimum-change dynamic programming.
- Lewis P.O. (2001). A likelihood approach to estimating phylogeny from discrete
  morphological character data. Systematic Biology 50:913–925.
  https://doi.org/10.1080/106351501753462876
  Discrete-character likelihood and variable-character ascertainment.
- Csűrös M. (2008). Malin: maximum likelihood analysis of intron evolution in
  eukaryotes. Bioinformatics 24:1538–1539.
  https://doi.org/10.1093/bioinformatics/btn226
  Corresponding intron positions and gain/loss reconstruction.
- Márquez Y., Mantica F., Cozzuto L. et al. (2021). ExOrthist: a tool to infer exon
  orthologies at any evolutionary distance. Genome Biology 22:239.
  https://doi.org/10.1186/s13059-021-02441-9
  Exon sequence and intron-position context for correspondence; boundary
  conservation filters are not transplanted as evidence that boundaries never change.
- Young N.D. and Healy J. (2003). GapCoder automates the use of indel characters
  in phylogenetic analysis. BMC Bioinformatics 4:6.
  https://doi.org/10.1186/1471-2105-4-6
  Boundary-based indel character coding and applicability; this is not implemented
  as an automatic mutation-interval reconstruction in v18.
- TOGA2 preprint (2026), version 1, posted 30 June 2026.
  https://doi.org/10.64898/2026.06.30.735536
  Protein-coding projection and reference-absent introns. A reference-relative
  difference alone does not identify the phylogenetic gain/loss direction.

For command-line conventions, consult the official IQ-TREE documentation:
https://www.iqtree.org/doc/Command-Reference . IntraPhy has explicit logging,
output ownership and backups; it does not claim IQ-TREE-style checkpoint resume.
