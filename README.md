# INSIPHY

**INSIPHY** means **IN**tragenic **SI**nteny **PHY**logenetics. It is a
Python method package for reconstructing gene-internal structural evolution
from genome sequence, genome annotation and fixed phylogenetic trees.

## 中文介绍

INSIPHY 面向近缘物种之间的单基因或小型重复基因家族比较。上游流程先给出
同源基因或候选同源拷贝集合，推荐来源是 OrthoFinder、OMA、OrthoDB 或人工
整理的 duplication clade。INSIPHY 从这些已给定的 gene/copy 开始，分析每个
基因内部的 segment 同源关系、局部顺序、邻接关系和系统发育结构事件。多拷贝
同源集合可以额外提供 `copy_tree.tsv` 或 `gene_tree.tsv`；EG presence、EG
role、EG adjacency 和 source mixture 会在这棵 gene/copy tree 上推断，copy
multiplicity 保留在物种树上解释。

本项目只使用 genome FASTA、GFF/GTF annotation、gene/copy manifest 和物种树。
已有注释可能不完整，因此软件会结合基因组序列、剪接边界、phase、局部顺序、
拷贝背景和同源 segment graph，补全可能的隐藏片段或注释冲突候选。RNA-seq、
表达量、pathway 富集和全基因组 orthogroup 推断不属于当前输入模型。

核心问题：

1. 在上游已确定同源关系的基因集合内，推断外显子、CDS、UTR、候选外显子化来源片段
   和相邻结构之间的同源性、保守性与局部 synteny。内含子默认作为间隔、
   splice boundary、phase 和 motif 背景证据处理。
2. 在系统发育框架下，量化每个 gene-internal structural character 是否支持
   branch-level 结构变化、变化最可能落在哪条分支、统计支持有多强。INSIPHY
   报告可观测结构变化和统计证据；gene duplication、source joining、
   exonization、splice-boundary shift、segment split/fusion、gene conversion
   等机制解释由用户结合基因树、基因组位置、重复序列、表达或实验资料完成。

## Method Frame

INSIPHY implements five linked stages:

1. **Annotation completion**: genome sequence is checked against annotation to
   identify hidden segments, shifted splice boundaries and joined-segment
   candidates.
2. **Gene-internal element correspondence**: exon-like segment sequence,
   coverage, splice motif, intron phase, strand, boundary class and local order
   are combined into internal evidence clusters and then promoted to
   user-facing EGs when they represent exons, CDS/UTR intervals or
   sequence-supported candidate exonized source intervals. User-facing event
   calls and figures are organized around exon-like structural elements,
   splice boundaries and adjacency.
3. **Tree-guided progressive correspondence**: pairwise support is summarized
   by tree distance, promoted to element-level correspondence summaries, and
   reported as observed copy paths plus ancestral coverage summaries.
4. **Phylogenetic structural inference**: EG presence, EG role state,
   adjacency and source mixture are reconstructed on a supplied copy/gene tree
   when available, with species-tree fallback for single-copy cases. Copy
   multiplicity is reconstructed on the fixed species tree.
5. **Simulation calibration**: optional simulated truth sets summarize false
   positive rate, power, precision/recall, branch placement accuracy and
   bootstrap behavior. Calibration is kept separate from real-data event calls.

The package is a CLI/library. Workflow orchestration, cluster scheduling and
large project execution records stay outside the package.

Terminology boundary:

- **EG / exon-like group**: the user-facing visual and biological correspondence
  unit for exons, CDS intervals, UTRs and candidate exonized source intervals.
- **Internal homology component (`HC_*`)**: an implementation-level graph
  component retained for reproducibility and debugging. It is not a displayed
  biological unit.
- **Intron/context span**: an intronic or non-exonic interval used as splice
  boundary, phase, motif or source-context evidence. It is drawn as background
  context unless an event table supports a role-shift interpretation.

## Real-Data Quick Start

A case manifest contains the gene/copy set supplied by upstream homology
analysis:

```text
case_id	species	family_id	gene_id	gene_copy_id	genome_fasta	annotation_file
```

Build a real case from genome FASTA and GFF/GTF:

```bash
PYTHONPATH=src python3 -m insiphy.cli build-case \
  --manifest examples/real_cases/jingwei/manifest.tsv \
  --species-tree examples/real_cases/jingwei/species_tree.tsv \
  --copy-tree examples/real_cases/jingwei/copy_tree.tsv \
  --output-dir work/jingwei_case \
  --aligner minimap2 \
  --threads 4
```

Run the phylogenetic structural model:

```bash
PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir work/jingwei_case \
  --output-dir results/jingwei \
  --bootstrap-replicates 200 \
  --stochastic-maps 200 \
  --foreground-branches examples/real_cases/jingwei/foreground_branches.tsv \
  --seed 7
```

Generate colorblind-friendly SVG figures:

```bash
PYTHONPATH=src python3 -m insiphy.cli visualize \
  --input-dir work/jingwei_case \
  --result-dir results/jingwei \
  --output-dir results/jingwei_figures
```

The default correspondence encoding is color. Use
`--correspondence-encoding pattern` for texture/line-style encoding when a
color-independent figure is needed.

Check available local alignment backends:

```bash
PYTHONPATH=src python3 -m insiphy.cli inspect-aligners
```

Summarize operating characteristics on simulated truth sets:

```bash
PYTHONPATH=src python3 -m insiphy.cli calibrate \
  --output-dir results/calibration \
  --scenario negative_control \
  --scenario exonization \
  --replicates 20 \
  --bootstrap-replicates 100 \
  --seed 101
```

## Main Outputs

- `case_summary.tsv`
- `annotation_completion_candidates.tsv`
- `element_correspondence.tsv`
- `element_phylogenetic_coverage.tsv`
- `internal_homology_assignments.tsv`
- `internal_homology_graph_edges.tsv`
- `segment_conservation.tsv`
- `segment_correspondence.tsv`
- `progressive_correspondence.tsv`
- `progressive_element_correspondence.tsv`
- `ancestral_element_graph.tsv`
- `ancestral_intragenic_paths.tsv`
- `internal_homology_phylogenetic_coverage.tsv`
- `transcript_paths.tsv`
- `intron_sites.tsv`
- `copy_relationships.tsv`
- `ancestral_state_probabilities.tsv`
- `branch_event_probabilities.tsv`
- `candidate_structural_events.tsv`
- `event_support_summary.tsv`
- `interpretation_hints.tsv`
- `character_model_scores.tsv`
- `model_fit.tsv`
- `hypothesis_tests.tsv`
- `hypothesis_bootstrap.tsv`
- `branch_history_posteriors.tsv`
- `foreground_tests.tsv`
- `phylogeny_scope.tsv`
- `baseline_comparison.tsv`
- `intragenic_graph_edges.tsv`
- `alignment_backend_report.tsv`
- `calibration_operating_characteristics.tsv` from `calibrate`
- `calibration_replicates.tsv` from `calibrate`

`hypothesis_tests.tsv` reports the invariant/no-change null model versus a
one-rate CTMC/Mk model on the active tree for that character, including
likelihoods, LRT statistic, p value, BH q value, fitted rate, AIC and BIC.
`phylogeny_scope.tsv` records whether each layer used `copy_tree.tsv`,
`gene_tree.tsv` or `species_tree.tsv`. The main layers are
`element_presence`, `element_role_state`, `element_adjacency_state`,
`source_mixture` and `copy_multiplicity`. `hypothesis_bootstrap.tsv`
contains empirical p values when bootstrap is requested.
`branch_history_posteriors.tsv` contains stochastic-map summaries for event
placement along branches. `candidate_structural_events.tsv` and
`event_support_summary.tsv` report `structural_change_type`,
`structural_pattern` and `call_scope`, so core structural events, copy-context
evidence and ambiguous paralogous-similarity evidence remain separable.
`interpretation_hints.tsv` is a non-statistical helper table with possible
biological readings and caveats; it is not used for p values, q values,
support tiers or benchmark scores.

`element_correspondence.tsv` and `element_phylogenetic_coverage.tsv` are the
primary biological correspondence tables. `progressive_element_correspondence.tsv`
summarizes near-species, within-clade and deep-tree support for each EG.
`internal_homology_*` tables expose implementation-level graph components for
reproducibility and debugging.

Visualization outputs:

- `intragenic_synteny.svg`: gene-internal exon-like structure by species/copy.
  EG labels mark exon-like correspondence groups, gray spans mark introns or
  other context, and links/ribbons connect corresponding blocks across tracks.
  The default encoding uses a colorblind-aware palette plus labels;
  `--correspondence-encoding pattern` switches to texture, line style and labels.
- `phylogenetic_event_map.svg`: active phylogenetic tree with structural-event
  markers and support summaries.
- `integrated_phylo_synteny.svg`: active phylogenetic tree and exon-like
  synteny tracks in one figure.
- `visualization_manifest.tsv`: figure inventory.

## Documentation

- `docs/input_format.md`: manifest and TSV input formats.
- `docs/method.md`: biological model and computational stages.
- `docs/statistical_model.md`: p values, q values, bootstrap and branch
  posterior interpretation.
- `docs/real_data_benchmark.md`: accession-level Drosophila benchmark plan.
- `docs/algorithm_engineering.md`: runtime, memory, multithreading and code
  quality policy.
- `docs/publication_gap.md`: remaining work before manuscript-scale claims.
