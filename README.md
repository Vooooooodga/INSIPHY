# INSIPHY

**INSIPHY** means **IN**tragenic **SI**nteny **PHY**logenetics. It is a
method package for reconstructing and comparing gene-internal structural
evolution from genome sequence and genome annotation.

## 中文介绍

INSIPHY 面向近缘物种之间的单基因或单基因家族比较。项目输入限定为基因组
序列和已有注释信息，例如 genome FASTA、GFF/GTF annotation 和物种树。方法
重点是基因内部结构的 synteny：在不同物种、不同拷贝之间识别可对应的外显子、
CDS、内含子来源片段和邻接关系，并在系统发育框架下解释这些结构如何变化。

当前版本围绕两个核心问题展开：

1. 在近缘物种之间进行基因内部结构的同源推断，判断片段之间的对应关系、
   保守性、拷贝背景和局部顺序是否支持同源。
2. 在系统发育框架下分析演化事件是否涉及基因内部结构变化，并用片段存在、
   外显子化、来源混合、物理邻接、分裂/融合和拷贝数变化来解释事件过程。

INSIPHY 的设计目标是补充传统 exon orthology 或 annotation-transfer 方法在
复杂外显子演化场景下的不足。已有注释可能不完整，因此方法会把序列相似性、
局部 synteny、剪接边界、phase 和拷贝上下文一起纳入证据，报告可能的隐藏片段
或注释冲突候选。表达量、RNA-seq、pathway 富集和转录组组装不属于当前输入模型。

## What It Does

INSIPHY implements three linked tasks:

1. sequence-supported annotation completion;
2. homologous segment group and correspondence scoring;
3. fixed-tree structural inference for segment presence, role, adjacency,
   source mixture, and copy multiplicity.

The first implementation focuses on duplicated and chimeric genes. The bundled
demos are curated method fixtures:

- `demos/jingwei`: Adh/yande source mixture and hidden annotation control.
- `demos/sdic`: AnxB10/sw source mixture, intron-derived segment role change,
  and copy ambiguity.

The demos exercise method behavior and output semantics. They are not complete
accession-level reanalyses.

## Quick Start

Run the bundled demos without installation:

```bash
PYTHONPATH=src python -m insiphy.cli run \
  --input-dir demos/jingwei \
  --output-dir demo_results/jingwei

PYTHONPATH=src python -m insiphy.cli run \
  --input-dir demos/sdic \
  --output-dir demo_results/sdic
```

Prepare tables from a genome FASTA and GFF/GTF annotation:

```bash
PYTHONPATH=src python -m insiphy.cli extract-gene \
  --genome genome.fa \
  --annotation annotation.gff3 \
  --gene-id GeneA \
  --family-id family_a \
  --species SpeciesA \
  --gene-copy-id SpeciesA_GeneA \
  --output-dir work/family_a

PYTHONPATH=src python -m insiphy.cli derive-tables \
  --input-dir work/family_a \
  --identity-threshold 0.7
```

Add `species_tree.tsv`, then run:

```bash
PYTHONPATH=src python -m insiphy.cli run \
  --input-dir work/family_a \
  --output-dir results/family_a
```

Run local tests:

```bash
PYTHONPATH=src python -m unittest discover -s tests
```

Run a small simulated benchmark:

```bash
PYTHONPATH=src python -m insiphy.cli simulate \
  --output-dir simulated/sim_gene \
  --seed 7

PYTHONPATH=src python -m insiphy.cli run \
  --input-dir simulated/sim_gene \
  --output-dir simulated/sim_gene_results

PYTHONPATH=src python -m insiphy.cli benchmark \
  --input-dir simulated/sim_gene \
  --output-dir simulated/sim_gene_results
```

## Main Outputs

- `annotation_completion_candidates.tsv`
- `hsg_assignments.tsv`
- `segment_conservation.tsv`
- `segment_correspondence.tsv`
- `ancestral_state_probabilities.tsv`
- `branch_event_probabilities.tsv`
- `candidate_structural_events.tsv`
- `character_model_scores.tsv`
- `model_comparison.tsv`
- `baseline_comparison.tsv`
- `intragenic_graph_edges.tsv`
- `demo_summary.tsv`

`candidate_structural_events.tsv` includes a biological `event_class` field,
so low-level state changes can be interpreted as exonization, source joining,
new adjacency, segment loss/gain or copy expansion candidates.

Biological sources and demo scope are documented in `docs/data_sources.md`.
Input formats and method details are documented in `docs/input_format.md` and
`docs/method.md`.
Reference-method notes and remaining publication gaps are documented in
`docs/literature_review.md` and `docs/publication_gap.md`.
