# INSIPHY

**INSIPHY** means **IN**tragenic **SI**nteny **PHY**logenetics. It is a
method package for reconstructing and comparing gene-internal structural
evolution from genome sequence and genome annotation.

## 中文介绍

INSIPHY 面向近缘物种之间的单基因或单基因家族比较。项目输入限定为基因组
序列和已有注释信息，例如 genome FASTA、GFF/GTF annotation 和物种树。方法
重点是基因内部结构的 synteny：在不同物种、不同拷贝之间识别可对应的外显子、
CDS、内含子来源片段和邻接关系，并在系统发育框架下解释这些结构如何变化。
同源基因或候选同源拷贝集合由上游方法提供，例如 OrthoFinder、OMA、
OrthoDB 或人工整理的 duplication clade；INSIPHY 不负责全基因组 orthogroup
推断。

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

Version 0.4 adds manifest-level `source_label`, `copy_role` and `role_hint`
support so source/background and derived copies can be carried into HSG source
mixture inference. Version 0.3 adds transcript-aware extraction, splice/frame-aware hidden segment
scans, graph-based homologous segment correspondence, copy relationship calls,
branch-length-aware CTMC/Mk model fitting and invariant-model LRT p values for
structural characters.

The distributed package is a Python CLI/library. It does not require Nextflow,
Snakemake or a workflow engine. On the R730 server, Nextflow+Slurm was used only
to create formal, auditable real-demo run records under the local project
standards.

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
PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir demos/jingwei \
  --output-dir demo_results/jingwei

PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir demos/sdic \
  --output-dir demo_results/sdic
```

Prepare tables from a genome FASTA and GFF/GTF annotation:

```bash
PYTHONPATH=src python3 -m insiphy.cli extract-gene \
  --genome genome.fa \
  --annotation annotation.gff3 \
  --gene-id GeneA \
  --family-id family_a \
  --species SpeciesA \
  --gene-copy-id SpeciesA_GeneA \
  --source-label source_A \
  --copy-role source \
  --output-dir work/family_a

PYTHONPATH=src python3 -m insiphy.cli derive-tables \
  --input-dir work/family_a \
  --identity-threshold 0.7
```

Add `species_tree.tsv`, then run:

```bash
PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir work/family_a \
  --output-dir results/family_a
```

Run local tests:

```bash
PYTHONPATH=src python3 -m unittest discover -s tests
```

Run a small simulated benchmark:

```bash
PYTHONPATH=src python3 -m insiphy.cli simulate \
  --output-dir simulated/sim_gene \
  --seed 7

PYTHONPATH=src python3 -m insiphy.cli run \
  --input-dir simulated/sim_gene \
  --output-dir simulated/sim_gene_results

PYTHONPATH=src python3 -m insiphy.cli benchmark \
  --input-dir simulated/sim_gene \
  --output-dir simulated/sim_gene_results
```

Prepare a real case from curated local genome and annotation files:

```bash
PYTHONPATH=src python3 -m insiphy.cli inspect-annotation \
  --annotation annotation.gff3 \
  --query jingwei \
  --output-dir work/inspect

PYTHONPATH=src python3 -m insiphy.cli build-case \
  --manifest examples/real_cases/jingwei/manifest.tsv \
  --species-tree examples/real_cases/jingwei/species_tree.tsv \
  --output-dir work/jingwei_case

PYTHONPATH=src python3 -m insiphy.cli scan-hidden-segments \
  --source-fasta source_segments.fa \
  --target-fasta target_gene_interval.fa \
  --output-dir work/hidden_scan
```

## Main Outputs

- `annotation_completion_candidates.tsv`
- `hsg_assignments.tsv`
- `hsg_graph_edges.tsv`
- `segment_conservation.tsv`
- `segment_correspondence.tsv`
- `transcript_paths.tsv`
- `intron_sites.tsv`
- `copy_relationships.tsv`
- `ancestral_state_probabilities.tsv`
- `branch_event_probabilities.tsv`
- `candidate_structural_events.tsv`
- `character_model_scores.tsv`
- `model_fit.tsv`
- `hypothesis_tests.tsv`
- `model_comparison.tsv`
- `baseline_comparison.tsv`
- `intragenic_graph_edges.tsv`
- `demo_summary.tsv`
- real-case preparation outputs: `gene_candidate_report.tsv`,
  `case_provenance.tsv`, `case_build_report.tsv`, `hidden_segment_scan.tsv`

`candidate_structural_events.tsv` includes a biological `event_class` field,
so low-level state changes can be interpreted as exonization, source joining,
new adjacency, segment loss/gain, copy expansion, retrocopy-like or annotation
artifact candidates. `hypothesis_tests.tsv` reports the explicit statistical
test for each structural character: invariant/no-change null model versus a
one-rate CTMC/Mk model on the species tree, with likelihoods, LRT statistic,
p value, fitted rate, AIC and BIC.

Biological sources and demo scope are documented in `docs/data_sources.md`.
Input formats and method details are documented in `docs/input_format.md` and
`docs/method.md`.
Reference-method notes and remaining publication gaps are documented in
`docs/literature_review.md` and `docs/publication_gap.md`.
