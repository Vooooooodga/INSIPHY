# INSIPHY

**INSIPHY** (**IN**tragenic **SI**nteny **PHY**logenetics) 是一个面向单拷贝直系同源基因的基因内部结构比较和系统发育事件定位软件包。它接收上游已经确定的同源基因集合、基因组序列、注释文件和一棵有根物种树，在固定系统发育框架下分析同源基因内部结构如何保持、获得、丢失、分裂或融合。

## 方法边界

当前正式范围是单拷贝直系同源基因。推荐上游使用 OrthoFinder 或人工整理得到每个物种一个目标基因的输入集合。输入仅使用基因组序列和已有注释，以及上游同源分组和树；不接收转录组数据。由基因组和注释提取的转录本、CDS 或蛋白 FASTA 可作为比对中间序列。INSIPHY 不重新做全基因组 orthology 推断，不做 pathway 富集，也不替用户判定转座、选择、基因转换等分子机制。

多拷贝代码保留在 `--analysis-scope experimental-multicopy`，用于后续开发。正式解释先围绕单拷贝基因。

核心输入：

- genome FASTA；
- GFF3/GTF annotation；
- 每个物种一个目标基因的 manifest，或 OrthoFinder single-copy orthogroup；
- 有根物种树；可选似然模型需要分支长度或单位分支设置。

## 核心生物学观察

INSIPHY 把基因内部结构拆成三类可观察对象。三类对象分开建模，避免把一个注释缺失或预测候选直接解释成确定事件。

1. `exon_presence`：同源序列单元是否存在。这里关注 DNA 序列层面的存在、明确缺失或未知。
2. `exon_role`：存在的同源序列是否具有外显子角色。仅有 DNA 命中或蛋白投影候选不能自动成为确认外显子。
3. `splice_junction`：两个同源结构单元之间是否存在剪接连接。该层用于描述外显子分裂、融合和边界变化相关的结构模式。

`EG_*` 是同源结构单元的稳定表格标识。生物学对象是完整外显子或由明确比对坐标限定的同源外显子区块。内含子和其他非外显子片段提供边界、间隔、上下文和序列存在证据；它们不会在图中被当作已确认同源外显子块。

## 四个操作阶段

1. **结构提取与注释保留**
   从 genome 和 GFF3/GTF 提取目标基因、所有被选择的转录本路径、外显子、CDS/UTR 属性、内含子、剪接边界和搜索区间。v0.14 的设计以 all-transcript repertoire 为默认分析思想：保留所有注释转录本路径，重复的相同基因组结构去重，真实不同的剪接边界保留。

2. **序列辅助观察补全**
   在目标基因和限定扩展区间内使用外部成熟比对工具和受控短序列比对，记录同源 DNA、蛋白投影、边界冲突和候选预测。`predicted_exon_candidate`、`inferred_role=predicted_CDS`、`predicted_role=CDS` 只表示预测证据，正式图和矩阵不会把它们当作确认 CDS 同源。

3. **基因内部同源对应**
   用比对覆盖、边界投影、阅读框、相位、局部顺序和树上近缘关系构建同源结构单元。分裂和融合允许一般的 `1:n` 或 `n:1` 有序互补投影；多个片段都覆盖同一参考位置时记录为重复或歧义证据。

4. **系统发育事件定位**
   默认使用最大简约法在固定有根树上定位结构状态变化。可选 ER/ARD 或 foreground CTMC 模型报告拟合参数与诊断；LRT P 值要求两个比较模型都满足检验条件，条件后验要求所选模型拟合有效。分支是系统发育树上的边；结构事件是某个结构位点状态在该分支两端发生变化的模型结果。

## 默认统计解释

默认命令使用 `--model parsimony`。`required` 表示所有最少变化历史都需要该分支变化；`possible` 表示部分最少变化历史允许该变化。二者是树、结构观察和简约准则条件下的分类；P 值只出现在满足条件的可选模型比较中。

可选似然模型提供：

- ER 与 ARD 的嵌套 LRT；
- homogeneous 与 foreground rate multiplier 的嵌套 LRT；
- 节点状态和分支端点转移的条件后验；
- 期望转移次数。

P 值只回答两个预设模型之间的拟合差异。它不回答某一次分裂、融合或获得事件是否“真实发生”。参数不可识别、优化失败、边界解或观测位点太少时，概率输出会标记为不可用；可视化不会把无效拟合画成事件概率。

## 安装

```bash
python3 -m venv .venv
.venv/bin/pip install -e .
```

INSIPHY 调用成熟比对库和外部工具：Biopython `PairwiseAligner` 用于受限短序列，MAFFT/minimap2/miniprot 用于相应的序列问题。软件包内部不包含 workflow engine。

## OrthoFinder 导入

资源表：

```text
species	genome_fasta	annotation_file	assembly	annotation	release
```

导入 single-copy orthogroup：

```bash
insiphy import-orthofinder \
  --orthofinder-dir path/to/OrthoFinder/Results \
  --orthogroup OG0001234 \
  --genome-manifest genomes.tsv \
  --species-tree SpeciesTree_rooted.txt \
  --output-dir prepared/OG0001234
```

正式单拷贝导入先通过精确 ID 将成员映射到注释中的基因座。同一基因座的多个转录本或蛋白成员合并计数，随后要求每个目标物种恰好一个不同的基因座。ID 无法解析、零基因座或多个基因座会进入排除记录。

## 基本使用

构建 case：

```bash
insiphy build-case \
  --manifest prepared/OG0001234/manifest.tsv \
  --species-tree prepared/OG0001234/species_tree.tsv \
  --output-dir work/OG0001234 \
  --aligner mafft \
  --context-aligner minimap2 \
  --flank 1000 \
  --max-extension 10000 \
  --threads 8
```

默认简约事件定位：

```bash
insiphy run \
  --input-dir work/OG0001234 \
  --output-dir results/OG0001234 \
  --analysis-scope single-copy \
  --model parsimony \
  --evidence-aligner miniprot \
  --threads 8
```

可选 ER/ARD，使用完整命令在独立结果目录生成对应关系、序列证据和补注释结果，再进行推断：

```bash
insiphy run \
  --input-dir work/OG0001234 \
  --output-dir results/OG0001234_erard \
  --analysis-scope single-copy \
  --model er-ard \
  --evidence-aligner miniprot \
  --threads 8
```

可视化：

```bash
insiphy visualize \
  --input-dir work/OG0001234 \
  --result-dir results/OG0001234 \
  --output-dir figures/OG0001234
```

默认图使用色盲友好颜色表示外显子同源组成员，保留完整注释盒子和各条转录本路径。带状连线仅覆盖直接比对支持的碱基区间；蛋白对应使用投影回 CDS 的区间，不延伸到 UTR。只有同源分组、没有直接比对区间时保留盒子颜色，不绘制碱基对应带。候选、预测和未知结构使用不同边框和浅色样式，不接入确认外显子的连线。`--correspondence-encoding pattern` 提供纹理和线型版本。

## 主要输出

结构与对应关系：

- `element_correspondence.tsv`：同源结构单元成员，区分 `exon_like`、`candidate_source` 和 `absent`；
- `splice_boundary_correspondence.tsv`：投影到同源坐标的剪接边界；
- `observed_element_tree_coverage.tsv`：现生观察到的同源结构单元树上覆盖；
- `observed_intragenic_paths.tsv`：现生转录本路径和元素路径摘要；
- `structural_site_matrix.tsv`：三类结构位点在叶节点的观察状态和证据；
- `raw_gene_features.tsv`：保留搜索区间内的原始注释类型、坐标、父子关系和属性，标明属于目标基因或重叠上下文。保留行为的定向回归测试已通过，首轮真实运行已生成该表；额外结构类型的同源性和演化模型需要分别建立。

默认简约结果：

- `node_structural_states.tsv`；
- `branch_structural_events.tsv`；
- `structural_site_summary.tsv`。

可选概率结果：

- `model_fits.tsv`；
- `model_tests.tsv`；
- `node_state_posteriors.tsv`；
- `branch_transition_posteriors.tsv`；
- `structural_changes.tsv`；
- `run_parameters.json`；
- `excluded_families.tsv`。

旧的 `ancestral_element_graph.tsv` 和 `ancestral_intragenic_paths.tsv` 已从 v0.14 接口中移除。当前输出只描述现生观察和模型条件下的节点/分支状态，不宣称完整祖先转录本图已经联合重建。

## v0.14 状态

软件包版本为 **0.14.0**。运行 `20260918_121100_insiphy` 已完成：**136 项选定正式回归测试通过**（Slurm `61625`，17.063 秒），五案例 × 1/16 线程的十个三模型推断及绘图任务完成，五组比较的六张核心表均一致。复用原有 case/evidence，未重复比对。

v0.14 是一次正确性修复版本。已实现的修复包括：多转录本 repertoire 语义、预测角色与确认角色分离、精确 ID 映射、同源候选可视化区分、无效统计拟合不画成事件概率、现生路径输出重命名、单拷贝简约默认路径、可选似然有效性守卫。

本轮结果及仍存在的限制：

- **RpL32**：20 个结构位点均无观测状态差异、最少变化数为零，仍有未知观察。
- **Hdac3**：精确恢复 62 bp 内含子的结构差异，四物种为 1 有 / 3 无。Dana 分支分裂与 melanogaster 亚群祖先分支融合具有相同简约代价。[分支结果](/data/projects/intragenic_structure/results/20260918_121100_insiphy/hdac3/threads_1/analysis/results_parsimony/branch_structural_events.tsv)
- **rec8：部分恢复。** 边界准确包围 39 bp 目标区间，物种状态为 1 有 / 1 无 / 2 未知；直接内含子 DNA 同源对应和外显子/内含子角色差异仍缺失。Cryo 分支融合与 Octo 分支分裂均为可能历史。[分支结果](/data/projects/intragenic_structure/results/20260918_121100_insiphy/spog_00055/threads_16/analysis/results_parsimony/branch_structural_events.tsv)
- **spo5**：55 bp 目标差异尚未恢复。[结构观察](/data/projects/intragenic_structure/results/20260918_121100_insiphy/spbc29a10_02/threads_1/analysis/results_parsimony/structural_site_matrix.tsv)
- **dsx**：EG_0029 的 Bter 冲突角色已改为未知，保留序列存在；剩余两个观测角色差异，各有五个物种未知。最终角色层 ER 不可识别，六项模型比较的 P/Q 均为 `NA`，节点和分支后验不输出。[最终模型诊断](/data/projects/intragenic_structure/results/20260918_121100_insiphy/dsx_bees/threads_16/analysis/results_er-ard/model_fits.tsv)

已完成 SVG 语义测试和指定连线端点的文件核对，未做实际渲染审阅。上述结果支持具体修复的正确性，尚不足以宣称全面生物学性能或发表就绪。132 项中间基线和 97 项首轮记录均作为历史保留，详见 [真实数据记录](docs/real_data_benchmark.md)。

## 文档

- `docs/method.md`：生物学对象、四个操作阶段和事件定位；
- `docs/input_format.md`：输入文件、ID 映射和输出接口；
- `docs/statistical_model.md`：简约、CTMC、LRT 和后验解释；
- `docs/algorithm_engineering.md`：算法复杂度、外部工具和代码策略；
- `docs/publication_gap.md`：v0.14 修复账本和发表前工作；
- `docs/publication_readiness_review.md`：当前可支持结论和待验收项目。
