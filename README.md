> **0.17.0：可判定范围与单目标组图。** 见[迁移说明](docs/MIGRATION_0.17.md)、[范围政策](docs/scope_policy.md)、[绘图说明](docs/target_views.md)及[本轮验证](validation/v017/README.md)。旧demo结果是历史记录；合成例子不是生物学准确性验证。

# IntraPhy

**IntraPhy** (*Phylogenetic inference of intragenic structure*；原项目名 **INSIPHY**) 是一个面向单拷贝直系同源基因的基因内部结构比较和系统发育事件定位软件包。它接收上游已经确定的同源基因集合、基因组序列、注释文件和一棵有根物种树，在固定系统发育框架下分析同源基因内部结构如何保持、获得、丢失、分裂或融合。

> **兼容性：** Python import namespace 仍为 `insiphy`；旧命令 `insiphy` 继续作为 `intraphy` 的兼容别名。

## 方法边界

当前正式范围是单拷贝直系同源基因。推荐上游使用 OrthoFinder 或人工整理得到每个物种一个目标基因的输入集合。输入仅使用基因组序列和已有注释，以及上游同源分组和有根树；不接收转录组数据。由基因组和注释提取的转录本、CDS 或蛋白 FASTA 可作为比对中间序列。IntraPhy 不重新做全基因组 orthology 推断，不做 pathway 富集，不推断组织或条件特异的转录本使用，不联合重建完整祖先转录本，也不替用户判定转座、选择、基因转换等分子机制。

多拷贝代码保留在 `--analysis-scope experimental-multicopy`，用于后续开发。正式解释先围绕单拷贝基因。

核心输入：

- genome FASTA；
- GFF3/GTF annotation；
- 每个物种一个目标基因的 manifest，或 OrthoFinder single-copy orthogroup；
- 有根物种树；可选似然模型需要分支长度或单位分支设置。

## 核心生物学观察

IntraPhy 把基因内部结构拆成三类可观察对象。三类对象分开建模，避免把一个注释缺失或预测候选直接解释成确定事件。

1. `exon_presence`：同源序列单元是否存在。这里关注 DNA 序列层面的存在、明确缺失或未知。
2. `exon_role`：存在的同源序列是否具有外显子角色。仅有 DNA 命中或蛋白投影候选不能自动成为确认外显子。
3. `splice_junction`：两个同源结构单元之间是否存在剪接连接。该层用于描述外显子分裂、融合和边界变化相关的结构模式。

`EG_*` 是同源结构单元的表格标识。它只在冻结的输入、参数和同源对应算法版本内保持稳定，不能跨数据修订或算法版本当作永久生物学编号。生物学对象是完整外显子或由明确比对坐标限定的同源外显子区块。内含子和其他非外显子片段提供边界、间隔、上下文和序列存在证据；它们不会在图中被当作已确认同源外显子块。

## 四个操作阶段

1. **结构提取与注释保留**
   从 genome 和 GFF3/GTF 提取目标基因、所有被选择的转录本路径、外显子、CDS/UTR 属性、内含子、剪接边界和搜索区间。默认使用 all-transcript repertoire：保留所有注释转录本路径，重复的相同基因组结构去重，真实不同的剪接边界保留。

2. **序列辅助观察补全**
   在目标基因和限定扩展区间内使用外部成熟比对工具和受控短序列比对，记录同源 DNA、蛋白投影、边界冲突和候选预测。短 DNA 的 feature-bounded 比对只产生候选；只有找到唯一、有序、位于同一转录本路径的双侧 flank 后，程序才提取两 flank 之间的基因组区间并重新比对。通过该锚定步骤的结果才可形成硬核酸 membership 或位置观察。`predicted_exon_candidate`、`inferred_role=predicted_CDS`、`predicted_role=CDS` 只表示预测证据，正式图和矩阵不会把它们当作确认 CDS 同源。

3. **基因内部同源对应**
   用比对覆盖、边界投影、阅读框、相位、局部顺序和树上近缘关系构建同源结构单元。确认的硬 membership 由 element、occurrence 和实际 matched subinterval 共同限定；一个 parent occurrence 可以贡献多个局部 membership。候选或未解析对应不创建硬观察。分裂和融合允许一般的 `1:n` 或 `n:1` 有序互补投影；多个片段都覆盖同一参考位置时记录为重复或歧义证据。

4. **系统发育事件定位**
   默认使用最大简约法在固定有根树上定位结构状态变化。`1↔2` 模式同时报告 junction 层的 intron gain/loss 和 `structural_relation=split/fusion`。`1↔n`（`n>=3`）只有在多个 cutpoint 于同一物种、同一 gene copy、同一 transcript path 共现时才总结为 `required` compound；缺少联合路径证据时标为 `possible_non_joint`。可选 ER/ARD 或 foreground CTMC 模型报告拟合参数与诊断；LRT P 值要求两个比较模型都满足检验条件，条件后验要求所选模型拟合有效。

## 默认统计解释

默认命令使用 `--model parsimony`。`required` 表示所有最少变化历史都需要该分支变化；`possible` 表示部分最少变化历史允许该变化。二者是树、结构观察和简约准则条件下的分类；P 值只出现在满足条件的可选模型比较中。

可选似然模型提供：

- ER 与 ARD 的嵌套 LRT；
- homogeneous 与 foreground rate multiplier 的嵌套 LRT；
- 节点状态和分支端点转移的条件后验；
- 期望转移次数。

P 值只回答两个预设模型之间的拟合差异。它不回答某一次分裂、融合或获得事件是否“真实发生”。CTMC 至少需要两个明确状态的 tips，且位点必须符合所选 ascertainment rule。一个 `linked_group_id` 含多个位点时，当前独立位点模型不计算该层 LRT，并报告 `correlated_linked_sites_not_modelled`。受检对比还需要有限的双侧 profile 区间。单基因少位点条件下的 P 值准确度尚未评估。

## 安装

```bash
python3 -m venv .venv
.venv/bin/pip install -e .
```

IntraPhy 调用成熟比对库和外部工具：Biopython `PairwiseAligner` 用于受限短序列，MAFFT/minimap2/miniprot 用于相应的序列问题。软件包内部不包含 workflow engine。

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
  --transcript-policy all \
  --coding-msa-mode linsi \
  --short-context-max-length 300 \
  --flank 1000 \
  --max-extension 10000 \
  --threads 8
```

先生成并冻结默认 repertoire 观察矩阵，同时完成默认简约事件定位：

```bash
insiphy run \
  --input-dir work/OG0001234 \
  --output-dir results/OG0001234 \
  --analysis-scope single-copy \
  --model parsimony \
  --annotation-view repertoire \
  --evidence-aligner miniprot \
  --threads 8
```

将 `results/OG0001234/structural_site_matrix.tsv` 作为本次分析的冻结 repertoire 矩阵。ER/ARD 和 foreground 必须通过 `--structural-site-matrix` 读取这一份文件，避免重新构造观察状态：

```bash
insiphy infer-phylogeny \
  --input-dir work/OG0001234 \
  --output-dir results/OG0001234_erard \
  --analysis-scope single-copy \
  --model er-ard \
  --annotation-view repertoire \
  --structural-site-matrix results/OG0001234/structural_site_matrix.tsv \
  --threads 8
```

```bash
insiphy infer-phylogeny \
  --input-dir work/OG0001234 \
  --output-dir results/OG0001234_foreground \
  --analysis-scope single-copy \
  --model foreground \
  --annotation-view repertoire \
  --foreground-branches foreground_branches.tsv \
  --structural-site-matrix results/OG0001234/structural_site_matrix.tsv \
  --threads 8
```

`canonical` 是独立的 view-specific 敏感性分析。它应在单独结果目录生成并冻结自己的矩阵；该矩阵不能与 repertoire 矩阵混用。

```bash
insiphy run \
  --input-dir work/OG0001234 \
  --output-dir results/OG0001234_canonical \
  --analysis-scope single-copy \
  --model parsimony \
  --annotation-view canonical \
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

默认每个目标单独输出结构、目标ribbon与条件历史三张图。所有图复用完整灰色原生注释背景；只高亮当前目标，其余区域保持灰色。不同target不同颜色、同一target三图同色并有图例。支持`--target FAMILY/layer/site`与`--target-manifest targets.tsv`。ribbon仅覆盖该目标实际匹配块，不从相同ID虚构连线。需要旧混合总览时使用`--layout legacy-overview`；`--correspondence-encoding pattern`仍适用于旧总览。


### 可判定范围（0.17）

默认不丢低覆盖字符；明确0与1均属于可判定观测。`--analysis-range high-coverage --min-callable-fraction 0.7`仅选择本次推断子集，不改变坐标或背景。完整矩阵、实际子集、逐物种原因及树上覆盖分别输出。

### 可运行合成示例

```bash
python tools/build_v017_example.py examples/synthetic_v017
# 解压后的浏览器打开 examples/synthetic_v017/figures_all/index.html
```

示例含负链、E3分段对应、E4有无与等简约历史、无转录本归属DNA、低覆盖、重叠／反链／UTR／ncRNA背景。对应块与状态是人工给定输入，检验范围、推断与绘图契约，不声称aligner独立恢复了这些真值。

## 主要输出

结构与对应关系：

- `element_correspondence.tsv`：同源结构单元成员，区分 `exon_like`、`candidate_source` 和 `absent`；
- `splice_boundary_correspondence.tsv`：投影到同源坐标的剪接边界；
- `observed_element_tree_coverage.tsv`：现生观察到的同源结构单元树上覆盖；
- `observed_intragenic_paths.tsv`：现生转录本路径和元素路径摘要；
- `structural_site_matrix.tsv`：三类结构位点在叶节点的观察状态和证据；
- `raw_gene_features.tsv`：保留搜索区间内的原始注释类型、坐标、父子关系和属性，标明属于目标基因或重叠上下文；额外结构类型只有在同源关系和状态定义明确时才进入正式模型；
- `alignment_backend_report.tsv`：记录比对 mode、backend 和执行来源；精确 executable version 由容器或执行 provenance 提供。

默认简约结果：

- `node_structural_states.tsv`；
- `branch_structural_events.tsv`；
- `structural_site_summary.tsv`；
- `compound_structural_events.tsv`：汇总同一路径上联合支持的多 cutpoint 结构模式及 `possible_non_joint` 情形；
- `phylogeny_scope.tsv`：记录每层总位点数、实际纳入位点数、annotation view 和不可用原因。

可选概率结果：

- `model_fits.tsv`；
- `model_tests.tsv`；
- `node_state_posteriors.tsv`；
- `branch_transition_posteriors.tsv`；
- `structural_changes.tsv`；
- `run_parameters.json`；
- `excluded_families.tsv`。

旧的 `ancestral_element_graph.tsv` 和 `ancestral_intragenic_paths.tsv` 已从当前接口中移除。当前输出只描述现生观察和模型条件下的节点/分支状态，不宣称完整祖先转录本图已经联合重建。

## 历史：v0.15 状态（非本次版本）

以下段落描述原0.15开发阶段；当前软件包为**0.17.0**。v0.15 的五案例正式运行尚未执行，结果状态记为 **pending**。在正式运行完成前，README 不提供 v0.15 的案例数值、恢复率、线程一致性或模型检验结论。

RpL32 预定作为定性保守对照；spo5、rec8 和 Hdac3 作为可观察结构差异的阳性 benchmark；dsx 只作描述性案例。三个阳性案例现有 manifests 使用人工整理的基因分组，尚无独立的上游全基因组单拷贝资格判定，因此正式结果必须在这一条件下解释。

运行 `20260918_121100_insiphy` 的全部数值属于 **v0.14 历史基线**，用于记录旧实现的行为，不能替代 v0.15 的正式评估。历史明细和 v0.15 待填写位置见 [真实数据记录](docs/real_data_benchmark.md)。

## 文档

- `docs/method.md`：生物学对象、四个操作阶段和事件定位；
- `docs/input_format.md`：输入文件、ID 映射和输出接口；
- `docs/statistical_model.md`：简约、CTMC、LRT 和后验解释；
- `docs/algorithm_engineering.md`：算法复杂度、外部工具和代码策略；
- `docs/publication_gap.md`：v0.14 历史修复账本和发表前工作；
- `docs/publication_readiness_review.md`：当前可支持结论和待验收项目。
