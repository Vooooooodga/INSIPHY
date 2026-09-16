# INSIPHY

**INSIPHY** (**IN**tragenic **SI**nteny **PHY**logenetics) 是一个从基因组序列、
基因注释和固定物种树推断基因内部结构演化的 Python 软件包。

## 项目定位

INSIPHY 接收上游已经确定的同源基因集合。推荐使用 OrthoFinder 给出的
single-copy orthogroup，也支持人工整理的单拷贝直系同源基因集合。软件不在全
基因组重新搜索基因，也不重新判定基因层级的 orthology。

当前正式分析范围为单拷贝直系同源基因。多拷贝代码保留在
`--analysis-scope experimental-multicopy`，暂不用于正式结论。

输入只包括：

- genome FASTA；
- GFF3/GTF annotation；
- 每个物种一个基因的 manifest；
- 固定的有根物种树；可选概率分析还需要分支长度或显式的单位分支设置。

RNA-seq、表达量、pathway 富集和机制判定均不在当前模型中。软件报告结构变化
及其统计不确定性；转座、选择、基因转换等机制需要使用者结合额外证据解释。

## 生物学问题

1. **注释补全**：在给定基因及可配置的扩展区间内，用序列相似性、剪接边界、阅读框和局部
   顺序识别漏注释或边界不完整的外显子候选。
2. **基因内部同源对应**：以完整外显子为观察单位，将 CDS、UTR、phase 和
   reading frame 保存为属性，并把外显子边界投影到同源序列坐标。`EG_*` 只是
   表格中的稳定标识符，生物学对象是同源外显子或由明确比对坐标限定的同源
   外显子区块。内含子只提供剪接边界和间隔信息。
3. **系统发育事件定位**：把可观察的基因结构编码为三类二状态位点，在固定有根树
   上寻找所需变化次数最少的历史，保留同等简约的替代解释。速率估计和概率重建
   保留为可选分析。

三类结构位点为：

- `exon_presence`：同源外显子序列明确缺失 / 存在；
- `exon_role`：同源序列存在时，处于非外显子状态 / 外显子状态；
- `splice_junction`：一对同源结构单元之间无剪接连接 / 有剪接连接。

未搜索到、覆盖不足、对应关系含混或注释冲突的观察记为 `unknown`。简约重建允许
未知观察对应两种状态，可选似然分析对未知状态求和。序列存在与外显子身份分别记录；
仅有同源 DNA 命中不能确认外显子身份或实际转录。

## 统计模型

默认 `--model parsimony` 使用等代价获得/丢失的最大简约重建。程序对整棵树的
联合最优历史计算节点状态和分支端点组合；`required` 表示所有最简约历史均要求
该方向的分支变化，`possible` 表示部分最简约历史允许该变化。两者均为条件于
树、结构观察与简约准则的分类，不转换为概率或 P 值。只定位到分支，不估计年代。

下面的连续时间模型仅在指定 `--model er-ard` 或 `--model foreground` 时使用。

每个结构层内的所有同源位点共享参数，计算采用 Felsenstein pruning 和二状态
连续时间 Markov 链。该设计对应分子进化软件中“多个序列位点共同估计少量模型
参数”的基本原则。

可选 ER/ARD 检验为：

- H0 `ER`：获得率与丢失率相等；
- H1 `ARD`：分别估计获得率和丢失率；
- 统计量：`2 * (lnL_ARD - lnL_ER)`；
- 参考分布：满足可估计条件时使用渐近 `chi-square(df=1)`；
- 多个基因或结构层的有效检验按检验类型进行 Benjamini-Hochberg 校正。

`--model foreground` 比较全树共享的 ARD 模型和指定前景分支具有倍率参数的
ARD 模型。零假设下前景倍率为 1。

自动构建的位点按“至少一个物种观察到状态 1”进行检出条件修正。根节点状态 1
频率默认独立估计。若参数不可识别、落在边界、前景与背景混淆或嵌套模型优化
异常，`p_value` 报告为 `NA`。祖先状态和分支概率属于条件于极大似然参数的
经验贝叶斯结果。数值搜索范围受限、参数缺乏约束和优化失败分别报告；未计算
联合参数不确定性时，相应敏感性范围为 `NA`。

## 安装

```bash
python3 -m venv .venv
.venv/bin/pip install -e .
```

外显子对应可选择 MAFFT 或 minimap2。miniprot 只接受重建后的蛋白质查询并映射
到已知同源基因区间。内部短序列比对调用 Biopython，不包含自行实现的动态规划
后备算法。软件包内部不包含工作流引擎。

## OrthoFinder 导入

准备一个资源表：

```text
species	genome_fasta	annotation_file	assembly	annotation	release
```

导入一个 single-copy orthogroup：

```bash
insiphy import-orthofinder \
  --orthofinder-dir path/to/OrthoFinder/Results \
  --orthogroup OG0001234 \
  --genome-manifest genomes.tsv \
  --species-tree SpeciesTree_rooted.txt \
  --output-dir prepared/OG0001234
```

该命令生成标准 `manifest.tsv` 和 `species_tree.tsv`。任一目标物种含零个或多个
基因时，正式单拷贝导入会停止并写出排除原因。

## 完整使用

从真实 genome/annotation 提取基因结构并计算对应关系：

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

运行默认简约事件定位：

```bash
insiphy run \
  --input-dir work/OG0001234 \
  --output-dir results/OG0001234 \
  --analysis-scope single-copy \
  --model parsimony \
  --evidence-aligner miniprot \
  --threads 8
```

`--flank` 与 `--max-extension` 的单位为 bp，用于限定可搜索范围。原注释不被覆盖；
`gene_loci.tsv` 分别记录注释范围和实际搜索范围。

可选前景分支检验：

```bash
insiphy infer-phylogeny \
  --input-dir work/OG0001234 \
  --output-dir results/OG0001234_foreground \
  --model foreground \
  --foreground-branches foreground.tsv \
  --threads 8
```

`foreground.tsv` 可提供 `parent_id` 与 `child_id`，或提供
`branch_scope`（格式为 `parent_label->child_label`）。

生成默认彩色和备选纹理图：

```bash
insiphy visualize \
  --input-dir work/OG0001234 \
  --result-dir results/OG0001234 \
  --output-dir figures/OG0001234

insiphy visualize \
  --input-dir work/OG0001234 \
  --result-dir results/OG0001234 \
  --output-dir figures/OG0001234_pattern \
  --correspondence-encoding pattern
```

默认颜色采用色盲友好调色板，同时保留 `EG_*` 标签和同源连接线。纹理模式使用
斜线、点、网格和线型。

## 正式输出

- `element_correspondence.tsv`：完整外显子或明确同源外显子区块的对应成员；
- `splice_boundary_correspondence.tsv`：投影到同源序列坐标的剪接边界；
- `structural_site_matrix.tsv`：三类结构位点在叶节点的观察状态与证据；
- `node_structural_states.tsv`：默认简约重建允许的祖先状态；
- `branch_structural_events.tsv`：分支上的必要变化及可能变化，保留方向歧义；
- `structural_site_summary.tsv`：每个位点的观察数、缺失观察和最少变化次数。

可选概率分析另外输出：

- `model_fits.tsv`：ER、ARD 或前景模型的参数、profile-likelihood 95% 区间、
  lnL 和 AIC；
- `model_tests.tsv`：零假设、备择假设、LRT、P 值、q 值和可估计状态；
- `node_state_posteriors.tsv`：每个内部节点的状态后验概率；
- `branch_transition_posteriors.tsv`：每条分支两个方向的端点转移后验及期望转移数；
- `structural_changes.tsv`：每个位点与分支的双向变化概率，不自动判定历史方向；
- `excluded_families.tsv`：不满足单拷贝要求的家族；
- `run_parameters.json`：树、模型和运行参数。

图形输出包括基因内部 synteny、系统发育结构变化图，以及树与结构条带的整合图。
默认简约图使用实心符号表示必要变化、空心符号表示可能变化。概率分析的符号
大小和透明度对应条件概率，两种图的图例分别标明含义。

## 真实 Demo

`examples/real_cases/rpl32_control` 包含五种果蝇的真实单拷贝 RpL32 数据。新版按
完整外显子重新分析该案例，不预设其必须完全保守；不同转录本注释造成的结构
差异会在结果中保留为观察或未知状态。

Jingwei 和 Sdic 的多拷贝材料保留为实验性开发案例，当前版本不将它们作为正式
单拷贝方法的性能证据。

v0.13.0 两组真实数据已完成默认简约分析、可选 ER/ARD 分析和默认彩色绘图。
RpL32 的三个外显子组和两个剪接边界在五个物种中均有保守对应。七种蜂的 dsx
仍有大量未知观察，当前无法可靠定位获得、丢失或分裂事件；注释补全报告一个
44 bp 的预测编码区和五处边界冲突候选。这些候选尚缺少独立生物学确认。

结果和 SVG 图位于服务器的
`/data/projects/intragenic_structure/results/20260916_202436_insiphy`；参数、
原始观察数量和推断限制见 [真实数据报告](docs/real_data_benchmark.md)。
两组数据均缺少可用于辨别获得率与丢失率的状态差异，可选检验的 P 值为 `NA`。
这一轮尚未提供已知结构变化阳性案例的灵敏度或准确率评估。

## 文档

- `docs/method.md`：生物学状态、对应关系和计算阶段；
- `docs/statistical_model.md`：似然、模型比较和后验量的数学定义；
- `docs/input_format.md`：输入文件约定；
- `docs/algorithm_engineering.md`：复杂度、并行和内存策略；
- `docs/literature_review.md`：方法所依据的系统发育与基因结构文献；
- `docs/publication_gap.md`：当前可支持的结论及发表前仍需完成的实证工作。
