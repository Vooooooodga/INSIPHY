# INSIPHY 单拷贝结构演化：下一轮实施交接规格

日期：2026-09-18。核对基线：v0.14.0，提交 `8428b2d`。

**文档性质：待实施规格。本文新增接口、参数和行为均为计划，不代表当前版本已经实现。**
本轮仅编写方案，没有修改算法、运行测试、重新分析 demo、下载数据或发布版本。
实施目标建议为 v0.15.0；由下一位总工程师完成集成后决定是否具备发布条件。

## 0. 总工程师首先要作出的决定

采用以下固定范围，避免各个执行者分别设计一套方法：

1. 输入为上游确定的单拷贝同源基因、基因组 FASTA、已有注释和用户提供的有根树。默认接收 OrthoFinder 结果；软件内部不搜索整基因同源关系。
2. 全部已提供的转录本保留。研究对象首先是基因座的已注释结构集合，另提供转录本路径信息。软件不推断组织、性别或表达量。
3. 默认分析为固定树上的等代价最大简约推断，回答哪些结构差异存在、哪些分支定位得到全部或部分最简历史支持。
4. 保留可选 CTMC 似然分析。速率参数服务于这一概率模型；速率差异检验不列为每个单基因 demo 的必做成功条件。
5. 首轮完成下面 F1-F5、输出语义和真实案例重跑。复杂重复外显子簇保留观测，局部歧义明确报告；整基因多拷贝代码保留，暂不扩展。
6. 原始外显子、内含子、CDS、UTR 和其他已提供序列要素都保留其生物学身份。比对所需的局部子区间只表示对应范围，不自动成为一种新生物学元件。
7. 本轮不开发新的序列比对核心、不推断分子机制、不运行模拟校准、不以显著 P 值或吻合文献方向作为强制验收目标。

### 0.1 五项主问题及交付物

| ID | 已见问题 | 必须交付的改变 | 直接验收案例 |
|---|---|---|---|
| F1 | 短非编码序列沿用长序列种子/分数条件，真实短区间可能完全没有命中 | 基因内部锚点限定的短序列比对；分开序列对应与剪接位置对应 | spo5 55 bp、rec8 短区间 |
| F2 | 整个 occurrence 的坐标和角色被用于解释其中很短的比对片段 | 所有对应携带实际块坐标、父元件和转录本；冲突按实际交集判定 | dsx EG0013、EG0055、EG0029 |
| F3 | 多命中信息在表转换中丢失，重复和分裂的对应关系混杂 | 候选全集在接口内保留；有序对应链；重复覆盖与互补覆盖分开 | dsx、后续 PKM；确定性坐标回归 |
| F4 | 整段蛋白 0.70 identity 门槛阻断较远近缘物种的局部对应 | 家族蛋白 MSA、公用列坐标、局部证据规则；取消单一 identity 作为唯一入口 | spo5、rec8、保守 RpL32 |
| F5 | 外显子身份、预测、非编码身份和注释缺失混用 | 转录本级角色、基因座级观测、候选注释三者分开；未覆盖区域保持未知 | dsx、末端截短注释案例 |

**依赖顺序：坐标与表契约 → F2/F3 基础接口 → F1/F4 比对 → F5 状态 → 系统发育输出 → 从原始数据重跑。**
F1/F4 的算法工作可在接口冻结后并行。F5 的最终接入必须等待实际块坐标和候选歧义可用。

### 0.2 已有修复必须保留

外部 reviewer 审查的 `2c3e74c` 与当前基线不同。以下问题已有修复或部分修复，执行者应先读当前函数，避免重复覆盖：

| 内容 | 本轮处置 |
|---|---|
| 已知非外显子同源片段、明确缺失进入结构观测 | 保留 `candidate_source/absent` 记录；补实际区间语义 |
| 同源 DNA 与预测外显子分离 | 保留；落实到全部转录本和局部子区间 |
| 真实缺失要求有序锚点及区间证据 | 保留；进一步排除中间有未知碱基、比对不唯一的情况 |
| 最长 CDS 使用 CDS 覆盖并集 | 保留，不能退回转录本基因组跨度 |
| 剪接位置按实际序列映射、负链和 phase 处理 | 保留；不恢复长度比例投影和按外显子长度合并边界 |
| 树的连通性、根、叶标签及枝长检查 | 保留在输入入口；零枝长转移矩阵严格为单位阵 |
| 拟合有效性、LRT 可用性和后验输出分开 | 保留；禁止无效拟合回退后产生正常外观的后验 |
| ascertainment 规则和输入筛选配套 | 保留；新增表字段必须穿透到条件似然 |
| 多转录本路径保留、现生结构摘要正确命名 | 保留，不将现生多数路径改称联合祖先结构 |
| 图上的部分区间连线、候选与已知结构区分 | 保留，增加本规格中的实际渲染验收 |
| dsx EG0029 的 44 bp 预测与目标 300 bp 对应区间相交 | 这是需要保留的局部冲突；F2 修复不能将该条直接恢复为确定非外显子 |
| 旧 experimental-multicopy 校准不代表正式单拷贝模型 | 保留范围说明；不启动旧校准命令 |

已有“136 项选定回归通过”仅描述上一轮选定命令。最后一次五案例运行主要重新计算系统发育，没有重新生成全部比对证据。上述记录均不能替代本轮前端改变后的真实案例分析。

## 1. 分析层级与代码归属

保持五层，前三层构成用户最关注的生物学主体：

```text
FASTA + GFF/GTF + 已知单拷贝关系 + 有根树
  1. 原始元件与注释路径                    preprocess / annotation
  2. 序列对应、局部补注释候选              alignment / coding_correspondence
  3. 元件对应与可用于分析的结构观测        correspondence / structural_sites
  4. 固定树上的历史推断及可选模型比较      parsimony / structural_phylogeny
  5. 生物学结果表、树与基因内共线性图      visualization / docs
```

一条对应证据经过各层时始终携带来源，不能只传一个混合分数。
所有正式推断共用一次生成的 `structural_site_matrix.tsv`；parsimony 与 CTMC 不再各自重新解释注释。

### 1.1 现有入口及修改位置

| 模块/函数 | 修改任务 |
|---|---|
| `case.build_case` | 提取后先建立公用编码坐标，再生成对应；保持原有 manifest 和树入口 |
| `preprocess.extract_gene` | 原始基因范围、实际搜索范围、feature 所属关系及转录本路径明确保存 |
| `preprocess.match_evidence`、`MATCH_FIELDS` | 不再丢弃 alternative hits、块坐标、分子类型、未完成搜索原因 |
| `preprocess._projection_record/_projection_compatibility` | 使用真实投影区间判断一对多；不靠 occurrence 长度推测覆盖 |
| `preprocess.cluster_segments/derive_tables` | 将候选生成、对应链选择、元件归组拆为可独立调用的步骤 |
| `alignment.AlignmentStats` 及现有适配器 | 外部比对结果统一；短区间使用现有成熟库；同分映射保留 |
| `coding_correspondence.CodingProjectionIndex` | 从转录本两两 MAFFT 缓存改为家族 MSA 投影索引 |
| `correspondence.infer_correspondence` | membership 记录实际匹配子区间和父元件；候选、确定对应分开 |
| `annotation.generate_sequence_evidence/complete_annotation` | DNA 存在、预测角色、已有注释角色独立传递 |
| `structural_sites._element_site_rows` | 根据局部角色和覆盖范围构造状态 |
| `structural_sites._exon_prediction_overlaps_observation` | 以匹配块和预测 CDS 块逐块求交，保留 EG0029 真重叠 |
| `structural_sites._junction_site_rows` | 精确切点、相位、跨越覆盖及转录本来源 |
| `parsimony._parsimony_tables/_event_type/infer_single_copy_parsimony` | 保留全体最简历史的边际状态集合；统一矩阵读写与输出语义 |
| `structural_phylogeny.infer_single_copy_phylogeny` | 直接消费冻结矩阵，保持既有概率内核；补充统计适用条件 |
| `cli.py` | 新参数集中解析一次；命令之间共享同一个配置对象 |

函数名称以当前源码为准；拆函数可以，改名后需要保留明确的调用迁移。不要平行创建第二套 `v2_pipeline`。

## 2. 数据与坐标契约

### 2.1 坐标只允许在明确边界转换

1. 现有面向用户的 GFF/TSV `start/end` 继续为 **1-based、闭区间、参考基因组正向坐标**，`start <= end`。
2. 新的比对内存对象使用 **0-based、半开区间**；字段统一以 `start0/end0` 命名。禁止把两种约定放入无说明的同一四元组。
3. 旧字段通过唯一适配函数转入新对象；所有 TSV 写出经过唯一逆转换。无需重写未涉及的多拷贝模块。
4. `gene_strand`、`alignment_relative_strand` 分开；基因负链与序列比对反向具有不同含义。
5. GFF 长度为 `end-start+1`；半开长度为 `end0-start0`。闭区间 `[a,b]` 对应半开 `[a-1,b)`。
6. 对基因组区间 `[L,R)` 的转录方向序列，局部 `[u,v)` 映回正链为 `[L+u,L+v)`，负链为 `[R-v,R-u)`。
7. 剪接切点位于碱基之间，单独存 `cut0`。不要用相邻外显子的任意一个碱基坐标代替切点。

把这些适配函数放在现有坐标处理最集中的模块中；若确需共享，仅增加一个小型 `coordinates.py`。函数处理区间转换、方向和块求交，不负责生物学判定。

### 2.2 原始对象与分析对象

| 对象 | 必需字段 | 语义 |
|---|---|---|
| GeneLocus | family/species/gene_copy、contig、strand、original_bounds、search_bounds | 一个已确定同源关系的基因座 |
| TranscriptPath | transcript_id、parent_gene_id、ordered_feature_ids、CDS 状态、原注释属性 | 一条提供的注释路径 |
| FeatureOccurrence | occurrence_id、source_feature_id、feature_type、parent IDs、interval、phase | 不覆盖原注释的结构实例 |
| AlignmentCandidate | candidate_id、query/target IDs、blocks、score_scheme、backend、替代候选 | 一种局部序列对应 |
| ElementMembership | element_id、occurrence_id、transcript_id、matched_blocks、parent_feature_id、repeat_instance_id | 元件实际对应范围 |
| StructuralObservation | site_id、species、layer、state、applicability、reason、evidence_ids、transcript_scope | 进入系统发育的观测 |

`element_id` 是追踪对应关系的标识，显示名称仍然使用外显子、内含子位置、UTR 子区间等生物学名称。
同一个原始外显子可有多个 membership；它们不会凭拆分次数自动增加独立事件数。

同一 DNA 区间允许多个注释类型叠加。例如某段在一个转录本内为 CDS，在另一条路径内为内含子。不要把基因强行划成互斥且覆盖每个碱基的功能类别。
调控区、嵌套 ncRNA 等仅在输入注释实际提供时保留其类型和所属基因；邻近或重叠基因的 feature 不自动归入目标基因。

### 2.3 AlignmentCandidate 的具体字段

在现有 `AlignmentStats` 周围增加适配后的候选集合，避免每层复制大量字典：

```python
@dataclass(frozen=True)
class AlignedBlock:
    query_start0: int
    query_end0: int
    target_start0: int
    target_end0: int

@dataclass(frozen=True)
class CandidateSet:
    candidates: tuple
    enumeration_complete: bool
    incomplete_reason: str
```

上述示意兼容项目 Python 3.9 基线。正式字段类型使用项目现有风格，不引入新的对象验证框架。
每个候选至少还需保存：

```text
candidate_id, query_occurrence_id, target_occurrence_id
query_transcript_id, target_transcript_id
sequence_kind = nucleotide | amino_acid
backend, backend_version, score_scheme, raw_score
nt_identity, aa_identity, known_aligned_pairs, unknown_aligned_pairs
query_covered_bases, target_covered_bases, query_length, target_length
aligned_blocks, gap_blocks, relative_strand
mapping_quality, is_secondary, hit_count, alternative_candidate_ids
left_anchor_id, right_anchor_id, search_interval
enumeration_complete, incomplete_reason
```

`nt_identity` 与 `aa_identity` 不互相填写。N/X 不计入已知匹配数，未知数量单独报告。
coverage 用**覆盖并集长度/相应对象长度**计算；外显子部分对应还要报告匹配子区间覆盖，不能以较短一方覆盖充分宣布整个长外显子全部对应。
短序列库没有 MAPQ 时输出 NA，不填写 0.99 之类伪概率。
同分最优比对的数量、不同基因组命中位置数量、近优对应链数量分别命名。

### 2.4 输出表迁移

优先扩展现有表。实际文件名在任务开始时按代码写出函数确认，避免创建同义表：

| 表 | 增加或统一的内容 |
|---|---|
| `segment_occurrences.tsv` | 原 feature、原基因边界、搜索边界、转录本来源、原注释 partial 属性 |
| 现有 match evidence 表 | 候选 ID、替代候选、真实块坐标、分子类型、score_scheme、搜索完整性 |
| `element_correspondence.tsv` | 父元件、matched_blocks、resolved/ambiguous/candidate、歧义原因 |
| `annotation_completion_candidates.tsv` | DNA 存在证据、预测角色、实际预测块、与原注释冲突的交集 |
| `structural_site_matrix.tsv` | 下述完整观测字段及发现规则 |
| 现有节点/分支结果表 | inference_method、support_kind、条件范围、定位集合及 unavailable_reason |

观测字段统一为：

```text
family_id, site_id, layer, species, state
applicability, observation_reason, transcript_scope
parent_feature_ids, member_interval_ids, evidence_ids
discovery_rule, discovery_species, observation_mask
linked_group_id, site_kind, observation_source, annotation_completeness
```

`state` 为 `0/1/?`；`applicability` 为 `applicable/inapplicable/undetermined`。
同一字符的 `observation_mask` 在模型计算入口构造一次；不要逐表重算成不同物种集合。
当前 parsimony 和 CTMC 写矩阵时存在较短列清单；必须共用一个 writer，保证新增元信息不被静默丢弃。
新表含 schema 版本；旧输入只经过一个兼容适配器。旧数据缺实际块坐标时保持 unavailable，不用父区间伪造块。
ID 由有序的来源字段和稳定编号生成；不增加摘要算法或文件指纹系统。

## 3. F4：编码区的共同坐标与同源证据

### 3.1 执行顺序

1. 对每条输入注释路径构建 CDS，保持转录顺序、负链方向和遗传密码表。
2. 保留 CDS 每个碱基到原 feature/基因组的映射。跨外显子密码子由拼接 CDS 的实际三个碱基组成。
3. GFF phase 用于检查阅读框连续性及处理已标注的部分起始 CDS；不能逐个 CDS 删除 phase 个碱基后再拼接。
4. 完整 CDS 末端终止密码子单独处理；内部终止、未知密码子和 frame 冲突记录来源。局部不可翻译不应使其 DNA 和原注释记录消失。
5. 同一基因座中完全相同的蛋白只送一条序列进 MSA，保存其全部 transcript 别名和不同基因组路径。仅序列相似不能合并转录本。
6. 一个家族建立一次蛋白 MSA；生成 `residue -> MSA column -> CDS bases -> genome` 投影索引。
7. 将同一 MSA 的局部覆盖用于全部 coding occurrence 对应。两两外显子 DNA 比对作为补充证据，不再先决定哪些蛋白能够被比较。

这里采用跨物种共同列、保留转录路径的思想。文献已明确讨论转录本边界细分与相同编码区的多路径表示；工程实现仍需保留原始外显子身份。[跨物种剪接图方法论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8327911/)

### 3.2 MAFFT 调用约定

新增 `--coding-msa-mode linsi|einsi`，默认 `linsi`。实际安装版本先读本地帮助，以下为官方文档对应的目标调用：

```text
linsi: mafft --amino --localpair --maxiterate 1000 --thread T --threadit 0 INPUT
einsi: mafft --amino --genafpair --ep 0 --maxiterate 1000 --thread T --threadit 0 INPUT
```

L-INS-i 面向较少序列的精细局部比较；E-INS-i 用于较长不可比区间。保留默认蛋白替换矩阵并记录版本，不自行混合 DNA 和蛋白评分。[MAFFT 官方参数说明](https://mafft.cbrc.jp/alignment/software/manual/manual.html)
`--threadit 0` 固定迭代阶段的并行行为；官方说明多线程迭代可能使重复运行不同。[MAFFT 多线程说明](https://mafft.cbrc.jp/alignment/software/multithreading.html)

MSA 使用稳定排序的输入 ID。MAFFT 的 guide tree 属于比对算法内部；历史推断仍使用用户提供的有根树。
借鉴 progressive alignment 的方式体现在家族公共坐标和局部逐步细化；本轮不运行全基因组 Cactus/HAL，不新增全基因组输入需求。
去重后超过约 200 条蛋白时报告序列规模；这一数字是精细模式的使用建议，不作为删除转录本的规则。仍按显式指定的模式运行并记录资源；极大重复簇按第 12 节报告局部分辨限制。本轮只承诺上述两个编码 MSA 模式，不暗中切换尚未接入的模式。

### 3.3 局部证据怎样评分

保留三类量，禁止混成一个“同源概率”：

1. **序列分数**：既有比对算法给出的 score，附 scoring scheme。
2. **对应证据等级**：以下确定性规则的结论，附每条规则是否满足。
3. **历史支持**：后续树上最简历史集合或条件概率，附模型。

编码区候选使用共同 MSA 列上的局部残基对应，记录 identity、替换矩阵得分、gap 比例和两侧可对应列。
取消 `MIN_AA_IDENTITY=0.70` 作为全体编码区的强制排除条件；保留 `--identity-threshold` 对明确使用它的 DNA 模式的含义，避免同名参数继续暗中控制蛋白。

默认对应等级的具体规则：

| 等级 | 必需条件 | 后续用途 |
|---|---|---|
| `resolved_local` | 共同 MSA 中有实际残基对应；同一基因内位置相容；不存在竞争重复位置；有双侧锚点，或完整末端加单侧锚点 | 仅实际覆盖子区间可赋观测 |
| `supported_unanchored` | 局部序列匹配但无足够位置约束 | 候选表与细节图，暂不硬赋祖先分析状态 |
| `ambiguous_mapping` | 两种可接受对应给出不同目标区间/切点 | 按受影响范围写未知，已一致覆盖区域仍可使用 |
| `uncovered` | 所问区间没有可用映射 | 未知并注明覆盖原因 |

MSA 中“两个残基排在同一列”本身不足以确认同源。`resolved_local` 还须通过独立局部锚点检查：每侧默认取至多 15 个可用残基，至少 8 个已知配对，BLOSUM62 配对分数和大于 0，且在该基因的候选位置中无并列竞争。长度不足时输出 `anchor_short`，仍保留剪接位置的其他证据。
这里的 15/8 和分数阈值是**预先声明的工程判据**，有效范围尚待真实数据比较；不能写成文献已验证的普适灵敏度。对微外显子不要求其自身有 8 个残基，锚点来自其两侧。
前后 phase 相容用于编码对应佐证。所研究的 phase/边界变化不能同时作为强制否决条件，否则会系统性过滤真实变化。

函数目标接口：

```text
build_family_coding_alignment(transcripts, alignment_options)
  -> FamilyCodingProjection
FamilyCodingProjection.project(occurrence_id)
  -> 实际 MSA 列与基因组块
FamilyCodingProjection.candidates(query_id, target_copy_id)
  -> CandidateSet
```

对 CDS/UTR 混合外显子，仅 CDS 部分获得蛋白证据；UTR 使用 DNA 证据或保留未覆盖。

## 4. F1：短片段、内含子和非编码对应

### 4.1 双侧锚点限定搜索

同源搜索始终在给定的基因及其声明的扩展范围内。这里的左右锚点来自目标基因内部可对应片段，不要求相邻基因共线性。

```text
确定左右同源锚点
  -> 在各物种转录方向上定位两锚点之间的实际区间
  -> 分开“有没有同源 DNA”和“此位置有没有剪接间隔”两个问题
  -> 短区间用无种子依赖的成熟配对比对
  -> 将实际块及替代 gap 位置送回 correspondence
```

锚点必须同 contig、方向相容、顺序相容，并实际夹住目标区间。只检查两个 occurrence ID 存在于同一基因不够。
未找到足够锚点时保留候选搜索结果；不赋 true absence。

### 4.2 后端与明确评分

新增 `--short-context-max-length 300`，定义短片段路由阈值；这只是计算配置。
查询长度不超过该值且已经有锚点限定时，使用现有 Biopython `PairwiseAligner` 适配器。
长区间继续使用现有 minimap2/LASTZ 适配器；读取官方预设和本地帮助后传参。不得将适合长片段的最低累计分数直接套到几十 bp 查询。[minimap2 参数说明](https://lh3.github.io/minimap2/minimap2.html)

短 DNA 命名评分 `nt_blastn_v1`：匹配 +2，错配 -3，gap open -7，gap extend -2，长度 n 的 gap 得分为 `-7-2*(n-1)`。
显式设置所有评分，避免 Biopython 版本默认值变化。使用支持模糊碱基的矩阵或明确的未知碱基处理；N 不构成已知匹配。[Biopython 配对比对文档](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html)
这些值用于局部序列优化，不提供外显子身份的概率。

两种已命名的比对范围：

| 范围 | 调用方式 | 证据含义 |
|---|---|---|
| 已锚定的整个间隔，双方均完整 | global，内外 gap 评分均显式设置 | 比较两个锚点之间完整间隔的对应与缺口 |
| 预测缺失注释的小片段对较长限定区间 | local | 只声明命中块；未命中部分保留未知 |

当前 `MAX_INTERNAL_DP_CELLS=250000` 暂时保留。超过该规模先利用锚点缩小区间；仍超出时使用声明的外部后端或记录 `bounded_alignment_unavailable`，不把大区间无声截断。
数组/短序列范围扩大应有实测时间和内存记录后再改参数。服务器资源充分允许合理并行，无需因此取消算法的规模约束。

### 4.3 替代比对与短序列接受条件

读取最优比对迭代器，保留不同目标位置或切点的最优解。默认最多物化 64 个不同块结构；到达上限必须标记 `enumeration_complete=false`，不能宣称唯一对应。
同一块结构的不同 CIGAR 序列可归并，但影响切点的 gap 位置必须保留。
该库输出的是最优解；“已枚举全部最优解”不代表已经搜索全部近优解。后端局部命中、其他已注释区间和外部比对提供的替代候选一并保留。

短 DNA `resolved_local` 需要：双侧基因内锚点、已知查询覆盖默认至少 0.80、非未知配对至少 12 nt、identity 默认至少 0.70、实际块在允许候选中无竞争位置。
这组门槛仅为首轮可重复的判定规则，运行前冻结；对结构更短或更分化的元件输出候选，仍可独立报告可定位的剪接切点。
不得以短片段获得正分直接证明同源。3-11 nt 微外显子仅凭自身序列通常缺少分辨力，两侧锚点及结构对应需要独立给出。[微外显子综述](https://pmc.ncbi.nlm.nih.gov/articles/PMC5863539/)

### 4.4 序列缺失必须单独判定

`sequence_absent` 仅在以下条件共同成立时赋值：

1. 同一同源区间双侧锚点位置已确定。
2. 目标组装在两锚点及中间区域连续，没有未确定碱基导致的信息缺口。
3. 序列比对明确支持该查询单元相对于两侧的删除缺口，而非仅没有局部命中。
4. 主要可接受替代比对对这一删除结论一致。

目标中存在无法对应的替代序列、锚点间无法比对、contig 终止、重复多命中均赋未知并保留原因。
“该同源剪接位置没有内含子”可由连续跨越切点的 CDS 支持，即使祖先内含子的 DNA 已经无法比较。该记录进入 junction 层；它不自动产生对整个未知非编码序列的 presence 判定。

## 5. F2/F3：局部对应、顺序和重复

### 5.1 先确定候选，再建立有序对应

`preprocess.match_evidence` 必须原样传递候选集合的关键信息。现有 `hit_count/alternative_hits/mapping_quality` 的中途丢失在这里修复。
候选去重键为来源对象、实际块、方向、后端和评分方案；两条转录本指向相同基因组片段时不能当作两份独立序列支持。

对一个基因对、一个锚点间隔构建有向无环图：

1. 节点为候选匹配块；只比较同一分子类型、同一评分方案。
2. 若两个块在 query/target 转录方向均按顺序且不重叠，则允许连接。
3. 固定锚点约束路径起止。对同一原始区间的重叠候选，先按投影端点细分候选覆盖，再禁止重复累计同一配对碱基。
4. 边界角色、内含子是否存在不加入匹配奖励，避免优先选择“更保守”的结构。
5. 顺序不相容的候选保留为 `noncollinear_candidate`；本轮不赋具体重排机制。

同一评分方案内路径得分为所覆盖非重叠块的序列分数之和。跨块 gap 延伸处理沿用明确后端语义；若无法恢复后端的跨块得分，仅把该指标命名为 `sum_of_block_scores`，不能与完整比对 score 混比。
DNA/蛋白分数不相加；二者给出相容约束或冲突。原有经验混合权重最多保留作诊断列，不再单独决定进入结构矩阵的资格。

### 5.2 对应链算法与歧义

令候选块分数为 w(v)，则：

```text
F(v) = w(v) + max(0, max F(u), u 可在 v 前面)
B(v) = w(v) + max(0, max B(z), z 可在 v 后面)
S*   = 最佳可行完整路径分数
经过 v 的最好路径分数 = F(v) + B(v) - w(v)
```

固定起止锚点时，未到达锚点的路径得分为不可行，不能用上述 0 重新开始。实现应分别区分局部模式和固定锚点模式。
第一版用排序加显式前驱，最坏 `O(C^2)`，C 为该局部间隔候选块数。先限定区间与去重，重复簇巨大时不枚举全部路径。
若真实案例显示这里耗时显著，再用坐标压缩与区间最大值索引优化单调二维链；不提前引入复杂数据结构。

保留所有满足 `best_through >= S* - delta` 的候选及其兼容关系；默认 `delta=max(2*match_reward, 0.05*abs(S*))` 仅用于**同一核酸评分方案**。
蛋白对应采用同一替换矩阵的命名分差策略，第一版只将最优并列及已报告竞争重复位置列为歧义，不借用上述核酸数值。
分差是敏感性阈值，无概率含义。有限候选集合中的唯一最优必须表述为“在已报告候选中唯一”。
`enumeration_complete=false` 且被截断范围可能影响该位点时，输出 `candidate_search_incomplete`；其他独立锚定区间仍可继续。

不需要输出指数规模的全部路径。保存候选子图及每个目标位点可能的映射集合即可。
链中的局部独立歧义采用区间级 unknown；不要因一个重复小块使整个基因所有元件失效。

### 5.3 分裂、融合和重复的判定

| 观测关系 | 计算要求 | 允许的结构描述 |
|---|---|---|
| 一个外显子对应多个后代外显子 | 多个后代片段在共同序列轴上互补、顺序一致；相关剪接切点可定位 | 一对多结构对应；树上再判断 split/fusion 方向 |
| 多个片段覆盖共同轴的相同位置 | 保留各 repeat_instance；候选竞争区间相同 | 内部重复对应或无法分辨的重复实例 |
| 仅部分片段覆盖 | 未覆盖部分不补齐；各切点检查是否实际跨越 | 部分一对多对应，未覆盖事件未知 |
| 同一物种不同转录本拆分方式不同 | 保留每条路径；基因座角色与路径状态分别计算 | 已注释替代路径；不直接写为物种间 gain/loss |

不限制 1→2。对于 1→n，按所有实际可定位切点生成 n-1 个 junction 观测，保留同一个 `linked_group_id`。
允许少量边界移动造成的不同覆盖，精确记录插入/删除块；不能用一个固定百分比强迫覆盖严格互补。
“n-1 个切点改变”与“发生了 n-1 次独立分子事件”分别报告，后者当前不推断。
文献有 1→3、3→1 等实际结构情形，可作为定义依据；这些情形不要求立即采用高维联合历史模型。[真菌内含子演化研究](https://academic.oup.com/mbe/article/38/10/4166/6192812)

### 5.4 局部角色冲突

将 predicted CDS 的块与 membership 的 matched_blocks 求交；只有相交部分产生冲突。
一个 40 kb 内含子其他位置存在预测 CDS，不会改变其远处已匹配 UTR 的角色。
对于 dsx EG0029，目标同源区间约为 `2476402-2476701`，已有 44 bp 预测约为 `2476398-2476441`，确实相交：DNA presence 保持 1，受影响的角色保留待定。
边界冲突可把一条 membership 分为受影响和未受影响的局部证据，但不把这些子段自动解释成新外显子。

## 6. F5：注释与结构状态的明确语义

### 6.1 注释补全及搜索范围

保留三个范围：原注释 gene bounds、该基因真实关联子 feature 的覆盖并集、允许的 genome search bounds。
原 gene 行 `1-2000`，关联 exon 延伸至 `3000` 时，提取时应接纳明确 parent 关系的 feature；记录 gene 行与子 feature 范围不一致。
末端 CDS 对应延伸到原 gene 行以外时，在声明的搜索范围内投影，生成补注释候选。
默认继续保留 `--flank 1000 --max-extension 10000`，并报告左右是否触及搜索上限；触及上限的末端未发现不得赋缺失。
邻近其他基因的注释仅作所属关系记录；跨入其区域不自动重新指定基因归属。

miniprot/CESAR 类编码投影可提示遗漏或边界冲突。已知参考结构的保守倾向可能影响投影，原注释、无角色偏好的序列对应和预测应独立保存。
CESAR 2.0 的边界变化及多外显子相连实例适合本轮比较；其默认物种参数不能无说明移植到真菌。[CESAR 2.0 论文](https://publications.mpi-cbg.de/Sharma_2017_6922.pdf)

### 6.2 转录本级角色

给每条路径、每个匹配子区间记录：

```text
path_role = exonic | intronic | uncovered | conflicting
coding_role = CDS | five_prime_UTR | three_prime_UTR | noncoding_exon | mixed | unknown
position_role = first | internal | last | single | unknown
annotation_source, original_attributes, partial_start, partial_end
```

UTR 长度变化、可变末端外显子和 CDS 分裂分别输出。DNA 和已有注释可支持“给定注释路径中的结构差异”；实际 poly(A) 使用、性别特异剪接及表达调控需要外部资料解释。[可变末端加工综述](https://pmc.ncbi.nlm.nih.gov/articles/PMC8795681/)
内含子中可存在其他已注释序列类型，数据结构允许叠加；软件不根据一个 motif 自动赋予调控功能。[内含子生物学综述](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2023.1150212/full)

### 6.3 基因座级状态表

| 实际证据 | DNA presence | exon role | junction |
|---|---:|---:|---|
| 同源 DNA 与至少一条已有路径外显子相交，实际覆盖充分 | 1 | 1，标记 repertoire-any | 对具体已注释连接另行计算 |
| 同源 DNA 存在，在所给全部覆盖该区间的完整路径中为内含子 | 1 | 0，标记 annotation-conditional | 可定位切点另行计算 |
| 同源 DNA 存在，缺少覆盖该处的可用路径 | 1 | ? | ? |
| 仅有新的编码投影、无原注释外显子 | 1，前提是序列对应成立 | ?，预测另列 | 预测连接只作为候选 |
| 已确认相对于锚点的序列删除 | 0 | ?，inapplicable | 依赖该 DNA 的连接为 inapplicable |
| 正确对齐的连续 CDS 跨过待检内含子切点 | 按序列对象另定 | 按已有路径另定 | 0，前提是覆盖切点两侧 |
| 对应位置有明确已注释内含子 | 按序列对象另定 | 按已有路径另定 | 1 |
| 原注释和局部预测冲突 | 可仍为 1 | ?，conflicting_annotation | 受影响连接为 ? |

`annotation-conditional=0` 表示提供的注释路径不将该区间作为外显子；它不能保证所有组织和条件下都无外显子使用。
缺失转录本、局部模型不足与真正非外显子仍可在 DNA+注释输入中不可分；结果必须保留这个条件。
同一基因座任一提供路径为 exonic 时，gene-repertoire 的 exon role 为 1；其他路径 intronic 记录为 alternative usage，不自动构成矛盾。

### 6.4 junction 识别

1. CDS 接头按共同 MSA 列加密码子内切点定位，phase 缺失时明确记录，不能补成 0。
2. 非编码接头按已确定核酸比对列定位。
3. 两处实际不同的切点即使相距 1-3 nt 也保留；同源位置只能依据列映射等价合并。
4. 替代 gap 放置产生切点范围时保存集合；不同候选不一致时该切点未知。
5. 没有该内含子须有同一注释路径连续跨越位置；仅仅“该物种只有一个外显子”不够。
6. canonical/noncanonical motif 保存为描述；现有 GT-AG、GC-AG、AT-AC 支持保留。不把未列 motif 一律写成不存在。

附近剪接位点变化与注释差异在真实研究中经常纠缠；修订采用可追踪坐标和独立注释记录，不照搬其他数据集的重注释代价。[ReSplicer 方法论文](https://academic.oup.com/gbe/article/8/8/2340/2198117)

## 7. 系统发育：本版具体建模范围

### 7.1 三个现有层保留，输出的生物学名称校正

| 内部层名 | 建议用户名称 | 0→1 / 1→0 的结构含义 |
|---|---|---|
| `exon_presence` | 同源序列单元存在 | sequence gain / sequence loss |
| `exon_role` | 给定注释集合中的外显子身份 | exon-role gain / exon-role loss |
| `splice_junction` | 同源剪接位置的内含子 | intron gain / intron loss；在跨同源外显子内部时附 split/fusion 关系 |

兼容已有层字段；修正 `parsimony._event_type` 将 DNA presence 直接命名为 exon gain/loss 的歧义。
目前不增加“转座子介导”“基因转换”“选择压力”“新功能化”等自动类别。
顺序/方向/已注释 UTR/ncRNA 的对应结果可以输出；只有已有明确同源位置与状态定义的对象进入上述正式模型。其余作为描述性结构结果，不丢弃原数据。

### 7.2 默认最大简约推断

对位点 j，叶节点状态集合 A_i 为 `{0}`、`{1}` 或 `{0,1}`。代价 `c(a,b)=0` 当相同，其他为 1。

```text
D_v(s) = sum over children u [ min_t (D_u(t) + c(s,t)) ]
D_tip(s) = 0 if s in A_tip else +infinity
minimum_changes = min_s D_root(s)
```

先后序计算 D，再用树外侧消息获取每个节点、每条分支在**全部最简历史中**允许的状态/端点状态对。
使用 `prefix/suffix` 子消息和计算外侧值，避免 `infinity - infinity`。保留多分叉，未知 tip 不任意填成 0。
令分支可行状态对为 E：若所有 E 中状态均 0→1，输出 `required_gain`；若部分为 0→1，输出 `possible_gain`。丢失同理。
定位支持为集合关系，不能用“最简历史有几条”直接制造概率。
全未知位点保留在观测表，推断标记 `no_observed_states`；全部已知状态一致时报告观测保守范围和未观测物种数量。

分裂/融合是对定位的结构转换作出的名称：一个外显子内部 n-1 个已定位接头改变，可形成 1↔n 描述。
同一 junction 层的多个接头在同一分支都 required 时可写 `supported_compound_pattern`；仅分别 possible 时输出候选组合，不能断言共同发生，更不能相乘为联合概率。跨层结果继续分别陈述，不使用该标签拼接联合历史。
树上的变化放置在一个分支区间内；没有时间标定不报告绝对发生时间。

### 7.3 跨层依赖的处理决定

本轮继续独立的三层分析，避免在修复前端同时改变全部历史模型。
每层输出条件明确；不得将三个边际最优状态拼成一条确定祖先转录本。
presence 为 0 时角色属于 inapplicable，模型入口不能将它当普通 role=0。
部分 tip 的 DNA 消失而其他 tip 有角色差异时，角色推断只表述为仍可观察基因座上的条件历史，不提供联合存在-角色概率。

若未来确需联合存在与角色，可另立三状态模型 `{DNA absent, DNA present/nonexonic, DNA present/exonic}`，为 DNA 存在但角色未知设置后两状态的允许集合。其转换代价、采样条件和与当前输出的区别需单独设计。
**该联合模型列为后续设计项，本轮执行者不得自行接入默认推断。**

### 7.4 可选 CTMC 保留的模型

继续使用固定树上的二状态连续时间马尔可夫模型：

```text
Q = [[-g, g], [l, -l]]
P(t) = exp(Q*t)
L_tip = (1,0), (0,1), or (1,1)
L_v(s) = product_u sum_t P_st(branch_u) * L_u(t)
L_site = sum_s pi_root(s) * L_root(s)
```

实现保持现有 log-space 剪枝、模式压缩、多初值优化与 profile 计算。无需重写成熟且已核对的核心公式。
g/l 的单位由输入枝长单位决定；单位枝长模型没有“每百万年”的意义。
同一层内多个结构位置贡献似然；物种提供各位置的系统发育观测，碱基数、软件运行次数、线程数都不充当结构位点样本量。

| 菜单 | H0 | H1 | 报告限制 |
|---|---|---|---|
| `parsimony` | 无假设检验 | 全部最简状态历史 | 无 P 值 |
| `er-ard` | g=l | g、l 独立 | 比较获得/丢失转换率参数化；不判断某一次分裂是否发生 |
| `foreground` | 预定前景 m=1 | 前景 g、l 同乘自由 m | 比较前景总体转换速率；前景必须分析前指定 |

根频率选项必须在嵌套模型间相同。`estimated` 为独立根参数；`stationary` 由 g/l 给定；`fixed` 用用户值。
维数差由实际自由参数计算；当前 ER/ARD 常规差 1，前景乘数常规差 1，不把这一数字硬套所有将来模型。
`LR=2*(logL1-logL0)`，数值上 H1 不应明显差于 H0；允许浮点容差，超出时报告拟合问题。
输入树枝长和转移率不同时任意自由缩放估计，以免引入新的不可识别性。

这种“固定数据、固定树、明确自由参数与嵌套关系”的组织方式与 PAML 的模型比较原则一致。PAML 也明确提醒额外分支参数可造成不可识别。[PAML BASEML 指南](https://github.com/abacus-gene/paml/wiki/BASEML)
菜单式入口只负责明确问题、模型及必要输出；不自动寻找解释最漂亮的模型。[HyPhy 方法入口](https://www.hyphy.org/methods/)

### 7.5 统计输出及不可识别情形

每个 family×layer 至少报告：总位点数、纳入数、可变位点数、已知 tip 数分布、相关位点组数、参数个数、根设置、枝长设置、发现规则、拟合状态。
分别保留 `fit_status`、`posterior_available`、`lrt_available`；不能用一个布尔值代替三者。
无观测对比、拟合失败、边界解、开放 profile 区间、病态信息矩阵均写具体原因。
profile 上限到达搜索范围与有限置信上限分别报告。开放区间时，祖先概率的参数敏感性范围写 unavailable；不能仅用 MLE 一点生成零宽范围。

常规卡方近似只有在正则条件下提供渐近解释。几个结构位置的单基因数据通常无法仅凭 Hessian 判断近似误差。
本轮不做模拟校准。可选检验输出增加 `reference_distribution=asymptotic_chi_square` 和 `small_sample_accuracy=unassessed`；资料不足时 P=NA。
不要自行设置“超过若干位点即可靠”的新阈值。采用已声明的数值可识别条件，诚实保留小样本适用范围。
多重比较 BH 仅作用于预先定义同类问题的有效 P 值，报告纳入检验数量；BH 不改善模型设定或小样本参考分布。

祖先状态概率与分支端点概率为给定模型、估计参数和观测条件下的概率；当前未对全部参数或比对不确定性积分，标记 `conditional_at_estimates`。
期望转换次数可大于分支两端不同的概率，两者字段不能合并。

### 7.6 发现规则与相关结构位点

在状态构造后、模型前统一处理：

| 选项 | 选择规则 | 条件似然 |
|---|---|---|
| `observed-at-least-one` | 已知 tip 中至少一个 1 | 按该位点的已知 tip mask 除以被发现概率 |
| `variable-only` | 已知 tip 同时有 0 和 1 | 条件于在已观察 tip 中可变 |
| `complete-universe` | 用户给出独立定义的完整元件目录，保留 all-zero | 不按出现 1 截断 |

发现概率与相同 root/branches/parameters 计算。记录被排除对象及原因，不能先按另一规则丢弃再切换选项。
真实注释缺失可能依赖元件状态；当前 mask 条件化没有建立完整注释漏检概率模型，需明确该限制。
同一删除可能影响多个相邻位置；相同状态模式可以压缩计算，但仍须保留位点数和 `linked_group_id`。
对同一物理切点、同一对应单元的重复记录去重；不同物理位点即使模式相同也不能随意当一个生物学观测。
不同层似然不相加作一个总检验，事件表不简单相加作独立事件总数。

## 8. 完整执行步骤与 subagent 分工

### 8.1 工作方式

下一位总工程师负责决定范围、冻结共享接口、整合补丁和解释结果。执行者一次领取边界明确的任务，不必由总工程师实时跟随每条命令。
每个任务结束提交：修改文件、接口变化、已解决问题、未解决问题、允许命令的执行结果。等待提交和阶段汇总即可。
只使用当前客户端实际提供的 agent/model 选项；GPT-5.5/5.6 是否可选由运行环境决定，不能在任务报告中虚构模型。
同一工作区避免两个执行者同时修改同一文件。共享文件的最终修改由指定集成人员负责。
执行者只运行总工程师明确列出的脚本/工作流，不获得无限制下载、安装、重试或新建工作流权限。

### 8.2 可直接派发的任务单

| 任务 | 写入范围 | 依赖 | 必须交付 |
|---|---|---|---|
| T0 契约集成 | 公共字段、坐标适配、`io.py` 必需部分、CLI 参数草案 | 无 | 第 2 节契约、唯一矩阵 writer，旧字段兼容；由一个执行者完成 |
| T1 短比对 | `alignment.py`、`test_protein_alignment.py` 中适配器相关部分或指定测试类 | T0 | F1、同分候选、明确评分和规模状态；不动结构状态规则 |
| T2 编码坐标 | `coding_correspondence.py`、`test_coding_correspondence.py` | T0；调用 T1 冻结接口 | 家族 MSA、codon/genome 映射、局部锚点规则 |
| T3 对应与链 | `preprocess.py` 对应部分、`correspondence.py`、指定 evidence 测试类 | T0；接入 T1/T2 | 候选无损传递、局部 membership、有序链、重复/互补区分 |
| T4 注释观测 | `annotation.py`、`structural_sites.py`、`test_observations.py` | T3 | 转录本角色、局部冲突、剪接状态、末端搜索状态 |
| T5 系统发育 | `parsimony.py`、`structural_phylogeny.py`、`test_statistics.py` | T0；最终需 T4 矩阵 | 共享矩阵、结果语义、条件范围；保留概率内核 |
| T6 真实案例资料 | 现有 `docs/real_positive_cases.md`、`docs/data_sources.md` 的指定段落 | 无 | FDPS/MAMSTR/PKM 的原始文献证据和资源清单；不默认执行下载 |
| T7 命令/结果集成 | `case.py`、`cli.py`、既有可视化模块、相关文档 | T1-T5 | 用户场景全链、图、版本迁移 |
| T8 服务器执行 | 既有工作流及其已批准参数/协议 | 集成冻结后 | 从 FASTA/GFF 重跑，正式日志和病例解释 |

T1 与 T2 若需同改 `alignment.py`，T2 只提交接口需求，由 T1 写入。T3 与 T4 对 `preprocess.py` 提取层的需求同理由 T3 集成。
T5 可先处理 writer/输出字段及数值边界，最终验收等待观测矩阵冻结。
测试文件由任务表指定所有者；同一文件包含其他人的未提交改动时不覆盖。

### 8.3 派发模板

```text
任务：Tn，目标问题 Fx。
基线：当前提交及相关未提交变化。
允许读取：本规格指定模块及相应文献/官方帮助。
允许修改：逐项列出的文件和函数。
禁止修改：其他代理拥有的文件；服务器标准；旧多拷贝算法。
必须遵守：坐标、候选集合、观测和统计契约。
允许执行：总工程师明确列出的脚本/工作流；未列出时仅阅读和编辑。
提交：改动摘要、接口、明确的未完成条件、执行状态；不自行发布。
失败：记录命令及首个可行动错误，停止该执行步骤，交总工程师安排具体修复。
```

### 8.4 集成顺序

1. 读取本规格与当前工作区改动；确认服务器标准和上游单拷贝范围。
2. 完成 T0，冻结坐标和公共字段；在代码中保留简短语义注释。
3. 并行 T1/T2/T6；T5 可做不依赖前端的输出修订。
4. T3 接入候选及编码坐标；先完成局部对应，再接入 chain。
5. T4 根据真实块构建观测；逐项核对 F1-F5 的因果链，不提前重跑旧矩阵演示。
6. T5 消费同一矩阵；T7 完成 CLI 和可视化。
7. 冻结算法/参数，登记预先声明的真实案例问题和预期可评估范围。
8. T8 通过现有工作流一次开展选定回归及真实分析；失败由总工程师定位后派发有范围的修复，再明确重跑受影响部分。
9. 按病例汇总改善和剩余失败原因；更新方法说明、当前限制和真实基准文档。
10. 完成代码和文档后再决定版本及 GitHub 发布；本规划轮不修改版本号或远端。

## 9. 用户命令、参数与服务器执行

### 9.1 新参数保持克制

本轮新增参数限制为：

```text
--coding-msa-mode linsi|einsi
--short-context-max-length 300
--annotation-view repertoire|canonical
```

`--annotation-view repertoire` 默认与当前 `--transcript-policy all` 配合；canonical 只在用户明确选择时使用，并在结果中记录选取规则。
现有 `--transcript-policy` 控制提取，`annotation-view` 控制汇总，不允许请求已被提取阶段舍弃的路径。
评分、锚点窗口、候选分差和候选物化上限先放入一个命名配置对象并在结果中写出；没有需求时不为每个内部常量增加 CLI 开关。
用户显式指定后端失败应报告；自动模式选择了哪个后端、为何选择应可读，不进行隐蔽 fallback。

### 9.2 目标用户操作

以下是**实现后**需要在现有服务器工作流中实际完成的命令形态。命令参数值由既有案例表提供，不要求用户编辑 Python。

```bash
insiphy import-orthofinder \
  --orthofinder-dir ORTHOFINDER_RESULTS \
  --orthogroup ORTHOGROUP \
  --genome-manifest GENOME_MANIFEST \
  --species-tree ROOTED_TREE \
  --output-dir analysis/upstream

insiphy build-case \
  --manifest analysis/upstream/manifest.tsv \
  --species-tree analysis/upstream/species_tree.tsv \
  --output-dir analysis/case \
  --transcript-policy all --coding-msa-mode linsi \
  --aligner mafft --context-aligner minimap2 --threads 16 \
  --flank 1000 --max-extension 10000

insiphy run --input-dir analysis/case \
  --output-dir analysis/results_parsimony \
  --analysis-scope single-copy --model parsimony \
  --annotation-view repertoire --evidence-aligner miniprot --threads 16

insiphy visualize --input-dir analysis/case \
  --result-dir analysis/results_parsimony \
  --output-dir analysis/figures --correspondence-encoding color
```

其它有明确既定单拷贝关系的案例可以直接提供 manifest，不能伪称已经全部经过 OrthoFinder。
可选 CTMC 由 `infer-phylogeny --model er-ard|foreground` 使用同一冻结观测矩阵。实现时增加内部 matrix 读入口即可；若新增 CLI 路径参数，明确命名并写入帮助。
当前工作流为不同结果目录复制对应表，下一轮应明确传递冻结矩阵及其观察设定；不同模型不得重新推断同源关系。

### 9.3 本服务器路径与约束

```text
源码：/data/projects/intragenic_structure/software/intragenic_structure_v1
现有工作流：/data/projects/intragenic_structure/scripts/20260916_145605_insiphy
案例表：/data/projects/intragenic_structure/inputs/insiphy_real_cases/structural_cases.tsv
最近完整前端结果：/data/projects/intragenic_structure/results/20260918_111929_insiphy
最近再推断结果：/data/projects/intragenic_structure/results/20260918_121100_insiphy
现有依赖容器：/data/containers/insiphy/0.12.4/image.sif
```

容器版本不代表挂载源码版本；执行记录需分别写明。复用前先读取既有 `container.yaml` 和调用方式，不重复安装系统软件。
现有 `params.yaml` 的 `parameters.stage=phylogeny` 只做再推断。本轮前端改变必须切回 `main.nf` 的完整 `RUN_INSIPHY_CASE` 分支，不能继续复用旧 prepared results 作为新前端成果。
按照服务器标准更新该工作流的 `params.yaml/protocol.md` 及运行目录；正式执行入口仅为该目录的 `bash run.sh`。
软件包依赖保持 Python 库及外部比对器；服务器的 Nextflow/Slurm 属于现有执行设施，留在包外。
执行前阅读 `STANDARDS_INDEX.md` 路由的文件；下载/新容器分别增加相应模块规范。规范缺失即停止该操作并询问用户。
命令默认直连，清除代理环境。下载代理需针对具体文件重新获批；不将历史 git 代理许可扩展到数据或容器。
本计划授权范围不包含额外完整性扫描、摘要计算、系统安装或新建通用执行框架。

## 10. 真实案例：逐例问题和接受结果

文献结论只交给评价和解释步骤，不输入同源判断、状态生成或树上方向选择。
每例先冻结输入 accession/version、基因 ID、转录本来源、树、焦点位置及文献证据等级，再开始运行。

| 案例 | 本轮提问 | 可以接受的结果 | 不应强制的结果 |
|---|---|---|---|
| spo5 / `spbc29a10_02` | 55 bp 区间及两侧位置是否能够对应；编码/内含子角色能否区分 | 比对块、角色、切点和受限原因明确；修复短后端完全无搜索机会的问题 | 必须得到某方向的显著事件 |
| rec8 / `spog_00055` | 39/43 bp 区间及负链切点的同源对应 | 正确方向/phase；DNA 对应与 junction 分开报告 | 找不到内含子 DNA 即判角色变化不存在 |
| Hdac3 | Dana 62 bp 内含子与其他物种连续编码结构 | 找到确切切点；保留等代价获得/丢失定位；需要时加入有资料的外群 | 在现有树下强制唯一获得 |
| RpL32 | 已覆盖单拷贝结构是否保守，未知是否有具体原因 | 每层覆盖率及已观察保守结构；无差异可以无可用 LRT | 以一个保守基因估计总体误报率 |
| dsx 七蜂 | 共同编码骨架、UTR/非编码对应、多路径及注释不完整 | 分开编码/UTR/非编码、已知与候选；局部冲突可追踪 | 所有蜂的全部末端外显子完全相同 |
| FDPS，待准备 | 文献 3 个外显子与 1 个连续外显子关系 | 两个切点及互补覆盖，检验 3↔1 支持 | 从外显子数量直接判定同源 |
| MAMSTR，待准备 | 已报道约 30 bp acceptor 变化能否被保留 | 两个精确切点/对应范围，排除旧宽容差合并 | 把相近位置自动称为同一位置 |
| PKM，待准备 | 互斥外显子 9/10 是否保留重复实例和路径 | 同源候选及各路径，歧义范围可读 | 将互斥路径直接标成物种间获得/丢失 |

FDPS/MAMSTR 由 CESAR 2.0 实例出发，重新核对现行组装与单拷贝资格；PKM 作为重复/路径压力案例，其资格和注释需逐物种确认。[CESAR 2.0 实例](https://publications.mpi-cbg.de/Sharma_2017_6922.pdf)、[PKM 转录本研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC7835739/)

### 10.1 当前坐标定位信息

这些坐标帮助定位旧失败，不作为算法强制命中清单：

```text
spo5: S. pombe NC_003423.3:2536843-2536897:+，55 bp
rec8: Cryo NW_013185626.1:118741-118779:-，39 bp coding
      Octo NW_013185621.1:171120-171162:-，43 bp intronic
dsx:  EG0013，Acer UTR 与 Amel 长内含子中的局部 DNA 对应
      EG0055，Bpas 非编码外显子与 Bter 局部内含子对应
      EG0029，Bign CDS/UTR 混合外显子与 Bter 局部对应及 44 bp 预测冲突
```

EG 编号可能随新对应改变；通过来源基因、feature 和坐标追踪，不能要求新算法维持偶然的旧编号。
spo5/rec8 的结构定义参考原始跨物种比较，重新定位时必须核对组装版本。[Zhu 与 Niu 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061683)
dsx 编码核心保守与某些蜂类末端/雌性相关路径有差异可以同时成立；文献 cDNA 用作独立背景证据，分析输入继续限定基因组与注释。[蜂类 dsx 原始研究](https://link.springer.com/article/10.1007/s13592-020-00735-8)

### 10.2 真实案例评价分开三项

1. **结构可评估覆盖**：文献定位点中，多少在当前输入和同源比对中可评估；未知按原因计数。
2. **结构判定**：可评估点的角色/剪接切点/互补对应是否与独立资料一致，给逐点表。
3. **历史定位**：在给定树和外群下，支持哪些分支、有哪些等价方向；结构对比正确但定位不唯一应单列。

只有具备明确已知正例和可定义负例范围时计算 precision/recall，并同时给分子与分母。
未全面审查的位点不能当作负例。几篇文献的事件目录不足以估计全基因组误报率。
文献候选、实验支持、现代注释支持分别标注。新论文的 candidate 清单不自动升为真实阳性。

### 10.3 后续资料执行者的具体任务

优先准备 FDPS 与 MAMSTR，随后 PKM。每例交付现有 docs 表中的：论文图/补充表位置、物种、组装与注释来源、基因别名、原始坐标、单拷贝证据、树/外群、预期可评估问题。
只在总工程师指定的数据下载任务中获得下载权限；沿用服务器现有下载规范与目录。
真菌公开研究的对齐数据可用于参考病例选择，但主要分析仍重新从 FASTA/GFF 生成，不把已知对齐真值传入被测方法。[真菌研究代码与资料](https://github.com/Brookesloci/fungi_intron_paper_2020/)
新的大规模内含子获得候选集可辅助找病例，需逐个核对，不能假定其所有候选已经证实，也不导入其整基因最佳命中规则替代 OrthoFinder 输入。[2025 年研究](https://academic.oup.com/gbe/article/17/6/evaf091/8134177)、[原始代码](https://github.com/celinehohzm/intron_gain)

## 11. 明确授权范围内的验收

本节描述下一轮实施任务的验收内容。本规划轮不执行这些命令。执行时将准确测试类/脚本列入已有工作流，避免默认扫描全部历史 toy 测试。

### 11.1 确定性回归

| 位置 | 必须覆盖 |
|---|---|
| 比对适配 | 短查询实际进入正确后端；未知碱基；同分 gap；候选截断状态；无命中与失败分别记录 |
| 编码映射 | 正负链、跨外显子密码子、部分 CDS、CDS/UTR 混合、重复蛋白别名 |
| 对应区间 | 互补 1→3、重复覆盖、仅覆盖半个父外显子、真实多命中穿透序列化 |
| 注释角色 | 同一位置不同转录本；DNA 存在角色未知；真实相交预测与父区间远处预测分开 |
| junction | 精确切点、相位、连续跨越才赋 0；近邻不同边界不合并 |
| 树及简约 | 多分叉、未知 tip、全未知、全部最简分支集合，不将可能事件写为必需 |
| CTMC | 既有枚举/闭式公式回归、零枝长、开放 profile、选择规则一致、非法拟合无后验 |
| 接口 | 原字段兼容、新元信息穿透、不同命令读同一矩阵 |

这些是代码公式和确定语义的回归断言，不作随机模拟实验，也不作为发表级统计误差率估计。
任何新的预期值须由区间/状态定义或独立公式确定，不能先运行现代码再把输出抄成答案。

### 11.2 真实用户功能覆盖

至少覆盖已有五例的原始数据完整流程；新增 FDPS/MAMSTR 以资料就绪为前提。
执行 OrthoFinder 导入、直接 manifest、all transcript、一个明确 canonical 对照、默认 parsimony、ER/ARD、预定 foreground、color 图及一份 pattern 图。
1/16 线程比较输入到观测和分支结论；浮点似然预先规定 `abs <= 1e-8` 或 `rel <= 1e-7` 容差，位点状态/坐标/事件集合需一致。
这不要求所有绘图字节或所有浮点打印完全相同；稳定排序后生物学表必须可比较。
在 spo5/rec8 上做预先声明的 L-INS-i/E-INS-i 敏感性对照，记录对应和切点变化；结果不稳定的局部不能隐藏。
短 DNA identity/coverage 规则仅对指定真实区间作声明好的参数敏感性分析：spo5/rec8 采用 identity `{0.60,0.70,0.80}` 与 coverage `{0.60,0.80}`，其余条件固定。复用已保存候选只重算接受规则，不重复运行同一比对；该对照只评估阈值敏感性，不评估后端搜索不确定性。不选择最能吻合文献的参数作为事后默认。

### 11.3 可视化验收

默认同源外显子/已确认片段用颜色，pattern 为替代模式。两者分别能读图，不同时叠加作为强制编码。
主图左侧有根树，右侧基因横向条带，实际对应区间之间有连续 ribbon/线条。
内含子画连接线，原始外显子画实体，CDS/UTR 用高度或边框区分；颜色表示对应，不能同时编码置信概率。
候选在独立轨道并给标签；主图不得用整段长内含子连线代替其内部小块。
同一物种多转录本分轨，原始基因长度较长时允许明确标注的断轴/压缩内含子，图例说明比例。
采用适合色觉差异的调色方案及稳定标签；重复颜色时仍可通过标签区分。复杂基因使用局部放大图。
每个重点案例实际打开渲染图，检查条带、连接、物种顺序、文字边界和候选图例。记录人工阅读结论；仅生成 SVG 文件不算完成这项验收。

## 12. 资源、复杂度与简洁工程要求

### 12.1 资源分配

设总线程预算 T。独立家族/局部任务并行数 W，每个外部工具线程 t，满足 `W*t <= T`。
默认一个家族蛋白 MSA 使用 T；短配对任务使用 W 个单线程 worker；CTMC 按 family×layer 分配，BLAS 各为 1。
不要每个 Python worker 内再启动 T 线程 MAFFT/minimap2。
有序结果按来源键统一排序后写出；worker 完成顺序不参与编号和图排列。

### 12.2 内存和复杂度目标

| 步骤 | 目标 |
|---|---|
| genome 提取 | 使用既有索引/受限范围读取；每个 worker 不复制整个基因组 |
| family MSA | 去重完全相同蛋白，仅保留一次序列和映射别名 |
| CDS/genome 投影 | 以连续块和紧凑数组保存；不为每个物种对建立大量 Python 碱基二元组 |
| 短 DP | `O(m*n)`，限定锚点区间；沿用已有单次规模约束 |
| 对应链 | 初版 `O(C^2)` 局部候选；仅有证据显示瓶颈再优化索引 |
| 简约 | U 个独特观测模式、N 树节点时约 `O(U*N*K^2)`，当前 K=2 |
| CTMC | 同上每次似然，转移矩阵按枝长/模型缓存；优化调用次数另计 |
| 图 | 根据实际 membership 块生成连接；不绘制所有候选的全连接网络 |

缓存仅限运行内、按明确对象 ID 和配置；不增加永久缓存数据库和哈希命名层。
较大重复簇可保持 raw feature 完整、分块计算可确定的外围。超过当前对应分辨能力时仅该重复区标记 unresolved，报告实例数、候选数、覆盖范围。
运行记录使用现有 Slurm/工作流时间和最大 RSS，禁止声称由五个案例推导出普适加速倍数。

### 12.3 避免过度防御式代码的具体规则

1. 输入边界检查一次：树、坐标、物种/基因关系、文件字段。模块内部消费已明确的类型契约。
2. 预期的生物学不确定性用状态和原因表达；外部程序失败用明确异常。二者不混用。
3. 禁止宽泛 `except Exception: return unknown`，禁止缺字段时用任意默认坐标继续推断。
4. 不为可能的未来后端引入插件总线、服务容器、五层配置工厂或通用工作流框架。
5. 新增函数围绕真实复杂度：坐标转换、候选映射、状态构造、树上计算；同一判据只实现一次。
6. 保留原始证据和有用状态；不每一步重新读所有文件做完整性扫描。
7. 文件写出沿用已有 API；手工代码改动用 `apply_patch`，不覆盖用户未提交变化。
8. 保留现有依赖优先。大改包结构、Python 基线或额外依赖需要列出明确必要性。

## 13. 完成条件与不能过度宣称的部分

### 13.1 工程完成

- F1-F5 均有明确调用路径、输出字段和相应回归。
- 默认一次 `run` 完成分析，模型菜单共用同一观测矩阵。
- 原始 feature、实际对应块、转录本、观测和树上结果能逐级追溯。
- 选定真实用户操作及 1/16 线程分析完成，失败/不可分析的部分有具体原因。
- 图实际查看，README、输入格式、统计模型和当前限制与实现一致。

### 13.2 生物学完成

- 基础场景包括保守外显子、同源切点有无、局部角色差异、1↔n 互补对应、边界变化。
- 已知短片段有合适的搜索机会；能否确认同源由证据决定。
- dsx 的结果能分开编码骨架、UTR/非编码和注释不完整，避免把整个基因概括成结构不保守。
- 每个真实焦点给出现生结构、序列对应、观测编码、树上定位四步说明。
- 不确定历史保留集合；未知注释不伪装缺失；参数不可估计不输出确定的统计支持。

### 13.3 当前方法适用范围

适合：已确定单拷贝关系、近缘物种、至少部分基因内部区域可可靠比对、用户接受基于已有注释集合的角色解释。
对基因其他构成序列，提供已有注释和可对应范围的结构比较；正式历史分析仅限状态定义和同源位置足够明确的对象。
暂有限制：极短无锚点片段、大规模近同重复外显子簇、严重不完整组装、缺少可用注释路径、无法定位的基因末端、复杂局部重排。
这些限制按区间报告，尽量保留同一基因中可分析的普通部分。
整基因多拷贝、完整联合祖先转录本图、实际转录使用、演化机制和功能推断留在本版范围之外。

发表准备度还需要独立真实病例、比较方法的公平任务划分，以及统计适用范围的论证。
完成本轮不能自动宣称所有 P 值具有已知有限样本错误率，也不能将全部未知减少视为准确率提高。

## 14. 总工程师最后应交付给用户的报告

### 14.1 遇到分歧时怎样决策

| 观察到的情况 | 总工程师的动作 |
|---|---|
| 原始 feature 正确、导出矩阵错误 | 修状态构造和来源投影，保留概率内核 |
| 已知短区间根本未进入合适的比对 | 修路由及搜索范围，再检查对应；暂不改树上模型 |
| 结构对比可靠、两个历史等代价 | 保留两种定位；有独立可用外群时扩充真实输入 |
| 多条路径给出不同角色 | 先确定 repertoire 或 canonical 分析对象，保留原路径 |
| 模型可拟合、P 值不显著 | 解释当前模型比较和参数区间，不为获得显著性改默认 |
| 更低阈值减少未知但新增矛盾 | 按独立证据比较覆盖和错误；保留候选，不能仅凭未知减少接受更低阈值 |
| 一个复杂重复区无法分辨 | 局部报告未解决，继续分析该基因其他确定区域 |
| 现有代码已经正确实现某项建议 | 保留该实现，补缺失接口或必要说明，避免重复重构 |

### 14.2 交付格式

报告以结果为主，逐例提供下表，不只给“demo 跑通”：

```text
案例/输入版本/物种和基因数
可评估元件数、切点数、未知数及原因
文献焦点 -> 当前对应坐标 -> 观测状态 -> 树上可能分支
相对上一版本改变了什么、直接涉及哪项修复
统计模型、参数/NA 原因、P 值具体比较的假设
颜色主图和需要的局部放大图
未解决的问题、是否限制本基因的分析范围
```

项目汇总再给 F1-F5 与已有修复的完成清单、实际执行命令范围、1/16 线程及资源记录、下一版仍未解决事项。
用户不需要从内部编号猜生物学含义；报告使用外显子、内含子位置、编码区、UTR 和注释路径的名称。

## 15. 阅读与来源交接

以下阅读任务服务于具体实现，不要求执行者遍历所有参考仓库：

| 主题 | 必读内容 | 应提取的实现信息 |
|---|---|---|
| 蛋白共同坐标与多路径 | 上述跨物种剪接图论文；原始 `subexons.py/alignment.py` | 边界细分、来源追踪、同列对应与重复限制；不照抄其全局阈值 |
| 精确内含子位置 | [GenePainter 官方帮助](https://genepainter.motorprotein.de/help_commandline) | MSA 列及 phase 的表示，避免比例坐标 |
| 短元件 | 微外显子综述与 Biopython 官方接口 | 短序列的分辨限制、无种子 DP、替代 gap 位置 |
| 边界与注释错误 | ReSplicer 原文与其位点处理代码 | 序列/注释冲突如何留下证据；研究特定代价不直接移植 |
| 编码投影 | CESAR 2.0 原文与官方说明 | 连续外显子和 splice shift 的真实例子、模型的参考保守倾向 |
| 单拷贝位置历史 | 真菌研究的方法部分与原始资料 | 位置对应、1↔n、生物学事件与统计字符数量的区别 |
| 统计菜单 | PAML 官方指南相关模型/祖先状态段落，HyPhy 方法说明 | 条件似然、根参数、嵌套模型、可识别性和问题驱动入口 |
| 真实案例 | Zhu/Niu、dsx、FDPS/MAMSTR、PKM 原始研究 | 图表级证据、输入版本、局限与单拷贝资格 |

原始实现入口：[跨物种剪接图代码](https://github.com/PhyloSofS-Team/thoraxe/tree/master/thoraxe/subexons)、[ReSplicer](https://github.com/csurosm/ReSplicer)、[CESAR 2.0](https://github.com/hillerlab/CESAR2.0)。
借鉴算法或数据组织时记录实际读到的文件/函数及版本。独立实现保留论文引用；如直接复制任何第三方代码，必须遵守其许可证和署名要求。
本文的数值判据、接口和实施顺序为项目设计建议；参考文献支撑生物学现象及已有方法原则，不为这些新默认参数提供未经测量的准确率保证。
