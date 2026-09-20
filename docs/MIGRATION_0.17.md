# 0.16 → 0.17：可判定范围与单目标组图

基线为本对话交付的完整0.16.0源码。没有将GitHub中更早的0.15替换进来；没有改远程仓库。

## 有意改变的行为

1. 较短的合法首末外显子／转录本不再因未到达gene envelope而自动变为partial。明确的partial属性仍有效；未标记不意味着所有转录本完整性已经验证。
2. 同一父外显子的两个互补、非重叠子区间，可以保持该父ID并按坐标次序组成1:n或n:1链。重复覆盖仍不能当作互补。
3. DNA物理对应和已注释转录本对应分别记录。pathless DNA候选不再仅因其他候选有路径而退出；它不因此成为已证实exon或剪接路径。
4. 短序列仍有默认12配对碱基要求。仅在既有路由已经获得独立双侧锚点、唯一且完整枚举的有界匹配时，3–11 nt候选可走更严格上下文路线（identity≥max(0.90,原阈值)，coverage≥0.95）。这是未校准的候选资格，不是同源概率；普通孤立9nt匹配仍拒绝。没有新splice predictor。
5. 一个路径中的中间外显子未能对应时，不再将两侧已对应外显子构造成新的剪接连接；输出未解析原因。
6. `visualize`默认改为`target-groups`，每个目标输出三张图；需要原0.16总览时显式使用`--layout legacy-overview`。旧draw_synteny等Python入口仍保留。

## 保留的科学契约

仍为固定树上的单拷贝分析。Sankoff、ER/ARD、前景过程及发现偏差规则不重写；unknown和inapplicable不转为0。低覆盖本身不是错误，阈值也不证明可靠性。正式CTMC仍条件于给定mask、对应和模型，未完成统计校准。

## 新增范围接口

`run`和`infer-phylogeny`增加：

```text
--analysis-range all|high-coverage      默认all
--min-callable-fraction FLOAT           默认0.70
```

- all：完整观测目录参与相应模型原有的可估计性规则；低覆盖对象仍可能输出局部或不可估结论。
- high-coverage：只让“明确0或1数 / 全部树tips数 ≥ 阈值”的字符进入本次拟合。
- 0.50、0.70、0.90覆盖摘要是同一矩阵的计数视图，不是重复拟合。
- 高覆盖过滤不改变原DNA、坐标、注释、比对或路径。
- `structural_site_matrix.tsv`始终为完整的有效观测目录；`analysis_structural_site_matrix.tsv`为本次实际子集。显式mask原值另存到逐观测审阅表。
- `analysis_scope.json`明确本次selection rule，与模型自身ascertainment分别解释。
- 现有schema-v3字段未改名。`ObservationMatrix.rows`是实际子集，`.full_rows`保留完整快照，`.selection`为只读范围信息。

## 新增绘图接口

```bash
insiphy visualize --input-dir INPUT --result-dir RESULT --output-dir FIGURES \
  --target 'FAMILY/exon_presence/EG_0003'

insiphy visualize --input-dir INPUT --result-dir RESULT --output-dir FIGURES \
  --target-manifest targets.tsv
```

每个目标的结构、ribbon、条件性历史三视图使用同色；不同目标不同色且有图例。不同层的同名site须使用完整ID以避免混淆。预测同源位置不等于真实功能。

图形以完整**提供的**原生注释为灰色背景，不声称注释已覆盖真实所有转录本。没有raw注释时明确标记仅有准备阶段坐标。每个locus保持自己的原生线性坐标，未对trim区拼接压缩。选择目标、改变范围或缺少该目标都不改变同family的背景。

新模块：`observations/eligibility.py`、`reporting/target_model.py`、`target_scene.py`、`target_views.py`。不存在运行时下载旧代码、额外哈希或新的外部依赖。

## 测试迁移

原有四个旧绘图测试显式选择legacy-overview以继续覆盖旧契约；另以新增测试覆盖默认组图。一个旧partial预期改为正确的非自动partial预期，属于上述有意bugfix，不是删除回归。

默认不需要重建现有全基因组数据。要应用修复到自己的真实数据，应在新输出目录重新执行相应准备／推断阶段，不能仅把旧图重命名为0.17结果。
