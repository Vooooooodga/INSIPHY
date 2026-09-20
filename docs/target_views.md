# 单目标、固定背景的组图

## 三张图

每个target目录中输出：

- `structure.svg`：完整灰色注释＋目标叠加，不画ribbon；
- `ribbons.svg`：相同背景＋只属于目标、且有实际匹配块支持的ribbon；
- `phylogeny.svg`：相同结构及目标ribbon，叠加该层该site的既有固定树事件。

一组只表达一个family/layer/site。纯RNA连接以两个切点显示，不把整条内含子DNA虚构成一条同源ribbon。没有实际匹配块时不画带，unknown或位置有歧义时采用目标色虚线边框。

每个family的背景几何只构建一次，并复用相同SVG组。提供的所有物种、注释区段、不同路径、CDS/UTR、反链／重叠及其他标签仍在灰色轨道中。未进入当前模型不等于没有注释、没有结构或未被研究。

ribbon只连接相邻面板物种的实际目标匹配，避免跨越未知物种而掩盖缺口；这不是祖先遗传路径或因果方向。树为固定输入的拓扑布局，不将横向画布长度解释为时间。实心事件标记表示所有等简约历史所需，空心表示可能或给定模型条件的结果。无效拟合不制造概率；不画虚构祖先外显子。

## 目标选择

没有选择参数时，使用完整结构目录的所有目标。很大的目录建议显式选择；没有隐藏的“显著才画”门槛。

```bash
insiphy visualize --input-dir input --result-dir results --output-dir figures \
  --target 'family/exon_presence/EG_0003' \
  --target 'family/exon_role/EG_0004'
```

只给site ID时必须唯一。E3同时存在于presence和role时，须给完整ID。

## targets.tsv

必需列：`target_id, family_id, layer, site_id`（TSV）。可选：`element_id, label, color`。

```text
target_id	family_id	layer	site_id	label
exon3	GENE	exon_presence	EG_3	Exon 3 homologous DNA
exon4	GENE	exon_presence	EG_4	Exon 4 homologous DNA
```

可以重复同一target行指定多个原生区间，此时另外必须提供`species,gene_copy_id,contig,start,end,strand`；start/end为1-based闭区间。建议提供`occurrence_id`，否则没有足够信息将该手动区间接到已存在的ribbon匹配。手动区间不能创造新的推断状态。

同一target_id不能混合层、site、标签和颜色。颜色可提供`#RRGGBB`，同一图集不能重复。自动色按完整目录分配后才执行selector，因此只画某个目标时不会改变其原色。大目录仍需依赖目标名称、图例和线型；不同hex颜色不是无限数量的知觉区分保证。

## 完整背景的含义

它是“完整提供的原生注释”，不是软件恢复了全部真实结构。优先读取`raw_gene_features.tsv`；缺少时回退到准备阶段occurrence/path记录，并在图中注明。纯DNA无转录本归属时使用DNA轨道，不伪造RNA路径。

不改变原坐标，不把trim区压缩掉；不同物种locus可使用各自线性尺度，图上明确标注，不能用盒子宽度直接推断跨物种长度比例。反链按5′→3′显示，坐标仍为原生坐标。

## 输出

`index.html`汇总目标；`visualization_manifest.tsv`记录target、颜色、视图及拟合资格；`target_intervals.tsv`记录实际高亮坐标；每family一份`native_background_*.json`记录几何以便检查。没有哈希要求。

## 旧图

```bash
insiphy visualize --input-dir input --result-dir results --output-dir old_figures \
  --layout legacy-overview --correspondence-encoding pattern
```

`pattern`只影响旧总览，新组图按显式目标颜色和可信性线型表达，不自动继承多目标混合色图。
