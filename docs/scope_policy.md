# 可判定范围，不是删去不保守序列

## 三种量

`state1_fraction_called`描述有结构的比例；`callable_fraction_panel`描述整个面板中能看清0或1的比例；`callable_fraction_applicable`仅作为角色层辅助分母。高覆盖筛选使用第二种。

10种中2有内含子、8种同源位置明确无内含子：100%可判定，保留。2有、8种未知：20%。role层2适用、8个DNA确实缺失：适用内可达100%，但全树仅20%；两种分母均输出。

## 单一owner

`observations/eligibility.py`消费已由观测层判定的结构矩阵。它不会根据英文说明猜角色、降低比对阈值或重新扫描DNA。明确mask的输入值在audit保留，在推断中有效值为unknown。

输出：

- `structural_observation_eligibility.tsv`：原值、有效值、mask、不适用及原原因；
- `structural_character_eligibility.tsv`：0/1/unknown分母、实际物种、可判定MRCA、根两侧覆盖、是否进入本次分析；
- `scope_coverage_views.tsv`：显示阈值下的计数，不附加模型；
- `analysis_scope.json`：实际选择规则与限制；
- 完整矩阵与实际分析矩阵分别保存。

已知0不因缺少对应序列而被当成未知，前提是原观测模块已获得明确缺失证据。程序无法将仅仅无命中升级为可靠0。一个层可判定，不要求其他层也成功。

## 坐标与统计单位

屏蔽仅发生在状态／字符选择层。A—未知B—C仍保留原始B范围，不能生成RNA A–C连接。原始外显子切为多个绘图块，不会因此创建额外统计字符；原linked_group继续传给数值模块。仍未实现任意相关多字符的完整联合事件模型。

## 条件性解释

70%是可调整的工程默认，不是准确性保证，也不是生物学阈值。近缘tips集中产生的高覆盖不保证深层方向可识别。固定mask、state-dependent漏检和候选发现过程的差别仍存在。原observed-at-least-one、variable-only、complete-universe规则保持独立，未以改名规避ascertainment。

## 边界

没有新增完整祖先转录本、复杂倒位、多拷贝、RNA活动或选择系数模型。提供的注释不是真实所有组织和时间的RNA repertoire。预测剪接信号分不代表外显子概率。
