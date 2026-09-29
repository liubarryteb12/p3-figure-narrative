# p3-figure-narrative — P3 出图/叙事层

> 手册 v1.2.8（`AI_AGENT_EXECUTION_MANUAL.md`）五层拓扑中的 **P3 层**独立仓库。
> 职责（手册 §2.1）：**出图 + 表 + figure_index.json** —— 把 P2 的分析数据
> 渲染成符合 SCI 投稿规范的图。

## 内容

| 文件 | 职责 |
|---|---|
| `scripts/plot/utils.R` | 画图基建：语义色板 PAL（up=#B2182B / down=#2166AC / ns=#BDBDBD…）、pal_diverging(RdBu)/pal_sequential(viridis)/pal_categorical、图幅常量 W_SINGLE=89 / W_ONE_HALF=136 / W_DOUBLE=183 mm、`theme_paper()`（图例框外右侧+双 vertical）、`save_png()`（PNG 300dpi cairo + 同名 PDF 矢量）、`fig_filename()` 五段图名、`fmt_p()`/`sig_stars()` |
| `scripts/plot/volcano.R` | 火山图：红 Up / 蓝 Down / 灰 NS、阈值虚线、副标题三区计数 |
| `scripts/plot/heatmap.R` | 热图：行 Z-score 夹 ±3、RdBu 发散色、hclust 排序、标签预算算术判据（决策落 label_decisions.csv） |
| `scripts/plot/enrichment.R` | GO/KEGG 三点式点图：size=Count、colour=-log10(p.adj)、x=GeneRatio |
| `scripts/runner.R` | 分析 runner 参考实现（P2 联调版）：CHIP-06/RNA-07 确定性模拟 + 出图段接线示例 |
| `scripts/preflight.R` | L3 前置检查参考实现 |

## 图规范（硬约束）

1. **图名五段坐标**：`<阶段>-<模块>-<图>-unit<N>-<名称>`（geo=01 / scrna=02 / spatial=03）
2. **图幅毫米制**：单栏 89mm / 1.5 栏 136mm / 双栏 183mm 封顶；300dpi + 同名 PDF
3. **图例**：一律框外右侧纵向（`legend.direction` + `legend.box` 均 vertical）
4. **配色一门禁**：语义色板全篇一致（红=高/实验、蓝=低/对照、灰=NS）
5. **统计三件套进图**：HR(95%CI)+p；p<0.001 显示 `< 0.001`；ns 显式标注

## 与其他层的关系

- 上游 **P2**（geo/scrna/spatial 三 skill 仓 + biopipeline-template 执行端）产出 `tables/*.csv`
- 下游 **P4**（p4-typesetting 仓）消费 `figures/` + `figure_index.json` 生成三线表/图注/SCI_HANDOFF.md
- MVP 完整闭环仍在 biopipeline-template 仓（P2–P4 串联）；本仓是 P3 层的独立载体，便于单独迭代出图能力

## 历史

本仓前身为 `MVP-bioworkflow-version2`（P1–P4 端到端验证体系第二版），2026-09-28
经用户拍板重命名为 P3 专职仓库，旧内容已替换（git 历史可回溯）。
