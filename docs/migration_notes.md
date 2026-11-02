# 迁移说明（原 43 个 R 文件 → 新代码库）

本文档记录整理过程，方便你（和审稿人）追溯：每个原始脚本去了哪里、哪些问题被修复、哪些地方需要你确认。

## 1. 总体情况

- 原始文件：43 个 `.R`，原本平铺在一个文件夹里，共约 7,500 行。
- 其中 1 组是完全重复文件（MD5 相同）：`all_neurons_analyse.R` 与 `GO_Figure_6_2 (2).R`。
- 另有若干“同一分析的不同版本”（`inhibit.R` → `inhibit_6_29.1.R`、`cell_chat_5.09.R` → `cell chat-6.29.R` 等），已合并为同一模块中的最终版本。
- 整理后：11 个按执行顺序编号的分析模块 + 3 个公共函数文件 + 配置 + 文档 + 测试。

原始文件完整保留在 `archive/original_scripts/`，**不参与**流水线执行，仅作溯源。

## 2. 文件对照表

| 原始脚本 | 归入模块 | 处理方式 |
|---|---|---|
| `single cell-20240607.R` | 01 | 10x 读取循环、MT/HB 质控、merge、Harmony 流程 → `read_cellranger_dirs()` / `apply_qc_filter()` / `integrate_atlas()` |
| `single cell-20250802.R` | 01 | 同上（另一批数据的同一流程），合并去重 |
| `single cell-20250801.R` | 01 | h5 读取 → `read_h5_matrix()` |
| `single cell(读取非标准文件).R` | 01 | MatrixMarket 三件套读取 → `read_mtx_triplet()` |
| `正式实验代码9.5.R` | 01 / 03 | 队列分组映射 → `data/metadata/cohort_mapping.csv`；分组配色 → `config/palettes.R` |
| `QC.R` | 01 / 02 | QC 小提琴图 → `plot_qc_violin()`；marker 热图 → `build_marker_heatmap()` |
| `质控.R` | 01 / 02 | 血红蛋白比例 QC，并入 `add_qc_metrics()` |
| `doubleFinder_28.1.R` | 01 / env | 安装步骤 → `env/install_packages.R`；marker DotPlot → `02_annotation.R` |
| `正式code908.R` | 02 / 03 | 细胞注释、分组 UMAP、点图、比例图 |
| `正式code908-2.R` | 02 | 与 `full_codes.R` 大幅重叠，取并集 |
| `full_codes.R` | 02 / 03 / 05 / 10 | 注释、UMAP 分面、比例图、抑制性亚群、PRMT8 热图 |
| `marker gene.R` | 02 / 10 | 指定顺序的 marker DotPlot → `plot_marker_dotplot()` |
| `inhibit.R` | 05 | 抑制性神经元亚聚类 + 亚型注释 |
| `inhibit_2.R` | 05 | PVALB 亚型 DEG + GO → `05_inhibitory_neurons.R` |
| `inhibit_3.R` | 05 | UMAP 配色、亚型比例、火山图 |
| `inhibit_4.R` | 04 / 05 | 全抑制性神经元火山图（标注 ASD 基因）、PV 亚型 DEG |
| `inhibit_5_28.1.R` | 04 / 05 | `run_deg_by_celltype()` 的雏形（保留）、Venn、DEG 计数热图 |
| `inhibit_6_29.1.R` | 03 / 05 | 细胞比例、marker 小提琴、HP_AUTISM 打分、热图 |
| `excit_2.R` | 06 | 兴奋性神经元亚聚类 + 皮层分层 + 亚型 DEG |
| `all_neurons_analyse.R` | 05 / 06 / 07 | 神经元再聚类、E/I 注释、兴奋性亚群、GO 条形图 |
| `GO_Figure_6_2 (2).R` | — | **与 `all_neurons_analyse.R` 完全重复（MD5 相同），已弃用** |
| `DEG_enrichment_network.R` | 07 | GO 富集 + 网络图（节点饼图、Jaccard 边） |
| `Fig.6_D.R` | 04 / 07 | 按细胞类型的 KOF vs WT DEG + GO 网络图 |
| `GO_Figure_6.R` | 07 | GO 条形图（条内嵌基因名） |
| `GO_Figure_6_2 (1).R` | 07 | 同上，但输入为 KEGG 结果 → 合并为同一函数 |
| `code_for_GO_enrichment.R` | 06 / 07 | 兴奋性 DEG 导出 + GO 条形图 |
| `最新的修改figure.R` | 04 / 07 | 按 celltype×group 组合身份做 DEG（已由 `run_deg_by_celltype()` 覆盖） |
| `cell chat.R` | 08 | CellChat 单组 + 多组、圆形图美化（统一标尺） |
| `cell_chat_5.09.R` | 08 | 神经元亚群注释 + 4 组 CellChat + 中心度 |
| `cell chat-6.29.R` | 08 | 上一版的重写版，作为最终流程 |
| `PRMT8KO_snRNA_seq_cellchat.R` | 08 | 突触相关通路气泡图、信息流、整体对比 |
| `code_for_cellchat.R` | 08 | 手工录入数值画的 outgoing/incoming 散点 → 改为从 CellChat 结果直接计算并导出 |
| `WGCNA.R` | 09 | hdWGCNA 全流程 + PRMT8 模块 + ggraph 重绘 |
| `Heatmap-last.R` | 10 | 离子通道基因热图 |
| `Heatmap-last-2.R` | 10 | PRMT8 FeaturePlot + 基因集热图 |
| `热图绘制0916.R` | 10 | GABA/谷氨酸受体等基因集热图 |
| `热图绘制0916-2.R` | 05 / 10 | 抑制性 cluster 再聚类、cluster 比例、GO 条形图 |
| `heatmap_correlation_region.R` | 10 | 脑区相关性热图（ComplexHeatmap） |
| `相关性热图code.R` | — | **模拟数据的配色模板，不产出任何论文图，仅存档** |
| `细胞比例.R` | 03 | 堆积比例柱状图 |
| `ASD_S1.R` | 11 | GSE203201 文件整理、读取、Harmony、D30 DEG、PRMT8 时序 |
| `single cell-ASD.R` | 11 | GSE203201 h5 读取、注释、limma + GSEA、热图 |
| `single cell-20250811.R` | 05 / 10 | AD vs NC DEG、PRMT8 高/低分组、四分位 DEG、GSEA |

## 3. 修复的问题

按严重程度排列。前 6 条是**会直接报错或产出错误结果**的问题。

| # | 位置 | 问题 | 处理 |
|---|---|---|---|
| 1 | `single cell-ASD.R:73-74` | 字符串跨行断裂（`"Neural` 换行 `crest/Mes"`），整个文件无法解析 | 重写为完整的细胞类型标签；注释改为标记法 |
| 2 | `excit_2.R:129-130` | `DEGs_excit_KOFvsWT <- write.csv()` 空参数调用，必然报错 | 删除，改为显式写出 CSV |
| 3 | `inhibit_4.R:7-20` | `FindMarkers(object = pv, ...)` 中 `pv` 在本文件未定义；随后引用未定义的 `deg.sig` | 在模块中先按亚型取子集再计算，变量名统一 |
| 4 | `inhibit_3.R:163-171` | `VIP <- subset(...)` 之后却对 `pv` 求 DEG，对象名与内容不一致 | 统一为按亚型循环 |
| 5 | `inhibit_4.R:7` | `subset(inh, subset = inhibit_celltype == "PV")`，但标签是 `"PVALB"` → 得到 0 个细胞 | 标签统一为 `PVALB`，并在配置中校验 |
| 6 | `full_codes.R:31` | `save.image(file = scRNAharmony_9.5)` 传入对象而非文件名；`:283` `ggsave("UMAP_total.png", )` 尾逗号 | 改为 `write_object()` / `save_figure()` |
| 7 | `full_codes.R:55-56` | 该代码块在 `df` 创建前引用 `df` 的范围（`df` 在 200 行后才定义） | 先取坐标与 metadata 再绘图 |
| 8 | `inhibit_6_29.1.R:136-165` | 200 个基因的向量硬编码在脚本里 | 移到 `data/metadata/hp_autism_genes.txt` |
| 9 | `code_for_cellchat.R:9-24` | 手抄的 outgoing/incoming 数值（无法溯源） | 改为从 CellChat 对象计算并导出表 |
| 10 | `full_codes.R:114` 等 | 组内细胞数（109037 等）硬编码 | 改为从对象统计；文档中注明出处 |
| 11 | 多个文件 | 绝对路径（`D:/scrnadata`、`E:/scdata`、`G:/sn_RNAseq_data`、`C:/Users/黄廉/Desktop/...`） | 全部集中到 `config/config.R`，支持环境变量覆盖 |
| 12 | `single cell-20250811.R:86` | `stat_compare_means()` 未加载 ggpubr，且未做多重检验 | 改为样本水平的 Wilcoxon，并做 BH 校正 |
| 13 | `single cell-20250801.R:29` | `ScaleData(features = rownames(obj))` 对全部基因标准化（极慢、易爆内存） | 只对高变基因 |
| 14 | `full_codes.R:314` | `theme(margin = margin())` 不是合法主题元素 | 改为 `plot.margin` |
| 15 | `inhibit_5_28.1.R` 等 | `dplyr::top_n()` 在 dplyr ≥1.1 已弃用，且对并列值不稳健 | 改为显式排序后取前 n |
| 16 | `热图绘制0916-2.R:111` | `rm(list = ls())` 会清空整个会话 | 删除 |
| 17 | 多个文件 | `p`、`pv`、`inh`、`newdata`、`scRNA_seurat` 等对象名在不同脚本里指向不同版本，无法判断哪个是最终数据 | 统一对象名（`human`/`mouse`/`inh`/`excit`）并落盘为 `.rds` |
| 18 | `cell chat.R` | 直接用 `@net$count`、`@netP$centr` 等 slot，跨 CellChat 版本不稳定 | 保留但固定版本（1.6.1），并在文档注明 |
| 19 | `QC.R:283` | `subset(scRNA_HKU, ...)` 重新赋值后没有保存，后续脚本用的仍是旧对象 | 在模块中显式保存过滤后对象 |
| 20 | `01` 相关脚本 | Seurat v5 分层数据未 `JoinLayers` 就做 `FindMarkers` | 统一在 `ensure_joined_layers()` 中处理 |

## 4. 新增内容（原脚本没有的）

- `run_all.R`：一条命令跑完整流程，或按模块前缀跑其中几步。
- `tests/smoke_test.R`：不依赖数据的自检（语法解析、配置一致性、配色覆盖、辅助函数），`Rscript tests/smoke_test.R` 即可运行。
- 样本水平的比例检验与 PRMT8 表达检验（原来是细胞水平直接比较，审稿人通常会质疑伪重复）。
- 集群映射的严格校验：聚类数与映射长度不一致时立刻报错，而不是静默错标。
- 所有图同时输出 `cairo_pdf`（可进 Illustrator）与 300 dpi PNG。
- 每次运行自动记录 `results/session_info.txt`。

## 5. 需要你确认的地方

1. **小鼠数据的存放**：`data/processed/mouse_ko_integrated.rds` 需要由你现有对象导出（脚本按你的 `scRNA_HKU` 结构编写：`group` 与 `celltype` 两列）。导出方式见 `README.md` 的“Preparing the input objects”。
2. **疾病队列的样本元数据**：`cohort_mapping.csv` 里的 54 个样本已按原脚本整理，请核对（尤其 `GSE161045` 同时出现在 AD 与 PD）。
3. **HP_AUTISM 基因集来源**：请在 Methods 里补上数据库/文献出处。
4. **脑区相关性矩阵**：`correlations_PRMT8_brain_region.xlsx` 的来源与相关系数类型（Pearson/Spearman）。
5. **图号**：`docs/figure_mapping.md` 里有一列 “Manuscript panel” 待填，填完就是投稿时的图-文件对照表。
