# PRMT8 抑制性神经元单核 RNA 测序分析代码

本仓库是论文 **PRMT8 在抑制性神经元中的作用（人类多队列 snRNA-seq + 小鼠条件敲除/化学遗传学）**
的完整可复现分析代码，从 Cell Ranger 原始输出一直到成图。

原始 43 个脚本的整理记录见 [`docs/migration_notes.md`](docs/migration_notes.md)。

---

## 分析内容

| 数据集 | 设计 | 在论文中的作用 |
|---|---|---|
| 人类多队列脑 snRNA-seq | AD、PD、MS、EPI 与对照，54 个样本，9 个 GEO 系列 | PRMT8 在疾病中的表达与共表达 |
| 小鼠 *Prmt8* 条件敲除 + DREADD | `WT`、`KOF`、`KOF_vector_CNO`、`KOF_Gq_CNO` | 因果验证：基因敲除与 Gq 激活 |
| 人脑类器官（GSE203201） | ASD vs 对照，D0–D100 | PRMT8 的发育轨迹 |

分析包括：质控与 Harmony 整合、细胞注释、细胞比例、按细胞类型的差异表达、
GO/KEGG 富集与通路网络、CellChat 细胞通讯、hdWGCNA 共表达模块、模块打分与全套出版级图。

---

## 目录结构

```
config/config.R      所有参数：路径、阈值、对比组、基因集
config/palettes.R    配色与因子顺序
R/00_setup.R         入口（定位仓库根目录、载入配置、加载包）
R/01 … R/11          按执行顺序编号的 11 个分析模块
R/utils/             io / analysis / plotting 三个公共函数文件
data/metadata/       队列映射、marker 表、基因集
data/raw/            输入数据（不入库）
data/processed/      流程产生的 Seurat 对象（不入库）
results/             图（PDF+PNG）与表（CSV）
env/                 安装脚本与版本清单
docs/                Methods、图-文件对照、数据可用性、迁移说明
tests/smoke_test.R   不需要数据的自检
archive/original_scripts/  原始 43 个脚本（仅存档，不参与运行）
```

---

## 三步跑起来

### 1. 安装依赖

```r
source("env/install_packages.R")
```

开发环境：R 4.3.0、Seurat 5.3.0、CellChat 1.6.1、clusterProfiler 4.10.0
（完整版本见 `env/package_versions.tsv`）。

### 2. 指定数据位置

两种方式，任选其一：

- 把数据放进 `data/raw/`（路径见 `config/config.R` 的默认值）；
- 用环境变量指向实际位置（推荐，数据通常不在 C 盘）：

```r
# 写入 ~/.Renviron
PRMT8_INPUT_DIR       = "D:/scrnadata"
PRMT8_MOUSE_ATLAS_RDS = "D:/scrnadata/scRNA_HKU.rds"
PRMT8_ASD_10X_DIR     = "D:/scrnadata/GSE203201"
PRMT8_BRAIN_COR_XLSX  = "D:/scrnadata/correlations_PRMT8_brain_region.xlsx"
```

### 3. 运行

```bash
Rscript run_all.R           # 全部模块
Rscript run_all.R 01 02 05  # 只跑这几步
Rscript run_all.R --list    # 查看模块说明
```

也可以在 RStudio 里逐个运行 `R/01_…` 到 `R/11_…`：每个模块都是独立脚本，
读取上一步写出的 `.rds`，因此可以随时中断、从任意一步继续。

第一次运行前建议先自检（不需要数据，几秒钟）：

```bash
Rscript tests/smoke_test.R        # 不需要数据、不需要 Seurat
Rscript tests/integration_test.R  # 需要 Seurat + ComplexHeatmap
```

前者解析所有脚本、校验配置与配色，并用玩具数据测试公共函数；后者用合成的
Seurat 对象把分析与绘图函数整条跑一遍，能在没有真实数据时发现接口错误。

---

## 准备输入对象

模块 02 之前需要两个对象：

1. **`data/processed/human_atlas_integrated.rds`** —— 由模块 01 从 Cell Ranger 目录
   自动构建；若整合已在别处完成，也可以直接提供。必需列：`orig.ident`、`group`、
   `seurat_clusters`，以及 `umap` 降维。
2. **`data/processed/mouse_ko_integrated.rds`** —— 小鼠敲除图谱，从你现有的对象导出：

```r
library(Seurat)
scRNA_HKU$group <- factor(scRNA_HKU$group,
                          levels = c("WT", "KOF", "KOF_vector_CNO", "KOF_Gq_CNO"))
saveRDS(scRNA_HKU, "data/processed/mouse_ko_integrated.rds")
```

必需列：`orig.ident`、`group`、`celltype`、`seurat_clusters`，以及 `umap` 降维。
若缺 `mt_percent`，模块 01 会自动计算。

---

## 输出

每张图同时输出两种格式：`cairo_pdf`（可进 Illustrator 精修、投稿用）与 300 dpi PNG（预览用）。
文件名前缀是模块编号，任何时候都能反查到脚本：

```
results/figures/05_prmt8_violin_mouse.pdf
results/tables/05_prmt8_test_mouse.csv
```

全部输出与对应图号的清单见 [`docs/figure_mapping.md`](docs/figure_mapping.md)
（其中的 “Manuscript panel” 一列需要你填）。

---

## 写 Methods 时可直接取用

[`docs/methods.md`](docs/methods.md) 逐节记录了：质控阈值、整合参数、聚类分辨率、
注释依据、统计检验方法（含样本水平检验）、富集参数、CellChat 与 hdWGCNA 设置、
软件版本。可以直接改写进 Methods。

---

## 上传到 GitHub

```bash
cd PRMT8-inhibitory-neuron-scRNAseq
git init
git add .
git commit -m "Initial commit: PRMT8 interneuron snRNA-seq analysis pipeline"

# 先在 GitHub 上建好空仓库，然后：
git branch -M main
git remote add origin https://github.com/MillerMore1/scRNA_codes.git
git push -u origin main

git tag -a v1.0.0 -m "Version accompanying the manuscript"
git push origin v1.0.0
```

上传前请替换这些占位内容：

| 位置 | 需要填什么 |
|---|---|
| `CITATION.cff` | 姓名、ORCID、单位、邮箱 |
| `LICENSE` | 版权人 |
| `README.md` / `README.zh-CN.md` | 联系邮箱 |
| `docs/data_availability.md` | 小鼠数据登录号、仓库地址 |

首次推送后建议在 Zenodo 归档打标签的版本，拿到 DOI 后填回
`docs/data_availability.md`。每次 push 会自动运行
`.github/workflows/ci.yml` 里的自检。

## 联系与许可

- 引用格式见 [`CITATION.cff`](CITATION.cff)（作者、单位、邮箱待你补全）。
- 代码许可：[MIT](LICENSE)。
- 原始数据集的登录号与使用条款见 [`docs/data_availability.md`](docs/data_availability.md)。

## 联系方式

William Huang — MillerMore1@users.noreply.github.com
