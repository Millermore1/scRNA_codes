# QC
library(Seurat)
library(ggplot2)
library(dplyr)

# 提取所需列
df <- scRNA_HKU@meta.data %>%
  select(nCount_RNA, nFeature_RNA, celltype)

# 按 UMI 中位数对细胞类型排序（和参考图一样，从大到小排列）
df$celltype <- factor(
  df$celltype,
  levels = df %>%
    group_by(celltype) %>%
    summarise(med = median(nCount_RNA, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(celltype)
)

p <- ggplot(df, aes(x = celltype, y = nCount_RNA)) +
  geom_boxplot(fill = "#a9cce3", color = "black",
               outlier.shape = NA) +        # 不显示离群点
  scale_y_continuous(
    breaks = c(0, 10000, 20000),
    labels = c("0", "10k", "20k")
  ) +
  labs(title = "Number of unique transcripts (UMIs) detected\nper cell across major cell types",
       x = NULL, y = "UMI counts") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15))

p

table(scRNA_HKU$celltype)


library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(viridis)
Idents(scRNA_HKU) <- "celltype"

# 1. 找每个细胞类型的 marker
markers <- FindAllMarkers(
  scRNA_HKU,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.5
)

# 2. 每种细胞类型取 top 8 marker
top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_max(order_by = avg_log2FC, n = 8) %>%
  ungroup()

# 3. 细胞太多时抽样，加快绘图并避免内存爆炸
set.seed(42)
cells_use <- unlist(lapply(levels(Idents(scRNA_HKU)), function(ct) {
  cells <- WhichCells(scRNA_HKU, idents = ct)
  sample(cells, min(length(cells), 300))
}))



# ---------- 1. 每个细胞类型的平均表达 ----------
avg_exp <- AverageExpression(
  scRNA_HKU,
  assays   = "RNA",
  features = top_markers$gene,
  group.by = "celltype",
  slot     = "data"
)$RNA

avg_exp <- as.matrix(avg_exp)

# ---------- 2. 按基因排序（同一细胞类型的 marker 排在一起）----------
top_markers <- top_markers %>%
  arrange(factor(cluster, levels = colnames(avg_exp)), desc(avg_log2FC))
avg_exp <- avg_exp[top_markers$gene, ]

# ---------- 3. 行 Z-score 标准化 ----------
mat <- t(scale(t(avg_exp)))
mat[is.na(mat)] <- 0

# ---------- 4. 颜色与注释 ----------
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#3B4CC0", "#F7F7F7", "#B40426")     # Nature 常用蓝-白-红
  # 若想保留紫黄风格，换成：viridis::viridis(3)
)

gene_group <- factor(
  top_markers$cluster,
  levels = colnames(avg_exp)
)

group_cols <- setNames(
  brewer.pal(n = length(levels(gene_group)), name = "Set2"),
  levels(gene_group)
)

col_anno <- HeatmapAnnotation(
  `Cell type` = colnames(mat),
  col  = list(`Cell type` = group_cols),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  simple_anno_size   = unit(4, "mm")
)

row_anno <- rowAnnotation(
  Marker = gene_group,
  col    = list(Marker = group_cols),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(3, "mm")
)

# ---------- 5. 绘制 ----------
ht <- Heatmap(
  mat,
  name = "Z-score",
  col  = col_fun,
  
  cluster_rows        = FALSE,       # 保持 marker 按细胞类型分组的顺序
  cluster_columns     = FALSE,       # 保持细胞类型顺序
  show_column_dend    = FALSE,
  show_row_dend       = FALSE,
  
  row_names_gp        = gpar(fontsize = 7, fontface = "italic"),
  column_names_gp     = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot    = 45,
  
  top_annotation      = col_anno,
  right_annotation    = row_anno,
  
  row_split           = gene_group,  # 按细胞类型分块
  row_title           = NULL,
  row_gap             = unit(1.5, "mm"),
  
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  ),
  
  border            = FALSE,
  use_raster        = TRUE,
  raster_quality    = 4
)

draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))

library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)

Idents(scRNA_HKU) <- "celltype"

# ---------- 1. 找 marker，每种细胞类型取 top 20（想要更多就改 n）----------
markers <- FindAllMarkers(
  scRNA_HKU,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.5
)

top_markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%   # <- 改这里控制基因数
  ungroup()

# ---------- 2. 每种细胞类型抽样 300 个细胞（想让图更密就调大）----------
set.seed(42)
n_cells <- 300                                   # <- 改这里控制列数

cells_use <- unlist(lapply(levels(Idents(scRNA_HKU)), function(ct) {
  cells <- WhichCells(scRNA_HKU, idents = ct)
  sample(cells, min(length(cells), n_cells))
}))

# ---------- 3. 提取 scale.data（对 marker 基因做缩放）----------
obj_sub <- subset(scRNA_HKU, cells = cells_use)
obj_sub <- ScaleData(obj_sub, features = top_markers$gene, verbose = FALSE)

mat <- GetAssayData(obj_sub, assay = "RNA", slot = "scale.data")[top_markers$gene, ]
mat <- as.matrix(mat)
mat[mat > 2.5]  <- 2.5    # 截断极值，防止个别离群细胞拉爆色阶
mat[mat < -2.5] <- -2.5

# ---------- 4. 细胞按类型排序、基因按所属类型排序 ----------
celltype_vec <- factor(obj_sub$celltype, levels = levels(Idents(scRNA_HKU)))
mat <- mat[, order(celltype_vec)]
celltype_vec <- celltype_vec[order(celltype_vec)]

top_markers <- top_markers %>%
  arrange(factor(cluster, levels = levels(celltype_vec)), desc(avg_log2FC))
mat <- mat[top_markers$gene, ]

# ---------- 5. 注释与配色 ----------
ct_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(nlevels(celltype_vec)),
  levels(celltype_vec)
)

col_anno <- HeatmapAnnotation(
  `Cell type` = celltype_vec,
  col  = list(`Cell type` = ct_cols),
  annotation_name_gp  = gpar(fontsize = 10, fontface = "bold"),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    ncol = 1
  ),
  simple_anno_size = unit(4, "mm")
)

gene_group <- factor(top_markers$cluster, levels = levels(celltype_vec))

col_fun <- colorRamp2(
  c(-2, 0, 2.5),
  c("#440154", "#3B528B", "#FDE725")     # 紫-黄（viridis 风格，与你示例图一致）
  # 想换 Nature 蓝白红就改成: c("#3B4CC0", "#F7F7F7", "#B40426")
)

# ---------- 6. 绘制 ----------
ht <- Heatmap(
  mat,
  name = "Scaled\nexpression",
  col  = col_fun,
  
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  show_column_names = FALSE,            # 细胞太多，不显示列名
  
  row_names_gp     = gpar(fontsize = 9, fontface = "italic"),
  
  top_annotation   = col_anno,
  
  column_split     = celltype_vec,      # 按细胞类型分列分块
  column_title     = levels(celltype_vec),
  column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
  column_title_rot = 45,
  column_gap       = unit(1.5, "mm"),
  
  row_split        = gene_group,        # 基因按所属细胞类型分块
  row_title        = NULL,
  row_gap          = unit(1, "mm"),
  
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  
  border         = FALSE,
  use_raster     = TRUE,                # 栅格化，几千列也很流畅
  raster_quality = 4
)

draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))

# ---------- 7. 导出 ----------
pdf("celltype_marker_heatmap_percell.pdf", width = 14, height = 14)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()


head(scRNA_HKU)
scRNA_HKU <- subset(scRNA_HKU, subset = mt_percent < 5)


library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ---------- 1. 设定分组顺序：WT 在最前 ----------
group_order <- c("WT", "KOF_vector_CNO", "KOF_Gq_CNO", "KOF")
scRNA_HKU$group <- factor(scRNA_HKU$group, levels = group_order)

# ---------- 2. 蓝色系配色（与顺序一一对应）----------
group_cols <- c(
  "WT"             = "#08306B",   # 深海军蓝（对照组，最醒目）
  "KOF_vector_CNO" = "#2171B5",   # 中蓝
  "KOF_Gq_CNO"     = "#6BAED6",   # 浅蓝
  "KOF"            = "#C6DBEF"    # 淡蓝
)

# ---------- 3. 绘图 ----------
qc_metrics <- c("nFeature_RNA", "nCount_RNA", "mt_percent", "HB_percent")
ylabs <- c("Detected genes", "UMI counts", "MT (%)", "HB (%)")

plots <- lapply(seq_along(qc_metrics), function(i) {
  VlnPlot(
    scRNA_HKU,
    features = qc_metrics[i],
    group.by = "group",
    pt.size  = 0,
    cols     = group_cols
  ) +
    geom_boxplot(width = 0.15, fill = "white",
                 outlier.size = 0, alpha = 0.9) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = NULL, y = ylabs[i], title = NULL) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1,
                                 face = "bold", colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.title.y = element_text(face = "bold")
    )
})

p_qc <- wrap_plots(plots, nrow = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18, face = "bold"))

p_qc

ggsave("QC_violin_by_group_blue.pdf", p_qc, width = 9, height = 7)











