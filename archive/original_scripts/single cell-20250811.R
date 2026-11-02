# 提取Exci.neuron亚群的细胞名
exci_cells <- WhichCells(newdata, expression = celltype == "Exci.neuron")

# 用subset函数提取该亚群，得到一个新Seurat对象
exci_seurat <- subset(newdata, cells = exci_cells)

# 标准化
exci_seurat <- NormalizeData(exci_seurat)
# 找高变基因
exci_seurat <- FindVariableFeatures(exci_seurat)
# 标准化
exci_seurat <- ScaleData(exci_seurat)
# PCA
exci_seurat <- RunPCA(exci_seurat)
# 最近邻
exci_seurat <- FindNeighbors(exci_seurat, dims = 1:10)
# 聚类
exci_seurat <- FindClusters(exci_seurat, resolution = 0.05)
# UMAP或tSNE降维
exci_seurat <- RunUMAP(exci_seurat, dims = 1:10)
# 或
# exci_seurat <- RunTSNE(exci_seurat, dims = 1:10)


FeaturePlot(exci_seurat, features = c("PRMT8"))

FeaturePlot(
  exci_seurat,
  features = c("PRMT8"),
  cols = c("lightgrey", "#2166AC", "#B2182B"),  # 更深的蓝和红
  min.cutoff = "q05",    # 可选，剔除极端低值
  max.cutoff = "q95"     # 可选，剔除极端高值，增强对比
)

head(exci_seurat)
table(exci_seurat@meta.data$orig.ident)


exci_seurat@meta.data$group <- ifelse(
  grepl("^AD", exci_seurat@meta.data$orig.ident), "AD",
  ifelse(
    grepl("^NC", exci_seurat@meta.data$orig.ident), "NC",
    ifelse(
      grepl("^PD", exci_seurat@meta.data$orig.ident), "PD",
      "Other"
    )
  )
)
table(exci_seurat@meta.data$group)

Idents(exci_seurat) <- exci_seurat@meta.data$group
de_markers <- FindMarkers(
  exci_seurat,
  ident.1 = "AD",
  ident.2 = "NC",
  min.pct = 0.25,         # 可根据需要调整
  logfc.threshold = 0.25  # 可根据需要调整
)

sig_degs <- subset(de_markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
head(sig_degs)
de_markers$significance <- "Not Sig"
de_markers$significance[de_markers$p_val_adj < 0.05 & de_markers$avg_log2FC > 0.5] <- "Up"
de_markers$significance[de_markers$p_val_adj < 0.05 & de_markers$avg_log2FC < -0.5] <- "Down"
library(ggplot2)

ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano plot: AD vs NC", x = "avg_log2FC", y = "-log10(adj P-value)") +
  theme_bw()

# 提取PRMT8的表达数据和组别信息
expr_data <- FetchData(exci_seurat, vars = c("PRMT8", "group"))

# 找出表达量大于0的细胞
cells_with_expr <- rownames(expr_data)[expr_data$PRMT8 > 0]

# 子集Seurat对象，只保留表达量大于0的细胞
exci_seurat_filtered <- subset(exci_seurat, cells = cells_with_expr)

# 绘制剔除表达量为0的PRMT8的小提琴图并进行t检验
VlnPlot(exci_seurat_filtered, features = "PRMT8", group.by = "group", pt.size = 0) +
  stat_compare_means(method = "t.test")

DotPlot(newdata, features = c("CLDN5","PECAM1","TMEM119","P2RY12","MAP2","GAD1","SLC1A3","GFAP","MOG","MBP","SLC17A7","PRMT8"))
library(ggplot2)
DotPlot(
  newdata,
  features = c("CLDN5","PECAM1","TMEM119","P2RY12","MAP2","GAD1","SLC1A3","GFAP","MOG","MBP","SLC17A7","PRMT8"),
  dot.scale = 15,
  col.min = 0,
  col.max = 1
) + RotatedAxis() + theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold")
  )

# 提取PRMT8表达数据
prmt8_expr <- FetchData(exci_seurat, vars = "PRMT8")

# 根据PRMT8表达水平进行分组（比如中位数分组）
median_expr <- median(prmt8_expr$PRMT8)
exci_seurat$PRMT8_group <- ifelse(prmt8_expr$PRMT8 > median_expr, "High", "Low")

# 查看分组情况
table(exci_seurat$PRMT8_group)

# 进行差异表达分析（在这两个组之间）
diff_genes <- FindMarkers(exci_seurat, ident.1 = "High", ident.2 = "Low", group.by = "PRMT8_group")

# 查看差异基因
head(diff_genes)

significant_genes <- subset(diff_genes, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
write.csv(significant_genes, file = "DEGs_PRMT8_high_vs_low.csv")


library(dplyr)
# 根据PRMT8表达水平，将细胞分为4组
exci_seurat$PRMT8_quartile <- ntile(prmt8_expr$PRMT8, 4)
# 转换为标签
exci_seurat$PRMT8_quartile <- factor(exci_seurat$PRMT8_quartile, labels = c("Q1_Lowest", "Q2_Lower", "Q3_Upper", "Q4_Highest"))

# 查看分组情况
table(exci_seurat$PRMT8_quartile)

# 设置分组信息
Idents(exci_seurat) <- "PRMT8_quartile"

# 进行差异表达分析，比较Q4和Q1
diff_genes <- FindMarkers(
  exci_seurat,
  ident.1 = "Q4_Highest",
  ident.2 = "Q1_Lowest",
  group.by = "PRMT8_quartile",
  min.pct = 0.1,            # 只考虑表达在至少10%的细胞中
  logfc.threshold = 0.25  # 只考虑对数FC大于0.25的基因
)

# 查看差异表达结果
head(diff_genes)

# 筛选显著差异基因（调整p值<0.05，绝对logFC>0.5）
significant_genes <- subset(diff_genes, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)

# 导出差异基因
write.csv(significant_genes, file = "DEGs_PRMT8_Q4_vs_Q1.csv")


# 读取数据
data <- read.csv("GOGP.csv", header = TRUE, stringsAsFactors = FALSE)

# 为避免横坐标标签重叠，建议横向条形图
barplot(
  height = data$Enrichment,
  names.arg = data$Channel,
  horiz = TRUE,
  las = 1, # 标签横向显示
  col = "skyblue",
  main = "Channel Enrichment 条形图",
  xlab = "Enrichment值"
)

# 可选：调整边距以防止标签被遮挡
par(mar = c(5,20,4,2)) # 下、左、上、右边距
barplot(
  height = data$Enrichment,
  names.arg = data$Channel,
  horiz = TRUE,
  las = 1,
  col = "skyblue",
  main = "Channel Enrichment 条形图",
  xlab = "Enrichment值"
)


rm(sig)

