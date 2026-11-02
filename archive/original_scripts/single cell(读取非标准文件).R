
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(harmony)
library(Matrix)

# 设置文件路径
path <- "C:/Users/黄廉/Desktop/GSE186537"

# 读取细胞元数据（假设文件名为 Human_cell_meta.txt）
cell_meta <- read.table(paste0(path, "/Human_cell_meta.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 读取基因信息（假设文件名为 Human_genes.txt）
genes <- readLines(paste0(path, "/Human_genes.txt"))

# 读取计数矩阵（假设文件名为 Human_counts.mtx）
library(Matrix)
counts <- readMM(paste0(path, "/Human_counts.mtx"))

# 设置行名和列名
rownames(counts) <- genes
colnames(counts) <- cell_meta$cell_id  # 根据实际列名调整

# 设置列名为 "cell_name" 列中的值
colnames(counts) <- cell_meta$cell_name

# 再次确认
head(colnames(counts))
#建立seurat对象
seurat_obj <- CreateSeuratObject(counts = counts, meta.data = cell_meta)


seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
DimPlot(seurat_obj)
library(ggplot2)

# 运行UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# 提取UMAP坐标
umap_coords <- Embeddings(seurat_obj, reduction = "umap")

# 查看结构
str(umap_coords)

# 获取类别信息
cluster_info <- seurat_obj@meta.data$seurat_clusters

# 创建数据框
umap_df <- data.frame(umap_coords, cluster = cluster_info)

# 绘图
library(ggplot2)
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size=0.5) +
  theme_minimal()

# 查看细胞总数
# 获取细胞名的数量，即细胞总数
n_cells <- ncol(seurat_obj)
print(n_cells)

