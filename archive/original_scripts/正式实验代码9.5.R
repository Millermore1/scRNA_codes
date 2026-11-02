# 定义每个组别对应的orig.ident列表
group_mapping <- list(
  "AD" = c("GSE157827_AD1", "GSE157827_AD10", "GSE157827_AD11", "GSE157827_AD12",
           "GSE157827_AD2","GSE157827_AD3","GSE157827_AD4","GSE157827_AD5","GSE157827_AD6","GSE157827_AD7","GSE157827_AD8","GSE157827_AD9",
           "GSE161045_AD13", "GSE161045_AD14", "GSE161045_AD15", "GSE161045_AD16","GSE161045_AD13","GSE161045_AD14","GSE161045_AD15","GSE161045_AD16"),
  "MS" = c("MS122", "MS121", "MS176", "MS242"),
  "PD" = c("GSE140231_PD_C1", "GSE140231_PD_C1B", "GSE140231_PD_C2B",
           "GSE140231_PD_N1B", "GSE140231_PD_N2B", "GSE140231_PD_N3",
           "GSE140231_PD_N4", "GSE140231_PD_N4B", "GSE140231_PD_N5", "GSE140231_PD_N5B",
           "GSE161045_PD1", "GSE161045_PD2", "GSE161045_PD3", "GSE161045_PD4"),
  "EPI" = c("GSE302285_Epi1", "GSE302285_Epi2", "GSE302285_Epi3",
            "GSE302285_Epi4", "GSE302285_Epi5", "GSE302285_Epi6"),
  "NC" = c("CO14","CO25","CO28","CO39","GSE157827_NC1","GSE157827_NC2","GSE157827_NC3","GSE157827_NC4","GSE157827_NC5","GSE157827_NC6","GSE157827_NC7",
           "GSE157827_NC8","GSE157827_NC9","SD48.16")
)

# 获取当前meta.data中的orig.ident
orig_ids <- scRNA_harmony@meta.data$orig.ident

# 初始化一个空向量存放每个细胞对应的组别
group_labels <- rep(NA, length(orig_ids))

# 遍历每个组别，给对应的细胞赋值
for (group_name in names(group_mapping)) {
  cells_in_group <- group_mapping[[group_name]]
  group_labels[orig_ids %in% cells_in_group] <- group_name
}

# 将新分组信息加入meta.data
scRNA_harmony@meta.data$group <- group_labels

# 查看每个组别的细胞数
table(scRNA_harmony@meta.data$group)


library(Seurat)
library(Matrix)

# 读取数据
expression_matrix <- read.delim("E:/scdata/MSCtr_snRNA_ExpressionMatrix_R.txt", header = TRUE, row.names = 1)

# 转换为稀疏矩阵
expr_mat <- as(as.matrix(expression_matrix), "dgCMatrix")

# 创建 Seurat 对象
seurat_obj <- CreateSeuratObject(counts = expr_mat)

# 查看对象信息
print(seurat_obj)
head(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")

cell_count <- ncol(scRNA_harmony)
print(paste("总细胞数为：", cell_count))

rm(expression_matrix)

# 定义每个组别对应的颜色（请根据实际组别调整）
custom_colors <- c(
  "AD" = "#1f78b4",    # 蓝色
  "MS" = "#7DC69B",    # 绿色
  "PD" = "#F2A1A7",    # 红色
  "EPI" = "#D8EEFB",
  "NC" = "#ED2225"# 橙色
  # 如果有其他组别，也可以继续添加
)

# 绘制UMAP并指定颜色
DimPlot(scRNA_harmony, reduction = "umap", group.by = "group") +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +  # 使用简洁主题
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  ggtitle("UMAP Clustering by Group")



