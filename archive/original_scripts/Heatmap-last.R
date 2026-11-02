library(Seurat)
library(ggplot2)
head(scRNA_seurat)
expr_matrix <- GetAssayData(scRNA_seurat, assay = "RNA", slot = "data")

genes_of_interest <- c(
  "CACNA1B", "CACNA1A", "CACNA1C", "CACNA1D", "CACNA2D1", 
  "HCN1", "HCN2", "HCN3", "HCN4"
)

# 确保基因在表达矩阵中
genes_present <- intersect(genes_of_interest, rownames(expr_matrix))
if(length(genes_present) == 0){
  stop("no interesting genes")
}

expr_subset <- expr_matrix[genes_present, ]

# 转置为样本为行，基因为列的数据框
expr_df <- as.data.frame(t(expr_subset))

# 添加分组信息（假设列名是 'group'）
expr_df$group <- scRNA_seurat@meta.data$group

# 计算每个组内这些基因的平均表达量
avg_expr_by_group <- expr_df %>%
  group_by(group) %>%
  summarise(across(all_of(genes_present), mean, na.rm = TRUE))

# 转换为矩阵，行是组，列是基因
heatmap_data <- as.data.frame(avg_expr_by_group)
rownames(heatmap_data) <- heatmap_data$group
heatmap_data$group <- NULL

# 转置，使得行是基因，列是组
heatmap_data_t <- t(heatmap_data)

# 进行Z-score标准化（对每个基因，即行进行标准化）
heatmap_data_scaled <- t(scale(t(heatmap_data_t)))

# 绘制热图
library(pheatmap)

pheatmap(
  heatmap_data_scaled,
  cluster_rows = TRUE,      # 可以选择聚类基因
  cluster_cols = FALSE,     # 不聚类组
  color = colorRampPalette(c("white", "#205EA7"))(20),
  display_numbers = FALSE,
  fontsize = 10,
  main = "Gene Expression in Different Groups",
  border_color = NA,
  cellwidth = 15,
  cellheight = 12
)

# 如果需要确保 "NC" 组在第一列，先获取所有列名
col_names <- colnames(heatmap_data_scaled)

# 将 "NC" 组移到第一列，其余的保持原顺序
if("NC" %in% col_names){
  col_names_new <- c("NC", setdiff(col_names, "NC"))
} else {
  stop("没有找到名为 'NC' 的组")
}

# 重新排列列
heatmap_data_scaled <- heatmap_data_scaled[, col_names_new]

# 重新绘制热图
pheatmap(
  heatmap_data_scaled,
  cluster_rows = TRUE,      # 可以选择聚类基因
  cluster_cols = FALSE,     # 不聚类组
  color = colorRampPalette(c("white", "#205EA7"))(20),
  display_numbers = FALSE,
  fontsize = 10,
  main = "Gene Expression in Different Groups",
  border_color = NA,
  cellwidth = 15,
  cellheight = 12
)

