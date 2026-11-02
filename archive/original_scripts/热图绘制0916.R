#heatmap的绘制
# 提取表达矩阵（归一化表达值，通常用于比较）
expr_matrix <- GetAssayData(p, assay = "RNA", slot = "data")

# 确保基因在表达矩阵中
genes_present <- intersect(gene6, rownames(expr_matrix))
if(length(genes_present) == 0){
  stop("没有找到感兴趣的基因，请确认基因名是否正确。")
}

expr_subset <- expr_matrix[genes_present, ]

# 转置为样本为行，基因为列的数据框
expr_df <- as.data.frame(t(expr_subset))

# 添加分组信息（假设列名是 'group'）
expr_df$group <- p@meta.data$group

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



genes1 <- c(
  # 离子通道相关
  "CACNA1A", "CACNA1B", "CACNA1E",
  "HCN1", "HCN2", "HCN3", "HCN4",
  "SCN1A", "SCN1B", "SCN8A",
  "VAMP2", "SLC32A1"
)



# 神经递质受体
gene2 <- c(
  "GABRA1", "GABRA2", "GABRB2", "GABRB3", "GABRD", "GABRG2",
  "GABBR1", "GABBR2",
  "GRIN1", "GRIN2A", "GRIN2B", "GRIN3A", "GRIN3B",
  "GRIA1", "GRIA2", "GRIA3", "GRIA4",
  "GRIP1", "GABBR1",
  "NTRK2"
)


gene4 <- c(
  "CAMK2A", "CAMK2B", "MAPK3", "MTOR",
  "PIK3CA", "PRKACA", "PRKCA",
  "CYFIP2", "RAC1",
  "LIMK1",
  "HOMER1", "PUM2", "PVALB",
  "CDC42", "CPEB1", "CFL1", "CTTN", "DIAPH3", "MARK1", "MARK2"
)

gene5 <- c(
  "DLG4", "SHANK1", "SHANK2", "SHANK3",
  "SAPAP1", "SAPAP2", "SYNGAP1",
  "VAMP2", "SNAP25",
  "SYN1", "SYN2",
  "WAVE1", "WASL", "ARC"
)

gene6 <- c("NLGN2", "NLGN3", "NRGN1", "PRMT8","BDNF","GAD67","PVALB","MYH10","ARID1B","CHD8","UBE3A")






































