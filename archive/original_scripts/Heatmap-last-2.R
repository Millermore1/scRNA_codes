library(Seurat)
library(ggplot2)
DimPlot(scRNA_seurat)
head(scRNA_seurat)
table(scRNA_seurat$group)

FeaturePlot(scRNA_seurat, features = "PRMT8")
#PRMT8 expression in all cells
FeaturePlot(
  scRNA_seurat,
  features = "PRMT8",
  min.cutoff = "q1",
  max.cutoff = "q20",
  pt.size = 1.2,           
  order = TRUE,
  cols = c("#E0E0E0", "#3A7CA5", "#0D3B66"),
  raster = TRUE,              
  raster.dpi = c(600, 600)   
) + 
  theme_minimal() +            
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  labs(title = "PRMT8 Expression",
       subtitle = "Single-cell RNA Sequencing") +
  guides(color = guide_colorbar(
    title = "Expression\nLevel",
    barwidth = 0.8,
    barheight = 10,
    frame.colour = "black",
    ticks.colour = "black"
  ))


#extract all neurons transcriptome data 

#heatmap of different gene sets
expr_matrix <- GetAssayData(scRNA_seurat, assay = "RNA", slot = "data")

genes1 <- c(
  # 离子通道相关
  "CACNA1A", "CACNA1B", "CACNA1E",
  "HCN1", "HCN2", "HCN3", "HCN4",
  "SCN1A", "SCN1B", "SCN8A"
)

# 确保基因在表达矩阵中
genes_present <- intersect(genes1, rownames(expr_matrix))
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


rm(scRNA_seurat)













