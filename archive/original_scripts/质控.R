DimPlot(scRNA_seurat)
head(scRNA_seurat)
grep("^MT", rownames(scRNA_seurat), value = TRUE)
scRNA_seurat$mt_percent <- NULL
scRNA_seurat$percent.hb <- PercentageFeatureSet(scRNA_seurat, features = c("HBB", "HBA2", "HBA1"))
VlnPlot(scRNA_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.hb"), ncol = 3)
my_colors <- brewer.pal(8, "Set2") 
# plot
VlnPlot(scRNA_seurat, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.hb"), 
        ncol = 3, 
        pt.size = 0,         # 先将 pt.size 设为 0，避免点太拥挤
        cols = my_colors,    # 使用自定义配色
        combine = TRUE) &    # 使用 & 符号对所有子图应用相同的 ggplot 主题
  geom_jitter(size = 0.1, show.legend = FALSE, alpha = 0.2) & # 重新叠加透明的小点，增加质感
  theme_bw() &               # 使用白色背景主题
  theme(
    axis.title.x = element_blank(),        # 隐藏 X 轴标题
    axis.text.x = element_text(angle = 45, hjust = 1), # X 轴标签旋转，防止重叠
    plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中加粗
    legend.position = "none"               # 隐藏右侧冗余图例
  )


head(scRNA_seurat)
