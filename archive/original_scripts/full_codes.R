#Full codes 
#doading pachages
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(viridis)
library(ggplot2)

# 
DotPlot(scRNAharmony_9.5, features = c("CLDN5","PECAM1","TMEM119","P2RY12","MAP2","GAD1","SLC1A3","GFAP","MOG","MBP","SLC17A7","PRMT8"))
p <- c("SLC17A7","SLC17A6","SLC17A8","GAD1","GAD2","SLC32A1")
DotPlot(scRNAharmony_9.5,features = p)
DimPlot(scRNAharmony_9.5)

DimPlot(scRNAharmony_9.5, reduction = "umap", 
        cols = viridis(8)) +
  theme(legend.position = "right")

#细胞注释
bfreaname.pbmc <- scRNAharmony_9.5
new.cluster.ids <- c("Excitatory.neurons", "Oligodendrocytes","Inhibitory.neurons", "Astrocytes", "Inhibitory.neurons", "Microglia",  
                     "Inhibitory.neurons", "Endothelial.cells")
names(new.cluster.ids) <- levels(scRNAharmony_9.5) 
scRNAharmony_9.5 <- RenameIdents(scRNAharmony_9.5, new.cluster.ids)
scRNAharmony_9.5$celltype <- Idents(scRNAharmony_9.5)
DimPlot(scRNAharmony_9.5)
head(scRNAharmony_9.5)
save.image(file = scRNAharmony_9.5)
#绘制总的umap图
color <- c("#919ac2","#ffac98","#70a4c8","#a5a9af","#63917d","#dbd1b4","#6e729a","#9ba4bd","#c5ae5f","#b9b8d6")
library(ggh4x)
p <- ggplot(scRNAharmony_9.5, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
  geom_point(size = 0.03,shape = 16, stroke = 0) +
  scale_color_manual(values = color) +
  facet_grid(~label) +
  theme_classic() + # 使用简洁主题
  theme( plot.background = element_blank(), # 移除背景
         panel.grid.major = element_blank(), # 移除主要网格线
         panel.grid.minor = element_blank(), # 移除次要网格线
         margin = margin(), # 增加顶部边距
         axis.title.x = element_blank(), # 移除x轴标题
         axis.title.y = element_blank(), # 移除y轴标题
         axis.text = element_blank(), # 移除坐标轴刻度标签
         axis.ticks = element_blank(), # 移除坐标轴刻度线
         axis.line = element_line(colour ="black", size = 0.3,
                                  arrow = arrow(length = unit(0.1,"cm"))), # 加粗坐标轴并添加箭头
         strip.background = element_rect(fill ='#e6bac5',color=NA),
         strip.placement ='outside',
         strip.text = element_text(size = 8),
         legend.position ="none",
         aspect.ratio = 1) +
  scale_x_continuous(limits = c(min(df$UMAP_1), max(df$UMAP_1)) ) + # 设置x轴范围
  scale_y_continuous(limits = c(min(df$UMAP_2), max(df$UMAP_2)) ) # 设置y轴范围
p

ggsave(filename ="BCells_Umap.png",height = 2,width = 2,plot = p,dpi = 300,bg="white")

DotPlot(scRNAharmony_9.5, features = c("CLDN5","PECAM1","ABCB1",
                                       "P2RY12","DOCK8","CSF1R",
                                       "SLC1A3","GFAP","AQP4",
                                       "GAD1","GAD2","PVALB",
                                       "MOG","MBP","OLIG1",
                                       "GRIN1","CAMK2A","CBLN2"))



DotPlot(scRNA_seurat, features = c(
  "CLDN5","PECAM1","ABCB1",
  "P2RY12","DOCK8","CSF1R",
  "GFAP","AQP4",
  "MOG","MBP","OLIG1",
  "GAD1","GAD2","PVALB",
  "GRIN1","CAMK2A","CBLN2"
), dot.scale=8) +
  scale_color_gradient(low = "lightblue", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    axis.text.y = element_text(size=10),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=10),
    plot.title = element_text(hjust=0.5, size=14, face="bold")
  ) +
  labs(
    title = "Gene Expression DotPlot",
    color = "Expression Level"
  )
# 查看当前分组类别
levels(Idents(scRNA_seurat))
# 设置分组的顺序
Idents(scRNA_seurat) <- factor(Idents(scRNA_seurat), levels = c("Excitatory.neurons", "Inhibitory.neurons", "Oligodendrocytes","Astrocytes","Microglia","Endothelial.cells"))

table(scRNA_seurat$group)
table(scRNA_seurat$group, useNA = "ifany")
# 查看NA细胞的orig.ident
na_cells_orig <- scRNA_seurat@meta.data$orig.ident[is.na(scRNA_seurat$group)]

# 根据orig.ident中的内容，重新赋值group
# 示例：
scRNA_seurat$group[is.na(scRNA_seurat$group) & grepl("PD", scRNA_seurat@meta.data$orig.ident)] <- "PD"
scRNA_seurat$group[is.na(scRNA_seurat$group) & grepl("NC", scRNA_seurat@meta.data$orig.ident)] <- "NC"

# 加载包
library(ggplot2)
library(dplyr)
library(scales)
library(plotrix)

# 定义数据
groups <- c("AD", "EPI", "MS", "NC", "PD")
counts <- c(109037, 34083, 11208, 100777, 31919)

# 创建数据框
data <- data.frame(
  group = groups,
  count = counts
)

# 计算比例
data <- data %>%
  mutate(
    fraction = count / sum(count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n = -1))
  )

# 定义颜色
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# 绘制环形比例图
ggplot(data, aes(ymin = ymin, ymax = ymax, xmin = 3, xmax = 4, fill = group)) +
  geom_rect(color = "white") +  # 添加白色边框
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +  # 控制环形厚度和中心空白
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "lines"),
    plot.title = element_text(hjust = 0.5, size=18)
  ) +
  scale_fill_manual(values = colors) +
  ggtitle("Sample Distribution by Group")



#计算细胞比例
library(dplyr)
# 载入必要的包
library(dplyr)
library(ggplot2)

# 提取元数据
meta_data <- scRNA_seurat@meta.data

# 计算每个细胞亚型的数量
celltype_counts <- meta_data %>%
  group_by(celltype) %>%
  summarise(count = n()) %>%
  arrange(desc(count))  # 按数量降序排列

# 创建渐变色标
color_palette <- colorRampPalette(c("#6BAED6", "#08306B"))(nrow(celltype_counts))

# 绘制精美柱状图
ggplot(celltype_counts, aes(x = reorder(celltype, -count), y = count, 
                            fill = reorder(celltype, -count))) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
  geom_text(aes(label = count), 
            vjust = -0.5, 
            size = 3.5,
            fontface = "bold",
            color = "black") +
  scale_fill_manual(values = color_palette) +
  labs(title = "Single-cell RNA Sequencing Cell Type Distribution",
       subtitle = "Cell counts by cell type",
       x = "Cell Type",
       y = "Cell Count",
       caption = "Data source: scRNA-seq experiment") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40", margin = margin(b = 20)),
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 1, margin = margin(t = 15)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # Expand y-axis for label display



#9.08_2
library(Seurat)
library(ggplot2)
library(ggh4x)

# 1. 提取 UMAP 坐标
df <- as.data.frame(Embeddings(scRNAharmony_9.5, reduction = "umap"))

# 2. 合并 meta.data
df <- cbind(df, scRNAharmony_9.5@meta.data)

# 查看前几行
head(df)

# 3. 自定义颜色
color <- c("#919ac2","#ffac98","#70a4c8","#a5a9af",
           "#63917d","#dbd1b4","#6e729a","#9ba4bd",
           "#c5ae5f","#b9b8d6")

# 改列名
colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.03, shape = 16, stroke = 0) +
  scale_color_manual(values = color) +
  facet_grid(~group) +
  theme_classic() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0), # ✅ 正确写法
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3,
                             arrow = arrow(length = unit(0.1, "cm"))),
    strip.background = element_rect(fill = '#e6bac5', color = NA),
    strip.placement = 'outside',
    strip.text = element_text(size = 8),
    legend.position = "none",
    aspect.ratio = 1
  ) +
  scale_x_continuous(limits = c(min(df$UMAP_1), max(df$UMAP_1))) +
  scale_y_continuous(limits = c(min(df$UMAP_2), max(df$UMAP_2)))

# 保存图片
ggsave("BCells_Umap.png", height = 2, width = 2, plot = p, dpi = 300, bg = "white")




#绘制总的umap
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.03, shape = 16, stroke = 0) +
  scale_color_manual(values = color) +
  theme_classic() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.3,
                             arrow = arrow(length = unit(0.1, "cm"))),
    strip.background = element_rect(fill = '#e6bac5', color = NA),
    strip.placement = 'outside',
    strip.text = element_text(size = 8),
    legend.position = "none",  # 如果想看图例可以改成 "right"
    aspect.ratio = 1
  ) +
  scale_x_continuous(limits = c(min(df$UMAP_1), max(df$UMAP_1))) +
  scale_y_continuous(limits = c(min(df$UMAP_2), max(df$UMAP_2)))

# 保存
ggsave("UMAP_total.png", height = 2, width = 2, plot = p, dpi = 300, bg = "white")
ggsave("UMAP_total.png", )

saveRDS(inhibitory_cells_qc, file = "inhibitory_cells_qc.rds")


#提取抑制性神经元细胞亚群
head(scRNAharmony_9.5)

# 1. 提取抑制性神经元亚群
inhibitory_cells <- subset(scRNAharmony_9.5, subset = celltype == "Inhibitory.neurons")
inhibitory_cells <- NormalizeData(inhibitory_cells)
inhibitory_cells <- FindVariableFeatures(inhibitory_cells)
inhibitory_cells <- ScaleData(inhibitory_cells)
inhibitory_cells <- RunPCA(inhibitory_cells, npcs = 20)
dims_to_use <- 1:10
inhibitory_cells <- FindNeighbors(inhibitory_cells, dims = dims_to_use)
inhibitory_cells <- FindClusters(inhibitory_cells, resolution = 0.05)
inhibitory_cells <- RunUMAP(inhibitory_cells, dims = dims_to_use)
DimPlot(inhibitory_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Re-clustered Inhibitory Neurons")

saveRDS(inhibitory_cells, file = "Inhibitory_Neurons_Reclustered.rds")

color2 <- c("")
p <- ggplot(inhibitory_cells, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
  geom_point(size = 0.03,shape = 16, stroke = 0) +
  scale_color_manual(values = color) +
  facet_grid(~label) +
  theme_classic() + # 使用简洁主题
  theme( plot.background = element_blank(), # 移除背景
         panel.grid.major = element_blank(), # 移除主要网格线
         panel.grid.minor = element_blank(), # 移除次要网格线
         margin = margin(), # 增加顶部边距
         axis.title.x = element_blank(), # 移除x轴标题
         axis.title.y = element_blank(), # 移除y轴标题
         axis.text = element_blank(), # 移除坐标轴刻度标签
         axis.ticks = element_blank(), # 移除坐标轴刻度线
         axis.line = element_line(colour ="black", size = 0.3,
                                  arrow = arrow(length = unit(0.1,"cm"))), # 加粗坐标轴并添加箭头
         strip.background = element_rect(fill ='#e6bac5',color=NA),
         strip.placement ='outside',
         strip.text = element_text(size = 8),
         legend.position ="none",
         aspect.ratio = 1) +
  scale_x_continuous(limits = c(min(df$UMAP_1), max(df$UMAP_1)) ) + # 设置x轴范围
  scale_y_continuous(limits = c(min(df$UMAP_2), max(df$UMAP_2)) ) # 设置y轴范围
p


markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

# 提取UMAP坐标和簇信息
umap_df <- Embeddings(inhibitory_cells, reduction = "umap") %>% 
  as.data.frame()
umap_df$cluster <- inhibitory_cells@active.ident
# 如果有label列，也加入
# umap_df$label <- inhibitory_cells$your_label_column
color2 <- c("#F7A6AC","#F7B7D2","#70a4c8","#EEF0A7","#B2DBB9","#B8E5FA","#6e729a","#9ba4bd","#c5ae5f","#b9b8d6")
p <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.03, shape = 16, stroke = 0) +
  scale_color_manual(values = color2) +
  theme_classic() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # margin = margin(),  # 删除这行或改为 margin(t=5, r=5, b=5, l=5)
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour ="black", size = 0.3,
                             arrow = arrow(length = unit(0.1,"cm"))),
    strip.background = element_rect(fill ='#e6bac5', color=NA),
    strip.placement ='outside',
    strip.text = element_text(size = 8),
    legend.position ="none",
    aspect.ratio = 1
  ) +
  coord_cartesian(xlim = range(umap_df$umap_1), ylim = range(umap_df$umap_2))
print(p)
DimPlot(inhibitory_cells)
DotPlot(inhibitory_cells, features = c("PRMT8", "PVALB", "SST", "VIP", "CR", "NPY", "GAD1", "GAD2"))
rm(dims_to_use)

inhibitory_markers <- FindAllMarkers(object = inhibitory_cells, test.use="wilcox" ,
                                     only.pos = TRUE,
                                     logfc.threshold = 0.25)
inhibitory_markers_all =inhibitory_markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
top100 = inhibitory_markers_all %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
View(top100)
write.csv(top100,"inhibitory_markers_all4.csv",row.names = T)


PV <- c("PVALB","LHX6","ERBB4","GAD1","GAD2")
SST <- c("SST","NPY","TAC1","LHX6")
VIP <- c("VIP","CCK","PROX1")
CR <- c("CALB2","NPY","CCK")
# 定义基因列表
PV_genes <- c("PVALB","LHX6","ERBB4","GAD1","GAD2","LAMP5","NRG1","PAX6")
SST_genes <- c("SST","NPY","TAC1","LHX6")
VIP_genes <- c("VIP","CCK","PROX1")
CR_genes  <- c("CALB2","NPY","CCK","PRMT8")

# 合并所有基因，确保显示所有标志物
all_genes <- unique(c(PV_genes, SST_genes, VIP_genes, CR_genes))

# 使用DotPlot
DotPlot(object = inhibitory_cells, features = all_genes) +
  # 增强颜色对比度
  scale_color_gradient(low = "lightblue", high = "red") +
  # 改善标签
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=10),
    plot.title = element_text(hjust=0.5, size=14, face="bold")
  ) +
  labs(
    title = "抑制性神经元亚型标志物表达",
    color = "表达水平"
  )


p <- DimPlot(inhibitory_cells, group.by = "RNA_snn_res.0.05", label = TRUE)
head(inhibitory_cells)

FeaturePlot(inhibitory_cells, features = "PRMT8", pt.size = 0.8,min.cutoff = 0, max.cutoff = 2) +
  scale_color_gradientn(colors = c("#EEF0A7", "#5086C4", "#882255"))



# 1. 计算线粒体基因比例
# 假设线粒体基因以"^MT-"开头
inhibitory_cells[["percent.mt"]] <- PercentageFeatureSet(inhibitory_cells, pattern = "^MT-")

# 2. （可选）如果有血红蛋白基因，也可以计算
# inhibitory_cells[["percent.HB"]] <- PercentageFeatureSet(inhibitory_cells, pattern = "^HB")

# 3. 根据筛选条件子集
inhibitory_cells_qc <- subset(inhibitory_cells,
                              subset = 
                                nFeature_RNA > 300 & nFeature_RNA < 5000 &  # 基因数筛选
                                percent.mt < 10 &                        # 线粒体比例筛选
                                # 如果有血红蛋白比例，也可以加入筛选
                                # percent.HB < 3 &
                                nCount_RNA < quantile(inhibitory_cells$nCount_RNA, 0.97) &  # 总UMI数筛选
                                nCount_RNA > 1000)

# 这样，你得到的 `inhibitory_cells_qc` 就是经过质控过滤后的细胞子集


DimPlot(inhibitory_cells_qc)
table(inhibitory_cells_qc$group)


# 1. 提取PRMT8的表达矩阵（以归一化值为例）
expr_matrix <- GetAssayData(object = inhibitory_cells_qc, assay = "RNA", slot = "data")
# 2. 提取PRMT8的表达值
prmt8_expr <- expr_matrix["PRMT8", ]  # 如果基因名不同，请调整

# 3. 将表达值加入meta.data
inhibitory_cells_qc$PRMT8_expr <- as.numeric(prmt8_expr)

# 4. 计算每个组的平均表达
library(dplyr)
expr_summary <- inhibitory_cells_qc@meta.data %>%
  group_by(group) %>%
  summarise(mean_PRMT8 = mean(PRMT8_expr, na.rm = TRUE),
            median_PRMT8 = median(PRMT8_expr, na.rm = TRUE),
            cell_count = n())

print(expr_summary)
library(ggplot2)
ggplot(inhibitory_cells_qc@meta.data, aes(x = group, y = PRMT8_expr)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("Group") +
  ylab("PRMT8 Expression")

table(inhibitory_cells_qc$group, useNA = "ifany")
na_cells_orig <- inhibitory_cells_qc@meta.data$orig.ident[is.na(inhibitory_cells_qc$group)]
table(na_cells_orig)

# 查看NA细胞的orig.ident
na_cells_orig <- inhibitory_cells_qc@meta.data$orig.ident[is.na(inhibitory_cells_qc$group)]

# 根据orig.ident中的内容，重新赋值group
# 示例：
inhibitory_cells_qc$group[is.na(inhibitory_cells_qc$group) & grepl("PD", inhibitory_cells_qc@meta.data$orig.ident)] <- "PD"
inhibitory_cells_qc$group[is.na(inhibitory_cells_qc$group) & grepl("NC", inhibitory_cells_qc@meta.data$orig.ident)] <- "NC"

DimPlot(inhibitory_cells_qc)



expr_data <- GetAssayData(inhibitory_cells_qc, assay="RNA", slot="data")  # 或 "counts" 根据你的数据类型


genes_of_interest <- c("ADNP", "DGCR2", "FAM19A1", "GABRA3", "GABRA5", "GABRB1", 
                       "GABRG1", "GABRG3", "GRIN2C", "GRIN2D", "ITGB3", "KCNMA1",
                       "KCNN3", "MECP2", "MED13L", "MEF2", "NEUROCAN", "NEUROG2",
                       "NR4A1", "PPP1R9A", "RHOA", "ROCK2", "RTN4", "S6K1",
                       "SAPAP3", "SAPAP4", "STX1A", "TBR1","PRMT8")

# 检查基因是否在表达矩阵中
genes_present <- intersect(genes_of_interest, rownames(expr_data))
genes_missing <- setdiff(genes_of_interest, rownames(expr_data))
if(length(genes_missing) > 0){
  message("以下基因在表达矩阵中未找到：", paste(genes_missing, collapse=", "))
}


library(dplyr)

# 构建一个数据框，列出细胞信息
meta_df <- inhibitory_cells_qc@meta.data

# 提取表达矩阵对应的基因
expr_sub <- expr_data[genes_present, ]

# 转换为数据框，行是基因，列是细胞
expr_df <- as.data.frame(t(expr_sub))
expr_df$cell <- rownames(expr_df)

# 结合meta数据
full_df <- left_join(meta_df, expr_df, by = c("cell" = "cell"))

# 计算每个组中每个基因的平均表达
avg_expr <- full_df %>%
  group_by(group) %>%
  summarise(across(all_of(genes_present), mean, na.rm=TRUE))

# 转换为矩阵，用于绘图
heatmap_mat <- as.data.frame(avg_expr)
rownames(heatmap_mat) <- heatmap_mat$group
heatmap_mat$group <- NULL
heatmap_matrix <- as.matrix(heatmap_mat)

