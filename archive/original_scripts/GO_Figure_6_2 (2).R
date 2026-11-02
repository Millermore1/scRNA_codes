library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
Idents(scRNA_seurat)
head(scRNA_seurat)
table(scRNA_seurat$celltype)

neurons <- subset(
  scRNA_seurat,
  subset = celltype %in% c("Excitatory.neurons", "Inhibitory.neurons")
)

DefaultAssay(neurons) <- "RNA"

neurons <- NormalizeData(neurons)
neurons <- FindVariableFeatures(neurons, nfeatures = 3000)
neurons <- ScaleData(neurons)
neurons <- RunPCA(neurons)

ElbowPlot(neurons, ndims = 50)

neurons <- FindNeighbors(neurons, dims = 1:30)
neurons <- FindClusters(neurons, resolution = 0.03)
neurons <- RunUMAP(neurons, dims = 1:30)

DimPlot(neurons, group.by = "seurat_clusters", label = TRUE)
DimPlot(neurons, group.by = "celltype", label = TRUE)

table(neurons$seurat_clusters, neurons$celltype)

neuron_markers <- FindAllMarkers(
  neurons,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)



write.csv(neuron_markers, "Neuron_subcluster_markers.csv", row.names = FALSE)
saveRDS(neurons, "Excitatory_Inhibitory_neurons_reclustered.rds")


library(Seurat)
library(ggplot2)

# 1. 定义 Marker 列表
PN_markers  <- c("SLC17A7", "CAMK2A", "SATB2", "CUX1", "CUX2", "RORB", "BCL11B", "TLE4", "FEZF2", "NEUROD6")
PV_markers  <- c("PVALB", "KCNC1", "GAD1", "GAD2", "SLC6A1")
SST_markers <- c("SST", "NPY", "NOS1", "GAD1", "GAD2", "SLC6A1")
VIP_markers <- c("VIP", "CALB2", "CRH", "HTR3A", "GAD1", "GAD2", "SLC6A1")

# 2. 合并并去重
neuron_markers <- unique(c(PN_markers, PV_markers, SST_markers, VIP_markers))

# 3. 检查并筛选对象中实际存在的基因
neuron_markers_use <- neuron_markers[neuron_markers %in% rownames(Excitatory_Inhibitory_neurons_reclustered)]

# 打印数据集中缺失的 Marker 基因
missing_genes <- setdiff(neuron_markers, neuron_markers_use)
if (length(missing_genes) > 0) {
  cat("未在对象中找到的基因:", paste(missing_genes, collapse = ", "), "\n")
}

# 4. 设置分群并作图
Idents(Excitatory_Inhibitory_neurons_reclustered) <- "seurat_clusters"

p1 <- DotPlot(
  Excitatory_Inhibitory_neurons_reclustered,
  features = neuron_markers_use,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  labs(
    title = "Neuron Subtype Marker DotPlot",
    x = "Marker Genes",
    y = "Neuron Clusters"
  )

p1

DimPlot(Excitatory_Inhibitory_neurons_reclustered,label = TRUE)

new.cluster.ids <- c("PN", "SST", "SST", "VIP", 
                     "PN", "PV", "PN", 
                     "PV", "PN", "PV", "SST","PN")

names(new.cluster.ids) <- levels(Excitatory_Inhibitory_neurons_reclustered)
All_E_I_neurons <- RenameIdents(Excitatory_Inhibitory_neurons_reclustered, new.cluster.ids)
DimPlot(All_E_I_neurons)
rm(All_E_I_neurons)







DimPlot(All_E_I_neurons_annotated)
head(All_E_I_neurons_annotated)
Idents(All_E_I_neurons_annotated) <- "RNA_snn_res.0.03"
DotPlot(All_E_I_neurons_annotated,
        features = "PRMT8")

# Seurat 中注释 PV 神经元的推荐 marker 组合
pv_markers <- c(
  "PVALB",           # 核心定义
  "GAD1", "GAD2",    # GABAergic
  "SLC32A1",         # VGAT
  "KCNC1", "KCNC2",  # fast-spiking
  "SOX6", "LHX6",    # MGE lineage
  "NXPH1", "NRIP3",  # PV enriched
  "PTHLH",           # chandelier
  "FGF13"            # chandelier
)

# DotPlot 可视化
DotPlot(All_E_I_neurons_annotated, features = pv_markers) + RotatedAxis()



neurons <- subset(
  scRNA_seurat,
  subset = celltype %in% c("Excitatory.neurons")
)

DefaultAssay(neurons) <- "RNA"

neurons <- NormalizeData(neurons)
neurons <- FindVariableFeatures(neurons, nfeatures = 3000)
neurons <- ScaleData(neurons)
neurons <- RunPCA(neurons)

ElbowPlot(neurons, ndims = 50)

neurons <- FindNeighbors(neurons, dims = 1:30)
neurons <- FindClusters(neurons, resolution = 0.03)
neurons <- RunUMAP(neurons, dims = 1:30)

DimPlot(neurons, group.by = "seurat_clusters", label = TRUE)
DimPlot(neurons, group.by = "celltype", label = TRUE)

table(neurons$seurat_clusters, neurons$celltype)

neuron_markers <- FindAllMarkers(
  neurons,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)


library(Seurat)
library(ggplot2)
library(dplyr)

# ========== 1. 删除 cluster 3 ==========
E_neurons_PRMT8 <- subset(E_neurons_PRMT8, subset = RNA_snn_res.0.03 != "3")

# ========== 2. 重新标准化与聚类 ==========
E_neurons_PRMT8 <- NormalizeData(E_neurons_PRMT8)
E_neurons_PRMT8 <- FindVariableFeatures(E_neurons_PRMT8, selection.method = "vst", nfeatures = 2000)
E_neurons_PRMT8 <- ScaleData(E_neurons_PRMT8)
E_neurons_PRMT8 <- RunPCA(E_neurons_PRMT8, features = VariableFeatures(object = E_neurons_PRMT8))
E_neurons_PRMT8 <- FindNeighbors(E_neurons_PRMT8, dims = 1:30)
E_neurons_PRMT8 <- FindClusters(E_neurons_PRMT8, resolution = 0.03)
E_neurons_PRMT8 <- RunUMAP(E_neurons_PRMT8, dims = 1:30)


# 2. 全部 marker 基因
markers <- c(
  "PRMT8",
  "SLC17A7", "SATB2", "SLC30A3",
  "CUX2", "LHX2", "OTOF",
  "RORB", "RSPO1", "SCNN1A", "RSPO2", "DEPTOR",
  "OSR1", "RXFP1",
  "CAR3",
  "FEZF2", "BCL6", "FAM84B", "POU3F1", "BCL11B",
  "TLE4", "FOXP2",
  "LY6G6E",
  "CPLX3", "CTGF"
)

# 3. 按指定 9 个基因的综合表达量排序聚类（高的在左）
sort_genes <- c("SATB2", "SLC17A7", "BCL11B", "TLE4", "FOXP2", "RXFP1", "RORB", "CUX2", "PRMT8")

p_sort <- DotPlot(E_neurons_PRMT8, features = sort_genes, group.by = "RNA_snn_res.0.03", dot.min = 0)

cluster_order <- p_sort$data %>%
  filter(features.plot %in% sort_genes) %>%
  group_by(id) %>%
  summarise(score = sum(avg.exp.scaled, na.rm = TRUE)) %>%
  arrange(desc(score)) %>%
  pull(id)

E_neurons_PRMT8$RNA_snn_res.0.03 <- factor(
  as.character(E_neurons_PRMT8$RNA_snn_res.0.03),
  levels = cluster_order
)

# 4. 绘制（点更大、颜色更深）
DotPlot(
  object = E_neurons_PRMT8,
  features = sort_genes,
  group.by = "RNA_snn_res.0.03",
  cols = c("#E8E8E8", "darkred"),   # 浅灰 → 深红
  dot.scale = 8,                    # 点更大（默认6）
  dot.min = 0.05
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "italic", color = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Cortical Excitatory Neuron Markers")

head(E_neurons_PRMT8)
table(E_neurons_PRMT8$RNA_snn_res.0.03)
Idents(E_neurons_PRMT8)
library(Seurat)
library(dplyr)
library(openxlsx)

# 1. 计算差异表达基因（每个 cluster vs 其余）
Idents(E_neurons_PRMT8) <- "RNA_snn_res.0.03"
all_markers <- FindAllMarkers(
  object = E_neurons_PRMT8,
  only.pos = FALSE,        # TRUE = 只保留该cluster高表达基因
  min.pct = 0.25,          # 基因必须在≥25%细胞中表达
  logfc.threshold = 0.25,  # log2FC cutoff
  test.use = "wilcox"      # Wilcoxon秩和检验
)

# 添加显著性标记
all_markers <- all_markers %>%
  mutate(significance = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up",
    p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
    TRUE ~ "NS"
  ))

# 2. 导出为 Excel（每个 cluster 一个 sheet）
wb <- createWorkbook()

clusters <- sort(unique(all_markers$cluster))

for (cl in clusters) {
  sheet_name <- paste0("Cluster_", cl)
  addWorksheet(wb, sheetName = sheet_name)
  
  df <- all_markers %>% 
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC))
  
  writeData(wb, sheet = sheet_name, x = df, rowNames = FALSE)
}

saveWorkbook(wb, file = "PRMT8_DEG_all_clusters.xlsx", overwrite = TRUE)

# 3. 特定 cluster 两两比较（示例：cluster 4 vs cluster 6）
# deg_4vs6 <- FindMarkers(
#   object = E_neurons_PRMT8,
#   ident.1 = "4",
#   ident.2 = "6",
#   min.pct = 0.25,
#   logfc.threshold = 0.25
# )
#write.xlsx(deg_4vs6, "PRMT8_cluster4_vs_6.xlsx", rowNames = TRUE)

getwd()


Goenrich <- read.csv("C:/Users/黄廉/Desktop/Ex_GO_enrichment.csv")
head(Goenrich)

library(ggplot2)
library(dplyr)
library(ggtext)

# 1. 按 Enrichment 升序排列并更新因子顺序（确保最高值位于图表最顶端）
Goenrich <- Goenrich %>%
  arrange(Enrichment) %>%
  mutate(Description = factor(Description, levels = Description))

# 2. 绘制按 Enrichment 排序后的富集条形图
ggplot(Goenrich, aes(x = Enrichment, y = Description, fill = Enrichment)) +
  geom_col(width = 0.85) +
  # 将 Description 嵌入条形图内部左侧
  geom_text(aes(label = Description, x = max(Enrichment) * 0.02),
            hjust = 0, 
            fontface = "bold.italic", 
            size = 4, 
            color = "#1A1A1A") +
  # 渐变配色
  scale_fill_gradient(low = "#C4E2E3", high = "#187C7C",
                      name = "Enrichment") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  # 标题与坐标轴
  labs(
    title = "GO Enrichment Analysis",
    y = "Biological Processes",
    x = "**Fold Enrichment**<br><br>*(Prmt8+)*<br>**excitatory neuron's biological function**"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.05, margin = margin(b = 10)),
    
    # 隐藏原生 Y 轴文本
    axis.text.y = element_blank(),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 10)),
    axis.ticks.y = element_blank(),
    
    # 解析 Markdown 格式的 X 轴标题
    axis.title.x = element_markdown(size = 13, hjust = 0.4, lineheight = 1.2, margin = margin(t = 5)),
    axis.text.x = element_text(size = 10, face = "bold", color = "#333333"),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11, margin = margin(b = 8)),
    legend.text = element_text(size = 10, face = "bold")
  ) +
  guides(fill = guide_colorbar(barwidth = 1.2, barheight = 5, ticks.linewidth = 1))



library(Seurat)
library(dplyr)

# 1. 过滤提取：剔除 RNA_snn_res.0.03 为 "3" 的细胞
E_neurons_sub <- subset(E_neurons_PRMT8, subset = RNA_snn_res.0.03 != "3")

# 验证提取结果（确认 Cluster 3 已被移除）
table(E_neurons_sub$RNA_snn_res.0.03)

# 2. 重新聚类流程 (Re-clustering)
# 剔除特定细胞群后，数据的基因变异度与主成分分布会发生变化，需重新寻找高变基因并降维

# 2.1 重新选择高变基因与标准缩放
E_neurons_sub <- FindVariableFeatures(E_neurons_sub, selection.method = "vst", nfeatures = 2000)
E_neurons_sub <- ScaleData(E_neurons_sub)

# 2.2 重新计算 PCA
E_neurons_sub <- RunPCA(E_neurons_sub, npcs = 30, verbose = FALSE)

# 【注意】：如果你的多队列数据之前使用了 Harmony 进行批次校正，需在此处重新运行 Harmony：
# E_neurons_sub <- RunHarmony(E_neurons_sub, group.by.vars = "dataset_id") 

# 2.3 重新构建近邻图与分群（可同时测试多个 resolution）
# dims 请根据 ElbowPlot(E_neurons_sub) 评估选择（通常 1:15 到 1:20）
E_neurons_sub <- FindNeighbors(E_neurons_sub, dims = 1:20)
E_neurons_sub <- FindClusters(E_neurons_sub, resolution = c(0.01, 0.03, 0.05, 0.1))

# 2.4 重新计算 UMAP 坐标
E_neurons_sub <- RunUMAP(E_neurons_sub, dims = 1:20)

# 3. 可视化检查新的 UMAP 结果
DimPlot(E_neurons_sub, group.by = "RNA_snn_res.0.03", label = TRUE, reduction = "umap")

library(Seurat)
library(dplyr)
library(ggplot2)

# 1. 定义需展示的兴奋性神经元 Marker 基因
sort_genes <- c("SATB2","BCL11B", "TLE4", "FOXP2","NEUROD6", "RXFP1", "RORB", "CUX2", "PRMT8")

# 2. 检查 Marker 基因是否全在重聚类数据的基因列表中
valid_genes <- intersect(sort_genes, rownames(E_neurons_sub))
if(length(valid_genes) < length(sort_genes)){
  message("注意：部分基因未在数据集中检测到: ", setdiff(sort_genes, valid_genes))
}

# 3. 计算各 Cluster 的表达得分并排序
p_sort <- DotPlot(
  object = E_neurons_sub, 
  features = valid_genes, 
  group.by = "RNA_snn_res.0.03", 
  dot.min = 0
)

cluster_order <- p_sort$data %>%
  filter(features.plot %in% valid_genes) %>%
  group_by(id) %>%
  summarise(score = sum(avg.exp.scaled, na.rm = TRUE)) %>%
  arrange(desc(score)) %>%
  pull(id)

# 4. 重新设定 Cluster 的 Factor Levels 顺序
E_neurons_sub$RNA_snn_res.0.03 <- factor(
  as.character(E_neurons_sub$RNA_snn_res.0.03),
  levels = cluster_order
)

# 5. 绘制美化后的 DotPlot
DotPlot(
  object = E_neurons_sub,
  features = valid_genes,
  group.by = "RNA_snn_res.0.03",
  cols = c("#E8E8E8", "darkred"),   # 浅灰 → 深红
  dot.scale = 10,                    # 点的比例（放大展示）
  dot.min = 0.01                    # 仅显示表达细胞比例 > 5% 的点
) + 
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(title = "Cortical Excitatory Neuron Subpopulations (Re-clustered)")








