library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)

head(scRNA_HKU)
excit_neurons <- subset(
  scRNA_HKU,
  idents = "Excitatory neurons"
)
table(scRNA_HKU$celltype)
excit_neurons <- NormalizeData(excit_neurons)
excit_neurons <- FindVariableFeatures(
  excit_neurons,
  selection.method = "vst",
  nfeatures = 2000
)

excit_neurons <- ScaleData(
  excit_neurons,
  features = VariableFeatures(excit_neurons)
)

excit_neurons <- RunPCA(
  excit_neurons,
  features = VariableFeatures(excit_neurons)
)

excit_neurons <- FindNeighbors(
  excit_neurons,
  dims = 1:20     
)

excit_neurons <- FindClusters(
  excit_neurons,
  resolution = 0.1  
)

# 7. 可视化（UMAP/tSNE）
excit_neurons <- RunUMAP(
  excit_neurons,
  dims = 1:20
)

DimPlot(
  excit_neurons,
  reduction = "umap",
  label = TRUE
)
saveRDS(excit_neurons, file = "excit_neurons.rds")

mycol <- c("#6D256B","#D8258D","#F599AF","#F58625","#0092BE","#0058A3","#263A83","#C6DEAF","#037F41","#9E86B5")


myc < c("#0092BE","#0058A3","#037F41","#263A83")

# 看核心 marker 的 DotPlot，方便一次比较所有 cluster
markers <- c("Slc17a7",      # excitatory
             "Cux1", "Cux2", "Calb1",      # L2/3 IT
             "Rorb",                       # L4/L5 IT
             "Bcl11b", "Fezf2",           # L5 ET/PT
             "Foxp1",                     # L6 IT
             "Tle4", "Foxp2", "Pcp4")     # L6 CT

DotPlot(excit_neurons, features = markers) + RotatedAxis()


markers_L5 <- c(
  "Rorb",          # L4/L5 IT
  "Bcl6", "Satb2", # L5 IT 相关
  "Bcl11b", "Fezf2", "Etv1", "Tbr1", "Crym"  # L5 ET/PT
)

DotPlot(excit_neurons, features = markers_L5) + RotatedAxis()



new.cluster.ids <- c("Exc_L2/3","Exc_L6","Exc_L4","Exc_L4","Exc_L5","Exc_L2/3","Exc_L5","Exc_L5","Exc_L5","Exc_L5","Exc_L6","Exc_L2/3","Exc_L2/3")

names(new.cluster.ids) <- levels(excit_neurons)
excit_neurons2 <- RenameIdents(excit_neurons, new.cluster.ids)
levels(excit_neurons2) #查看是否已改名
DimPlot(excit_neurons2, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()

# 将当前身份（Exc_L2/3, Exc_L4, Exc_L5, Exc_L6 等）写入 meta.data
excit_neurons2$subtype <- Idents(excit_neurons2)
saveRDS(excit_neurons2,file = "excit2.rds")
# 查看 meta.data 是否多了一列
head(excit_neurons2@meta.data)

table(excit_neurons2$group)

myc < c("#0092BE","#0058A3","#037F41","#263A83")
DimPlot(excit_neurons2, reduction = 'umap', 
        label = TRUE, pt.size = 0.5) + NoLegend()
# 自定义颜色（注意是 <- 不是 <）
myc <- c("#0092BE", "#0058A3", "#037F41", "#263A83")

p <- DimPlot(
  excit_neurons2,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.7,
  cols = myc        # 使用自定义颜色
) + 
  NoLegend() +
  theme_classic() +
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )

p

Idents(excit_neurons2) <- "group"
table(Idents(excit_neurons2))
deg_KOF_vs_WT <- FindMarkers(
  excit_neurons2,
  ident.1 = "KOF",
  ident.2 = "WT",
  logfc.threshold = 0.25,   # 可按需要调
  min.pct = 0.1,
  test.use = "wilcox"
)

DEGs_excit_KOFvsWT <- 
write.csv()
head(deg_KOF_vs_WT)

# 把基因名加成一列（方便导出）
deg_KOF_vs_WT$gene <- rownames(deg_KOF_vs_WT)

# 设置显著性筛选条件（可按需要调整）
deg_KOF_vs_WT_sig <- subset(
  deg_KOF_vs_WT,
  p_val_adj < 0.05 & abs(avg_log2FC) > 0.25
)

# 导出全部结果
write.csv(deg_KOF_vs_WT,
          file = "DEG_KOF_vs_WT_all.csv",
          row.names = FALSE)

# 导出筛选后的显著DEGs
write.csv(deg_KOF_vs_WT_sig,
          file = "DEG_KOF_vs_WT_excit_sig.csv",
          row.names = FALSE)
head(excit_neurons2)

library(dplyr)
library(ggplot2)
library(scales)  # 用于百分比坐标轴（可选）

meta <- excit_neurons2@meta.data

# 统计每个 group 中各 subtype 的细胞数，并计算比例
df_prop <- meta %>%
  count(group, subtype) %>%                    # 相当于 table(group, subtype)
  group_by(group) %>%
  mutate(prop = n / sum(n)) %>%               # 在每个 group 内求比例
  ungroup()

# 画堆积条形图（按比例）
ggplot(df_prop, aes(x = group, y = prop, fill = subtype)) +
  geom_bar(stat = "identity", position = "fill") +  # position = "fill" 即按比例堆积
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Group", y = "Cell proportion", fill = "Subtype") +
  theme_classic()
library(RColorBrewer)

ggplot(df_prop, aes(x = group, y = prop, fill = subtype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Blues") +  # 蓝色系配色
  labs(x = "Group", y = "Cell proportion", fill = "Subtype") +
  theme_classic()

mouse_layer_markers <- list(
  L2_3 = c("Cux1", "Cux2", "Reln"),
  L4   = c("Rorb", "Rprm"),
  L5   = c("Bcl11b", "Fezf2", "Tbr1"),
  L6   = c("Tle4", "Foxp2", "Pcp4")
)



library(Seurat)
library(ggplot2)

mouse_layer_markers <- list(
  L2_3 = c("Cux1", "Cux2", "Reln"),
  L6   = c("Tle4", "Foxp2", "Pcp4"),
  L4   = c("Rorb", "Rprm"),
  L5   = c("Bcl11b", "Fezf2", "Tbr1")
)

DotPlot(
  excit_neurons2,
  features = mouse_layer_markers,
  group.by = "subtype"
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightblue", high = "#0058A3") +
  theme_bw()


head(excit_neurons2)


# 看一下 group 里都有哪些分组
table(excit_neurons2$group)

# 用 group 作为当前身份
Idents(excit_neurons2) <- "group"

deg_KOF_vs_Ctrl <- FindMarkers(
  excit_neurons2,
  ident.1 = "KOF",
  ident.2 = "WT",   # 换成你实际的对照组名字
  logfc.threshold = 0.25,
  min.pct = 0.1,
  test.use = "wilcox"   # 默认即可
)

# 按校正后的 p 值排序
deg_KOF_vs_Ctrl <- deg_KOF_vs_Ctrl[order(deg_KOF_vs_Ctrl$p_val_adj), ]
head(deg_KOF_vs_Ctrl)


# 选 top 20 上调和下调基因
up_genes   <- rownames(deg_KOF_vs_Ctrl[deg_KOF_vs_Ctrl$avg_log2FC > 0, ])[1:20]
down_genes <- rownames(deg_KOF_vs_Ctrl[deg_KOF_vs_Ctrl$avg_log2FC < 0, ])[1:20]
features_use <- c(up_genes, down_genes)

DoHeatmap(
  excit_neurons2,
  features = features_use,
  group.by = "group"
) + scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))


library(ggplot2)

deg_KOF_vs_Ctrl$gene <- rownames(deg_KOF_vs_Ctrl)

# 选 p_adj 小的 top 20
top20_deg <- deg_KOF_vs_Ctrl[1:20, ]

ggplot(top20_deg, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC,
                      fill = avg_log2FC > 0)) +
  geom_col() +
  coord_flip() +
  labs(x = "Gene", y = "avg_log2FC (KOF vs Ctrl)", fill = "Up in KOF") +
  theme_bw()

head(scRNA_HKU)