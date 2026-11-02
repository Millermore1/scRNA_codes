library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(harmony)
setwd("D:/scrnadata/")
load("scRNA_HKU.Rdata")
save(inh,file = "inhibit_new")
load("inhibit_new")
head(scRNA_HKU)
#抑制性神经元的分析
Idents(scRNA_HKU) <- "celltype"


# 2. 提取抑制性神经元
inh <- subset(scRNA_HKU, idents = "Inhibitory neurons")
table(inh$celltype)
inh <- NormalizeData(inh)
inh <- FindVariableFeatures(inh)
inh <- ScaleData(inh)
inh <- RunPCA(inh)
inh <- FindNeighbors(inh, dims = 1:20)
inh <- FindClusters(inh, resolution = 0.1)  # 
inh <- RunUMAP(inh, dims = 1:20)
DimPlot(inh,label = TRUE)

inh_markers <- FindAllMarkers(object = inh, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

#
inh_markers_1 =inh_markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
top200 = inh_markers_1 %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.csv(top200,"12cluster_top200.csv",row.names = T)

inhibitory.markers <- list(
  VIP   = c("Vip", "Crh", "Igf1", "Reln", "Vipr2"),
  SST   = c("Sst", "Gpr149", "Thsd7a", "Chodl"),
  PV    = c("Pvalb", "Gabra1", "Gabra4", "Kcnc1"),
  TH    = c("Th", "Slc6a3", "Ddc", "Slc18a2"),
  LAMP5 = c("Lamp5", "Reln", "Lhx1", "Lhx2", "Nos1"),
  PAX6  = c("Pax6", "Vip", "Reln", "Igf1")
)


# 按 VIP → SST → PV → TH → LAMP5 → PAX6 的顺序拼接
features.use <- c(
  inhibitory.markers$VIP,
  inhibitory.markers$SST,
  inhibitory.markers$PV,
  inhibitory.markers$TH,
  inhibitory.markers$LAMP5,
  inhibitory.markers$PAX6
)

# 去掉重复基因（比如 Reln / Vip 在多个 list 里）
features.use <- unique(features.use)

# 基础 DotPlot
p <- DotPlot(
  inh,
  features = features.use,
  group.by = "seurat_clusters"   # 或者 "group" / "celltype" 等
)

# 美化：颜色、字体、轴标签等
p + 
  scale_color_gradient(
    low  = "lightgrey",
    high = "firebrick"           # 可换成 "dodgerblue4" 等
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y  = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 9)
  ) +
  labs(
    x     = "Seurat clusters",
    y     = "Marker genes",
    color = "Average expression",
    size  = "% expressing"
  )


bfreaname.pbmc <- inh
new.cluster.ids <- c("PVALB", "SST","PVALB", "VIP", "PVALB", "PVALB", "VIP", 
                     "PVALB", "VIP","SST","SST")
names(new.cluster.ids) <- levels(inh) 
inh$inhibit_celltype <- Idents(inh)
inh <- RenameIdents(inh, new.cluster.ids)
save(inh,file = "inhibit_new.Rdata")
DimPlot(inh)
head(inh)
inhibit_markers <- c()



