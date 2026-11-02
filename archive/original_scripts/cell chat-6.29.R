#rewrited codes
# ==============================================================================
# 第一部分：环境准备与 R 包加载
# ==============================================================================
# 加载必要的 R 包
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(viridis)
library(ggplot2)
library(CellChat)

# ==============================================================================
# 第二部分：Seurat 数据预处理与降维聚类
# ==============================================================================
# 预览原始数据（假设 scRNA_HKU 已在环境中）
head(scRNA_HKU)
table(scRNA_HKU$celltype)

# 提取所有的神经元（兴奋性与抑制性）进行后续的细胞通讯分析
scRNA_neurons <- subset(scRNA_HKU, celltype %in% c("Excitatory neurons", "Inhibitory neurons"))
table(scRNA_neurons$celltype)

# 释放内存，移除庞大的原始对象
rm(scRNA_HKU)

# 标准 Seurat 预处理流程
scRNA_neurons <- NormalizeData(scRNA_neurons)
scRNA_neurons <- FindVariableFeatures(scRNA_neurons, selection.method = "vst", nfeatures = 2000)
scRNA_neurons <- ScaleData(scRNA_neurons, features = VariableFeatures(scRNA_neurons))
scRNA_neurons <- RunPCA(scRNA_neurons, features = VariableFeatures(scRNA_neurons))

# 聚类与 UMAP 降维
scRNA_neurons <- FindNeighbors(scRNA_neurons, dims = 1:30)
scRNA_neurons <- FindClusters(scRNA_neurons, resolution = 0.05)
scRNA_neurons <- RunUMAP(scRNA_neurons, dims = 1:30)

# 基础可视化检查
DimPlot(scRNA_neurons, label = TRUE)
DimPlot(scRNA_neurons, group.by = "celltype")

# ==============================================================================
# 第三部分：神经元亚群精细注释
# ==============================================================================
# 定义 PV+ 等抑制性中间神经元的特征基因集
inhibitory.markers_unique <- list(
  VIP   = c("Vip", "Crh", "Reln", "Vipr2"), 
  SST   = c("Sst", "Gpr149", "Thsd7a", "Chodl"),
  PV    = c("Pvalb", "Gabra1", "Gabra4", "Kcnc1"),
  TH    = c("Th", "Slc6a3", "Ddc", "Slc18a2"),
  LAMP5 = c("Lamp5", "Lhx1", "Lhx2", "Nos1"),
  PAX6  = c("Pax6", "Igf1")
)

# 绘制气泡图检查 Marker 基因表达情况
DotPlot(scRNA_neurons, features = inhibitory.markers_unique) + 
  RotatedAxis() + 
  ggtitle("Inhibitory Neuron Markers")

# 根据聚类结果重命名细胞亚群（注意：这里的长度需与你的 clusters 数量严格一致）
new.cluster.ids <- c("PN", "PVALB", "PN", "PN", "SST", "PVALB", "PN", "VIP", 
                     "SST", "PN", "PN", "PN", "PN", "SST", "PN")
names(new.cluster.ids) <- levels(scRNA_neurons) 

# 保存原始聚类信息并重命名 Idents
scRNA_neurons$original_clusters <- Idents(scRNA_neurons)
scRNA_neurons <- RenameIdents(scRNA_neurons, new.cluster.ids)
scRNA_neurons$neuron_celltype <- Idents(scRNA_neurons)

# 可视化注释后的 UMAP 图
DimPlot(scRNA_neurons, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Neuron Subtypes Annotation")

# 检查分组信息（如 WT, KOF_Gq_CNO 等）
table(scRNA_neurons$group)

# ==============================================================================
# 第四部分：CellChat 细胞通讯分析 (单组独立计算)
# ==============================================================================
Idents(scRNA_neurons) <- "neuron_celltype"

# 1. 按 group 将 Seurat 对象拆分为列表
sc_list <- SplitObject(scRNA_neurons, split.by = "group")

# 2. 准备小鼠配受体数据库
CellChatDB <- CellChatDB.mouse
cellchat.list <- list()

# 3. 循环处理每个组的通讯推断
for (g in names(sc_list)) {
  print(paste0(">>> 正在推断分组通讯网络: ", g))
  
  # 提取数据与元信息
  data.input <- GetAssayData(sc_list[[g]], assay = "RNA", slot = "data")
  meta <- sc_list[[g]]@meta.data
  
  # 创建并初始化 CellChat 对象
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "neuron_celltype")
  cellchat@DB <- CellChatDB
  
  # 标准推断流程
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat.list[[g]] <- cellchat
}

# ==============================================================================
# 第五部分：网络中心度计算与多组合并
# ==============================================================================
# 必须先为 list 中的每个独立对象计算中心度，然后再合并，否则后续画图会报错
for (i in seq_along(cellchat.list)) {
  group_name <- names(cellchat.list)[i]
  print(paste0(">>> 正在为分组 [", group_name, "] 计算网络中心度..."))
  
  # 计算并覆写回列表
  cellchat.list[[i]] <- netAnalysis_computeCentrality(cellchat.list[[i]], slot.name = "netP")
}

# 验证第一个对象是否成功写入中心度
if (is.null(cellchat.list[[1]]@netP$centr)) {
  stop("❌ 中心度计算失败！请检查上方报错信息！")
} else {
  print("✅ 所有独立分组中心度计算完毕！准备合并...")
}

# 合并包含中心度的 CellChat 对象
cellchat.merged <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

# ==============================================================================
# 第六部分：高级可视化输出 (PDF)
# ==============================================================================

# 1. 独立分组信号角色散点图 (使用 Patchwork 拼接 2x2 排版)
print(">>> 正在生成 4 组独立信号角色对比图...")
plot_list <- list()
for (i in seq_along(cellchat.list)) {
  group_name <- names(cellchat.list)[i]
  plot_list[[i]] <- netAnalysis_signalingRole_scatter(cellchat.list[[i]], title = group_name)
}
combined_plot <- wrap_plots(plots = plot_list, ncol = 2)

pdf(file = "CellChat_Signaling_Role_Scatter_Patchwork.pdf", width = 16, height = 12)
print(combined_plot)
dev.off()

# 2. 气泡图：比较特定细胞 (SST) 给其它细胞 (PN, PVALB) 发送的信号在 4 组间的差异
pdf("CellChat_BubblePlot_SST_to_Others.pdf", width = 10, height = 8)
netVisual_bubble(cellchat.merged, 
                 sources.use = "SST",             # 发送方
                 targets.use = c("PN", "PVALB"),  # 接收方
                 comparison = c(1, 2, 3, 4),      # 对比所有 4 个组别
                 angle.x = 45)                    
dev.off()

# 3. 排名图：比较各组信号通路总信息流 (Information Flow) 的变化
pdf("CellChat_Information_Flow_Rank.pdf", width = 8, height = 10)
rankNet(cellchat.merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
dev.off()

# 4. 差异交互图：比较组 3 (KOF_Gq) 和 组 2 的通讯强度差异
# 注：这里的 comparison = c(3, 2) 意为以组2为基准，观察组3的变化
pdf("CellChat_Diff_Circle_Gq.pdf", width = 8, height = 8)
netVisual_diffInteraction(cellchat.merged, 
                          weight.scale = TRUE, 
                          measure = "weight", 
                          comparison = c(3, 2), 
                          title.name = "Differential interaction strength (KOF vs gq)")
dev.off()

# ==============================================================================
# 第七部分：持久化保存数据
# ==============================================================================
# 保存中间计算的大型对象，方便日后直接载入复现
print(">>> 正在保存 RDS 文件，请耐心等待...")
saveRDS(scRNA_neurons, file = "scRNA_neuron.rds")
saveRDS(cellchat.merged, file = "CellChat_Merged_4groups.rds")
saveRDS(cellchat.list, file = "CellChat_List_4groups.rds")
print("🎉 所有分析和数据保存大功告成！")