library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(viridis)
library(ggplot2)
head(scRNA_HKU)
table(scRNA_HKU$celltype)
rm(scRNA_HKU)
#extract all neurons to do cell communication calculate
scRNA_neurons <- subset(scRNA_HKU, celltype %in% c("Excitatory neurons", "Inhibitory neurons"))
table(scRNA_neurons$celltype)
scRNA_neurons <- NormalizeData(scRNA_neurons)
scRNA_neurons <- FindVariableFeatures(scRNA_neurons, selection.method = "vst", nfeatures = 2000)
scRNA_neurons <- ScaleData(scRNA_neurons, features = VariableFeatures(scRNA_neurons))
scRNA_neurons <- RunPCA(scRNA_neurons, features = VariableFeatures(scRNA_neurons))
scRNA_neurons <- FindNeighbors(scRNA_neurons, dims = 1:30)
scRNA_neurons <- FindClusters(scRNA_neurons, resolution = 0.05)
scRNA_neurons <- RunUMAP(scRNA_neurons, dims = 1:30)
DimPlot(scRNA_neurons, label = "true")
DimPlot(scRNA_neurons,group.by = "celltype")

#define the PV+ interneurons
inhibitory.markers_unique <- list(
  VIP   = c("Vip", "Crh", "Reln", "Vipr2"), 
  SST   = c("Sst", "Gpr149", "Thsd7a", "Chodl"),
  PV    = c("Pvalb", "Gabra1", "Gabra4", "Kcnc1"),
  TH    = c("Th", "Slc6a3", "Ddc", "Slc18a2"),
  LAMP5 = c("Lamp5", "Lhx1", "Lhx2", "Nos1"),
  PAX6  = c("Pax6", "Igf1")
)


rm(scRNA_HKU)
DotPlot(scRNA_neurons, features = inhibitory.markers_unique) + 
  RotatedAxis() + 
  ggtitle("Inhibitory Neuron Markers")

new.cluster.ids <- c("PN", "PVALB", "PN", "PN", "SST", "PVALB", "PN", "VIP", 
                     "SST", "PN", "PN", "PN", "PN", "SST", "PN")
names(new.cluster.ids) <- levels(scRNA_neurons) 
scRNA_neurons$original_clusters <- Idents(scRNA_neurons)
scRNA_neurons <- RenameIdents(scRNA_neurons, new.cluster.ids)
scRNA_neurons$neuron_celltype <- Idents(scRNA_neurons)
DimPlot(scRNA_neurons, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Neuron Subtypes Annotation")

head(scRNA_neurons)

table(scRNA_neurons$group)


library(CellChat)
library(Seurat)
library(patchwork)

Idents(scRNA_neurons) <- "neuron_celltype"

# 1. 按 group 拆分 Seurat 对象
sc_list <- SplitObject(scRNA_neurons, split.by = "group")

# 2. 准备小鼠配受体数据库
CellChatDB <- CellChatDB.mouse

# 3. 创建一个空的 list 来存放跑完的 CellChat 对象
cellchat.list <- list()

# 4. 循环处理每个组
for (g in names(sc_list)) {
  print(paste0("正在计算分组: ", g))
  
  # 提取该组的表达矩阵和 metadata
  data.input <- GetAssayData(sc_list[[g]], assay = "RNA", slot = "data")
  meta <- sc_list[[g]]@meta.data
  
  # 创建 CellChat 对象
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "neuron_celltype")
  cellchat@DB <- CellChatDB
  
  # 预处理与网络推断（标准流程）
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # 将算好的对象存入 list
  cellchat.list[[g]] <- cellchat
}

# 检查一下你的 list 名称和顺序 (通常是 KOF, KOF_Gq_CNO, KOF_vector_CNO, WT)
names(cellchat.list)

# 提取你希望比较的组的名称，并进行合并
# add.names 参数会把你的 group 名字打在合并后的对象上
cellchat.merged <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))



# 1. 为独立的 4 个组分别计算中心度 (这一步最花时间)
print("步骤 1/4: 正在为 cellchat.list 中的每个组计算中心度...")
cellchat.list <- lapply(cellchat.list, function(x) {
  # 计算每个小对象的中心度
  netAnalysis_computeCentrality(x, slot.name = "netP")
})

# 2. 重新合并对象，并强制覆盖掉旧的 cellchat.merged
print("步骤 2/4: 正在重新合并带有中心度信息的对象...")
cellchat.merged <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

# 3. 检查中心度是否成功写入合并后的对象
if (is.null(cellchat.merged@netP$centr)) {
  stop("合并失败：中心度信息依然为空，请检查 cellchat.list 是否正确！")
} else {
  print("步骤 3/4: 中心度合并成功！准备绘制散点图...")
}

# 4. 输出到 PDF 以防止画图超时
print("步骤 4/4: 正在生成 PDF 文件...")
pdf(file = "CellChat_Signaling_Role_Scatter_Fixed.pdf", width = 16, height = 12)

# 执行画图
netAnalysis_signalingRole_scatter(cellchat.merged, title = "Signaling roles across conditions")

# 保存并关闭文件
dev.off()

print("全部完成！请去工作目录查看 CellChat_Signaling_Role_Scatter_Fixed.pdf")




# ================= 请从这里开始完整复制 =================

# 1. 强制为 list 中的每个组逐一计算中心度（带进度提示）
for (i in seq_along(cellchat.list)) {
  group_name <- names(cellchat.list)[i]
  print(paste0(">>> 正在为分组 [", group_name, "] 计算中心度，请稍候..."))
  
  # 覆盖写入带有中心度信息的新对象
  cellchat.list[[i]] <- netAnalysis_computeCentrality(cellchat.list[[i]], slot.name = "netP")
}

# 2. 合并前硬核自检：验证第一个对象是否真的有了中心度
if (is.null(cellchat.list[[1]]@netP$centr)) {
  stop("❌ 步骤1计算失败！请检查上方是否有红字报错阻断了运行！")
} else {
  print("✅ 所有独立分组的中心度计算完毕！准备合并...")
}

# 3. 重新合并对象 (这个时候由于上面算好了，合并进去的就带中心度了)
cellchat.merged <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

# 4. 生成 PDF 并画图
print(">>> 正在生成 PDF 图像文件...")
pdf(file = "CellChat_Signaling_Role_Scatter_Success.pdf", width = 16, height = 12)

# 绘制散点图
netAnalysis_signalingRole_scatter(cellchat.merged, title = "Signaling roles across conditions")

# 保存并关闭文件
dev.off()

print("🎉 终极合并与画图全部完成！快去你的工作目录里寻找 PDF 文件吧！")

# ================= 复制到此结束 =================








# 确保加载了拼接图形的包
library(patchwork)

# 1. 创建一个空的列表，用来存放 4 张图
plot_list <- list()

# 2. 循环遍历 cellchat.list，为每个组单独画散点图
for (i in seq_along(cellchat.list)) {
  group_name <- names(cellchat.list)[i]
  
  # 直接对独立的包含中心度的对象画图
  p <- netAnalysis_signalingRole_scatter(cellchat.list[[i]], title = group_name)
  
  # 将画好的图存入列表
  plot_list[[i]] <- p
}

# 3. 将 4 张图拼接成 2行2列 的大图
combined_plot <- wrap_plots(plots = plot_list, ncol = 2)

# 4. 输出到 PDF (这次绝对不会报错了！)
print(">>> 正在生成最终的 4 组对比散点图 PDF...")
pdf(file = "CellChat_Signaling_Role_Scatter_Final.pdf", width = 16, height = 12)

# 打印拼接好的图形
print(combined_plot)

dev.off()

print("🎉 大功告成！请去工作目录查看 CellChat_Signaling_Role_Scatter_Final.pdf")



# 假设你想看 "SST" 发送给 "PN" 和 "PVALB" 的所有分子相互作用，在 4 个组中的区别
# comparison = c(1,2,3,4) 表示按你的 cellchat.list 里的顺序展示这 4 个组
pdf("CellChat_BubblePlot_SST_to_Others.pdf", width = 10, height = 8)

netVisual_bubble(cellchat.merged, 
                 sources.use = "SST",             # 发送方
                 targets.use = c("PN", "PVALB"),  # 接收方
                 comparison = c(1, 2, 3, 4),      # 对比的组别
                 angle.x = 45)                    # X轴标签倾斜度

dev.off()




pdf("CellChat_Information_Flow_Rank.pdf", width = 8, height = 10)

# 比较各组中所有信号通路的总信息流 (Information flow = 通讯概率之和)
rankNet(cellchat.merged, mode = "comparison", stacked = T, do.stat = TRUE)

dev.off()









# 假设你想比较 WT 和 KOF_Gq_CNO 的差异
# 需要先确认这两个组在你的 cellchat.list 里的索引号
# 假设 KOF_Gq_CNO 是第3个，WT 是第4个 (请根据 names(cellchat.list) 的实际顺序修改)

pdf("CellChat_Diff_Circle_Gq.pdf", width = 8, height = 8)

# comparison = c(4, 3) 意味着以组4(WT)为基准，看组3(KOF_Gq)的变化
head(cellchat.list)
netVisual_diffInteraction(cellchat.merged, 
                          weight.scale = T, 
                          measure = "weight", 
                          comparison = c(3, 2), 
                          title.name = "Differential interaction strength (KOF vs gq)")
dev.off()

saveRDS(scRNA_neurons, file = "scRNA_neuron.rds")
saveRDS(cellchat.merged, file = "CellChat_Merged_4groups.rds")
saveRDS(cellchat.list, file = "CellChat_List_4groups.rds")








