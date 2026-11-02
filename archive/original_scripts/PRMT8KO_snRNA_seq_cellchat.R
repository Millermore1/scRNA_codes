library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)

cellchat.merged <- readRDS("CellChat_Merged_4groups.rds")
# 1. 从合并对象的所有组别中提取通路，并取去重后的并集
all_pathways <- unique(unlist(lapply(cellchat.merged@netP, function(x) x$pathways)))

print("你的数据中包含的所有通路：")
print(all_pathways)

# 2. 定义你感兴趣的突触相关通路候选名单 (必须是大写)
synapse_candidates <- c("GLUTAMATE", "GABA", "NRXN", "NLGN", "NCAM", 
                        "CADM", "EPHA", "EPHB", "NRG", "PTN", "MK")

# 3. 取交集，找出你的数据中真正存在的突触通路
synapse_pathways_in_data <- intersect(synapse_candidates, all_pathways)

print("--------------------------------------------------")
print("在你的数据中实际有显著通讯的突触通路：")
print(synapse_pathways_in_data)


print("正在生成突触通路气泡图，请稍候...")
pdf("CellChat_Synapse_Bubble_Merged.pdf", width = 12, height = 8)

# 绘制多组对比的气泡图
netVisual_bubble(cellchat.merged, 
                 targets.use = c("PN", "PVALB"),       # 你想重点看的接收细胞
                 signaling = synapse_pathways_in_data, # 刚刚精准提取的 7 条突触通路
                 comparison = c(1, 2, 3, 4),           # 你的 4 个组别顺序
                 angle.x = 45,
                 title.name = "Synapse-related signaling to PN and PVALB")

dev.off()
print("✅ 气泡图已完成！请在工作目录查看 CellChat_Synapse_Bubble_Merged.pdf")



print("正在生成突触通路信息流条形图，请稍候...")
pdf("CellChat_Synapse_Information_Flow.pdf", width = 8, height = 6)

# 绘制多组信息流 (Information flow) 堆叠图
rankNet(cellchat.merged, 
        signaling = synapse_pathways_in_data, # 刚刚精准提取的 7 条突触通路
        mode = "comparison", 
        stacked = TRUE, 
        do.stat = TRUE)

dev.off()
print("✅ 信息流对比图已完成！请在工作目录查看 CellChat_Synapse_Information_Flow.pdf")



print("正在生成整体通讯数量与强度对比图...")
pdf("CellChat_Compare_Overall_Interactions.pdf", width = 10, height = 5)

# 将画布一分为二，左边画数量，右边画强度
par(mfrow = c(1,2), xpd=TRUE) 

# 比较总通讯数量
compareInteractions(cellchat.merged, show.legend = F, group = c(1,2,3,4)) 
# 比较总通讯强度 (weight)
compareInteractions(cellchat.merged, show.legend = F, group = c(1,2,3,4), measure = "weight")

dev.off()
print("✅ 总体对比图已保存为 CellChat_Compare_Overall_Interactions.pdf")



print("正在生成两组间通路变化散点图...")
pdf("CellChat_Signaling_Changes_Scatter_WT_vs_KOF.pdf", width = 8, height = 6)

# id.x = 4 (你的WT组), id.y = 1 (你的KOF实验组)
# 请根据你的实际组别数字编号进行修改！
netAnalysis_signalingChanges_scatter(cellchat.merged, id.x = 4, id.y = 1)

dev.off()
print("✅ 通路差异散点图已保存为 CellChat_Signaling_Changes_Scatter_WT_vs_KOF.pdf")



