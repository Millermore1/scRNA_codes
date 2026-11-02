library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)

species <- "mouse"  # 如果是人类数据改为 "human"

if (species == "human") {
  CellChatDB <- CellChatDB.human
} else if (species == "mouse") {
  CellChatDB <- CellChatDB.mouse
} else {
  stop("species 只能是 'human' 或 'mouse'")
}

# 1.2 确保有 celltype 注释
if (is.null(scRNA_HKU$celltype)) {
  stop("scRNA_HKU 对象中没有 $celltype 注释，请先完成细胞类型注释。")
}

# 可选：去掉 celltype 为 NA 的细胞
scRNA_HKU <- subset(scRNA_HKU, cells = colnames(scRNA_HKU)[!is.na(scRNA_HKU$celltype)])

# 1.3 确认使用的 assay 和数据类型
DefaultAssay(scRNA_HKU) <- "RNA"
# 如果还没 Normalize，可以先做一次
# scRNA_HKU <- NormalizeData(scRNA_HKU)

###############################
## 2. 单一条件的 CellChat 分析
##   （如果你只分析整体，不区分 group，直接做这一部分）
###############################

## 2.1 从 Seurat 提取表达矩阵和元数据
data.input <- GetAssayData(scRNA_HKU, assay = "RNA", slot = "data")

meta <- data.frame(
  celltype = scRNA_HKU$celltype,
  row.names = colnames(scRNA_HKU)
)

# 检查一下
print(table(meta$celltype))

## 2.2 创建 CellChat 对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

# 赋予数据库
cellchat@DB <- CellChatDB

# 如只看分泌信号，可用：
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
# cellchat@DB <- CellChatDB.use

## 2.3 预处理：筛选高表达基因和配体-受体对
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 可选：平滑表达数据
# cellchat <- projectData(cellchat)

## 2.4 计算通讯概率和路径
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)

# 过滤细胞太少的配对（默认最少10个细胞）
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

# 汇总网络（通路层面、配体-受体层面）
cellchat <- aggregateNet(cellchat)

## 2.5 保存结果对象
saveRDS(cellchat, file = "cellchat_single_condition.rds")

###############################
## 3. 导出结果表
###############################

# 3.1 所有配体-受体通讯
df.net.all <- subsetCommunication(cellchat)
write.csv(df.net.all, "CellChat_all_interactions.csv", row.names = FALSE)

# 3.2 示例：只看 Astrocyte → Microglia
if ("Astrocyte" %in% meta$celltype & "Microglia" %in% meta$celltype) {
  df.Astro_Micro <- subsetCommunication(
    cellchat,
    sources.use = "Astrocyte",
    targets.use = "Microglia"
  )
  write.csv(df.Astro_Micro, "CellChat_Astro_to_Microglia.csv", row.names = FALSE)
}

###############################
## 4. 可视化单条件网络
###############################

# 4.1 画整体细胞群间通讯的圈图
groupSize <- as.numeric(table(cellchat@idents))

pdf("CellChat_circle_overall_count.pdf", width = 7, height = 7)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
dev.off()

pdf("CellChat_circle_overall_weight.pdf", width = 7, height = 7)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights"
)
dev.off()

# 4.2 某个具体通路的网络图
pathways.show <- cellchat@netP$pathways
print(head(pathways.show))

# 举例：如果有 "TNF" 通路
if ("NRXN" %in% pathways.show) {
  pdf("CellChat_NRXN_circle.pdf", width = 6, height = 6)
  netVisual_aggregate(
    cellchat,
    signaling = "NRXN",
    layout = "circle"
  )
  dev.off()
  
  pdf("CellChat_NRXN_heatmap.pdf", width = 6, height = 6)
  netVisual_heatmap(
    cellchat,
    signaling = "NRXN",
    color.heatmap = "Reds"
  )
  dev.off()
}

# 4.3 信号发送/接收能力的气泡图（所有通路）
pdf("CellChat_overall_bubble.pdf", width = 9, height = 7)
netAnalysis_signalingRole_network(cellchat, signaling = NULL, width = 8, height = 2.5)
dev.off()

###############################
## 5. 多条件（group）差异通讯分析
##   （如果你有多个样本组，比如 KOF vs WT，用这一部分）
###############################

# 确认有 group 信息
if (!is.null(scRNA_HKU$group)) {
  
  print(table(scRNA_HKU$group))
  
  # 分组运行 CellChat：每个 group 独立建一个 CellChat 对象
  cellchat_list <- list()
  
  for (grp in unique(scRNA_HKU$group)) {
    message("Processing group: ", grp)
    
    obj.sub <- subset(scRNA_HKU, subset = group == grp)
    
    # 确保还有 celltype 信息
    obj.sub <- subset(obj.sub, cells = colnames(obj.sub)[!is.na(obj.sub$celltype)])
    
    data.input.sub <- GetAssayData(obj.sub, assay = "RNA", slot = "data")
    meta.sub <- data.frame(
      celltype = obj.sub$celltype,
      row.names = colnames(obj.sub)
    )
    
    cellchat.sub <- createCellChat(object = data.input.sub, meta = meta.sub, group.by = "celltype")
    cellchat.sub@DB <- CellChatDB
    
    cellchat.sub <- subsetData(cellchat.sub)
    cellchat.sub <- identifyOverExpressedGenes(cellchat.sub)
    cellchat.sub <- identifyOverExpressedInteractions(cellchat.sub)
    cellchat.sub <- computeCommunProb(cellchat.sub, raw.use = TRUE)
    cellchat.sub <- filterCommunication(cellchat.sub, min.cells = 10)
    cellchat.sub <- computeCommunProbPathway(cellchat.sub)
    cellchat.sub <- aggregateNet(cellchat.sub)
    
    cellchat_list[[grp]] <- cellchat.sub
  }
  
  # 保存每个组
  saveRDS(cellchat_list, file = "cellchat_list_by_group.rds")
  
  ###############################
  ## 6. 合并并比较不同条件
  ###############################
  
  # 合并 CellChat 对象（按 group 标记）
  cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
  
  # 按功能相似性计算网络相似度
  cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
  cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
  cellchat_merged <- netClustering(cellchat_merged, type = "functional")
  
  # UMAP 可视化不同组的信号网络
  pdf("CellChat_multi_condition_functional_embedding.pdf", width = 7, height = 6)
  netVisual_embeddingPairwise(cellchat_merged, type = "functional", label.size = 3.5)
  dev.off()
  
  # 也可以对结构相似性做类似分析（可选）
  # cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
  # cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
  # cellchat_merged <- netClustering(cellchat_merged, type = "structural")
  # pdf("CellChat_multi_condition_structural_embedding.pdf", width = 7, height = 6)
  # netVisual_embeddingPairwise(cellchat_merged, type = "structural", label.size = 3.5)
  # dev.off()
  
  # 单个通路在不同组之间的对比（举例：WNT）
  if ("WNT" %in% cellchat_merged@netP$pathways) {
    pdf("CellChat_WNT_compare_between_groups.pdf", width = 10, height = 5)
    netVisual_bubble(
      cellchat_merged,
      signaling = "WNT",
      comparison = c(1, 2),   # 第一组 vs 第二组，如 KOF vs WT
      angle.x = 45,
      remove.isolate = TRUE
    )
    dev.off()
  }
  
  # 保存合并对象
  saveRDS(cellchat_merged, file = "cellchat_merged_groups.rds")
}

###############################
## 脚本结束
###############################








#绘图
library(CellChat)
library(patchwork)
library(ggplot2)

# 读取合并后的多组比较结果
cellchat_merged <- readRDS("CellChat_Merged_4groups.rds")

# 继续进行多组比较的可视化，例如：
# 比较两组之间的总通讯数量和强度的差异
gg1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2))
gg2 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "weight")
gg1 + gg2

# 或者查看特定通路在两组间的气泡图差异
netVisual_bubble(cellchat_merged, signaling = "CXCL", comparison = c(1, 2), angle.x = 45)


# 读取包含各个分组独立对象的列表
cellchat_list <- readRDS("cellchat_list_by_group.rds")

# 例如，你想分别画出 WT 和 KOF 组内 "CXCL" 通路的网络图进行肉眼对比：
weight.max <- getMaxWeight(cellchat_list, attribute = c("idents","count"))

par(mfrow = c(1,4), xpd=TRUE)
for (i in 1:length(cellchat_list)) {
  netVisual_circle(cellchat_list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cellchat_list)[i]))
}



#绘图2
library(CellChat)
library(scales)

cellchat_list <- readRDS("cellchat_list_by_group.rds")

## 保证所有图细胞类型顺序一致
cell.order <- levels(cellchat_list[[1]]@idents)

## 统一配色：色盲友好 Okabe-Ito 风格
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)

color.use <- setNames(
  colorRampPalette(okabe_ito)(length(cell.order)),
  cell.order
)

## 提取并重排网络矩阵，保证行列顺序一致
mat.list <- lapply(cellchat_list, function(obj) {
  mat <- obj@net$count
  mat <- mat[cell.order, cell.order, drop = FALSE]
  return(mat)
})

## 统一边宽尺度，方便组间肉眼比较
edge.max <- max(unlist(lapply(mat.list, max)), na.rm = TRUE)

## 统一节点大小尺度
group.size.list <- lapply(cellchat_list, function(obj) {
  gs <- as.numeric(table(obj@idents)[cell.order])
  gs[is.na(gs)] <- 0
  return(gs)
})

vertex.max <- max(unlist(group.size.list), na.rm = TRUE)

## 输出矢量 PDF，适合投稿后期编辑
grDevices::cairo_pdf(
  filename = "CellChat_number_interactions_circle_publication.pdf",
  width = 3.8 * length(cellchat_list),
  height = 4.2,
  family = "Arial"
)

par(
  mfrow = c(1, length(cellchat_list)),
  mar = c(0.5, 0.5, 2.2, 0.5),
  oma = c(0, 0, 0, 0),
  xpd = NA
)

for (i in seq_along(cellchat_list)) {
  
  netVisual_circle(
    mat.list[[i]],
    color.use = color.use,
    vertex.weight = group.size.list[[i]],
    vertex.weight.max = vertex.max,
    vertex.size.max = 18,
    weight.scale = TRUE,
    edge.weight.max = edge.max,
    edge.width.max = 8,
    alpha.edge = 0.55,
    edge.curved = 0.20,
    label.edge = FALSE,
    vertex.label.cex = 0.85,
    vertex.label.color = "black",
    margin = 0.25,
    title.name = paste0("Number of interactions\n", names(cellchat_list)[i])
  )
}

dev.off()
# 绘图修改

library(CellChat)

cellchat_list <- readRDS("cellchat_list_by_group.rds")

## 1. 统一细胞类型配色（顺序必须与 idents 的 levels 一致！）
celltype_colors <- c(
  "Excitatory neurons" = "#E64B35",
  "Inhibitory neurons" = "#F39B7F",
  "Oligodendrocyte"    = "#3C5488",
  "OPC"                = "#4DBBD5",
  "Astrocyte"          = "#00A087",
  "Microglia"          = "#B09FC1",
  "Endothelial cells"  = "#E8C559"
)
ident_levels <- levels(cellchat_list[[1]]@idents)
color.use <- celltype_colors[ident_levels]
if (any(is.na(color.use))) stop("有细胞类型不在配色表中，请检查 ident_levels")

## 2. 统一标尺：边粗细 + 顶点大小，保证组间可比
weight.max      <- getMaxWeight(cellchat_list, attribute = c("idents", "count"))
vertex.size.max <- max(unlist(lapply(cellchat_list,
                                     function(x) max(as.numeric(table(x@idents))))))

## 3. 绘图：cairo_pdf 字体更清晰，组名大字居中
n_grp <- length(cellchat_list)
cairo_pdf("CellChat_circle_beautified.pdf",
          width = 4.5 * n_grp, height = 5.2, family = "Arial")
par(mfrow = c(1, n_grp), xpd = TRUE, mar = c(0.5, 0.5, 0.5, 0.5))

for (i in seq_along(cellchat_list)) {
  netVisual_circle(
    cellchat_list[[i]]@net$count,
    color.use        = color.use,          # 节点/边颜色
    vertex.weight    = as.numeric(table(cellchat_list[[i]]@idents)),
    vertex.size.max  = 14,                 # 统一最大节点
    weight.scale     = TRUE,
    edge.weight.max  = weight.max[2],      # 统一边粗细标尺（关键！）
    edge.width.max   = 10,
    alpha.edge       = 0.55,               # 边半透明，减少杂乱感
    edge.curved      = 0.25,               # 边弧度
    label.edge       = FALSE,
    vertex.label.cex = 0.85,
    vertex.label.color = "grey20",
    title.name       = NULL,               # 不要默认标题
    margin           = 0.25
  )
  ## 组名大字居中（仿参考图样式）
  text(0, 0, names(cellchat_list)[i], cex = 2.4, font = 2, col = "grey15")
}
dev.off()


#
library(CellChat)
library(scales)
library(RColorBrewer)

cellchat_list <- readRDS("cellchat_list_by_group.rds")

## 保证所有图细胞类型顺序一致
cell.order <- levels(cellchat_list[[1]]@idents)

## 统一配色：ColorBrewer Set2 风格（柔和期刊风，与 count 图区分）
set2_base <- brewer.pal(8, "Set2")

color.use <- setNames(
  colorRampPalette(set2_base)(length(cell.order)),
  cell.order
)

## 提取并重排互作强度矩阵，保证行列顺序一致
mat.list <- lapply(cellchat_list, function(obj) {
  mat <- obj@net$weight
  mat <- mat[cell.order, cell.order, drop = FALSE]
  return(mat)
})

## 统一边宽尺度，方便组间肉眼比较
edge.max <- max(unlist(lapply(mat.list, max)), na.rm = TRUE)

## 统一节点大小尺度
group.size.list <- lapply(cellchat_list, function(obj) {
  gs <- as.numeric(table(obj@idents)[cell.order])
  gs[is.na(gs)] <- 0
  return(gs)
})

vertex.max <- max(unlist(group.size.list), na.rm = TRUE)

## 输出矢量 PDF，适合投稿后期编辑
grDevices::cairo_pdf(
  filename = "CellChat_interaction_strength_circle_publication.pdf",
  width = 3.8 * length(cellchat_list),
  height = 4.2,
  family = "Arial"
)

par(
  mfrow = c(1, length(cellchat_list)),
  mar = c(0.5, 0.5, 2.2, 0.5),
  oma = c(0, 0, 0, 0),
  xpd = NA
)

for (i in seq_along(cellchat_list)) {
  
  netVisual_circle(
    mat.list[[i]],
    color.use = color.use,
    vertex.weight = group.size.list[[i]],
    vertex.weight.max = vertex.max,
    vertex.size.max = 18,
    weight.scale = TRUE,
    edge.weight.max = edge.max,
    edge.width.max = 8,
    alpha.edge = 0.55,
    edge.curved = 0.20,
    label.edge = FALSE,
    vertex.label.cex = 0.85,
    vertex.label.color = "black",
    margin = 0.25,
    title.name = paste0("Interaction strength\n", names(cellchat_list)[i])
  )
}

dev.off()