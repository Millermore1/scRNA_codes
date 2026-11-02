library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(harmony)
folder_path <- "C:/Users/黄廉/Desktop/GSE203201"
h5_files <- list.files(folder_path, pattern = "\\.h5$", full.names = TRUE)

# 读取所有 h5 文件并创建 Seurat 对象列表
seurat_list <- lapply(h5_files, function(f) {
  counts <- Read10X_h5(f)
  obj <- CreateSeuratObject(counts, project = basename(f))
  obj$sample <- basename(f) # 给每个细胞加上样本标签
  obj
})

# 合并所有 Seurat 对象
merged_seurat <- merge(x = seurat_list[[1]], y = seurat_list[-1])

# 1. 质控（QC）
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 根据需要设置阈值进行细胞过滤（示例参数，请根据实际数据调整）
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# 2. 标准化
merged_seurat <- NormalizeData(merged_seurat)

# 3. 识别高变基因
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

# 4. 数据缩放
merged_seurat <- ScaleData(merged_seurat)

# 5. 主成分分析（PCA）
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))
ElbowPlot(merged_seurat)

# 6. 最近邻查找和聚类
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.1)

# 7. UMAP可视化
merged_seurat <- RunUMAP(merged_seurat, dims = 1:20)
DimPlot(merged_seurat, reduction = "umap",label = TRUE)
head(merged_seurat)
# 如果是Seurat对象
merged_seurat$group <- ifelse(grepl("ASD", merged_seurat$orig.ident), "ASD", "Ctrl")
table(merged_seurat$group)

merged_seurat <- JoinLayers(merged_seurat, assay = "RNA")

cellmarkers <- FindAllMarkers(object = merged_seurat, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)


#2.对计算好的每cluster的marker基因进行筛选
all.markers = cellmarkers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
#筛选出P<0.05的marker基因
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
# top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
View(top10)
write.csv(top10,"12cluster_top10.csv",row.names = T)
DotPlot(merged_seurat, features = c("DES","DCN","PDGFRA"))




new.cluster.ids <- c("Neural Progenitors", "Neural Progenitors", "Neural Progenitors", "Neural Progenitors", "neuroepithelial", "Neurons","Neural
crest/Mes" )
names(new.cluster.ids) <- levels(merged_seurat) 
merged_seurat <- RenameIdents(merged_seurat, new.cluster.ids)
DimPlot(merged_seurat)
head(merged_seurat)
head(deg_ASD_vs_Ctrl)


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# 1. 提取上调基因
deg_up <- deg_ASD_vs_Ctrl %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)
gene_list <- rownames(deg_up)

# 2. SYMBOL转ENTREZID
gene_df <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_ids <- gene_df$ENTREZID

# 3. GO富集
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

# 4. 可视化
barplot(ego, showCategory=20, font.size=12, title="GO Enrichment")
dotplot(ego, showCategory=20, font.size=12, title="GO Enrichment")
head(merged_seurat)



library(Seurat)
library(pheatmap)

# 假设你的基因列表如下
gene_list <- c(
  "LAMP5", "PTN", "SIX3", "ATP1A2", "ID2", "HES5", "HES4", "CRABP1", "SEMA6A", "MEIS2",
  "CCND1", "EFNB2", "FLNA", "NREP", "NCALD", "PLPPR3", "RBFOX2", "KIF5C", "EPS8L1", "MDK",
  "IGFBP5", "GLI3", "CAMK2N1", "SRGAP3", "CHURC1", "NECTIN2", "NECTIN3", "RBM10", "MAPK12",
  "SPARC", "CALM1"
)

# 提取表达矩阵（以RNA为例，log-normalized）
expr_mat <- GetAssayData(merged_seurat, slot = "data")[gene_list, ]

# 获取细胞分组（orig.ident）
group_info <- merged_seurat@meta.data$orig.ident
names(group_info) <- colnames(merged_seurat)

# 可选：按组取均值，做热图
library(Matrix)
library(dplyr)

# 构建每组的均值表达矩阵
avg_expr <- sapply(unique(group_info), function(grp) {
  cells <- names(group_info[group_info == grp])
  if(length(cells) == 1) {
    expr_mat[,cells]
  } else {
    Matrix::rowMeans(expr_mat[,cells,drop=FALSE])
  }
})


colnames(avg_expr) <- unique(group_info)
rownames(avg_expr) <- gene_list

pheatmap(avg_expr, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, show_colnames = TRUE)
col_order <- c("ASD_1.h5", "ASD_2.h5", "ASD_3.h5", "Ctrl_1.h5", "Ctrl_2.h5", "Ctrl_3.h5")
avg_expr <- avg_expr[, col_order]

pheatmap(avg_expr, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = TRUE, show_colnames = TRUE)

colnames(avg_expr)
library(pheatmap)
my_palette <- colorRampPalette(c("#95BEEC", "white", "#611A44"))(100)
pheatmap(
  avg_expr,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = my_palette
)


heatmap

library(pheatmap)

# 自定义颜色渐变
my_palette <- colorRampPalette(c("#AED6D5", "white", "#B86078"))(100)

# 画热图
pheatmap(
  avg_expr,
  scale = "row",                # 行标准化
  cluster_rows = TRUE,          # 对行聚类
  cluster_cols = FALSE,         # 不对列聚类
  show_rownames = TRUE,         # 显示行名
  show_colnames = TRUE,         # 显示列名
  color = my_palette,           # 自定义配色
  border_color = NA,            # 去掉单元格边框（更美观）
  fontsize_row = 10,            # 行名字体大小
  fontsize_col = 10,            # 列名字体大小
  cellwidth = 20,               # 可选：单元格宽度
  cellheight = 15               # 可选：单元格高度
)

library(enrichplot)
desc <- gsea_result@result$Description
idx <- grep("protein-RNA complex assembly", desc, ignore.case = TRUE)
if(length(idx)>0){
  gseaplot2(gsea_result, geneSetID = idx[1], title = desc[idx[1]])
} else {
  print("没有找到该通路，请检查通路名称或基因集。")
}

DEG_results <- FindMarkers(
  merged_seurat, 
  ident.1 = "ASD", 
  ident.2 = "Ctrl", 
  group.by = "group",
  min.pct = 0.25,  # 最少表达百分比筛选
  logfc.threshold = 0.25 # 最小对数差异
)

head(merged_seurat)
library(Seurat)
library(limma)
library(clusterProfiler)


# 2. 按PRMT8表达量分组
prmt8_expr <- FetchData(merged_seurat, vars = "PRMT8")[,1]
prmt8_group <- ifelse(prmt8_expr > median(prmt8_expr), "High", "Low")
merged_seurat$PRMT8_group <- prmt8_group
table(merged_seurat$PRMT8_group)

# 3. 提取表达矩阵和分组信息
expr_mat <- as.matrix(GetAssayData(merged_seurat, slot = "data"))  # log-normalized表达量
group <- merged_seurat$PRMT8_group

# 4. 差异表达分析
design <- model.matrix(~ group)
fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)
deg <- topTable(fit, coef=2, number=Inf)

# 5. 构建GSEA排名向量
gene_list <- deg$logFC
names(gene_list) <- rownames(deg)
gene_list <- sort(gene_list, decreasing = TRUE)

# 6. 读取基因集和GSEA分析
gene_sets <- read.gmt("your_path/to/geneset.gmt") # 请替换为你的gmt文件路径
gsea_result <- GSEA(geneList = gene_list, TERM2GENE = gene_sets)


library(org.Hs.eg.db)
gsea_result <- gseGO(
  geneList = gene_list, 
  OrgDb = org.Hs.eg.db, 
  ont = "BP", 
  keyType = "SYMBOL"
)
dotplot(gsea_result)

ridgeplot(gsea_result)

head(gsea_result)


# 输出所有通路的名称（中文或英文，取决于你的.gmt文件）
all_pathways <- gsea_result@result$Description
print(all_pathways)


DimPlot(merged_seurat)
head(merged_seurat)
table(merged_seurat$group)
