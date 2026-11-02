library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)
library(harmony)
## 将文件夹中的文件名存到"dir_name"中
dir_name <- list.files("E:/scdata")
dir_name
####Read10X(data.dir绝对路径读取####
dir_name=list.files("E:/scdata/")
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("E:/scdata/", dir_name[i], sep = ""))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = dir_name[i],min.cells = 3, min.features = 300)
}
# a=Read10X("数据呢/T2//")

for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]#获取scRNAlist中的第i个Seurat对象
  # 计算线粒体比例
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")#计算以"MT-"开头的基因比例
  # 计算红细胞比例
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") #定义红细胞基因列表
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))#在Seurat对象的RNA数据中查找红细胞基因的索引位置
  HB_genes <- rownames(sc@assays$RNA)[HB_m]  # 获取匹配到的红细胞基因的行名
  HB_genes <- HB_genes[!is.na(HB_genes)]  # 删除NA值（未匹配到的基因）
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features=HB_genes)   #计算红细胞基因比例，将结果存储在名为"HB_percent"的新列中
  # 将sc赋值给scRNAlist[[i]]
  scRNAlist[[i]] <- sc
  # 删除sc
  rm(sc)
}
# 查看当前默认的分析对象
DefaultAssay(scRNAlist[[1]])
# 查看“RNA” assay中的数据层
slotNames(scRNAlist[[1]]@assays$RNA)

violin_before <- list()
for(i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                                pt.size = 0.01, 
                                ncol = 4) 
}



scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & 
                mt_percent < 10 & 
                HB_percent < 3 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & 
                nCount_RNA > 1000)})

View(scRNAlist[[1]]@meta.data)

scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])
## 统计细胞数
table(scRNAlist[[]]$orig.ident)

violin_after <- VlnPlot(scRNAlist,
                        features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                        pt.size = 0.01,
                        ncol = 4)
# 将图片输出到画板上 
violin_after




# harmony整合是基于PCA降维结果进行的。
scRNAlist <- NormalizeData(scRNAlist) %>% # 数据归一化处理
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% #筛选高变基因
  ScaleData() %>% #数据标准化
  RunPCA(npcs = 30, verbose = T)#npcs：计算和存储的PC数（默认为 50）
a=DimPlot(scRNAlist,reduction = "pca",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
a


scRNAlist <- JoinLayers(scRNAlist)#连接layers层的count数据


scRNA_harmony <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
#PCA图看到还有一点的批次效应（融合得比较好批次就弱）
b
# 合并图片
pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated
#去批次前后对比


# 后续都是基于Harmony矫正之后的数据，不是基因表达数据和直接的PCA降维数据。
# 设置reduction = 'harmony'，后续分析是基于Harmony矫正之后的数据。
####聚类、umap/tsne降维降维####
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution = 0.05)#分辨率可以自己调
##umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:10)



umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
# 合并图片
umap_integrated1
umap_integrated2
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)


markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
#筛选出P<0.05的marker基因
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
# top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
View(top50)
write.csv(top50,"12cluster_top10.csv",row.names = T)
getwd()

FeaturePlot(scRNA_harmony, features = c("SLC17A6", "GAD1"), reduction = "umap")

FeaturePlot(scRNA_harmony, features = c("SLC17A7", "GAD1", "GFAP", "MBP", "CX3CR1", "CLDN5"))
FeaturePlot(scRNA_harmony, features = c("PRMT8","PSD95","G3BP1","FYN","SYN1","LCK"))

new.cluster.ids <- c("Exci.neuron", "Oligodendrocytes","Astrocyte", "Inhibitory.Neurons", "Inhibitory.Neurons", "Microglia", "Oligodendrocytes", 
                     "Inhibitory.Neurons", "Endothelial")
names(new.cluster.ids) <- levels(scRNA_harmony)
p <- RenameIdents(scRNA_harmony, new.cluster.ids)
levels(p) #查看是否已改名
DimPlot(scRNA_harmony)

library(Seurat)
library(ggplot2)
library(RColorBrewer)

# 用深色调色板
mycolors <- brewer.pal(length(unique(p$seurat_clusters)), "Dark2")

DimPlot(
  p, reduction = "umap", label = TRUE, pt.size = 0.5, 
  cols = mycolors, label.size = 6
) + 
  NoLegend() +
  theme_void() + # 去除背景和边框
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    panel.border = element_blank() # 明确去除边框
  )

save(p, file = "2025.804.RData")
getwd()



newdata$celltype <- Idents(newdata)
head(newdata@meta.data)

head(newdata)


