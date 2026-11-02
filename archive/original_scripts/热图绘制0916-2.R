#
DotPlot(scRNA_seurat, features = c("PVALB", "SST", "VIP", "NPY","CALB2","RELN","GAD1","GAD2"),group.by = "RNA_snn_res.0.05")
DimPlot(scRNA_seurat,group.by = "RNA_snn_res.0.05",label = TRUE)
head(scRNA_seurat)

cells_subset <- subset(scRNA_seurat, subset = RNA_snn_res.0.05 == 2)

# 2. 进行数据预处理（可选，根据需要）
# 例如：归一化
cells_subset <- NormalizeData(cells_subset)
# 识别变量基因
cells_subset <- FindVariableFeatures(cells_subset)
# 数据缩放
cells_subset <- ScaleData(cells_subset)
# 进行PCA
cells_subset <- RunPCA(cells_subset)

# 3. 邻接图和聚类（参数可以调整）
cells_subset <- FindNeighbors(cells_subset, dims = 1:10)
cells_subset <- FindClusters(cells_subset, resolution = 0.01)

# 4. 结果可视化（如t-SNE或UMAP）
cells_subset <- RunUMAP(cells_subset, dims = 1:10)
DimPlot(cells_subset, reduction = "umap", label = TRUE)

DotPlot(p, features = c("GAD1", "GAD2", "SST", "PVALB", "VIP", "NPY", "CALB2", "RELN","PRMT8"),group.by = "RNA_snn_res.0.05")

rm(scRNA_seurat)
new.cluster.ids <- c("VIP+", "PV+","SST+","PV+", "PV+")
names(new.cluster.ids) <- levels(cells_subset)
p <- RenameIdents(cells_subset, new.cluster.ids)
levels(p) #查看是否已改名
DimPlot(cells_subset)


head(p)
DotPlot(p, features = c("GAD1","PVALB", "SST",  "VIP",  "CALB2", "RELN","PRMT8"),dot.scale = 15,group.by = "RNA_snn_res.0.05") +
  RotatedAxis() +  # 让基因名标签旋转，便于阅读
  scale_color_gradient(low = "blue", high = "red")  # 调整颜色渐变

Idents(p) <- p$RNA_snn_res.0.05

inhibitory_markers <- FindAllMarkers(object = p, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

all_inhibitory_markers.markers =inhibitory_markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)
#筛选出P<0.05的marker基因
top301 = all_inhibitory_markers.markers %>% group_by(cluster) %>% top_n(n = 300, wt = avg_log2FC) #将每个cluster lgFC排在前10的marker基因挑选出来
# top15 = all.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

write.csv(top301,"top301.csv",row.names = T)
getwd()


# 读取本地CSV文件
GOdata <- read.csv("C:/Users/黄廉/Desktop/GO.csv")
head(GOdata)
# 加载必要的包
library(ggplot2)

# 假设GOdata已存在
# 如果还没有安装ggplot2，可以先安装
# install.packages("ggplot2")

# 绘制条形图，颜色深浅反映Z.score
ggplot(GOdata, aes(x = reorder(GO.Biological.Processes, Z.score), y = Z.score, fill = Z.score)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # 横向条形图，便于阅读
  scale_fill_gradient(low = "#CCE5E5", high = "#0E8585") +
  labs(title = "GO Enrichment Analysis",
       x = "Biological Processes",
       y = "Z-score",
       fill = "Enrichment") +
  theme_minimal()
rm(inhibitory_markers)
head(p)

library(dplyr)


# 提取元数据
meta <- p@meta.data

# 筛选出cluster=11的细胞
cluster11_meta <- meta %>% filter(RNA_snn_res.0.05 == 11)

# 获取细胞ID（行名）
cells_in_cluster11 <- rownames(cluster11_meta)

# 如果需要在Seurat对象中筛选这些细胞
cells_in_cluster11_obj <- subset(p, cells = cells_in_cluster11)

library(dplyr)

proportions <- meta %>%
  filter(RNA_snn_res.0.05 == 11) %>%
  group_by(group) %>%
  summarise(cluster11_count = n())

total_cells <- meta %>%
  group_by(group) %>%
  summarise(total = n())

# 计算比例
proportions <- left_join(proportions, total_cells, by = "group") %>%
  mutate(proportion = cluster11_count / total)
print(proportions)
library(Seurat)
DimPlot(cells_subset)
rm(list = ls())
packageVersion("ggplot2")
packageVersion("Seurat")

install.packages("Seurat")
library(Seurat)





