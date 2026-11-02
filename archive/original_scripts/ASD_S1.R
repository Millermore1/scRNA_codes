library(Seurat)
library(dplyr)

# 1. 设置你的工作路径
data_dir <- "F:/sn_RNAseq_data/public_database/R_meta_data/ASD"
setwd(data_dir)

# 2. 获取目录下所有 .gz 压缩文件
files <- list.files(path = data_dir, pattern = "\\.gz$")

# 3. 通用正则匹配：同时兼容 '.' 和 '_' 连接的文件名后缀
sample_names <- unique(gsub("[\\._](barcodes|genes|features|matrix).*", "", files))

cat("共检测到", length(sample_names), "个样本：\n")
print(sample_names)

# 4. 自动创建文件夹并重命名移动文件为 Seurat 标准格式
for (s in sample_names) {
  sample_path <- file.path(data_dir, s)
  dir.create(sample_path, showWarnings = FALSE)
  
  # 匹配属于当前样本的文件
  s_files <- list.files(path = data_dir, pattern = paste0("^", gsub("\\+", "\\\\+", s), "[\\._]"))
  
  for (f in s_files) {
    old_file <- file.path(data_dir, f)
    
    # 根据关键字判定文件类型并归类
    if (grepl("barcodes", f)) {
      file.rename(from = old_file, to = file.path(sample_path, "barcodes.tsv.gz"))
    } else if (grepl("matrix", f)) {
      file.rename(from = old_file, to = file.path(sample_path, "matrix.mtx.gz"))
    } else if (grepl("features|genes", f)) {
      file.rename(from = old_file, to = file.path(sample_path, "features.tsv.gz"))
    }
  }
}

message("\n文件夹整理完成！开始批量读取并构建 Seurat 对象...")

# 5. 批量读取并自动解析所有样本的元数据 (Metadata)
seurat_list <- list()

for (s in sample_names) {
  message("正在读取: ", s)
  sample_path <- file.path(data_dir, s)
  
  # 读取 10x 矩阵
  cts <- Read10X(data.dir = sample_path)
  
  # 创建 Seurat 对象
  obj <- CreateSeuratObject(counts = cts, project = s, min.cells = 3, min.features = 200)
  
  # 智能提取分组 (Group: ASD / IMR / NC)
  if (grepl("ASD", s)) {
    obj$Group <- "ASD"
  } else if (grepl("IMR", s)) {
    obj$Group <- "IMR"
  } else if (grepl("NC", s)) {
    obj$Group <- "NC"
  } else {
    obj$Group <- "Other"
  }
  
  # 智能提取时间点 (Timepoint: D0, D6, D12, D30, D60, D100)
  tp <- regmatches(s, regexpr("D\\d+", s))
  if (length(tp) > 0) {
    obj$Timepoint <- tp
  } else {
    obj$Timepoint <- "Unknown"
  }
  
  seurat_list[[s]] <- obj
}

# 6. 一键合并所有样本
seu_asd <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = sample_names
)

# 7. 保存整合后的原始 Seurat 对象到本地
saveRDS(seu_asd, file = file.path(data_dir, "GSE203201_raw_seurat.rds"))

message("\n全部处理完毕！结果预览：")
print(seu_asd)

# 查看各个分组与时间点的细胞分布统计表
print(table(seu_asd$Timepoint, seu_asd$Group))

# 调整时间点因子顺序
seu_asd$Timepoint <- factor(
  seu_asd$Timepoint, 
  levels = c("D0", "D1", "D2", "D6", "D12", "D30", "D60", "D100")
)

# 再次查看表格，时间顺序就恢复正常了
table(seu_asd$Timepoint, seu_asd$Group)

# 计算线粒体基因百分比
seu_asd[["percent.mt"]] <- PercentageFeatureSet(seu_asd, pattern = "^MT-")

# 小鼠或人核基因组数据可同步查看核糖体基因 (可选)
seu_asd[["percent.rp"]] <- PercentageFeatureSet(seu_asd, pattern = "^RP[SL]")

# 绘制质控小提琴图
VlnPlot(seu_asd, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Group")

library(harmony)

# 1. 标准化与高变基因筛选
seu_asd <- NormalizeData(seu_asd)
seu_asd <- FindVariableFeatures(seu_asd, selection.method = "vst", nfeatures = 2000)
seu_asd <- ScaleData(seu_asd)
seu_asd <- RunPCA(seu_asd, npcs = 30)

# 2. 按样本批次 (orig.ident) 运行 Harmony 去批次
seu_asd <- RunHarmony(seu_asd, group.by.vars = "orig.ident")

# 3. 基于 Harmony 降维空间进行 UMAP 聚类
seu_asd <- RunUMAP(seu_asd, reduction = "harmony", dims = 1:30)
seu_asd <- FindNeighbors(seu_asd, reduction = "harmony", dims = 1:30)
seu_asd <- FindClusters(seu_asd, resolution = 0.1)

# 4. 可视化查看按组别和时间点的细胞分布
DimPlot(seu_asd, reduction = "umap", group.by = "Group_broad")
DimPlot(seu_asd, reduction = "umap", group.by = "Timepoint")
seu_asd$Group_broad <- ifelse(seu_asd$Group %in% c("IMR", "NC"), "Control", "ASD")
head(seu_asd)
table(seu_asd$Timepoint)
seu_asd$Timepoint <- sub(".*[\\._](D\\d+).*", "\\1", seu_asd$orig.ident)

# 2. 重新设置为正确的时间点因子顺序（去除了不存在的 D1 和 D2）
seu_asd$Timepoint <- factor(
  seu_asd$Timepoint, 
  levels = c("D0", "D6", "D12", "D30", "D60", "D100")
)

seu_D30 <- subset(seu_asd, subset = Timepoint == "D30")

table(seu_D30$Group_broad)



library(Seurat)
library(dplyr)
library(ggplot2)
library(Seurat)

# 1. 合并 Seurat v5 的数据图层 (Layers)
seu_D30 <- JoinLayers(seu_D30)

# 2. 重新设置默认分组
Idents(seu_D30) <- "Group_broad"

# 3. 运行差异基因分析
deg_D30_global <- FindMarkers(
  object = seu_D30,
  ident.1 = "ASD",
  ident.2 = "Control",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

# 查看前几行结果
head(deg_D30_global)
# 1. 把行名的基因名单独保存为 gene 列
deg_D30_global$gene <- rownames(deg_D30_global)

# 2. 导出为 .csv 文件（不保留默认数字行名）
write.csv(deg_D30_global, file = "DEGs_ASD_DAY30.csv", row.names = FALSE)


# 在 UMAP 上并排展示 ASD vs Control 的 PRMT8 表达情况
FeaturePlot(
  object = seu_asd, 
  features = "PRMT8", 
  split.by = "Group_broad",
  cols = c("lightgrey", "red")  # 表达量从低到高的颜色映射
)

library(Seurat)
library(ggplot2)

# 在 D30 数据集中按 Group_broad 绘制 PRMT8 小提琴图
VlnPlot(
  object = seu_D30,
  features = "PRMT8",
  group.by = "Group_broad",
  pt.size = 0.1,                      # 散点大小；如果觉得太密集，可设为 0 隐藏散点
  cols = c("ASD" = "#E64B35", "Control" = "#4DBBD5") # 自定义组别颜色
) + 
  labs(title = "PRMT8 Expression at D30", y = "Log-normalized Expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 全数据集：按时间点分组，按 Group_broad 拆分
VlnPlot(
  object = seu_asd,
  features = "PRMT8",
  group.by = "Timepoint",
  split.by = "Group_broad",
  pt.size = 0,                        # 细胞数量较多（10万+），建议设为 0 以防绘图卡顿
  cols = c("ASD" = "#E64B35", "Control" = "#4DBBD5")
) + 
  labs(title = "PRMT8 Expression Dynamics Across Timepoints", x = "Timepoint", y = "Expression") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))




# 1. 计算每个时间点 PRMT8 的平均表达量及标准误 (SEM)
prmt8_data <- FetchData(seu_asd, vars = c("PRMT8", "Timepoint"))

trend_df <- prmt8_data %>%
  group_by(Timepoint) %>%
  summarise(
    mean_exp = mean(PRMT8),
    sd_exp = sd(PRMT8),
    n = n(),
    se_exp = sd_exp / sqrt(n)
  )

# 2. 绘制表达量变化趋势折线图
ggplot(trend_df, aes(x = Timepoint, y = mean_exp, group = 1)) +
  geom_line(color = "#2b5c8f", size = 1.2) +
  geom_point(color = "#2b5c8f", size = 3) +
  geom_errorbar(aes(ymin = mean_exp - se_exp, ymax = mean_exp + se_exp), width = 0.15, color = "#2b5c8f") +
  labs(
    title = "PRMT8 Mean Expression Trajectory",
    x = "Timepoint",
    y = "Average Expression ± SEM"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 11, color = "black")
  )




library(Seurat)
library(dplyr)
library(ggplot2)

# 1. 提取 PRMT8 表达量以及对应的 Timepoint 和 Group_broad 信息
prmt8_data <- FetchData(seu_asd, vars = c("PRMT8", "Timepoint", "Group_broad"))

# 2. 按时间点 (Timepoint) 和组别 (Group_broad) 分组计算均值与标准误 (SEM)
trend_df_group <- prmt8_data %>%
  group_by(Timepoint, Group_broad) %>%
  summarise(
    mean_exp = mean(PRMT8),
    sd_exp   = sd(PRMT8),
    n        = n(),
    se_exp   = sd_exp / sqrt(n),
    .groups  = "drop"
  )

# 3. 绘制分组对比折线趋势图
p_line <- ggplot(trend_df_group, aes(x = Timepoint, y = mean_exp, color = Group_broad, group = Group_broad)) +
  geom_line(linewidth = 1.2) +                                # 表达量走势折线
  geom_point(size = 3) +                                      # 时间点数据点
  geom_errorbar(aes(ymin = mean_exp - se_exp, ymax = mean_exp + se_exp), width = 0.15, linewidth = 0.6) + # 误差棒 (SEM)
  scale_color_manual(values = c("ASD" = "#E64B35", "Control" = "#4DBBD5")) + # 自定义经典配色 (ASD: 红, Control: 蓝)
  labs(
    title = "PRMT8 Expression Dynamics: ASD vs Control",
    x = "Timepoint",
    y = "Average Expression ± SEM",
    color = "Group"
  ) +
  theme_classic() +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text    = element_text(size = 11, color = "black"),
    axis.title   = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_line)
head(seu_asd)
saveRDS(seu_asd, file = "ASD_GSE203201.rds")
Idents(seu_asd) <- "RNA_snn_res.0.1"
DimPlot(seu_asd)
