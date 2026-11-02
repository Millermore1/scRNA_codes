library(ggplot2)
library(Seurat)
library(dplyr)
library(scales)

DimPlot(excit2)
table(excit2$group)
# 1. 确保设置正确的默认分析层 (Assay)
DefaultAssay(excit2) <- "RNA"
# 2. 计算差异表达基因
# 注意：min.pct 指在至少某个比例的细胞中表达，可根据需要调整
kof_vs_wt <- FindMarkers(
  excit2, 
  ident.1 = "KOF", 
  ident.2 = "WT", 
  group.by = "group",
  test.use = "wilcox", 
  min.pct = 0.1, 
  logfc.threshold = 0.25
)

table(excit2$group)
kof_vector_vs_Gq <- FindMarkers(
  excit2, 
  ident.1 = "KOF_vector_CNO", 
  ident.2 = "KOF_Gq_CNO", 
  group.by = "group",
  test.use = "wilcox", 
  min.pct = 0.1, 
  logfc.threshold = 0.25
)
significant_genes_vector_Gq <- kof_vector_vs_Gq[kof_vector_vs_Gq$p_val_adj < 0.05,]


# 3. 筛选显著基因 (通常使用 padj < 0.05)
significant_genes_vector_vs_Gq <- kof_vs_wt[kof_vs_wt$p_val_adj < 0.05, ]


library(writexl)

# 2. 同样需要把行名（基因名）转换为单独的一列
significant_genes_df_Gq <- data.frame(
  Gene = rownames(significant_genes_vector_Gq),
  significant_genes_vector_Gq,
  stringsAsFactors = FALSE
)

# 3. 导出为 Excel 文件
write_xlsx(significant_genes_df_Gq, path = "Excit_significant_genes_df_Gq.xlsx")





# 4. 可视化查看结果
library(EnhancedVolcano)
EnhancedVolcano(kof_vector_vs_Gq,
                lab = rownames(kof_vector_vs_Gq),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'vector vs Gq: Differentially Expressed Genes',
                pCutoff = 0.05,
                FCcutoff = 0.5)




significant_genes_df <- data.frame(
  Gene = rownames(significant_genes),
  significant_genes,
  stringsAsFactors = FALSE
)


write_xlsx(significant_genes_vector_Gq, path = "Excit_vector_vs_Gq_significant_genes.xlsx")


subset(excit2,)
table(excit2$subtype)


exc_l5 <- subset(excit2, subset = subtype == "Exc_L5")


table(exc_l5$subtype)   # 应该全是 Exc_L5
head(exc_l5)
# 设定分组身份（把 group 改成你元数据中的列名，比如 condition / orig.ident）
Idents(exc_l5) <- "group"
table(Idents(exc_l5))

deg_l5_group <- FindMarkers(
  object = exc_l5,
  ident.1 = "KOF",
  ident.2 = "WT",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

head(deg_l5_group)

write.csv(deg_l5_group, "DEG_ExcL5_KOF_vs_WT.csv")
getwd()


# 1. 加载必要的包
library(dplyr)
library(stringr)
library(ggplot2)

# 2. 读取你的数据
rectome <- read.csv("F:/sn_RNAseq_data/PRMT8KO_snRNA_seq/Excit_neurons_vectorvsGq_r1_enrichment_CO .csv")

# 3. 数据清洗：把 | 替换成逗号加空格
plot_bar_data <- rectome %>%
  mutate(Clean_Genes = str_replace_all(Genes, "\\|", ", "))

# 4. 绘制纯蓝、不加粗的经典条形图
ggplot(plot_bar_data, aes(x = Z.score, y = reorder(Description, Z.score))) +
  # 将 fill 统一写在 geom_col 内部，实现纯色填充（这里用了很好看的 steelblue 钢蓝色）
  geom_col(width = 0.68, fill = "#9DDCAF") +
  # 右侧基因标签保持学术规范的斜体（非加粗）
  geom_text(aes(label = Clean_Genes), 
            hjust = 0, 
            nudge_x = 0.15, 
            size = 3.2, 
            fontface = "italic", 
            color = "grey20") +
  # 预留右侧空间防止文字被切掉
  expand_limits(x = max(plot_bar_data$Z.score) * 1.6) +
  # 极简经典主题
  theme_classic(base_size = 12) +
  labs(
    x = "Activation Z-score",
    y = NULL,
    title = "GO Pathway Enrichment & Driver Genes",
    subtitle = "Cell Type: Excitatory Neurons (PRMT8 KO vs WT)"
  ) +
  theme(
    # 所有字体的 face 全部调整为 "plain"（普通不加粗）
    plot.title = element_text(face = "plain", size = 13, margin = margin(b=5)),
    plot.subtitle = element_text(face = "plain", size = 10, color = "grey30"),
    axis.text.x = element_text(face = "plain", color = "black"),
    axis.text.y = element_text(face = "plain", color = "black", size = 10),
    axis.title.x = element_text(face = "plain", size = 11),
    # 辅助线和坐标轴线
    axis.line = element_line(color = "grey50"),
    panel.grid.major.x = element_line(color = "grey95", linetype = "dashed")
  )
getwd()