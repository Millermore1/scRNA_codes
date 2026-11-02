#1.29 整体分析
table(scRNA_HKU$group)
library(dplyr)
library(ggplot2)
meta_df <- scRNA_HKU@meta.data
colnames(meta_df)
# cell propotion
ggplot(meta_df, aes(x = group, fill = celltype)) +
  geom_bar(position = "fill") +    
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Group",
    y = "Cell proportion",
    fill = "Cell type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
# cell marker

library(Seurat)
library(ggplot2)

# 确认身份设置为 celltype
Idents(scRNA_HKU) <- "celltype"

# 选定要画的小鼠 marker 列表
markers_mouse <- c(
  # Excitatory neurons
  "Slc17a7", "Slc17a6", "Camk2a",
  # Inhibitory neurons
  "Gad1", "Gad2", "Slc6a1",
  # Oligodendrocyte
  "Mog", "Mbp", "Plp1",
  # OPC
  "Pdgfra", "Cspg4", "Olig2",
  # Astrocyte
  "Aqp4", "Gfap", "Aldh1l1",
  # Microglia
  "P2ry12", "Cx3cr1", "Tmem119",
  # Endothelial cells
  "Pecam1", "Cldn5", "Kdr"
)

# 画合并小提琴图
p <- VlnPlot(
  scRNA_HKU,
  features = markers_mouse,
  group.by = "celltype",  # 按细胞类型分组
  pt.size = 0,            # 不显示单细胞散点，图更干净
  ncol = 4                # 每行 4 个基因 panel，自行调节
) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1)
)

p


# 创建一个更简洁的版本，按细胞类型使用不同颜色
library(ggplot2)
library(dplyr)
library(tidyr)

# 提取表达数据
exp_data <- FetchData(
  scRNA_HKU, 
  vars = c(markers_mouse, "celltype")
)

# 转换为长格式
exp_long <- exp_data %>%
  pivot_longer(
    cols = all_of(markers_mouse),
    names_to = "gene",
    values_to = "expression"
  )

# 创建自定义分组
exp_long <- exp_long %>%
  mutate(
    cell_group = case_when(
      gene %in% markers_mouse[1:3] ~ "Excitatory",
      gene %in% markers_mouse[4:6] ~ "Inhibitory",
      gene %in% markers_mouse[7:9] ~ "Oligodendrocyte",
      gene %in% markers_mouse[10:12] ~ "OPC",
      gene %in% markers_mouse[13:15] ~ "Astrocyte",
      gene %in% markers_mouse[16:18] ~ "Microglia",
      gene %in% markers_mouse[19:21] ~ "Endothelial",
      TRUE ~ "Other"
    )
  )

# 绘制分组小提琴图
ggplot(exp_long, aes(x = celltype, y = expression, fill = celltype)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~ gene + cell_group, ncol = 3, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Cell Type",
    y = "Expression Level"
  ) +
  scale_fill_brewer(palette = "Set3")











## ===== 0. 加载包 =====
library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

## 假定你的对象名为 scRNA_HKU，meta.data 里有列：
##   group    （分组：KOF, KOF_Gq_CNO, KOF_vector_CNO, WT 等）
##   celltype （细胞类型：Excitatory neurons, Microglia 等）


## ===== 1. 定义 HP_AUTISM 基因集（MSigDB 原始人类基因名）=====

hp_autism_human <- c(
  "AARS1","ABCB7","ACTL6B","ADNP","ADSL","AHDC1","AKT1","ALAD","ALDH5A1",
  "ALG13","ALG14","ANKRD11","AP3B2","ARCN1","ARL6","ARNT2","ARVCF","ASCC3",
  "ASXL3","ATP1A2","ATP1A3","ATP6V1A","ATRX","BAZ1B","BBIP1","BBS1","BBS10",
  "BBS12","BBS2","BBS4","BBS5","BBS7","BBS9","BCKDK","BRD4","BUD23","CACNA1A",
  "CACNA1B","CACNA1C","CACNA2D1","CASZ1","CDK13","CDK19","CELF2","CEP19",
  "CEP290","CFAP418","CHD1","CHD7","CHD8","CHRNA7","CLIP2","CLTC","CNKSR2",
  "CNTNAP2","COMT","CREBBP","CYFIP2","CYP27A1","DALRD3","DEAF1","DHCR7",
  "DHDDS","DICER1","DNAJC30","DNM1","DOCK7","DPAGT1","DPYD","DPYS","DYRK1A",
  "EEF1A2","EHMT1","EIF4A2","EIF4H","ELN","EP300","ERI1","EXT2","FBXO28",
  "FGF12","FGFR1","FKBP6","FLCN","FLG","FMR1","FOXG1","FTCD","FTSJ1","FZR1",
  "GABBR2","GABRA2","GABRA5","GABRB2","GABRD","GABRG2","GAMT","GATM","GJA5",
  "GJA8","GLRA2","GNB1","GP1BB","GRIA3","GRIN2D","GTF2I","GTF2IRD1","GTF2IRD2",
  "H4C5","HCN1","HDAC4","HDAC8","HECTD4","HERC2","HESX1","HIRA","HNF1B",
  "HSPG2","IFNG","IFT172","IFT27","IFT74","IL1RAPL1","IQSEC2","JMJD1C",
  "KCNA2","KCNAB2","KCNB1","KCNC2","KDM5C","KLLN","KMT2A","LHX1","LIMK1",
  "LUZP1","LZTFL1","MAGEL2","MAN1B1","MAOA","MECP2","MED12","MED13L",
  "MEIS2","METTL27","MKKS","MKRN3","MKS1","MMP23B","NAA20","NAGA","NCF1",
  "NDN","NDP","NECAP1","NFIB","NHS","NIPA1","NIPA2","NIPBL","NLGN3","NONO",
  "NPAP1","NPHP1","NR2F1","NSUN2","NSUN6","NTRK2","NUS1","OCA2","OTX2",
  "PACS2","PAH","PAK2","PARS2","PDE4D","PDPN","PIGA","PIGL","PIK3CA",
  "PLXND1","PPP3CA","PRDM16","PRKCZ","PROKR2","PTCHD1","PTEN","PWAR1",
  "PWRN1","RAD21","RAI1","RERE","REV3L","RFC2","RPL10","RREB1","SATB1",
  "SATB2","SCAPER","SCLT1","SCN1A","SCN2A","SCN3A","SCN8A","SDCCAG8",
  "SDHB","SDHC","SDHD","SEC23B","SEC24C","SEMA3E","SETD2","SH2B1","SHANK3",
  "SIN3A","SKI","SLC13A5","SLC1A2","SLC35C1","SLC38A3","SLC9A6","SMAD4",
  "SMC1A","SMC3","SNORD115-1","SNORD116-1","SNRPN","SOX2","SOX3","SPEN",
  "SPTBN1","STS","STX1A","SUPT16H","SYNGAP1","SYNJ1","SZT2","TAF4","TAF6",
  "TAOK1","TBCK","TBL2","TBR1","TBX1","TCF12","TMEM270","TMLHE","TRAK1",
  "TRIM32","TSC1","TSC2","TTC8","TUBG1","UBA5","UBE3A","UBE4B","UFD1",
  "USF3","VPS37D","WDPCP","WFS1","WWOX","YWHAG","YY1","ZBTB20","ZFX"
)

## ===== 2. 将基因名自动匹配到对象（人/鼠风格都兼容） =====
## 如果你的表达矩阵是小鼠风格（首字母大写其余小写），则用 hp_autism_mouse；
## 如果是人类风格（全大写），则用 hp_autism_human。
## 下面自动选择重叠最多的那一种。

hp_autism_mouse <- str_to_title(tolower(hp_autism_human))  # 变成 Mouse 风格
all_genes <- rownames(scRNA_HKU)

genes_mouse <- intersect(all_genes, hp_autism_mouse)
genes_human <- intersect(all_genes, hp_autism_human)

cat("在对象中找到的 Mouse 风格基因数：", length(genes_mouse), "\n")
cat("在对象中找到的 Human 风格基因数：", length(genes_human), "\n")

genes_use <- if (length(genes_mouse) >= length(genes_human)) genes_mouse else genes_human

cat("用于打分的 HP_AUTISM 基因数：", length(genes_use), "\n")

if (length(genes_use) < 10) {
  warning("用于打分的基因少于 10 个，检查一下对象基因名风格是否正确！")
}

## ===== 3. AddModuleScore 计算每个细胞的 HP_AUTISM 打分 =====

Idents(scRNA_HKU) <- "celltype"  # 以后按细胞类型分层画图

scRNA_HKU <- AddModuleScore(
  scRNA_HKU,
  features = list(genes_use),
  name = "HP_AUTISM"
)
## 会生成 meta.data 列：HP_AUTISM1
head(scRNA_HKU@meta.data[, c("group", "celltype", "HP_AUTISM1")])

## ===== 4. 所有细胞整体：按 group 比较 HP_AUTISM 分数 =====

p_all <- VlnPlot(
  scRNA_HKU,
  features = "HP_AUTISM1",
  group.by = "group",
  pt.size = 0.1
) +
  labs(
    x = "Group",
    y = "HP_AUTISM score"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_all)

## ===== 5. 按细胞类型分面：每种 celltype 内的 group 比较 =====

df_score <- FetchData(
  scRNA_HKU,
  vars = c("HP_AUTISM1", "group", "celltype")
)

p_by_celltype <- ggplot(df_score, aes(x = group, y = HP_AUTISM1, fill = group)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.2, color = "black") +
  facet_wrap(~ celltype, scales = "free_y") +
  labs(
    x = "Group",
    y = "HP_AUTISM score"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_by_celltype)

## ===== 6.（可选）WT vs KOF 统计检验（按 celltype）=====

stats_WT_KOF <- df_score %>%
  filter(group %in% c("WT", "KOF")) %>%
  group_by(celltype) %>%
  summarize(
    p_value = wilcox.test(
      HP_AUTISM1[group == "WT"],
      HP_AUTISM1[group == "KOF"]
    )$p.value,
    .groups = "drop"
  ) %>%
  arrange(p_value)

print(stats_WT_KOF)









df_score <- FetchData(
  scRNA_HKU,
  vars = c("HP_AUTISM1", "group", "celltype")
)

# 按 celltype 和 group 计算平均得分
df_mean <- df_score %>%
  group_by(celltype, group) %>%
  summarise(
    mean_score = mean(HP_AUTISM1, na.rm = TRUE),
    .groups = "drop"
  )

library(pheatmap)
# 转成矩阵：行 = celltype，列 = group
mat_score <- df_mean %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = mean_score
  )

# 第一列是 celltype，作为行名
mat_score <- as.data.frame(mat_score)
rownames(mat_score) <- mat_score$celltype
mat_score$celltype <- NULL

mat_score <- as.matrix(mat_score)

# 可选：按行做 Z-score 标准化（便于比较不同细胞类型之间）
# mat_score_scaled <- t(scale(t(mat_score)))

# 画热图（不标准化版本）
pheatmap(
  mat_score,
  scale = "none",          # 如果上面用的是 mat_score_scaled，改成 "none" 一样
  clustering_method = "complete",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE, # 如果想在格子里显示具体数值，改 TRUE
  fontsize_row = 8,
  fontsize_col = 10,
  main = "HP_AUTISM score (mean per celltype × group)"
)





library(ComplexHeatmap)
help(package= "ComplexHeatmap")


GO1 <- read.csv("PV_Gq_vector_GO_UP.csv")
head(GO1)



# 1. 读入数据
GO1 <- read.csv("PV_Gq_vector_GO_UP.csv", stringsAsFactors = FALSE)

# 2. 把 LogP 变成正的 -LogP
GO1$minusLogP <- -GO1$LogP

# 3. 取前 50 个（按 -LogP 从大到小）
GO_plot <- GO1[order(GO1$minusLogP, decreasing = TRUE), ]
GO_plot <- head(GO_plot, 50)

library(ggplot2)

ggplot(GO_plot, aes(x = reorder(Description, minusLogP), 
                    y = minusLogP, 
                    color = Category)) +
  geom_point(size = 3) +
  xlab("GO/通路描述") +
  ylab("-LogP") +
  theme_bw() +
  coord_flip() +  # 横着放，文字更好看
  theme(
    axis.text.y = element_text(size = 8)
  )
ggplot(GO_plot, aes(x = Description, y = minusLogP, color = Category)) +
  geom_point(size = 3) +
  xlab("GO/通路描述") +
  ylab("-LogP") +
  theme_bw() +
  coord_flip





DimPlot(inh)
head(inh)
# 设定当前身份为 group
Idents(inh) <- "group"
table(inh$group)
# KOF vs Ctrl 差异分析（名字自己改成你实际的 group 名）
deg_KOF_vector_CNO_vs_KOF_Gq_CNO <- FindMarkers(
  object = inh,
  ident.1 = "KOF_vector_CNO",
  ident.2 = "KOF_Gq_CNO",   # 如果你的对照叫 WT 或其他，就改成对应名字
  logfc.threshold = 0.25,  # 
  min.pct = 0.1,
  test.use = "wilcox"
)

library(dplyr)
library(openxlsx) 

deg_KOF_vector_CNO_vs_KOF_Gq_CNO$gene <- rownames(deg_KOF_vector_CNO_vs_KOF_Gq_CNO)
up_deg <- deg_KOF_vector_CNO_vs_KOF_Gq_CNO %>%
  filter(
    avg_log2FC > 0,      # 
    p_val_adj < 0.05     # 
  )
write.xlsx(up_deg, file = "up_deg_KOF_vector_CNO_vs_KOF_Gq_CNO.xlsx", rowNames = FALSE)
setwd("D:/scrnadata")
getwd()






