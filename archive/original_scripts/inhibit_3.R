library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
#美化图片
Subclass_colors<-c("#69AED4","#B7B44B","#D96B7C","#69BE54","#939325",
                   "#E5E576","#127845","#761A55","#6A32A0","#77471F",
                   "#292774","#AA7947",
                   "#FEEE75",
                   "#61DDF2","#71ACD8","#BD80B8",
                   "#7472B5",
                   '#4B96D1', #Inh VIP
                   '#329CFC', # Inh LAMP5
                   "#8BCEA0", # Inh SST
                   "#D6569F",# Inh PVALB
                   "#FFB531", '#FDD3BA',# Inh PAX6
                   "#D62729","#FFBA7B","#B25828", # Ast Oli OPC 
                   "#CFA2D2","#9A6BBF", #Mic
                   "#E283B5","#F2C0D5","#6EC5A6","#A45627","#F38C64","#8BA0CC")
Subtype_colors<-c("#69AED4","#B7B44B","#D96B7C","#69BE54","#939325",
                  "#E5E576","#127845","#761A55","#6A32A0","#77471F",
                  "#292774","#AA7947",
                  "#DFAC7B","#7DC783","#0A8CED","#FEEE75","#E3AD28",
                  "#5576A6",
                  "#61DDF2","#4888C7","#71ACD8","#BD80B8",
                  "#34A07C","#AC4689",
                  "#7472B5",
                  '#FBFAD2','#CBE5C5',"#E83D36", "#62B8BF","#B0CFE0",
                  "#B6D88F","#FFDCAC","#DBDCDA","#F9DDE9", #Inh VIP
                  '#2FA147','#8466AF',"#EFE3CD", # Inh LAMP5
                  "#F8B6B4","#90D2D7","#B85E2C","#3581BA", # Inh SST
                  "#F5999D","#A27833","#FBBE6F","#BD96EA",# Inh PVALB
                  "#DFCCE4", '#FDD3BA',# Inh PAX6
                  "#FFDEDE","#FFB0BD","#D7517C","#991343", #Ast
                  "#F9CFA5","#f79c5e","#b55a00", # Oli
                  "#FFE2CF","#b38766","#75421f",# OPC
                  "#E1CAE5","#573B88", #Mic
                  "#7F0009","#6D0775", #Mac, T
                  "#E283B5","#F2C0D5","#6EC5A6","#A45627","#F38C64","#8BA0CC")


mycolors <- c("#0A8CED", "#D62729", "#D7517C")
head(inh)
DimPlot(
  inh,
  reduction = "umap",
  group.by = "inhibit_celltype",
  cols = mycolors,
  label = TRUE,
  label.size = 5,          # 调整标签大小
  repel = TRUE,
  pt.size = 0.5,          # 调整点的大小
  order = TRUE,           # 按顺序绘制，避免重叠
  raster = FALSE          # 矢量图，更清晰
) +
  theme_minimal() +       # 简洁背景
  theme(
    panel.grid = element_blank(),  # 移除网格线
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",     # 图例位置
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    title = "UMAP Visualization of Inhibitory Neuron Subtypes"
  )


# 假设你的细胞类型列叫 celltype，其中有一类是 "PVALB"
head(inh)
# 1. 新建一个只有两类的分组：PVALB vs Other
inh$PVALB_group <- ifelse(inh$inhibit_celltype == "SST", "SST", "Other")
inh$PVALB_group <- factor(inh$PVALB_group, levels = c("Other", "SST"))
# 一定要把 levels 顺序指定为：Other 在前，PVALB 在后

# 2. 画 UMAP，按这个分组上色
DimPlot(
  inh,
  reduction = "umap",
  group.by  = "PVALB_group",
  cols      = c("grey80", "#D62729"),  # Other = 灰, PVALB = 蓝
  pt.size   = 0.2
)

# 1. 设定分组列名：按哪个分组来算比例
group_col <- "group"    # 如果你的是 sample / condition，就改成对应的名字

# 2. 取出分组和细胞类型信息
meta <- inh@meta.data

# 确认这两列存在
table(meta[[group_col]])        # 看一下分组
table(meta$inhibit_celltype)    # 看一下细胞类型

# 3. 统计细胞数和比例
df.prop <- meta %>%
  count(!!sym(group_col), inhibit_celltype) %>%   # 每个分组、每个细胞类型的 cell 数
  group_by(!!sym(group_col)) %>%
  mutate(prop = n / sum(n)) %>%                   # 转成比例
  ungroup()

# 4. 画堆积比例柱状图
ggplot(df.prop, aes(x = !!sym(group_col), y = prop, fill = inhibit_celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = group_col, y = "Cell proportion", fill = "Cell type") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

# 如果想固定细胞类型的显示顺序，可以在这里自定义：
df.prop$inhibit_celltype <- factor(
  df.prop$inhibit_celltype,
  levels = c("PVALB", "SST", "VIP")  # 按你的实际类型改
)

# 自定义一套好看的配色（示例，可按需要改）
cell.cols <- c(
  "PVALB"   = "#0A8CED",
  "SST"     = "#D62729",
  "VIP"     = "#D7517C")

ggplot(df.prop, aes(x = !!sym(group_col), y = prop, fill = inhibit_celltype)) +
  geom_bar(stat = "identity", position = "fill", color = NA) +   # 去掉边框
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = cell.cols) +
  labs(
    x    = NULL,
    y    = "Cell proportion (%)",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y  = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 5)),
    legend.position   = "right",
    legend.title      = element_text(face = "bold"),
    legend.key.size   = unit(0.4, "cm"),
    legend.text       = element_text(size = 9),
    axis.ticks.length = unit(0.15, "cm")
  )

getwd()


library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(patchwork)
library(tinyarray)
library(dplyr)
library(tibble)
#caculate DEGs
VIP <- subset(inh, subset = inhibit_celltype == "VIP")
VIP.DEG_KOF_Gq_vs_vector <- FindMarkers(
  object = pv,
  ident.1 = "WT",
  ident.2 = "KOF",
  min.pct = 0.1,
  logfc.threshold = 0.25,   #
  test.use = "wilcox"
)
vip_deg.sig <- VIP.DEG_KOF_Gq_vs_vector[
  VIP.DEG_KOF_Gq_vs_vector$p_val_adj < 0.05 & abs(VIP.DEG_KOF_Gq_vs_vector$avg_log2FC) > 0.25, 
]
vip_deg.sig <- cbind(gene = rownames(vip_deg.sig), vip_deg.sig)
write.csv(
  deg.sig,
  file = "vip_deg.sig_WT_KOF_DEGs.csv",
  row.names = FALSE
)
vip_deg.sig.up.KOF <- subset(vip_deg.sig, avg_log2FC > 0)
vip_deg.sig.up.WT  <- subset(vip_deg.sig, avg_log2FC < 0)

write.csv(vip_deg.sig.up.KOF, "vip_deg.sig.up.KOF.csv", row.names = FALSE)
write.csv(vip_deg.sig.up.WT,  "vip_deg.sig.up.WT.csv",  row.names = FALSE)

deg <- pv.DEG %>%
  rownames_to_column("gene")    # 

deg <- deg %>%
  mutate(
    group = case_when(
      p_val_adj < 0.05 & avg_log2FC >  0.25 ~ "Up",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down",
      TRUE                                  ~ "NS"
    )
  )



deg_data <- pv.DEG_KOF_Gq_vs_vector %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%         # 基因名这一列叫 gene
  mutate(
    log2FoldChange = avg_log2FC,        # 为了后面统一用这个名字
    padj           = p_val_adj
  ) %>%
  filter(!is.na(padj))                  # 去掉 NA

## 2. 设置显著性分组（你可以根据自己阈值修改）
volcano_df <- deg_data %>%
  mutate(
    negLog10Padj = -log10(padj),
    sig = case_when(
      padj < 0.05 & log2FoldChange >= 1  ~ "Up",
      padj < 0.05 & log2FoldChange <= -1 ~ "Down",
      TRUE                               ~ "NS"
    )
  )

## 3. 画普通火山图
p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = negLog10Padj, color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +   # log2FC 阈值
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +# padj 阈值
  theme_classic(base_size = 14) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj. p-value)",
    color = "Regulation"
  )












