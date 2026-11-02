#1.28 analyse
library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)
library(dplyr)
DEG_KOF_Gq_vs_vector_1 <- FindMarkers(
  object = inh,
  ident.1 = "KOF_Gq_CNO",
  ident.2 = "KOF_vector_CNO",
  min.pct = 0.1,
  logfc.threshold = 0,   #
  test.use = "wilcox"
)

DEG_KOF_Gq_vs_vector_1_filtered <- DEG_KOF_Gq_vs_vector_1 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%                  # 
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%  # 
  arrange(desc(avg_log2FC))                       #


write.csv(DEG_KOF_Gq_vs_vector_1_filtered,
          file = "DEG_KOF_Gq_vs_vector_1.csv",
          row.names = FALSE)
#venn
# 安装一次即可
# install.packages("ggVennDiagram")

library(ggVennDiagram)
library(ggplot2)

# 提取基因向量
g1 <- unique(DEG_KOF_Gq_vs_vector_1_filtered$gene)
g2 <- unique(deg_filtered$gene)

gene_list <- list(
  KOF_Gq_vs_vector = g1,
  KOF_vs_WT        = g2
)

# Venn
p2 <- ggVennDiagram(
  gene_list,
  label_alpha = 0,
  label = "both", # 显示count和percent
  category.names = c("KOF_Gq vs Vector", "KOF vs WT"),
  edge_size = 0.2,
  set_size = 3.5,
  label_size = 5.5
) +
  # 经典学术配色
  scale_fill_gradient(
    low = "#FEF1F3",
    high = "#F9CCE0",
    name = "Gene count" # 添加图例名称
  ) +
  scale_color_manual(values = c("#FEF1F3", "#FEF1F3")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15),
    text = element_text(family = "Times New Roman") # 学术字体
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 0.5,
      barheight = 8,
      title.position = "top",
      title.hjust = 0.5
    )
  )

# all DEGs analyse
run_deg_by_celltype <- function(
    seu_obj,
    group.col    = "group",     # 分组列名
    celltype.col = "celltype",  # 细胞类型列名
    ident1,                     # 比较组1，如 "KOF"
    ident2,                     # 比较组2，如 "WT"
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox"
) {
  # 用 celltype 作为 Idents，方便用 subset.ident 限定细胞类型
  Idents(seu_obj) <- seu_obj[[celltype.col, drop = TRUE]]
  celltypes <- levels(Idents(seu_obj))
  
  deg_list <- list()
  
  for (ct in celltypes) {
    message("Processing cell type: ", ct)
    
    # 这个 celltype 内是否同时含有 ident1 和 ident2 两个 group
    meta_ct <- seu_obj@meta.data[Idents(seu_obj) == ct, ]
    groups_in_ct <- unique(meta_ct[[group.col]])
    if (!all(c(ident1, ident2) %in% groups_in_ct)) {
      message("  Skip: one of groups not present in ", ct)
      next
    }
    
    # 在该 celltype 内做 ident1 vs ident2
    deg_ct <- FindMarkers(
      object       = seu_obj,
      ident.1      = ident1,
      ident.2      = ident2,
      group.by     = group.col,
      subset.ident = ct,
      min.pct      = min.pct,
      logfc.threshold = logfc.threshold,
      test.use     = test.use
    )
    
    if (nrow(deg_ct) == 0) {
      message("  No DEGs found in ", ct)
      next
    }
    
    deg_ct$gene       <- rownames(deg_ct)
    deg_ct$celltype   <- ct
    deg_ct$comparison <- paste0(ident1, "_vs_", ident2)
    
    deg_list[[ct]] <- deg_ct
  }
  
  if (length(deg_list) == 0) {
    warning("No DEGs found for this comparison in any cell type.")
    return(NULL)
  }
  
  deg_all <- bind_rows(deg_list)
  rownames(deg_all) <- NULL
  return(deg_all)
}





DEG_KOF_vs_WT_allCelltypes <- run_deg_by_celltype(
  seu_obj = scRNA_HKU,
  group.col = "group",
  celltype.col = "celltype",
  ident1 = "KOF",
  ident2 = "WT"
)


DEG_Gq_vs_vector_allCelltypes <- run_deg_by_celltype(
  seu_obj = scRNA_HKU,
  group.col = "group",
  celltype.col = "celltype",
  ident1 = "KOF_Gq_CNO",
  ident2 = "KOF_vector_CNO"
)

head(DEG_Gq_vs_vector_allCelltypes)
head(DEG_KOF_vs_WT_allCelltypes)
table(DEG_KOF_vs_WT_allCelltypes$celltype)



# 合并所有对比
deg_all <- bind_rows(DEG_KOF_vs_WT_allCelltypes, DEG_Gq_vs_vector_allCelltypes)
# 注意列名：Seurat 版本不同可能叫 avg_log2FC / log2FC / avg_logFC
colnames(deg_all)

deg_filtered <- deg_all %>%
  filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25   # 如列名不同请改成你的实际列名
  )

# 按 comparison × celltype 统计 DEG 数量
deg_count <- deg_filtered %>%
  group_by(comparison, celltype) %>%
  summarise(n_deg = n(), .groups = "drop")

deg_count
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


library(ggplot2)
# 为了图更整齐，可以给 celltype 排个顺序（按总 DEG 数量排）
deg_count <- deg_count %>%
  group_by(celltype) %>%
  mutate(total_deg_celltype = sum(n_deg)) %>%
  ungroup()

deg_count$celltype <- factor(
  deg_count$celltype,
  levels = deg_count %>%
    arrange(desc(total_deg_celltype)) %>%
    pull(celltype) %>%
    unique()
)

ggplot(deg_count, aes(x = comparison, y = celltype, fill = n_deg)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(
    x = "Group comparison",
    y = "Cell type",
    fill = "Number of DEGs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

library(dplyr)

# 看一下列名，确认 logFC 那列叫什么
colnames(deg_all)

# 1）先按阈值过滤：显著且变化幅度超过阈值
deg_filtered_dir <- deg_all %>%
  filter(
    p_val_adj < 0.05,          # FDR 阈值，根据需要可调
    abs(avg_log2FC) > 0.25     # log2FC 阈值，按你自己的分析标准改
  ) %>%
  mutate(
    direction = ifelse(avg_log2FC > 0, "Up", "Down")  # >0 为上调，<0 为下调
  )

# 2）按「比较 × 细胞类型 × 方向」统计数量
deg_count_dir <- deg_filtered_dir %>%
  group_by(comparison, celltype, direction) %>%
  summarise(n_deg = n(), .groups = "drop")

deg_count_dir



library(ggplot2)

ggplot(deg_count_dir, aes(x = celltype, y = n_deg, fill = direction)) +
  geom_col(position = "dodge") +      # Up / Down 并列
  coord_flip() +                      # 细胞类型在 y 轴，更好读
  facet_wrap(~ comparison, ncol = 1) +# 每个 comparison 一行图
  scale_fill_manual(values = c(Up = "firebrick", Down = "steelblue")) +
  labs(
    x = "Cell type",
    y = "Number of DEGs",
    fill = "Direction"
  ) +
  theme_bw()


ggplot(deg_count_dir,
       aes(x = comparison, y = celltype, fill = n_deg)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient(low = "white", high = "darkred") +
  facet_wrap(~ direction) +   # 上调 / 下调 各一块
  labs(
    x = "Group comparison",
    y = "Cell type",
    fill = "Number of DEGs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

library(ggplot2)

ggplot(deg_count_dir,
       aes(x = comparison, y = celltype, fill = n_deg)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = n_deg), color = "black", size = 3) +  # 在格子上写数字
  scale_fill_gradient(low = "white", high = "#E83D36") +
  facet_wrap(~ direction) +   # Up / Down 分面
  labs(
    x = "Group comparison",
    y = "Cell type",
    fill = "Number of DEGs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )



