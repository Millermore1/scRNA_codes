library(dplyr)
library(ggplot2)
library(ggrepel)

#PV_DEGs caculate
head(inh)
PV <- subset(inh, subset = inhibit_celltype == "PV")
PV.DEG_KOF_Gq_vs_vector <- FindMarkers(
  object = pv,
  ident.1 = "WT",
  ident.2 = "KOF",
  min.pct = 0.1,
  logfc.threshold = 0,   #
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

# all inhibitory neurons
# 查看 group 列有哪些分组
table(inh$group)
# 用 group 作为当前身份（Idents）
Idents(inh) <- inh$group

deg_KOF_vs_CON_all <- FindMarkers(
  inh,
  ident.1 = "KOF",
  ident.2 = "WT",
  logfc.threshold = 0,   # 画火山图一般先不要太高阈值
  min.pct = 0.1
)
head(deg_KOF_vs_CON_all)

## 1. 把行名（基因名）变成一列
deg_df <- deg_KOF_vs_CON_all %>%
  tibble::rownames_to_column(var = "gene")

## 2. 设置阈值并标记上调/下调/不显著
logfc_cutoff <- 0.25   # log2FC 阈值，可根据需要调整
p_cutoff     <- 0.05   # FDR 阈值（用 p_val_adj）

deg_df <- deg_df %>%
  mutate(
    neg_log10_p = -log10(p_val_adj + 1e-300),  # 防止出现 log(0)
    sig = case_when(
      p_val_adj < p_cutoff & avg_log2FC >=  logfc_cutoff  ~ "Up",
      p_val_adj < p_cutoff & avg_log2FC <= -logfc_cutoff  ~ "Down",
      TRUE                                                ~ "NS"
    )
  )

## 3. p_volcano
p_volcano <- ggplot(deg_df, aes(x = avg_log2FC, y = neg_log10_p, color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
             linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff),
             linetype = "dashed", color = "black") +
  labs(x = "log2 Fold Change (KOF vs Control)",
       y = "-log10(adj p-value)",
       color = "Regulation",
       title = "Volcano plot: KOF vs Control") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(xlim = c(-3, 3))    # 这里调整 X 轴范围

genes_to_label <- c(
  "Auts2",   # 经典 ASD 基因
  "Tcf4",
  "Ctnnb1",
  "Chd7",
  "Zeb2",
  "Pcdh10",
  "Cdh10",
  "Dscam",
  "Dscaml1",
  "Kirrel3",
  "Ncam2",
  "Ntng1",
  "Cacna1c",
  "Cacna1d",
  "Scn3a",
  "Gabrg2",
  "Gad1",
  "Plp1",
  "Slc4a10",
  "Pgap1",
  "Dlx1",
  "Dclk2",
  "Cntnap4",
  "Sema5a",
  "Sema6a",
  "Pcdh7",
  "Cdh4",
  "Cpa6"     # 
)




library(tibble)

deg_out <- deg_KOF_vs_CON_all %>%
  rownames_to_column(var = "gene")

# 导出到当前工作目录
write.csv(deg_out,
          file = "deg_KOF_vs_CON_all.csv",
          row.names = FALSE)

p_volcano <- ggplot(deg_df, aes(x = avg_log2FC, y = neg_log10_p)) +
  # 1. 绘制散点 - 调整大小和透明度
  geom_point(aes(color = sig, fill = sig), 
             alpha = 0.7, 
             size = 2.5,
             shape = 21,           # 带边框的圆点
             stroke = 0.3) +       # 边框粗细
  
  # 2. 颜色方案 - 更专业的配色
  scale_color_manual(
    values = c("Up" = "#E64B35",     # 更柔和的红色
               "Down" = "#4DBBD5",   # 更柔和的蓝色
               "NS" = "#BDBDBD"),    # 中灰色
    name = "Expression"
  ) +
  scale_fill_manual(
    values = c("Up" = "#E64B35",
               "Down" = "#4DBBD5", 
               "NS" = "#BDBDBD"),
    name = "Expression"
  ) +
  
  # 3. 阈值线 - 更细的虚线
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
             linetype = "dashed", 
             color = "black",
             linewidth = 0.4,        # 更细的线
             alpha = 0.7) +
  geom_hline(yintercept = -log10(p_cutoff),
             linetype = "dashed", 
             color = "black",
             linewidth = 0.4,
             alpha = 0.7) +
  
  # 4. 添加显著基因标签（可选）
  # geom_text_repel(
  #   data = subset(deg_df, sig != "NS" & abs(avg_log2FC) > 1 & neg_log10_p > 10),
  #   aes(label = gene_symbol),      # 假设有基因名列
  #   size = 3,
  #   max.overlaps = 20,             # 最大标签数
  #   box.padding = 0.5,
  #   segment.color = "grey50",
  #   segment.size = 0.2,
  #   min.segment.length = 0
# ) +

# 5. 坐标轴和标签
labs(
  x = expression(log[2]("Fold Change")),  # 使用数学表达式
  y = expression(-log[10]("adjusted p-value")),
  title = "Volcano Plot: KOF vs Control",
  subtitle = paste("FC cutoff:", logfc_cutoff, "| p-value cutoff:", p_cutoff)
) +
  
  # 6. 主题美化
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",                    # 图例放在顶部
    legend.title = element_text(face = "bold"),
    legend.justification = "center",
    panel.grid.minor = element_blank(),         # 去掉次要网格线
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  
  # 7. 坐标轴范围
  coord_cartesian(
    xlim = c(-3, 3),
    ylim = c(0, max(deg_df$neg_log10_p, na.rm = TRUE) * 1.05)  # 自动调整y轴
  ) +
  
  # 8. 添加阈值区域阴影（可选）
  annotate("rect",
           xmin = -Inf, xmax = -logfc_cutoff,
           ymin = -log10(p_cutoff), ymax = Inf,
           fill = "#4DBBD5", alpha = 0.03) +
  annotate("rect",
           xmin = logfc_cutoff, xmax = Inf,
           ymin = -log10(p_cutoff), ymax = Inf,
           fill = "#E64B35", alpha = 0.03)

# 显示图形
print(p_volcano)

# 保存高质量图片
ggsave("volcano_plot_enhanced.pdf", 
       plot = p_volcano,
       width = 8, 
       height = 6,
       dpi = 300)


p_volcano <- ggplot(deg_df, aes(x = avg_log2FC, y = neg_log10_p)) +
  # 1. 散点
  geom_point(aes(color = sig, fill = sig), 
             alpha = 0.7, 
             size = 2.5,
             shape = 21,
             stroke = 0.3) +
  
  # 2. 颜色
  scale_color_manual(
    values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "#BDBDBD"),
    name = "Expression"
  ) +
  scale_fill_manual(
    values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "#BDBDBD"),
    name = "Expression"
  ) +
  
  # 3. 阈值线
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
             linetype = "dashed", color = "black",
             linewidth = 0.4, alpha = 0.7) +
  geom_hline(yintercept = -log10(p_cutoff),
             linetype = "dashed", color = "black",
             linewidth = 0.4, alpha = 0.7) +
  
  # 4. 只给指定基因加标签
  geom_text_repel(
    data = subset(deg_df, gene %in% genes_to_label),   # 注意这里用 gene
    aes(label = gene),                                 # 这里也用 gene
    size = 3,
    max.overlaps = 50,
    box.padding = 0.5,
    segment.color = "grey50",
    segment.size = 0.2,
    min.segment.length = 0
  ) +
  # 如果只想标显著差异的三个基因，用这一行替换上面那句 data：
  # data = subset(deg_df, gene %in% genes_to_label & sig != "NS"),
  
  # 5. 标签
  labs(
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("adjusted p-value")),
    title = "Volcano Plot: KOF vs Control",
    subtitle = paste("FC cutoff:", logfc_cutoff, "| p-value cutoff:", p_cutoff)
  ) +
  
  # 6. 主题
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.justification = "center",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  
  # 7. 坐标轴范围
  coord_cartesian(
    xlim = c(-3, 3),
    ylim = c(0, max(deg_df$neg_log10_p, na.rm = TRUE) * 1.05)
  ) +
  
  # 8. 阈值阴影
  annotate("rect",
           xmin = -Inf, xmax = -logfc_cutoff,
           ymin = -log10(p_cutoff), ymax = Inf,
           fill = "#4DBBD5", alpha = 0.03) +
  annotate("rect",
           xmin = logfc_cutoff, xmax = Inf,
           ymin = -log10(p_cutoff), ymax = Inf,
           fill = "#E64B35", alpha = 0.03)

print(p_volcano)


library(dplyr)
library(tibble)

deg_filtered <- deg_KOF_vs_CON_all %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%                  # 
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%  # 
  arrange(desc(avg_log2FC))                       #


write.csv(deg_filtered,
          file = "deg_KOF_vs_CON_all_filtered_sorted.csv",
          row.names = FALSE)
head(DEG_KOF_Gq_vs_vector_1_filtered)
head(deg_filtered)




