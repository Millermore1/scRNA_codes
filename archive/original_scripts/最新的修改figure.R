library(Seurat)
library(dplyr)

# 1. 设置对比的组别名称
group_gq <- "KOF_Gq_CNO"
group_vec <- "KOF_vector_CNO"

# 2. 创建或确认联合列 (Cell Type + Group)
scRNA_HKU$celltype.group <- paste(scRNA_HKU$celltype, scRNA_HKU$group, sep = "_")
Idents(scRNA_HKU) <- "celltype.group"

# 3. 获取所有唯一的细胞类型
celltypes <- unique(scRNA_HKU$celltype)

# 4. 循环计算每一类细胞中的 Gq vs Vector 差异
gq_deg_list <- list()

for (ct in celltypes) {
  ident_1 <- paste(ct, group_gq, sep = "_")
  ident_2 <- paste(ct, group_vec, sep = "_")
  
  # 检查数据中是否存在这两组且细胞数是否足够
  if (ident_1 %in% Idents(scRNA_HKU) && ident_2 %in% Idents(scRNA_HKU)) {
    cells_1 <- sum(Idents(scRNA_HKU) == ident_1)
    cells_2 <- sum(Idents(scRNA_HKU) == ident_2)
    
    if (cells_1 >= 5 & cells_2 >= 5) {
      message(paste0("Computing DEGs for: ", ct))
      
      # 计算差异
      # ident.1 是 Gq 组，ident.2 是 Vector 组
      degs <- FindMarkers(scRNA_HKU, 
                          ident.1 = ident_1, 
                          ident.2 = ident_2, 
                          logfc.threshold = 0.1, # 建议设低一点，捕捉细微变化
                          min.pct = 0.1)
      
      if (nrow(degs) > 0) {
        degs$gene <- rownames(degs)
        degs$celltype <- ct
        gq_deg_list[[ct]] <- degs
      }
    }
  }
}

# 5. 合并并过滤显著基因 (p_adj < 0.05)
all_gq_degs <- bind_rows(gq_deg_list) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(celltype, desc(avg_log2FC))

# 6. 保存结果
write.csv(all_gq_degs, "DEGs_KOF_Gq_vs_Vector.csv", row.names = FALSE)

read.csv("C:/Users/黄廉/Desktop/GO_all_group.csv")
library(ggplot2)
library(dplyr) # 用于数据排序

# 1. 读取数据 (使用你提供的路径和内容)
# 注意：你的路径包含中文，建议在 R 中确认编码，或者将文件移至英文路径。
# 这里我直接使用你提供的数据框内容创建一个示例，方便你运行。
go_data <- data.frame(
  ID = c("nervous system development", "synapse organization", "synapse", 
         "neuron differentation", "cell junction", "brain development", 
         "MHC class I peptide loading complex", "generation of neurons"),
  P.value = c(5.120, 2.069, 4.920, 3.550, 3.419, 2.840, 2.779, 3.203)
)

go_data <- read.csv("C:/Users/黄廉/Desktop/GO_all_group.csv", encoding = "UTF-8")

# 2. **关键步骤**：按 P.value 对 ID 进行重新排序，保证柱子整齐
go_data <- go_data %>%
  mutate(ID = reorder(ID, P.value)) # 按 P.value 从小到大排序 ID

# 3. 绘制优美的横版柱状图
ggplot(go_data, aes(x = ID, y = P.value)) +
  # 绘制柱子，使用更高级、不刺眼的专业配色 (e.g., 深蓝色)
  geom_bar(stat = "identity", fill = "#F19ABC", width = 0.7) + 
  
  # 交换 X 轴和 Y 轴，使图表变为横版
  coord_flip() + 
  
  # 使用极简主题，去掉背景网格
  theme_classic() + 
  
  # 在柱子顶端添加 P 值数值标签，方便直接读取 (可选)
  geom_text(aes(label = sprintf("%.2f", P.value)), # 保留两位小数
            hjust = -0.2, # 标签微调到柱子外面
            size = 4, color = "black") + 
  
  # 设置坐标轴标签和标题
  labs(x = "", # 纵轴（排序后的ID）通常不需要标题，ID本身很清晰
       y = "-log10(P-value)", # 假设你的P.value已经是负对数处理过的
       title = "GO Enrichment: Top 8 Biological Processes") +
  
  # 进一步精细化主题设置
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # 标题居中加粗
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 14, color = "black", face = "plain"), # ID标签清晰
    axis.title.x = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.5), # 强化坐标轴线
    # 如果 ID 标签太长被截断，可以取消下面这一行的注释来扩大左侧边距
    # plot.margin = margin(l = 20, r = 20, t = 10, b = 10)
  ) +
  
  # 适当扩大 Y 轴（现在的横轴）范围，防止文字标签超出画布
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))