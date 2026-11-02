# -------------------------------------------------------------------------
# 1. 加载必要的绘图包
# -------------------------------------------------------------------------
library(ggplot2)

# -------------------------------------------------------------------------
# 2. 构建原始数据集
# -------------------------------------------------------------------------
cellchat_data <- data.frame(
  Group = rep(c("KOF", "KOF_vector_CNO", "KOF_Gq_CNO", "WT"), each = 4),
  CellType = rep(c("PN", "PVALB", "VIP", "SST"), times = 4),
  Outgoing = c(
    5.14, 4.40, 4.46, 3.92,  # KOF 组
    4.66, 4.53, 4.24, 3.70,  # KOF_vector_CNO 组
    5.20, 4.42, 4.42, 3.82,  # KOF_Gq_CNO 组
    5.66, 4.80, 4.95, 4.80   # WT 组
  ),
  Incoming = c(
    4.98, 4.21, 4.46, 4.16,  # KOF 组
    4.30, 4.77, 4.11, 3.95,  # KOF_vector_CNO 组
    4.41, 4.80, 4.56, 3.88,  # KOF_Gq_CNO 组
    5.02, 5.35, 5.25, 4.35   # WT 组
  )
)

# 核心修改：创建一个组合文本标签，形如 "PN (WT)"，防止混淆
cellchat_data$Label <- paste0(cellchat_data$CellType, " (", cellchat_data$Group, ")")

# -------------------------------------------------------------------------
# 3. 定义专属映射（颜色配细胞，形状配分组）
# -------------------------------------------------------------------------
# 颜色保持原样（不加粗、高颜值学术色）
cell_colors <- c(
  "PN"    = "#DE3533",  # 红色
  "PVALB" = "#00A0E9",  # 蓝色
  "VIP"   = "#8957A1",  # 紫色
  "SST"   = "#78B360"   # 绿色
)

# 为 4 个组分配 4 种经典散点形状（圆形、三角形、正方形、菱形）
group_shapes <- c(
  "KOF"            = 16,  # 实心圆
  "KOF_vector_CNO" = 17,  # 实心正三角
  "KOF_Gq_CNO"     = 15,  # 实心正方形
  "WT"             = 18   # 实心菱起
)

# -------------------------------------------------------------------------
# 4. 开始绘制同坐标系叠加散点图
# -------------------------------------------------------------------------
ggplot(cellchat_data, aes(x = Outgoing, y = Incoming, color = CellType, shape = Group)) +
  # 绘制散点
  geom_point(size = 4.5, alpha = 0.85) +
  # 在点上方标注组合标签，字体设为普通（不加粗）
  geom_text(aes(label = Label), 
            vjust = -1.2, 
            hjust = 0.5, 
            size = 2.8, 
            color = "grey20",
            show.legend = FALSE) +
  # 映射颜色和形状
  scale_color_manual(values = cell_colors, name = "Cell Type") +
  scale_shape_manual(values = group_shapes, name = "Experimental Group") +
  # 设置坐标轴标签
  labs(
    x = "Outgoing interaction strength",
    y = "Incoming interaction strength",
    title = "CellChat Signaling Role Overlaid Scatter Plot"
  ) +
  # 
  expand_limits(y = max(cellchat_data$Incoming) * 1.05, x = max(cellchat_data$Outgoing) * 1.02) +
  # 极简经典边框主题
  theme_bw(base_size = 12) +
  # 样式微调：所有字体保持普通（不加粗）
  theme(
    plot.title = element_text(face = "plain", hjust = 0.5, size = 13, margin = margin(b=15)),
    axis.title.x = element_text(face = "plain", size = 11),
    axis.title.y = element_text(face = "plain", size = 11),
    axis.text = element_text(face = "plain", color = "black"),
    legend.text = element_text(face = "plain", size = 10),
    legend.title = element_text(face = "plain", size = 10),
    panel.grid.minor = element_blank() # 去除次要网格线
  )