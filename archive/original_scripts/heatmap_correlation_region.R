library(ComplexHeatmap)
library(circlize)
library(tibble)
library(dplyr)

# === 1. 数据准备 (保持不变) ===
cor_matrix <- prmt8_cor_data %>%
  column_to_rownames(var = "genes") %>%
  as.matrix()

# === 2. 调配现代单细胞流行色 (孔雀蓝 -> 纯白 -> 宝石红) ===
col_fun <- colorRamp2(
  c(-1, 0, 1), 
  c("#008080", "white", "#DC143C")
)

# === 3. 增加高级边缘注释 (在热图右侧画一个箱线图，展示该基因在所有脑区的相关性分布) ===
row_ha <- rowAnnotation(
  Distribution = anno_boxplot(
    cor_matrix, 
    outline = FALSE,            # 不画离群点
    box_width = 0.6,
    gp = gpar(fill = "#E0E0E0", col = "grey40") # 优雅的灰色箱线图
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "italic")
)

# === 4. 绘制满血版高定热图 ===
ht <- Heatmap(
  cor_matrix,
  name = "Correlation", 
  col = col_fun,
  
  # -- 核心升级 1：K-means 智能区块切割 --
  row_km = 2,                # 将基因切分为 2 个功能簇
  column_km = 2,             # 将脑区切分为 2 个特征簇
  row_title = "Gene Clusters",
  column_title = "Brain Regions",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # -- 聚类树 --
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  # -- 字体与边框 --
  row_names_gp = gpar(fontsize = 12, fontface = "italic", col = "black"), 
  column_names_gp = gpar(fontsize = 12, fontface = "plain", col = "black"),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 2), # 加粗白线，让格子更有立体感
  
  # -- 核心升级 2：优雅地显示内部数值 --
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 如果数值绝对值 > 0.5，用白色字体（在深色背景上），否则用黑色字体
    text_col <- ifelse(abs(cor_matrix[i, j]) > 0.6, "white", "black")
    grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, 
              gp = gpar(fontsize = 9, col = text_col))
  },
  
  # -- 绑定刚才做的边缘注释 --
  right_annotation = row_ha
)

# === 5. 渲染出图 ===
draw(ht, padding = unit(c(4, 4, 4, 4), "mm"))