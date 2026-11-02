
#好看的相关性热图
library(ComplexHeatmap)
library(circlize) # 用于生成连续颜色映射
# ---------------- 1. 模拟数据 ----------------
# 假设我们有一个 20x20 的相关性矩阵
set.seed(123)
mat <- matrix(rnorm(400, mean = 0, sd = 0.3), nrow = 20)
rownames(mat) <- paste0("Sample", 1:20)
colnames(mat) <- paste0("Sample", 1:20)
# 模拟注释数据框 (对应图左侧的 Patients, Assays, Type)
meta_df <- data.frame(
  Type = c(rep("MES-Hyp", 10), rep("MES-Ast", 10)),
  Assays = sample(c("10x 3' v2", "10x 3' v3", "SMART-Seq2"), 20, replace = TRUE),
  Patients = sample(c("P62", "P92", "BT338", "MGH105"), 20, replace = TRUE)
)
rownames(meta_df) <- rownames(mat)

# ---------------- 2. 定义颜色 ----------------
# 2.1 定义热图主体的连续颜色 (对应右上角的 PCC 图例: 蓝-白-红)
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("#3B75AF", "white", "#B72A2E")) 

# 2.2 定义注释条的离散颜色 (按需自定义)
ann_colors <- list(
  Type = c("MES-Hyp" = "#D51F8A", "MES-Ast" = "#517BCC"),
  Assays = c("10x 3' v2" = "#2C85B2", "10x 3' v3" = "#F28E2B", "SMART-Seq2" = "#59A14F")
  # Patients 颜色太多这里省略，实际中需要一一对应
)

# ---------------- 3. 创建注释对象 ----------------
# 创建左侧的行注释对象
left_ann <- rowAnnotation(
  df = meta_df,
  col = ann_colors,
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

# ---------------- 4. 绘制并分割热图 ----------------
# 注意图中有明显的白色十字分割线，这是通过 split 参数实现的
ht <- Heatmap(
  mat,
  name = "PCC",                  # 图例名称
  col = col_fun,                 # 颜色映射
  left_annotation = left_ann,    # 添加左侧注释
  row_split = meta_df$Type,      # 按 Type 将行切分为两块 (MES-Hyp 和 MES-Ast)
  column_split = meta_df$Type,   # 按 Type 将列也切分为两块
  cluster_rows = TRUE,           # 开启行聚类 (左侧树状图)
  cluster_columns = TRUE,        # 开启列聚类 (顶部树状图)
  show_row_names = FALSE,        # 隐藏行名 (图里没有样本名)
  show_column_names = FALSE,     # 隐藏列名
  row_title = NULL,
  column_title = NULL
)

# 绘制
draw(ht, merge_legend = TRUE)