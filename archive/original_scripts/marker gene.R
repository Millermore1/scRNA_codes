
# 1. 生成 DotPlot 对象
dp <- DotPlot(scRNA_seurat, features = c(
  "CLDN5","PECAM1","ABCB1",
  "P2RY12","DOCK8","CSF1R",
  "GFAP","AQP4",
  "GAD1","GAD2","PVALB",
  "MOG","MBP","OLIG1",
  "GRIN1","CAMK2A","CBLN2","PRMT8"
), dot.scale=8)

# 2. 查看数据中的列名
colnames(dp$data)
# 通常会看到类似 'features.plot' 和 'id' 这样的列

# 3. 设置簇（id）列为因子，并定义顺序
desired_order <- c("Excitatory.neurons", "Inhibitory.neurons", "Oligodendrocytes",
                   "Astrocytes", "Microglia", "Endothelial.cells")
dp$data$id <- factor(dp$data$id, levels = desired_order)

# 4. 绘制图形，保持其他参数不变
library(ggplot2)
ggplot(dp$data, aes(x = features.plot, y = id, fill = pct.exp)) +
  geom_point(shape = 21, size = 8) +
  scale_fill_gradient(low = "lightblue", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=10),
    axis.text.y = element_text(size=10),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=10),
    plot.title = element_text(hjust=0.5, size=14, face="bold")
  ) +
  labs(
    title = "Gene Expression DotPlot",
    color = "Expression Level"
  )

