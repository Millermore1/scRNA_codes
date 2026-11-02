# 载入必要的包
library(ggplot2)
library(dplyr)

# 计算每个group中每个celltype的数量
celltype_counts <- scRNA_seurat@meta.data %>%
  group_by(group, celltype) %>%
  summarise(count = n()) %>%
  ungroup()

# 计算每个group的总细胞数
group_totals <- celltype_counts %>%
  group_by(group) %>%
  summarise(total = sum(count))

# 合并，计算比例
celltype_proportions <- left_join(celltype_counts, group_totals, by = "group") %>%
  mutate(proportion = count / total)

# 查看比例数据
print(celltype_proportions)

# 绘制堆积比例柱状图
ggplot(celltype_proportions, aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "不同group中各细胞类型的比例",
       x = "Group",
       y = "比例",
       fill = "细胞类型") +
  theme_minimal()

color <- c("#919ac2","#ffac98","#70a4c8","#a5a9af","#63917d","#dbd1b4","#6e729a","#9ba4bd","#c5ae5f","#b9b8d6")

ggplot(celltype_proportions, aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = color) +
  labs(title = "cell.proportions",
       x = "Group",
       y = "proportions",
       fill = "cell_type") +
  theme_minimal()



ggplot(celltype_proportions, aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +  # 添加黑色边框
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = color) +
  labs(title = "cell.proportions",
       x = "Group",
       y = "proportions",
       fill = "cell_type") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # 取消背景网格线
    panel.border = element_blank() # 可以取消边框线（如果有的话）
  )



