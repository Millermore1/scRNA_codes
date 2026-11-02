#WGCNA
library(Seurat)
library(hdWGCNA)
library(WGCNA)
# 1.Inhibitory.neurons
seurat_inh <- subset(scRNA_seurat, celltype == "Inhibitory.neurons")

# 2. 初始化 hdWGCNA 对象
seurat_inh <- SetupForWGCNA(
  seurat_inh,
  gene_select = "fraction", # 过滤掉在极少数细胞中表达的基因 (Dropout 严重的基因)
  fraction = 0.05,          # 要求基因至少在 5% 的抑制性神经元中表达
  wgcna_name = "Inhibitory_Network"
)

# 3. 构建元细胞 (Metacells) - 已修复
seurat_inh <- MetacellsByGroups(
  seurat_obj = seurat_inh,
  group.by = c("celltype", "orig.ident"), # 同时输入细胞类型和样本来源
  ident.group = "celltype",               # 指定细胞类型标签是哪一列
  k = 25,             
  max_shared = 10     
)

# 4. 归一化元细胞表达矩阵
seurat_inh <- NormalizeMetacells(seurat_inh)

# 提取构建网络所需的数据
seurat_inh <- SetDatExpr(
  seurat_inh,
  group_name = "Inhibitory.neurons", 
  group.by = "celltype"
)

# 测试不同的软阈值
seurat_inh <- TestSoftPowers(seurat_inh, networkType = 'signed')

# 绘制软阈值图 (观察 R-squared 拐点)
PlotSoftPowers(seurat_inh)


# 假设你在上图中观察到最佳软阈值是 9 (⚠️ 请必须将 9 替换为你实际观察到的数值)
seurat_inh <- ConstructNetwork(
  seurat_inh, 
  softPower = 7, 
  setDatExpr = FALSE,
  tom_name = 'Inhib_TOM' 
)

# 绘制模块树状图 (看看基因被分成了多少个颜色模块)
PlotDendrogram(seurat_inh, main = 'Inhibitory Neurons Module Dendrogram')

# 计算每个模块的特征基因 (Module Eigengenes)
seurat_inh <- ModuleEigengenes(seurat_inh)

# 获取每个模块排名前 100 的核心基因
hub_df <- GetHubGenes(seurat_inh, n_hubs = 100)

# 搜索你的目标基因
target_gene <- "PRMT8" # 把这里换成你的基因，比如 "GAD1" 或 "PV" 相关的基因
target_info <- subset(hub_df, gene_name == target_gene)
print(target_info)

# hdWGCNA 需要数值型的表型数据。我们把 group 列变成 0 和 1
# 假设你的 group 列里疾病叫 "PD"，对照组叫 "Control" (请根据实际对照组名称调整)
traits_df <- data.frame(
  PD_Status = ifelse(seurat_inh$group == "PD", 1, 0) 
)
# 将行名与细胞对齐
rownames(traits_df) <- colnames(seurat_inh)
# 将数值型表型加入 metadata
seurat_inh <- AddMetaData(seurat_inh, metadata = traits_df)

# 计算模块与 PD 状态的相关性
seurat_inh <- ModuleTraitCorrelation(
  seurat_inh,
  traits = "PD_Status", 
  group.by = "celltype"
)

# 绘制最终的模块-表型相关性热图
PlotModuleTraitCorrelation(seurat_inh)

# 1. 计算模块的特征基因 (Module Eigengenes)
seurat_inh <- ModuleEigengenes(seurat_inh)

# 2. 计算每个基因与所在模块的连接度得分 (kME) —— 这是解决报错的关键！
seurat_inh <- ModuleConnectivity(
  seurat_inh,
  group.by = "celltype", 
  group_name = "Inhibitory.neurons"
)

# ==========================================
# 3. 再次运行你的画图代码
# ==========================================
target_gene <- "PRMT8" 
target_module <- "#F09BA0" 

ModuleNetworkPlot(
  seurat_inh,
  outdir = 'Hub_Network_Plots', 
  mods = target_module,         
  n_inner = 10,                 
  n_outer = 20,                 
  n_conns = 300,                
  label_center = TRUE,
  vertex.label.cex = 1.2        
)

ModuleNetworkPlot(
  seurat_inh,
  outdir = 'Hub_Network_Plots', 
  mods = target_module,         
  n_inner = 10,                 
  n_outer = 20,                 
  n_conns = 300,                
  
  # === 美化参数 ===
  label_center = FALSE,          # 强烈建议关闭！中间那个模块名字往往会挡住最核心的基因
  vertex.label.cex = 1.0,        # 稍微调小一点字体，防止 30 个基因互相重叠
  vertex.label.font = 3,         # [核心美化] 将字体设置为斜体 3，符合生物学论文中基因名的书写规范
  vertex.label.color = "grey10", # 将纯黑字体改为深灰，看起来更柔和
  vertex.size = 2,               # 设置基础节点的大小
  edge.alpha = 0.3               # [核心美化] 将连线设置为 30% 透明度，解决“黑乎乎一坨毛线”的问题
)


# 1. 自动提取 PRMT8 所在的准确模块名称 (不再手动输入)
modules_df <- GetModules(seurat_inh)
correct_module <- as.character(modules_df$module[modules_df$gene_name == "PRMT8"])

# 在控制台打印一下，看看它到底叫什么名字
print(paste("PRMT8 所在的真正模块名是:", correct_module))

# 2. 运行美化版的网络拓扑图
ModuleNetworkPlot(
  seurat_inh,
  outdir = 'Hub_Network_Plots', 
  mods = correct_module,         # <--- 直接使用自动提取的准确名称
  n_inner = 10,                 
  n_outer = 20,                 
  n_conns = 300,                
  
  # === 顺手加上的美化参数 ===
  label_center = FALSE,          # 关闭中心模块标签，防止挡住核心基因
  vertex.label.cex = 1.0,        # 字体大小适中
  vertex.label.font = 3,         # 基因名设置为符合发表规范的斜体
  vertex.label.color = "grey10", # 柔和的深灰字体
  edge.alpha = 0.3               # 连线半透明，避免黑糊糊一团
)

library(ggraph)
library(tidygraph)
library(igraph)

# === 1. 提取 PRMT8 朋友圈的底层连线数据 ===
# 提取蓝模块的前30个核心基因
hub_df <- GetHubGenes(seurat_inh, n_hubs = 30)
mod_hubs <- hub_df$gene_name[hub_df$module == correct_module]

# 确保主角 PRMT8 一定在名单里
if(!("PRMT8" %in% mod_hubs)) mod_hubs <- c("PRMT8", mod_hubs)[1:30]

# 从底层矩阵切出这 30 个基因的连接权重
TOM <- GetTOM(seurat_inh)
subTOM <- as.matrix(TOM[mod_hubs, mod_hubs])

# === 2. 构建画图所需的网络对象 ===
net_graph <- graph_from_adjacency_matrix(subTOM, mode = "undirected", weighted = TRUE, diag = FALSE)

tg <- as_tbl_graph(net_graph) %>%
  activate(edges) %>% 
  filter(weight > 0.05) %>% # 砍掉太弱的连线，保持画面清爽
  activate(nodes) %>%
  mutate(
    node_degree = centrality_degree(weights = weight),
    gene_group = ifelse(name == "PRMT8", "PRMT8", "Other") # 单独标记 PRMT8
  )

# === 3. 开始画图 ===
p <- ggraph(tg, layout = 'stress') + 
  geom_edge_link(aes(edge_width = weight), edge_colour = "grey80", edge_alpha = 0.5) +
  geom_node_point(aes(size = node_degree, color = gene_group), alpha = 0.9) +
  scale_size_continuous(range = c(3, 10), guide = "none") + 
  scale_color_manual(values = c("PRMT8" = "#E64B35", "Other" = "#4DBBD5"), guide = "none") +
  geom_node_text(aes(label = name, size = node_degree), repel = TRUE, max.overlaps = 50) +
  theme_void()

# 这一句是魔法！它会绕过 PDF 强制输出，直接把图打在你的 RStudio 右下角！
print(p)

# 1. 提取所有模块排名前 50 的 Hub 基因 (你可以把 50 改成你想要的数字)
hub_df <- GetHubGenes(seurat_inh, n_hubs = 50)

# 2. 专门筛选出 "blue" 模块的基因信息
blue_hubs_info <- subset(hub_df, module == "blue")

# 3. 按 kME (模块内连接度得分) 从高到低严谨排序
blue_hubs_info <- blue_hubs_info[order(blue_hubs_info$kME, decreasing = TRUE), ]

# 4. 纯净版：只提取基因名称并打印在控制台
blue_gene_names <- blue_hubs_info$gene_name
print("=== Blue 模块 Top 50 核心基因名单 ===")
print(blue_gene_names)

# 5. 实用版：将包含基因名和 kME 跑分的完整表格导出到你的电脑上
write.csv(blue_hubs_info, file = "PRMT8_Blue_Module_HubGenes.csv", row.names = FALSE)
print("✅ 完整表格已保存为: PRMT8_Blue_Module_HubGenes.csv (请在当前工作目录下查看)")



# 1. 安装并加载专门读取 Excel 的包 (如果还没装过的话)
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

# 2. 设定文件路径 (注意：我已经帮你把路径里的 \ 全部改成了 /)
file_path <- "F:/sn_RNAseq_data/public_database/correlations_PRMT8_brain_region.xlsx"

# 3. 读取 Excel 数据
# read_excel 默认会读取表格中的第一个 Sheet (工作表)
prmt8_cor_data <- read_excel(file_path)

# 4. 查看数据的前 6 行，检查是否读取成功
head(prmt8_cor_data)
















