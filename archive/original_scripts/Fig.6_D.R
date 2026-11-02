# 载入需要的包
library(Seurat)
library(dplyr)
library(tibble)

# 1. 确保当前对象的默认聚类标签是你的细胞类型
Idents(scRNA_HKU) <- "celltype"

# 2. 提取所有唯一的细胞类型名称
cell_types <- unique(scRNA_HKU$celltype)
# 如果细胞类型中有 NA，需要剔除
cell_types <- cell_types[!is.na(cell_types)]

# 3. 创建一个空的列表，用来存放每种细胞类型的计算结果
deg_list <- list()

# 4. 循环遍历每种细胞类型
for (ct in cell_types) {
  message(paste0("正在计算: ", ct, " 中的 KOF vs WT 差异基因..."))
  
  # 使用 tryCatch 防止某个亚群因细胞数过少而报错中断循环
  tryCatch({
    # 计算差异基因
    # ident.1 = "KOF", ident.2 = "WT" 意味着：
    # log2FC > 0 表示在 KOF 组中表达上调，log2FC < 0 表示在 KOF 组中表达下调
    markers <- FindMarkers(scRNA_HKU, 
                           ident.1 = "KOF", 
                           ident.2 = "WT", 
                           group.by = "group",     # 指定分组依据的 metadata 列名
                           subset.ident = ct,      # 限定在当前的细胞类型中
                           logfc.threshold = 0.25, # log2FC 阈值，默认 0.25，可改为 0 以获取更多基因
                           min.pct = 0.1)          # 基因至少要在 10% 的细胞中表达
    
    # 5. 筛选具有统计学意义的基因 (p_val_adj < 0.05) 并标注细胞类型
    sig_markers <- markers %>% 
      rownames_to_column(var = "gene") %>% 
      filter(p_val_adj < 0.05) %>% 
      mutate(celltype = ct) %>%                    # 添加一列标注当前的细胞类型
      arrange(desc(avg_log2FC))                    # 按照 log2FC 降序排列，方便查看
    
    # 将结果存入列表
    deg_list[[ct]] <- sig_markers
    
  }, error = function(e) {
    message(paste0("跳过 ", ct, "，原因: ", e$message))
  })
}

# 6. 将列表中所有细胞类型的结果合并成一个大的数据框
all_sig_degs <- bind_rows(deg_list)

# 7. 查看结果的前几行
head(all_sig_degs)

# 8. 将结果保存为本地 CSV 文件，方便你在 Excel 中查阅或后续用于富集分析
write.csv(all_sig_degs, file = "KOF_vs_WT_Significant_DEGs_All_Celltypes.csv", row.names = FALSE)

message("所有细胞类型的差异基因计算完毕，并已保存为 CSV 文件！")
getwd()

#!/usr/bin/env Rscript

# GO/pathway network visualization for KOF vs WT single-cell DEGs.
# Input columns expected:
# gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, celltype

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(igraph)
  library(ggraph)
  library(scatterpie)
  library(ggforce)
  library(patchwork)
})

deg_file <- "D:/scrnadata/KOF_vs_WT_Significant_DEGs_All_Celltypes.csv"
out_pdf <- "C:/Users/\u9EC4\u5EC9/Documents/Codex/2026-08-02/new-chat/outputs/KOF_WT_GO_pathway_network.pdf"
out_png <- "C:/Users/\u9EC4\u5EC9/Documents/Codex/2026-08-02/new-chat/outputs/KOF_WT_GO_pathway_network.png"

deg <- readr::read_csv(deg_file, show_col_types = FALSE) %>%
  mutate(
    direction = if_else(avg_log2FC > 0, "Upregulated", "Downregulated"),
    gene = as.character(gene),
    celltype = as.character(celltype)
  )

celltype_levels <- c(
  "Excitatory neurons", "Inhibitory neurons", "Oligodendrocyte",
  "OPC", "Microglia", "Astrocyte", "Endothelial cells"
)

celltype_cols <- c(
  "Excitatory neurons" = "#3B73B9",
  "Inhibitory neurons" = "#D95F5F",
  "Oligodendrocyte" = "#5AB4AC",
  "OPC" = "#8C6BB1",
  "Microglia" = "#C49A45",
  "Astrocyte" = "#6A994E",
  "Endothelial cells" = "#E9C46A"
)

deg <- deg %>%
  mutate(celltype = factor(celltype, levels = celltype_levels))

run_go <- function(genes) {
  entrez <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = unique(genes),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  ) %>%
    unname() %>%
    discard(is.na)
  
  if (length(entrez) < 10) {
    return(NULL)
  }
  
  enrichGO(
    gene = entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
}

make_enrichment_table <- function(direction_name, top_n_per_celltype = 8) {
  deg %>%
    filter(direction == direction_name) %>%
    group_by(celltype) %>%
    summarise(genes = list(unique(gene)), .groups = "drop") %>%
    mutate(enrich = map(genes, run_go)) %>%
    mutate(result = map(enrich, ~ if (is.null(.x)) tibble() else as_tibble(.x))) %>%
    dplyr::select(celltype, result) %>%
    unnest(result) %>%
    group_by(celltype) %>%
    slice_min(p.adjust, n = top_n_per_celltype, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      term_id = ID,
      term = Description,
      gene_symbols = str_split(geneID, "/")
    )
}

term_gene_sets <- function(enrich_tbl) {
  enrich_tbl %>%
    dplyr::select(term_id, gene_symbols) %>%
    distinct(term_id, .keep_all = TRUE) %>%
    deframe()
}

make_edges <- function(enrich_tbl, min_jaccard = 0.18) {
  sets <- term_gene_sets(enrich_tbl)
  ids <- names(sets)
  if (length(ids) < 2) return(tibble(from = character(), to = character(), weight = numeric()))
  
  combn(ids, 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      a <- sets[[pair[1]]]
      b <- sets[[pair[2]]]
      j <- length(intersect(a, b)) / length(union(a, b))
      if (is.na(j) || j < min_jaccard) return(NULL)
      tibble(from = pair[1], to = pair[2], weight = j)
    })
}

make_nodes <- function(enrich_tbl) {
  cell_counts <- enrich_tbl %>%
    group_by(term_id, term, celltype) %>%
    summarise(n_genes = max(Count), min_padj = min(p.adjust), .groups = "drop")
  
  nodes <- cell_counts %>%
    group_by(term_id, term) %>%
    summarise(
      gene_count = max(n_genes),
      min_padj = min(min_padj),
      .groups = "drop"
    )
  
  pies <- cell_counts %>%
    mutate(celltype = as.character(celltype)) %>%
    dplyr::select(term_id, celltype, n_genes) %>%
    pivot_wider(names_from = celltype, values_from = n_genes, values_fill = 0)
  
  missing_celltypes <- setdiff(celltype_levels, colnames(pies))
  for (ct in missing_celltypes) {
    pies[[ct]] <- 0
  }
  
  left_join(nodes, pies, by = "term_id")
}

plot_direction_network <- function(enrich_tbl, title, panel_fill = "#F2D7D7") {
  nodes <- make_nodes(enrich_tbl)
  edges <- make_edges(enrich_tbl)
  
  if (nrow(nodes) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No enriched terms", size = 5) +
        labs(title = title) +
        theme_void(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
    )
  }
  
  vertex_tbl <- tibble(name = nodes$term_id)
  g <- graph_from_data_frame(edges, vertices = vertex_tbl, directed = FALSE)
  
  if (ecount(g) == 0) {
    g <- make_empty_graph(n = nrow(nodes), directed = FALSE)
    V(g)$name <- nodes$term_id
  }
  
  layout <- create_layout(g, layout = "fr")
  node_df <- as_tibble(layout) %>%
    left_join(nodes, by = c("name" = "term_id"))
  
  if (max(node_df$gene_count, na.rm = TRUE) == min(node_df$gene_count, na.rm = TRUE)) {
    node_df$radius <- 0.075
  } else {
    node_df$radius <- scales::rescale(node_df$gene_count, to = c(0.045, 0.105))
  }
  
  n_labels <- min(8, nrow(node_df))
  
  label_df <- node_df %>%
    arrange(min_padj) %>%
    slice_head(n = n_labels)
  
  p <- ggraph(layout)
  
  if (ecount(g) > 0) {
    p <- p +
      geom_edge_link(aes(width = weight), color = "grey72", alpha = 0.55, show.legend = FALSE)
  }
  
  p <- p +
    ggforce::geom_mark_ellipse(
      data = node_df,
      aes(x = x, y = y),
      expand = unit(5, "mm"),
      radius = unit(5, "mm"),
      fill = panel_fill,
      alpha = 0.23,
      color = NA,
      show.legend = FALSE
    ) +
    scatterpie::geom_scatterpie(
      data = node_df,
      aes(x = x, y = y, r = radius),
      cols = celltype_levels,
      color = "white",
      linewidth = 0.25
    ) +
    geom_text(
      data = label_df,
      aes(x = x, y = y + radius + 0.055, label = term),
      size = 3.1,
      lineheight = 0.92,
      color = "grey12"
    ) +
    scale_fill_manual(values = celltype_cols, drop = FALSE) +
    scale_edge_width(range = c(0.25, 1.2)) +
    coord_equal(clip = "off") +
    labs(title = title, fill = "Cell type") +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "bottom",
      plot.margin = margin(10, 18, 10, 18)
    )
  
  p
}

up_tbl <- make_enrichment_table("Upregulated", top_n_per_celltype = 8)
down_tbl <- make_enrichment_table("Downregulated", top_n_per_celltype = 8)

readr::write_csv(up_tbl, "C:/Users/\u9EC4\u5EC9/Documents/Codex/2026-08-02/new-chat/outputs/KOF_WT_GO_upregulated_terms.csv")
readr::write_csv(down_tbl, "C:/Users/\u9EC4\u5EC9/Documents/Codex/2026-08-02/new-chat/outputs/KOF_WT_GO_downregulated_terms.csv")

p_up <- plot_direction_network(
  up_tbl,
  "Enriched GO terms/pathways of upregulated DEGs (KOF/WT)",
  panel_fill = "#EFA3A5"
)

p_down <- plot_direction_network(
  down_tbl,
  "Enriched GO terms/pathways of downregulated DEGs (KOF/WT)",
  panel_fill = "#B7D3EA"
)

combined <- p_up + p_down + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(out_pdf, combined, width = 13.5, height = 8.2, device = cairo_pdf)
ggsave(out_png, combined, width = 13.5, height = 8.2, dpi = 320)












table(scRNA_HKU$group)
library(Seurat)
library(dplyr)
library(tibble)
Idents(scRNA_HKU) <- "celltype"
cell_types <- unique(scRNA_HKU$celltype)
cell_types <- cell_types[!is.na(cell_types)]
# 3. 创建一个空的列表，用来存放每种细胞类型的计算结果
deg_list <- list()

# 4. 循环遍历每种细胞类型
for (ct in cell_types) {
  message(paste0("正在计算: ", ct, " 中的 KOF vs WT 差异基因..."))
  
  # 使用 tryCatch 防止某个亚群因细胞数过少而报错中断循环
  tryCatch({
    # 计算差异基因
    # ident.1 = "KOF", ident.2 = "WT" 意味着：
    # log2FC > 0 表示在 KOF 组中表达上调，log2FC < 0 表示在 KOF 组中表达下调
    markers <- FindMarkers(scRNA_HKU, 
                           ident.1 = "KOF_vector_CNO", 
                           ident.2 = "KOF_Gq_CNO", 
                           group.by = "group",     # 指定分组依据的 metadata 列名
                           subset.ident = ct,      # 限定在当前的细胞类型中
                           logfc.threshold = 0.25, # log2FC 阈值，默认 0.25，可改为 0 以获取更多基因
                           min.pct = 0.1)          # 基因至少要在 10% 的细胞中表达
    
    # 5. 筛选具有统计学意义的基因 (p_val_adj < 0.05) 并标注细胞类型
    sig_markers <- markers %>% 
      rownames_to_column(var = "gene") %>% 
      filter(p_val_adj < 0.05) %>% 
      mutate(celltype = ct) %>%                    # 添加一列标注当前的细胞类型
      arrange(desc(avg_log2FC))                    # 按照 log2FC 降序排列，方便查看
    
    # 将结果存入列表
    deg_list[[ct]] <- sig_markers
    
  }, error = function(e) {
    message(paste0("跳过 ", ct, "，原因: ", e$message))
  })
}

# 6. 将列表中所有细胞类型的结果合并成一个大的数据框
all_sig_degs <- bind_rows(deg_list)

# 7. 查看结果的前几行
head(all_sig_degs)

# 8. 将结果保存为本地 CSV 文件，方便你在 Excel 中查阅或后续用于富集分析
write.csv(all_sig_degs, file = "KOF_vector_CNO_vs_KOF_Gq_CNO_Significant_DEGs_All_Celltypes.csv", row.names = FALSE)

message("所有细胞类型的差异基因计算完毕，并已保存为 CSV 文件！")
getwd()