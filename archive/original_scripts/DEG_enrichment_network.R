######################################################################
## DEG 富集分析网络图（复现参考图样式）
## - 上调 / 下调 DEGs 分别做 GO BP 富集（clusterProfiler, mouse）
## - 网络节点 = 富集通路，节点大小 = 基因数（Gene hits）
## - 节点饼图 = 该通路命中基因在各细胞类型中的构成
## - 边 = 通路间共享基因（Jaccard 相似度）
######################################################################

## ---- 0. 依赖包（首次运行先安装） -----------------------------------
## install.packages(c("dplyr","tidyr","stringr","igraph","ggplot2",
##                    "ggrepel","ggforce","ggnewscale","patchwork","scatterpie"))
## if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## BiocManager::install(c("clusterProfiler", "org.Mm.eg.db"))

## 注意加载顺序：AnnotationDbi（clusterProfiler 依赖）会遮蔽 dplyr::select，
## BiocGenerics 会遮蔽 intersect/union —— 所以 dplyr 必须最后加载。
suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(ggnewscale)
  library(patchwork)
  library(scatterpie)
  library(stringr)
  library(tidyr)
  library(dplyr)
})

## ---- 1. 参数 --------------------------------------------------------

## Windows 路径推荐使用 /，不要直接用单个反斜杠 \
deg_file <- "D:/scrnadata/KOF_vector_CNO_vs_KOF_Gq_CNO_Significant_DEGs_All_Celltypes.csv"

top_n      <- 10
edge_jac   <- 0.05
min_shared <- 1
seed       <- 42

node_r     <- c(0.035, 0.085)

up_col     <- "#A65628"
down_col   <- "#1B9E8E"

out_prefix <- "KOF_vector_CNO_vs_KOF_Gq_CNO_DEG_enrichment_network"

## 细胞类型配色（按数据中的实际名称对应）
celltype_cols <- c(
  "Excitatory neurons" = "#E64B35",
  "Inhibitory neurons" = "#4DBBD5",
  "Oligodendrocyte"    = "#00A087",
  "OPC"                = "#3C5488",
  "Microglia"          = "#F39B7F",
  "Astrocyte"          = "#8491B4",
  "Endothelial cells"  = "#91D1C2"
)

## ---- 2. 读取 DEG 数据 ------------------------------------------------

deg <- read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)

## 检查必要列是否存在
stopifnot(all(c("gene", "avg_log2FC", "celltype") %in% colnames(deg)))

## 根据 avg_log2FC 判断上下调
deg$direction <- ifelse(deg$avg_log2FC > 0, "up", "down")

## ---- 3. GO BP 富集 --------------------------------------------------

run_enrich <- function(genes) {
  
  eg <- enrichGO(
    gene          = unique(genes),
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  if (is.null(eg) || nrow(as.data.frame(eg)) == 0) {
    return(NULL)
  }
  
  ## 去冗余；如果失败则跳过
  eg <- tryCatch(
    simplify(
      eg,
      cutoff = 0.6,
      by = "p.adjust",
      select_fun = min
    ),
    error = function(e) eg
  )
  
  as.data.frame(eg)
}

## ---- 4. 构建节点 / 边 ------------------------------------------------

norm01 <- function(v, to = c(0.1, 0.9)) {
  r <- range(v, na.rm = TRUE)
  if (diff(r) == 0) {
    return(rep(mean(to), length(v)))
  }
  (v - r[1]) / diff(r) * diff(to) + to[1]
}

build_graph_data <- function(enrich_df, deg_sub) {
  
  top <- enrich_df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n)
  
  term_genes <- strsplit(top$geneID, "/")
  names(term_genes) <- top$ID
  
  ## 节点：每个通路命中的基因在各细胞类型的计数（饼图构成）
  gene_ct <- deg_sub %>%
    distinct(gene, celltype)
  
  ct_names <- names(celltype_cols)
  
  nodes <- lapply(seq_len(nrow(top)), function(i) {
    
    g <- term_genes[[i]]
    
    tab <- gene_ct %>%
      filter(gene %in% g) %>%
      count(celltype) %>%
      tibble::deframe()
    
    row <- as.list(setNames(integer(length(ct_names)), ct_names))
    
    for (ct in names(tab)) {
      if (ct %in% ct_names) {
        row[[ct]] <- tab[[ct]]
      }
    }
    
    tibble(
      ID       = top$ID[i],
      Term     = top$Description[i],
      hits     = length(unique(g)),
      p.adjust = top$p.adjust[i],
      !!!row
    )
  }) %>%
    bind_rows()
  
  ## 边：共享基因 + Jaccard 相似度
  ids <- top$ID
  
  edges <- do.call(
    rbind,
    lapply(combn(ids, 2, simplify = FALSE), function(pr) {
      
      a <- term_genes[[pr[1]]]
      b <- term_genes[[pr[2]]]
      
      shared <- length(base::intersect(a, b))
      jac <- shared / length(base::union(a, b))
      
      if (shared >= min_shared && jac >= edge_jac) {
        data.frame(
          from   = pr[1],
          to     = pr[2],
          weight = jac
        )
      } else {
        NULL
      }
    })
  )
  
  list(
    nodes      = nodes,
    edges      = edges,
    term_genes = term_genes
  )
}

## ---- 5. 单个面板绘图 -------------------------------------------------

plot_panel <- function(gdata, direction_label) {
  
  nodes <- gdata$nodes
  
  edge_df <- if (is.null(gdata$edges)) {
    data.frame(from = character(0), to = character(0))
  } else {
    gdata$edges
  }
  
  g <- graph_from_data_frame(
    d = edge_df,
    vertices = nodes$ID,
    directed = FALSE
  )
  
  set.seed(seed)
  
  xy <- if (length(E(g)) > 0) {
    layout_with_fr(g, weights = E(g)$weight)
  } else {
    layout_in_circle(g)
  }
  
  nodes$x <- norm01(xy[, 1])
  nodes$y <- norm01(xy[, 2])
  nodes$r <- norm01(sqrt(nodes$hits), to = node_r)
  nodes$label <- str_wrap(nodes$Term, width = 20)
  
  ed <- NULL
  
  if (!is.null(gdata$edges)) {
    ed <- gdata$edges %>%
      left_join(
        nodes %>% dplyr::select(ID, x, y),
        by = c("from" = "ID")
      ) %>%
      left_join(
        nodes %>% dplyr::select(ID, xend = x, yend = y),
        by = c("to" = "ID")
      )
  }
  
  ## Gene hits 尺寸图例
  h_leg <- c(max(nodes$hits), min(nodes$hits))
  r_leg <- norm01(sqrt(h_leg), to = node_r)
  
  xL <- 0.02
  
  leg <- tibble(
    h  = h_leg,
    r  = r_leg,
    y  = 0.02 + r_leg,
    x0 = xL,
    x1 = xL + r_leg + 0.055
  )
  
  ## 图例专用虚拟点
  leg_dummy <- tibble(
    celltype = factor(names(celltype_cols), levels = names(celltype_cols)),
    x = -5,
    y = -5
  )
  
  ggplot() +
    {
      if (!is.null(ed)) {
        geom_segment(
          data = ed,
          aes(x = x, y = y, xend = xend, yend = yend),
          color = "grey65",
          linewidth = 0.45,
          alpha = 0.8
        )
      }
    } +
    geom_scatterpie(
      data = nodes,
      aes(x = x, y = y, r = r),
      cols = names(celltype_cols),
      color = NA,
      alpha = 0.95
    ) +
    scale_fill_manual(values = celltype_cols, guide = "none") +
    new_scale_fill() +
    geom_point(
      data = leg_dummy,
      aes(x = x, y = y, fill = celltype),
      shape = 21,
      size = 3,
      color = "white",
      stroke = 0.3
    ) +
    scale_fill_manual(
      name = "Cell type",
      values = celltype_cols
    ) +
    guides(
      fill = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 3)
      )
    ) +
    ## Gene hits 图例
    geom_circle(
      data = leg,
      aes(x0 = x0, y0 = y, r = r),
      color = "grey20",
      linewidth = 0.5
    ) +
    geom_segment(
      data = leg,
      aes(x = x0 + r, y = y, xend = x1, yend = y),
      color = "grey20",
      linewidth = 0.5
    ) +
    geom_text(
      data = leg,
      aes(x = x1 + 0.008, y = y, label = h),
      hjust = 0,
      size = 3.2
    ) +
    annotate(
      "text",
      x = xL,
      y = max(leg$y + leg$r) + 0.045,
      label = "Gene hits",
      hjust = 0,
      size = 3.4,
      fontface = "bold"
    ) +
    ## 通路标签
    geom_text_repel(
      data = nodes,
      aes(x = x, y = y, label = label),
      size = 3.2,
      fontface = "bold",
      lineheight = 0.95,
      segment.color = "grey50",
      segment.size = 0.3,
      min.segment.length = 0,
      max.overlaps = Inf,
      seed = seed
    ) +
    coord_fixed(
      xlim = c(-0.03, 1.03),
      ylim = c(-0.03, 1.03),
      expand = FALSE,
      clip = "off"
    ) +
    labs(title = direction_label) +
    theme_void(base_family = "sans") +
    theme(
      panel.background = element_rect(
        fill = "#F4F6F8",
        color = "grey35",
        linewidth = 0.6
      ),
      plot.title       = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      plot.margin      = margin(10, 10, 6, 10)
    )
}

## 标题色条
banner <- function(txt, col) {
  
  ggplot() +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      fill = col
    ) +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = txt,
      color = "white",
      fontface = "bold",
      size = 4.8
    ) +
    coord_cartesian(
      xlim = c(-1, 1),
      ylim = c(-1, 1),
      expand = FALSE
    ) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

## ---- 6. 运行：上调 / 下调 -------------------------------------------

res <- list()

for (d in c("up", "down")) {
  
  sub <- deg %>%
    filter(direction == d)
  
  message(sprintf("[%s] unique genes: %d", d, length(unique(sub$gene))))
  
  en <- run_enrich(sub$gene)
  
  if (is.null(en)) {
    message(sprintf("[%s] no significant terms", d))
    next
  }
  
  write.csv(
    en,
    sprintf("%s_%s_GO_BP.csv", out_prefix, d),
    row.names = FALSE
  )
  
  gd <- build_graph_data(en, sub)
  
  res[[d]] <- plot_panel(gd, d)
}

stopifnot(length(res) > 0)

mk <- function(d, col, label) {
  banner(label, col) / res[[d]] + plot_layout(heights = c(1, 13))
}

p_all <- NULL

if ("up" %in% names(res)) {
  p_all <- mk(
    "up",
    up_col,
    "Enrichment analysis of upregulated DEGs"
  )
}

if ("down" %in% names(res)) {
  
  pd <- mk(
    "down",
    down_col,
    "Enrichment analysis of downregulated DEGs"
  )
  
  p_all <- if (is.null(p_all)) {
    pd
  } else {
    p_all | pd
  }
}

p_all <- p_all +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

## ---- 7. 保存 ---------------------------------------------------------

ggsave(
  filename = sprintf("%s.png", out_prefix),
  plot = p_all,
  width = 13,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = sprintf("%s.pdf", out_prefix),
  plot = p_all,
  width = 13,
  height = 7,
  device = cairo_pdf,
  bg = "white"
)

message("Done: ", out_prefix, ".png / .pdf")
