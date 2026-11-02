# =============================================================================
# R/02_annotation.R
#
# Cell-type annotation and the panels that document it:
#   * cluster marker detection (saved as a table for the supplement),
#   * manual cluster -> cell type assignment (mappings in config/config.R),
#   * marker dot plots for both species,
#   * UMAP of cell types, overall and split by experimental group,
#   * marker heat map of the top markers per cell type.
#
# Outputs
#   data/processed/human_atlas_annotated.rds
#   data/processed/mouse_ko_annotated.rds
#   results/tables/02_markers_<dataset>.csv
#   results/figures/02_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "dplyr", "ggplot2", "patchwork", "ComplexHeatmap", "circlize",
     "RColorBrewer")

# -----------------------------------------------------------------------------
# 1. Cluster markers
# -----------------------------------------------------------------------------

#' Find positive markers for every cluster and keep the strongest per cluster.
find_cluster_markers <- function(obj, ident_col = "seurat_clusters",
                                 only_pos = TRUE, min_pct = 0.25,
                                 logfc = 0.25, top_n = 50) {
  obj <- ensure_joined_layers(obj)
  Seurat::Idents(obj) <- ident_col

  markers <- Seurat::FindAllMarkers(
    obj, only.pos = only_pos, min.pct = min_pct,
    logfc.threshold = logfc, test.use = "wilcox", verbose = FALSE
  )
  if (nrow(markers) == 0) return(markers)

  markers$gene <- rownames(markers)
  markers <- markers[order(markers$cluster, -markers$avg_log2FC), , drop = FALSE]

  keep <- stats::ave(seq_len(nrow(markers)), markers$cluster,
                     FUN = function(i) seq_along(i) <= top_n)
  markers[keep == 1, , drop = FALSE]
}

#' Save a marker table for manual inspection before annotation.
write_marker_table <- function(markers, dataset) {
  save_table(markers, paste0("02_markers_", dataset), digits = 4)
}

# -----------------------------------------------------------------------------
# 2. Annotation
# -----------------------------------------------------------------------------

#' Annotate the human atlas with the manual cluster mapping.
annotate_human_atlas <- function(obj,
                                 mapping = cfg$design$human_cluster_mapping,
                                 cluster_col = NULL,
                                 new_col = "celltype") {
  if (is.null(cluster_col)) {
    cluster_col <- if (all(grepl("^[0-9]+$", Seurat::Idents(obj)))) {
      NULL
    } else {
      "seurat_clusters"
    }
  }

  obj <- relabel_clusters(obj, mapping, new_col = new_col,
                          cluster_col = cluster_col)

  unknown <- setdiff(unique(as.character(obj[[new_col, drop = TRUE]])),
                     cfg$design$human_celltypes)
  if (length(unknown) > 0) {
    warning("Labels not in cfg$design$human_celltypes: ",
            paste(unknown, collapse = ", "), call. = FALSE)
  }
  obj$celltype <- as.character(obj[[new_col, drop = TRUE]])
  obj$celltype <- ordered_factor(obj$celltype, cfg$design$human_celltypes)
  obj
}

#' Validate the mouse annotation that arrives with the processed object.
validate_mouse_annotation <- function(obj) {
  if (!"celltype" %in% colnames(obj[[]])) {
    stop("Mouse object has no 'celltype' column. Annotate it first or supply ",
         "an object that already carries the annotation.", call. = FALSE)
  }
  obj$celltype <- ordered_factor(obj$celltype, cfg$design$mouse_celltypes)
  if (any(is.na(obj$celltype))) {
    warning(sum(is.na(obj$celltype)), " cells have an unknown cell type.",
            call. = FALSE)
  }
  obj
}

# -----------------------------------------------------------------------------
# 3. Figures
# -----------------------------------------------------------------------------

#' Marker dot plot with a fixed cell-type order.
plot_marker_dotplot <- function(obj, features, celltype_col = "celltype",
                                levels = NULL, palette = NULL,
                                title = NULL) {
  obj <- ensure_joined_layers(obj)
  Seurat::Idents(obj) <- celltype_col

  df <- obj[[]]
  if (!is.null(levels)) {
    keep <- levels[levels %in% unique(as.character(df[[celltype_col]]))]
    obj <- subset(obj, cells = colnames(obj)[
      as.character(df[[celltype_col]]) %in% keep
    ])
    Seurat::Idents(obj) <- celltype_col
    obj[[celltype_col]] <- ordered_factor(obj[[celltype_col, drop = TRUE]],
                                          levels)
    Seurat::Idents(obj) <- celltype_col
  }

  features <- intersect(features, rownames(obj))
  Seurat::DotPlot(obj, features = features, dot.scale = 7) +
    ggplot2::scale_colour_gradient(low = "#D9E4EC", high = "#B40426",
                                   name = "Average\nexpression") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                          colour = "black", face = "italic"),
      axis.text.y = ggplot2::element_text(colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "grey92",
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0)
    )
}

#' Draw and save the marker heat map (average expression, row Z-score).
save_marker_heatmap <- function(obj, markers, dataset,
                                celltype_col = "celltype",
                                top_n = cfg$params$annotation$marker_heatmap_top_n,
                                palette = NULL,
                                width = 6, height = 7) {
  sig <- markers[
    markers$p_val_adj < cfg$params$de$padj_cutoff &
      markers$avg_log2FC > cfg$params$annotation$marker_heatmap_min_fc,
    , drop = FALSE
  ]
  if (nrow(sig) == 0) {
    warning("No markers passed the heat-map thresholds; skipping.", call. = FALSE)
    return(invisible(NULL))
  }

  top <- do.call(rbind, lapply(split(sig, sig$cluster), function(d) {
    utils::head(d[order(-d$avg_log2FC), , drop = FALSE], top_n)
  }))
  top$cluster <- as.character(top$cluster)

  ht <- build_marker_heatmap(
    obj,
    data.frame(gene = top$gene, cluster = top$cluster,
               avg_log2FC = top$avg_log2FC),
    celltype_col = celltype_col,
    palette = palette
  )

  save_graphics(paste0("02_marker_heatmap_", dataset), width, height,
                ComplexHeatmap::draw(ht, padding = grid::unit(c(2, 2, 2, 2), "mm")))
  invisible(ht)
}

# -----------------------------------------------------------------------------
# 4. Run: human atlas
# -----------------------------------------------------------------------------

step_header("Human atlas annotation")

human <- read_seurat(cfg$inputs$human_atlas_rds, "human atlas")

human_markers <- find_cluster_markers(human, "seurat_clusters")
write_marker_table(human_markers, "human")

human <- annotate_human_atlas(human)
write_object(human, "human_atlas_annotated")

save_figure(
  plot_umap(human, "celltype", pal_celltype_human,
            title = "Human brain snRNA-seq atlas"),
  "02_umap_celltype_human", width = 5.2, height = 4.4
)

save_figure(
  plot_umap_facets(human, "celltype", "group", pal_celltype_human, ncol = 5),
  "02_umap_celltype_by_group_human", width = 11, height = 2.8
)

save_figure(
  plot_marker_dotplot(human, cfg$params$annotation$human_markers,
                      "celltype", cfg$design$human_celltypes,
                      title = "Marker genes of major cell classes"),
  "02_dotplot_markers_human", width = 6.5, height = 3.2
)

save_figure(
  plot_feature(human, "PRMT8", limits = c(0, 2), point_size = 0.05,
               title = "PRMT8 expression"),
  "02_umap_prmt8_human", width = 4.6, height = 4
)

save_marker_heatmap(human, human_markers, "human",
                    palette = pal_celltype_human,
                    width = 6, height = 6.5)

# -----------------------------------------------------------------------------
# 5. Run: mouse atlas
# -----------------------------------------------------------------------------

step_header("Mouse atlas annotation")

mouse <- read_seurat(cfg$inputs$mouse_atlas_rds, "mouse atlas")
mouse <- validate_mouse_annotation(mouse)
write_object(mouse, "mouse_ko_annotated")

save_figure(
  plot_umap(mouse, "celltype", pal_celltype_mouse,
            title = "Prmt8 cKO mouse snRNA-seq atlas"),
  "02_umap_celltype_mouse", width = 5.2, height = 4.4
)

save_figure(
  plot_umap_facets(mouse, "celltype", "group", pal_celltype_mouse, ncol = 4),
  "02_umap_celltype_by_group_mouse", width = 10, height = 2.8
)

save_figure(
  plot_marker_dotplot(mouse, cfg$params$annotation$mouse_markers,
                      "celltype", cfg$design$mouse_celltypes,
                      title = "Marker genes of major cell classes (mouse)"),
  "02_dotplot_markers_mouse", width = 7, height = 3.4
)

mouse_markers <- find_cluster_markers(mouse, "celltype", top_n = 30)
write_marker_table(mouse_markers, "mouse")

note("done: annotation objects and figures written to results/")
