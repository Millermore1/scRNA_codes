# =============================================================================
# R/utils/plotting.R
#
# Figure-building helpers shared by the pipeline modules. These wrap the
# ggplot2 patterns that were copy-pasted across the original scripts (UMAP
# panels, ordered dot plots, volcano plots, composition bars, marker heat maps)
# so that every figure in the paper uses the same theme, point sizes and fonts.
# =============================================================================

# -----------------------------------------------------------------------------
# Theme
# -----------------------------------------------------------------------------

#' Minimal publication theme used for all ggplot figures.
theme_pub <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.text       = ggplot2::element_text(colour = "black"),
      axis.title      = ggplot2::element_text(colour = "black"),
      axis.line       = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      axis.ticks      = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text      = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.key.size = grid::unit(0.35, "cm"),
      plot.title      = ggplot2::element_text(face = "bold", hjust = 0,
                                              size = base_size + 1),
      plot.subtitle   = ggplot2::element_text(colour = "grey30",
                                              size = base_size - 1)
    )
}

#' Theme for the axis-less UMAP panels.
theme_umap <- function(...) {
  ggplot2::theme_void(...) +
    ggplot2::theme(
      axis.line            = ggplot2::element_line(
        colour = "black", linewidth = 0.3,
        arrow = grid::arrow(length = grid::unit(0.08, "cm"))
      ),
      plot.margin          = ggplot2::margin(2, 2, 2, 2),
      legend.position      = "none",
      panel.background     = ggplot2::element_blank()
    )
}

# -----------------------------------------------------------------------------
# Dimensional reduction
# -----------------------------------------------------------------------------

#' Extract UMAP coordinates + metadata as a plain data frame.
umap_data <- function(obj, reduction = "umap") {
  emb <- Seurat::Embeddings(obj, reduction = reduction)[, 1:2, drop = FALSE]
  df <- data.frame(emb, obj[[]], check.names = FALSE)
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df
}

#' Single UMAP coloured by a metadata column.
plot_umap <- function(obj, colour_by, palette = NULL, title = NULL,
                      label = FALSE, point_size = 0.1, reduction = "umap") {
  df <- umap_data(obj, reduction)
  levels_present <- unique(as.character(df[[colour_by]]))

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$UMAP_1, y = .data$UMAP_2,
    colour = .data[[colour_by]]
  )) +
    ggplot2::geom_point(size = point_size, stroke = 0) +
    ggplot2::ggtitle(title) +
    theme_umap()

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_colour_manual(
      values = pal_for(levels_present, palette), name = colour_by
    )
  }
  p
}

#' UMAP split into one facet per group, all panels sharing one colour scale.
#'
#' This is the "cell type composition across conditions" panel: one facet per
#' experimental group, identical embedding coordinates and colours so that
#' groups can be compared by eye.
plot_umap_facets <- function(obj, colour_by, facet_by,
                             palette = NULL,
                             point_size = 0.03,
                             facet_fill = "#e6bac5",
                             ncol = NULL,
                             reduction = "umap") {
  df <- umap_data(obj, reduction)
  levels_present <- unique(as.character(df[[colour_by]]))

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data$UMAP_1, y = .data$UMAP_2,
    colour = .data[[colour_by]]
  )) +
    ggplot2::geom_point(size = point_size, shape = 16, stroke = 0) +
    ggplot2::facet_wrap(
      stats::as.formula(paste("~", facet_by)),
      ncol = ncol, scales = "fixed"
    ) +
    ggplot2::coord_fixed() +
    theme_umap() +
    ggplot2::theme(
      panel.border  = ggplot2::element_rect(fill = NA, colour = "grey80",
                                            linewidth = 0.3),
      strip.text    = ggplot2::element_text(size = 9, face = "bold"),
      axis.line     = ggplot2::element_blank(),
      axis.ticks    = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_colour_manual(
      values = pal_for(levels_present, palette), name = colour_by
    )
  }
  p
}

#' Expression of one gene on the UMAP, with a sequential colour scale.
plot_feature <- function(obj, gene, palette = c("#E0E0E0", "#3A7CA5", "#0D3B66"),
                         limits = NULL, point_size = 0.3, title = NULL,
                         reduction = "umap") {
  args <- list(object = obj, features = gene, reduction = reduction,
               pt.size = point_size, order = TRUE)
  if (!is.null(limits)) {
    args$min.cutoff <- limits[1]
    args$max.cutoff <- limits[2]
  }

  p <- do.call(Seurat::FeaturePlot, args) +
    ggplot2::scale_colour_gradientn(colours = palette, name = gene) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
  p
}

# -----------------------------------------------------------------------------
# Dot plots
# -----------------------------------------------------------------------------

#' Dot plot with clusters ordered by a summary expression score.
#'
#' Reproduces the ordering trick used for the excitatory neuron dot plots: the
#' clusters are sorted by their summed scaled expression across a panel of
#' marker genes, so the gradient from layer 2/3 to layer 6 reads left to right.
dotplot_ordered <- function(obj, features, group.by,
                            order_by = features,
                            high_at = c("top", "bottom"),
                            cols = c("#EDEDED", "#B40426"),
                            dot.scale = 8, dot.min = 0.05,
                            angle = 45, cluster_col = NULL) {
  high_at <- match.arg(high_at)

  p <- Seurat::DotPlot(obj, features = features, group.by = group.by,
                       dot.min = 0, dot.scale = dot.scale)

  dat <- p$data
  score <- dat[dat$features.plot %in% order_by, , drop = FALSE]
  score <- stats::aggregate(avg.exp.scaled ~ id, data = score,
                            FUN = sum, na.rm = TRUE)
  lv <- score$id[order(score$avg.exp.scaled, decreasing = TRUE)]
  if (high_at == "top") lv <- rev(lv)   # last factor level is drawn at the top

  dat$id <- factor(as.character(dat$id), levels = lv)

  ggplot2::ggplot(dat, ggplot2::aes(
    x = .data$features.plot, y = .data$id,
    size = .data$pct.exp, colour = .data$avg.exp.scaled
  )) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_gradient(low = cols[1], high = cols[2],
                                   name = "Average\nexpression") +
    ggplot2::scale_size(range = c(0, 6), name = "% expressing") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = angle, hjust = 1,
                                           colour = "black", face = "italic"),
      axis.text.y  = ggplot2::element_text(colour = "black"),
      panel.grid.major = ggplot2::element_line(colour = "grey92",
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# -----------------------------------------------------------------------------
# Volcano plot
# -----------------------------------------------------------------------------

#' Volcano plot of one DE comparison.
#'
#' @param deg A data frame containing `gene`, `avg_log2FC` and `p_val_adj`.
#' @param label_genes Genes to label (only those present in `deg` are drawn).
plot_volcano <- function(deg, label_genes = NULL, title = NULL,
                         palette = pal_direction,
                         fc_cutoff = cfg$params$de$fc_cutoff,
                         padj_cutoff = cfg$params$de$padj_cutoff,
                         xlim = c(-3, 3)) {
  df <- deg
  df$neg_log10_p <- -log10(df$p_val_adj + 1e-300)
  df$sig <- ifelse(
    df$p_val_adj < padj_cutoff & df$avg_log2FC >=  fc_cutoff, "Up",
    ifelse(df$p_val_adj < padj_cutoff & df$avg_log2FC <= -fc_cutoff,
           "Down", "NS")
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$avg_log2FC,
                                        y = .data$neg_log10_p)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$sig), size = 1.1,
                        alpha = 0.75) +
    ggplot2::scale_colour_manual(values = palette, name = "Regulation") +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
                        linetype = "dashed", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff),
                        linetype = "dashed", linewidth = 0.35) +
    ggplot2::labs(x = expression(log[2]~"fold change"),
                  y = expression(-log[10]~"adjusted p-value"),
                  title = title) +
    ggplot2::coord_cartesian(xlim = xlim) +
    theme_pub()

  if (!is.null(label_genes)) {
    lab <- df[df$gene %in% label_genes, , drop = FALSE]
    if (nrow(lab) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = lab, ggplot2::aes(label = .data$gene),
        size = 2.8, max.overlaps = 50, min.segment.length = 0,
        segment.colour = "grey60", segment.size = 0.2, seed = cfg$params$seed
      )
    }
  }
  p
}

# -----------------------------------------------------------------------------
# Composition and summary panels
# -----------------------------------------------------------------------------

#' Stacked cell-type proportion bar plot, one bar per group.
plot_composition <- function(comp, palette = NULL, levels = NULL,
                             ylab = "Cell proportion", xlab = NULL,
                             y_percent = TRUE) {
  df <- comp
  if (!is.null(levels)) {
    df$group <- ordered_factor(df$group, levels)
  }
  ct_levels <- if (!is.null(palette)) names(palette) else unique(df$celltype)
  df$celltype <- ordered_factor(df$celltype, ct_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$proportion,
                                        fill = .data$celltype)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::labs(x = xlab, y = ylab, fill = "Cell type") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (y_percent) {
    p <- p + ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0)
    )
  }
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(
      values = pal_for(unique(as.character(df$celltype)), palette)
    )
  }
  p
}

#' Tile plot of the number of DEGs per comparison and cell type.
plot_deg_tile <- function(counts, fill_high = "#E83D36",
                          facet_by_direction = TRUE) {
  df <- counts
  df$celltype <- ordered_factor(df$celltype, unique(df$celltype))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$comparison,
                                        y = .data$celltype,
                                        fill = .data$n_deg)) +
    ggplot2::geom_tile(colour = "grey85") +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_deg), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = fill_high,
                                 name = "DEGs") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_pub() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid  = ggplot2::element_blank()
    )

  if (facet_by_direction && length(unique(df$direction)) > 1) {
    p <- p + ggplot2::facet_wrap(~ direction)
  }
  p
}

# -----------------------------------------------------------------------------
# Marker heat map
# -----------------------------------------------------------------------------

#' ComplexHeatmap of cell-type marker genes (averaged expression, Z-scored).
#'
#' @param obj Seurat object with `celltype` in the metadata.
#' @param top_markers Data frame with columns `cluster`, `gene`, `avg_log2FC`.
#' @return A ComplexHeatmap `Heatmap` object (call `draw()` on it).
build_marker_heatmap <- function(obj, top_markers, celltype_col = "celltype",
                                 palette = NULL,
                                 col_fun = circlize::colorRamp2(
                                   c(-2, 0, 2), c(col_zscore[["low"]],
                                                  col_zscore[["mid"]],
                                                  col_zscore[["high"]])
                                 )) {
  avg <- Seurat::AverageExpression(
    obj, assays = "RNA", features = unique(top_markers$gene),
    group.by = celltype_col, slot = "data"
  )$RNA
  avg <- as.matrix(avg)
  avg <- avg[unique(top_markers$gene), , drop = FALSE]

  mat <- t(scale(t(avg)))
  mat[is.na(mat)] <- 0

  celltypes <- colnames(mat)
  group_cols <- if (is.null(palette)) {
    stats::setNames(recycle_colours(length(celltypes)), celltypes)
  } else {
    pal_for(celltypes, palette)
  }

  row_group <- factor(top_markers$cluster[match(rownames(mat), top_markers$gene)],
                      levels = celltypes)

  ComplexHeatmap::Heatmap(
    mat,
    name = "Z-score",
    col  = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    row_split = row_group, row_title = NULL, row_gap = grid::unit(1, "mm"),
    column_split = factor(celltypes, levels = celltypes),
    row_names_gp = grid::gpar(fontsize = 6, fontface = "italic"),
    column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    column_names_rot = 45,
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Cell type` = celltypes,
      col = list(`Cell type` = group_cols),
      annotation_name_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      simple_anno_size = grid::unit(3, "mm")
    ),
    border = FALSE, use_raster = TRUE, raster_quality = 3
  )
}

# -----------------------------------------------------------------------------
# Enrichment bar plot
# -----------------------------------------------------------------------------

#' GO/KEGG bar figure with the term name and member genes written inside.
#'
#' Expects the layout produced by `enrich_to_bar_df()`: one row per term with
#' `score`, `Description`, `gene_label`, `y`, `term_y` and `gene_y`.
plot_go_bar_with_genes <- function(bar_df, fill = "#FFDE7B",
                                   x_label = "Fold enrichment",
                                   base_family = "sans") {
  if (is.null(bar_df) || nrow(bar_df) == 0) return(NULL)
  x_max <- max(bar_df$score) * 1.08

  ggplot2::ggplot(bar_df) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = .data$score,
                                    ymin = .data$y + 0.03,
                                    ymax = .data$y + 0.47),
                       fill = fill, colour = NA) +
    ggplot2::geom_text(ggplot2::aes(x = 0.05, y = .data$term_y,
                                    label = .data$Description),
                       hjust = 0, size = 3.1, fontface = "bold",
                       family = base_family, colour = "black") +
    ggplot2::geom_text(ggplot2::aes(x = 0.05, y = .data$gene_y,
                                    label = .data$gene_label),
                       hjust = 0, size = 2.9, fontface = "italic",
                       family = base_family, colour = "grey20") +
    ggplot2::scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.5, max(bar_df$y) + 0.65),
                                expand = c(0, 0)) +
    ggplot2::labs(x = x_label, y = NULL) +
    ggplot2::theme_classic(base_family = base_family, base_size = 10) +
    ggplot2::theme(
      axis.line.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.line.x  = ggplot2::element_line(linewidth = 0.4),
      axis.text.x  = ggplot2::element_text(colour = "black")
    )
}

# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

#' Panel label for a figure (a, b, c ...).
panel_tag <- function(x) {
  ggplot2::labs(tag = x) +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))
}
