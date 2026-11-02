# =============================================================================
# R/10_figures_heatmaps.R
#
# Expression heat maps and the PRMT8 correlation panel:
#   * mean expression of curated gene sets (ion channels, GABA/glutamate
#     receptors, synaptic genes) per cell class and per disease cohort,
#   * a cell class x cohort heat map of PRMT8 expression,
#   * the PRMT8 brain-region correlation matrix (public dataset),
#   * PRMT8 feature plots with quantile-based colour limits.
#
# Outputs
#   results/tables/10_*.csv
#   results/figures/10_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "dplyr", "ggplot2", "pheatmap", "ComplexHeatmap", "circlize",
     "readxl", "RColorBrewer")

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------

#' Mean expression matrix: genes (rows) x groups (columns).
mean_expression_matrix <- function(obj, genes, group_col = "group",
                                   celltype_col = NULL, celltype = NULL) {
  if (!is.null(celltype) && !is.null(celltype_col)) {
    obj <- subset(obj, cells = colnames(obj)[
      as.character(obj[[celltype_col, drop = TRUE]]) == celltype
    ])
  }
  genes <- intersect(genes, rownames(obj))
  if (length(genes) == 0) return(NULL)

  mat <- Seurat::AverageExpression(obj, assays = "RNA", features = genes,
                                   group.by = group_col, slot = "data")$RNA
  as.matrix(mat)
}

#' Row Z-score, clipped at +/- 2.5 to keep the colour scale readable.
zscore_rows <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z[z > 2.5]  <- 2.5
  z[z < -2.5] <- -2.5
  z
}

#' Heat map of a gene set across groups, drawn with pheatmap.
draw_gene_set_heatmap <- function(mat, zscore = TRUE,
                                  colours = c("#2166AC", "white", "#B2182B"),
                                  title = NULL, cluster_rows = TRUE,
                                  main_cex = 0.9) {
  if (is.null(mat) || nrow(mat) == 0) return(invisible(NULL))

  m <- if (zscore) zscore_rows(mat) else as.matrix(mat)
  pheatmap::pheatmap(
    m,
    cluster_rows = cluster_rows,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(colours)(100),
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 9,
    main = title,
    border_color = NA,
    cellwidth = 26,
    cellheight = 12,
    silent = FALSE
  )
}

# -----------------------------------------------------------------------------
# 2. Gene-set heat maps (human atlas)
# -----------------------------------------------------------------------------

step_header("Gene-set heat maps")

human <- read_seurat(p_processed("human_atlas_annotated.rds"),
                     "annotated human atlas")

gene_sets <- list(
  ion_channels   = cfg$genes$ion_channels,
  gaba_receptors = cfg$genes$gaba_receptors,
  synaptic_prmt8 = c("NLGN2", "NLGN3", "NRXN1", "PRMT8", "BDNF", "GAD1",
                     "PVALB", "MYH10", "ARID1B", "CHD8", "UBE3A")
)

for (gs in names(gene_sets)) {
  mat <- mean_expression_matrix(human, gene_sets[[gs]], "group")
  if (is.null(mat)) {
    warning("No genes of '", gs, "' present in the human object.", call. = FALSE)
    next
  }
  save_table(as.data.frame(as.table(mat), responseName = "mean_expression"),
             paste0("10_heatmap_", gs, "_human"))

  save_graphics(paste0("10_heatmap_", gs, "_human"), 4.2, 4.8, {
    draw_gene_set_heatmap(mat, title = gs)
  })
}

# -----------------------------------------------------------------------------
# 3. PRMT8 per cell class and cohort
# -----------------------------------------------------------------------------

step_header("PRMT8 heat map")

prmt8_mat <- tryCatch({
  genes <- intersect("PRMT8", rownames(human))
  if (length(genes) == 0) stop("PRMT8 absent")

  avg <- Seurat::AverageExpression(human, assays = "RNA", features = genes,
                                   group.by = c("celltype", "group"),
                                   slot = "data")$RNA
  as.matrix(avg)
}, error = function(e) {
  note("could not compute the PRMT8 matrix: ", conditionMessage(e))
  NULL
})

if (!is.null(prmt8_mat)) {
  save_table(as.data.frame(as.table(prmt8_mat), responseName = "mean_expression"),
             "10_prmt8_by_celltype_group_human", digits = 4)
  save_graphics("10_prmt8_heatmap_human", 4, 2.4, {
    pheatmap::pheatmap(prmt8_mat, cluster_rows = FALSE, cluster_cols = FALSE,
                       display_numbers = TRUE, number_format = "%.2f",
                       fontsize = 8, cellwidth = 26, cellheight = 16,
                       color = grDevices::colorRampPalette(c("white", "#205EA7"))(60),
                       border_color = NA,
                       main = "PRMT8 mean expression")
  })
}

# -----------------------------------------------------------------------------
# 4. PRMT8 feature plots
# -----------------------------------------------------------------------------

if ("PRMT8" %in% rownames(human)) {
  save_figure(
    plot_feature(human, "PRMT8", limits = c("q1", "q20"), point_size = 0.05,
                 title = "PRMT8 expression across the human atlas"),
    "10_prmt8_featureplot_human", width = 4.6, height = 4
  )
}

# -----------------------------------------------------------------------------
# 5. PRMT8 brain-region correlation matrix
# -----------------------------------------------------------------------------

step_header("PRMT8 brain-region correlation")

cor_file <- cfg$inputs$brain_cor_xlsx
if (file.exists(cor_file)) {
  cor_df <- readxl::read_excel(cor_file)
  cor_df <- as.data.frame(cor_df)

  if (!"genes" %in% colnames(cor_df)) {
    warning("Expected a 'genes' column in ", basename(cor_file),
            " -- columns found: ", paste(colnames(cor_df), collapse = ", "),
            call. = FALSE)
  } else {
    cor_matrix <- cor_df |>
      tibble::column_to_rownames(var = "genes") |>
      as.matrix()

    save_table(cor_df, "10_prmt8_brain_region_correlation", digits = 4)

    col_fun <- circlize::colorRamp2(c(-1, 0, 1),
                                    c("#008080", "white", "#DC143C"))

    row_ha <- ComplexHeatmap::rowAnnotation(
      Distribution = ComplexHeatmap::anno_boxplot(
        cor_matrix, outline = FALSE, box_width = 0.6,
        gp = grid::gpar(fill = "#E0E0E0", col = "grey40")
      ),
      annotation_name_gp = grid::gpar(fontsize = 9, fontface = "italic")
    )

    ht <- ComplexHeatmap::Heatmap(
      cor_matrix,
      name = "Correlation",
      col  = col_fun,
      row_km = 2, column_km = 2,
      row_title = "Genes", column_title = "Brain regions",
      cluster_rows = TRUE, cluster_columns = TRUE,
      row_names_gp = grid::gpar(fontsize = 9, fontface = "italic"),
      column_names_gp = grid::gpar(fontsize = 9),
      column_names_rot = 45,
      rect_gp = grid::gpar(col = "white", lwd = 1.2),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(
          sprintf("%.2f", cor_matrix[i, j]), x, y,
          gp = grid::gpar(
            fontsize = 7,
            col = ifelse(abs(cor_matrix[i, j]) > 0.6, "white", "black")
          )
        )
      },
      right_annotation = row_ha
    )

    save_graphics("10_prmt8_brain_region_correlation", 7, 6,
                  ComplexHeatmap::draw(
                    ht, padding = grid::unit(c(3, 3, 3, 3), "mm")
                  ))
  }
} else {
  note("brain-region correlation file not found -- skipping (",
       cor_file, ")")
}

# -----------------------------------------------------------------------------
# 6. Ordered marker dot plot (human atlas)
# -----------------------------------------------------------------------------

step_header("Ordered marker dot plot")

dot_order <- c("Excitatory.neurons", "Inhibitory.neurons", "Oligodendrocytes",
               "Astrocytes", "Microglia", "Endothelial.cells")
human$celltype <- ordered_factor(human$celltype, dot_order)

save_figure(
  plot_marker_dotplot(human, cfg$params$annotation$human_markers, "celltype",
                      levels = dot_order,
                      title = "Major cell-class markers (ordered)"),
  "10_dotplot_markers_ordered_human", width = 6.5, height = 3.2
)

note("done: heat maps and correlation panel complete")
