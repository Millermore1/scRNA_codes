# =============================================================================
# R/06_excitatory_neurons.R
#
# Excitatory-neuron subclustering, used for the layer-resolved figures:
#   * re-cluster excitatory neurons and annotate cortical layers
#     (L2/3, L4, L5, L6) from the layer marker panels,
#   * order clusters by their combined layer-marker expression,
#   * PRMT8 distribution across clusters,
#   * layer composition across the experimental groups,
#   * layer-specific differential expression.
#
# Outputs
#   data/processed/mouse_excitatory.rds
#   results/tables/06_*.csv
#   results/figures/06_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "dplyr", "ggplot2", "patchwork", "openxlsx")

# -----------------------------------------------------------------------------
# 1. Subset and re-cluster
# -----------------------------------------------------------------------------

step_header("Mouse excitatory neurons: subclustering")

mouse <- read_seurat(p_processed("mouse_ko_annotated.rds"), "annotated mouse atlas")

excit <- subset(mouse, cells = colnames(mouse)[
  as.character(mouse$celltype) == "Excitatory neurons"
])
note("excitatory neurons: ", ncol(excit), " cells")

excit <- run_seurat_standard(
  excit,
  nfeatures  = cfg$params$excitatory$nfeatures,
  npcs       = cfg$params$excitatory$npcs,
  dims       = cfg$params$excitatory$dims,
  resolution = cfg$params$excitatory$resolution
)

cluster_col <- paste0("RNA_snn_res.", cfg$params$excitatory$resolution)
save_figure(
  plot_umap(excit, cluster_col,
            palette = setNames(recycle_colours(30), as.character(0:29)),
            title = "Excitatory neuron subclusters"),
  "06_umap_excitatory_clusters_mouse", width = 4.6, height = 4
)

# -----------------------------------------------------------------------------
# 2. Layer markers and annotation
# -----------------------------------------------------------------------------

layer_panel <- unique(unlist(cfg$genes$excit_layers))
layer_panel <- intersect(layer_panel, rownames(excit))

save_figure(
  plot_marker_dotplot(excit, layer_panel, cluster_col,
                      title = "Cortical layer markers"),
  "06_dotplot_layer_markers_mouse", width = 5.5, height = 3.2
)

step_header("Annotate cortical layers")

excit <- relabel_clusters(excit, cfg$design$mouse_excitatory_mapping,
                          new_col = "subtype", cluster_col = cluster_col)
excit$subtype <- ordered_factor(excit$subtype, names(pal_excit_subtype))

save_figure(
  plot_umap(excit, "subtype", pal_excit_subtype,
            title = "Excitatory neuron layers"),
  "06_umap_excitatory_subtypes_mouse", width = 4.4, height = 3.8
)

# The excitatory dot plot is ordered by combined layer-marker expression so the
# gradient L2/3 -> L6 reads along one axis.
save_figure(
  dotplot_ordered(excit, intersect(cfg$genes$excit_ordering, rownames(excit)),
                  group.by = "subtype",
                  order_by = intersect(cfg$genes$excit_ordering,
                                       rownames(excit)),
                  high_at = "top"),
  "06_dotplot_excitatory_ordered_mouse", width = 5.4, height = 3.2
)

# -----------------------------------------------------------------------------
# 3. PRMT8 across layers
# -----------------------------------------------------------------------------

if ("Prmt8" %in% rownames(excit)) {
  prmt8_df <- Seurat::FetchData(excit, vars = c("Prmt8", "subtype", "group"))
  names(prmt8_df)[1] <- "expr"

  save_figure(
    ggplot2::ggplot(prmt8_df, ggplot2::aes(x = .data$subtype, y = .data$expr,
                                           fill = .data$group)) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
      ggplot2::geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1,
                            linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = pal_group_mouse) +
      ggplot2::labs(x = NULL, y = "Prmt8 expression") +
      theme_pub(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
    "06_prmt8_by_layer_mouse", width = 5, height = 3.4
  )
}

# -----------------------------------------------------------------------------
# 4. Layer composition
# -----------------------------------------------------------------------------

comp <- celltype_composition(excit[[]], "group", "subtype",
                             cfg$design$mouse_group_order)
save_table(comp, "06_layer_composition_mouse")

save_figure(
  plot_composition(comp, pal_excit_subtype,
                   levels = cfg$design$mouse_group_order),
  "06_layer_composition_mouse", width = 4, height = 3.6
)

# -----------------------------------------------------------------------------
# 5. Layer-specific DEGs
# -----------------------------------------------------------------------------

step_header("Layer-specific differential expression")

layer_deg <- list()
for (st in unique(as.character(excit$subtype))) {
  cells <- colnames(excit)[as.character(excit$subtype) == st]
  if (length(cells) < 50) next

  sub <- subset(excit, cells = cells)
  for (i in seq_len(nrow(cfg$design$mouse_contrasts))) {
    res <- run_deg_pair(sub,
                        ident1 = cfg$design$mouse_contrasts$ident1[i],
                        ident2 = cfg$design$mouse_contrasts$ident2[i],
                        group_col = "group")
    if (is.null(res)) next
    res$subtype <- st
    layer_deg[[length(layer_deg) + 1]] <- res
  }
}

if (length(layer_deg) > 0) {
  layer_deg <- do.call(rbind, layer_deg)
  save_table(layer_deg, "06_deg_excitatory_layers_mouse", digits = 4)
  save_table(filter_deg(layer_deg), "06_deg_excitatory_layers_significant_mouse",
             digits = 4)
}

# -----------------------------------------------------------------------------
# 6. Cluster markers in a multi-sheet workbook (supplementary table)
# -----------------------------------------------------------------------------

markers <- find_cluster_markers(excit, cluster_col, top_n = 100)

wb <- openxlsx::createWorkbook()
for (cl in sort(unique(markers$cluster))) {
  sheet <- paste0("Cluster_", cl)
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet = sheet,
                      x = markers[markers$cluster == cl, ],
                      rowNames = FALSE)
}
ensure_dir(cfg$dirs$tables)
openxlsx::saveWorkbook(wb, p_tab("06_excitatory_cluster_markers.xlsx"),
                       overwrite = TRUE)
note("saved table: 06_excitatory_cluster_markers.xlsx (",
     length(unique(markers$cluster)), " sheets)")

write_object(excit, "mouse_excitatory")

note("done: excitatory neuron module complete")
