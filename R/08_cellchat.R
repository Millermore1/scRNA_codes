# =============================================================================
# R/08_cellchat.R
#
# Cell-cell communication analysis with CellChat (mouse ligand-receptor
# database), performed on the neuron compartment (PN, PVALB, SST, VIP) of the
# four experimental groups.
#
# The module
#   1. re-clusters all neurons and annotates PN / PVALB / SST / VIP,
#   2. infers one communication network per group,
#   3. merges the networks and compares them (overall strength, differential
#      interactions, information flow, signalling roles),
#   4. exports every inferred interaction as a table and the main comparison
#      figures.
#
# Outputs
#   data/processed/mouse_neurons.rds
#   results/objects/cellchat_*.rds
#   results/tables/08_cellchat_interactions_<group>.csv
#   results/figures/08_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "CellChat", "ggplot2", "patchwork", "dplyr", "scales", "ggrepel")

# -----------------------------------------------------------------------------
# 1. Neuron compartment and subtype annotation
# -----------------------------------------------------------------------------

step_header("Neuron compartment for CellChat")

mouse <- read_seurat(p_processed("mouse_ko_annotated.rds"), "annotated mouse atlas")

neurons <- subset(mouse, cells = colnames(mouse)[
  as.character(mouse$celltype) %in% c("Excitatory neurons", "Inhibitory neurons")
])
note("neurons: ", ncol(neurons), " cells")

neurons <- run_seurat_standard(
  neurons,
  nfeatures  = cfg$params$neuron_subclustering$nfeatures,
  npcs       = cfg$params$neuron_subclustering$npcs,
  dims       = cfg$params$neuron_subclustering$dims,
  resolution = 0.05
)

cluster_col <- "RNA_snn_res.0.05"
save_figure(
  plot_umap(neurons, cluster_col,
            palette = setNames(recycle_colours(30), as.character(0:29)),
            title = "Neuron subclusters"),
  "08_umap_neuron_clusters", width = 4.6, height = 4
)

# Interneuron marker panel used to decide the PN / PVALB / SST / VIP labels.
panel <- unique(unlist(cfg$genes$inhib_subtypes))
panel <- intersect(panel, rownames(neurons))
save_figure(
  plot_marker_dotplot(neurons, panel, cluster_col,
                      title = "Interneuron markers across neuron clusters"),
  "08_dotplot_neuron_clusters", width = 7, height = 3.4
)

neurons <- relabel_clusters(neurons, cfg$design$mouse_neuron_mapping,
                            new_col = "neuron_celltype",
                            cluster_col = cluster_col)
neurons$neuron_celltype <- ordered_factor(neurons$neuron_celltype,
                                          names(pal_neuron_subtype))

save_figure(
  plot_umap(neurons, "neuron_celltype", pal_neuron_subtype,
            title = "Neuronal subtypes"),
  "08_umap_neuron_subtypes", width = 4.4, height = 3.8
)

write_object(neurons, "mouse_neurons")

# -----------------------------------------------------------------------------
# 2. One CellChat network per group
# -----------------------------------------------------------------------------

step_header("CellChat inference")

Seurat::Idents(neurons) <- "neuron_celltype"
sc_list <- Seurat::SplitObject(neurons, split.by = "group")
sc_list <- sc_list[vapply(sc_list, ncol, integer(1)) >= cfg$params$cellchat$min_cells]

db <- switch(cfg$params$cellchat$db,
             mouse  = CellChat::CellChatDB.mouse,
             human  = CellChat::CellChatDB.human,
             stop("cfg$params$cellchat$db must be 'mouse' or 'human'"))

cellchat_list <- list()

for (g in names(sc_list)) {
  step_header(paste0("group: ", g, " (", ncol(sc_list[[g]]), " cells)"))

  data_input <- Seurat::GetAssayData(sc_list[[g]], assay = "RNA", slot = "data")
  meta <- data.frame(neuron_celltype = sc_list[[g]]$neuron_celltype,
                     row.names = colnames(sc_list[[g]]))

  cc <- CellChat::createCellChat(object = data_input, meta = meta,
                                 group.by = "neuron_celltype")
  cc@DB <- db

  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = cfg$params$cellchat$type)
  cc <- CellChat::filterCommunication(cc,
                                      min.cells = cfg$params$cellchat$min_cells)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  cc <- CellChat::netAnalysis_computeCentrality(cc, slot.name = "netP")

  cellchat_list[[g]] <- cc

  interactions <- CellChat::subsetCommunication(cc)
  interactions$group <- g
  save_table(interactions, paste0("08_cellchat_interactions_", g), digits = 5)
}

write_object(cellchat_list, "cellchat_list")

# -----------------------------------------------------------------------------
# 3. Merge and compare
# -----------------------------------------------------------------------------

step_header("Compare communication networks")

ordered_groups <- intersect(cfg$design$mouse_group_order, names(cellchat_list))
cellchat_list <- cellchat_list[ordered_groups]

merged <- CellChat::mergeCellChat(cellchat_list, add.names = names(cellchat_list))
write_object(merged, "cellchat_merged")

# 3a. Overall number and strength of interactions.
save_graphics("08_cellchat_overall_interactions", 8, 3.6, {
  graphics::par(mfrow = c(1, 2), xpd = TRUE)
  CellChat::compareInteractions(merged, show.legend = FALSE,
                                group = seq_along(cellchat_list))
  CellChat::compareInteractions(merged, show.legend = FALSE,
                                group = seq_along(cellchat_list),
                                measure = "weight")
})

# 3b. Circle plots: count and weight, one panel per group, shared scales.
cell_order <- levels(cellchat_list[[1]]@idents)
color_use  <- pal_for(cell_order, pal_neuron_subtype)

count_mats <- lapply(cellchat_list, function(o) o@net$count[cell_order, cell_order])
weight_mats <- lapply(cellchat_list, function(o) o@net$weight[cell_order, cell_order])
vertex_sizes <- lapply(cellchat_list, function(o) {
  as.numeric(table(o@idents)[cell_order])
})

draw_circle_panels <- function(mats, sizes, title_prefix, edge_max) {
  graphics::par(mfrow = c(1, length(mats)), mar = c(0.5, 0.5, 2.2, 0.5),
                oma = c(0, 0, 0, 0), xpd = NA)
  for (i in seq_along(mats)) {
    CellChat::netVisual_circle(
      mats[[i]],
      color.use = color_use,
      vertex.weight = sizes[[i]],
      vertex.weight.max = max(unlist(sizes)),
      vertex.size.max = 16,
      weight.scale = TRUE,
      edge.weight.max = edge_max,
      edge.width.max = 8,
      alpha.edge = 0.55,
      edge.curved = 0.2,
      label.edge = FALSE,
      vertex.label.cex = 0.8,
      margin = 0.25,
      title.name = paste0(title_prefix, "\n", names(mats)[i])
    )
  }
}

save_graphics("08_cellchat_circle_interactions", 4.2 * length(cellchat_list), 4.4, {
  draw_circle_panels(count_mats, vertex_sizes, "Number of interactions",
                     max(unlist(count_mats), na.rm = TRUE))
})

save_graphics("08_cellchat_circle_strength", 4.2 * length(cellchat_list), 4.4, {
  draw_circle_panels(weight_mats, vertex_sizes, "Interaction strength",
                     max(unlist(weight_mats), na.rm = TRUE))
})

# 3c. Differential interaction between KOF and WT.
if (all(c("KOF", "WT") %in% names(cellchat_list))) {
  idx <- match(c("KOF", "WT"), names(cellchat_list))
  save_graphics("08_cellchat_diff_KOF_vs_WT", 5.2, 5.2, {
    CellChat::netVisual_diffInteraction(
      merged, weight.scale = TRUE, measure = "weight",
      comparison = idx,
      title.name = "Differential interaction strength (KOF vs WT)"
    )
  })
}

# 3d. Information flow across groups.
save_graphics("08_cellchat_information_flow", 8, 6, {
  CellChat::rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
})

# 3e. Signalling roles (outgoing vs incoming), one panel per group.
save_graphics("08_cellchat_signaling_roles", 10, 4.6, {
  plots <- lapply(seq_along(cellchat_list), function(i) {
    CellChat::netAnalysis_signalingRole_scatter(cellchat_list[[i]],
                                                title = names(cellchat_list)[i])
  })
  graphics::par(mfrow = c(1, length(plots)))
  for (p in plots) print(p)
})

# 3f. Synapse-related pathways, compared across all groups.
step_header("Synapse-related pathways")

all_pathways <- unique(unlist(lapply(cellchat_list, function(o) o@netP$pathways)))
synapse_pathways <- intersect(cfg$params$cellchat$synapse_pathways, all_pathways)
note("synapse pathways found: ", paste(synapse_pathways, collapse = ", "))

if (length(synapse_pathways) > 0) {
  save_graphics("08_cellchat_synapse_bubble", 10, 7, {
    CellChat::netVisual_bubble(
      merged,
      targets.use = c("PN", "PVALB"),
      signaling   = synapse_pathways,
      comparison  = seq_along(cellchat_list),
      angle.x     = 45,
      title.name  = "Synapse-related signalling to PN and PVALB"
    )
  })

  save_graphics("08_cellchat_synapse_infoflow", 7, 5.5, {
    CellChat::rankNet(merged, signaling = synapse_pathways,
                      mode = "comparison", stacked = TRUE, do.stat = TRUE)
  })

  # Per-pathway circle plots for the strongest pathways.
  for (pw in utils::head(synapse_pathways, 3)) {
    save_graphics(paste0("08_cellchat_pathway_", pw), 6, 5.5, {
      CellChat::netVisual_aggregate(merged, signaling = pw, layout = "circle")
    })
  }
}

# 3g. Outgoing vs incoming signalling strength, one point per cell group and
#     condition. The values are read from CellChat itself instead of being
#     transcribed by hand (which is what the original script did).
extract_signaling_roles <- function(cc, group_name) {
  p <- tryCatch(CellChat::netAnalysis_signalingRole_scatter(cc),
                error = function(e) NULL)
  if (is.null(p) || is.null(p$data)) return(NULL)

  d <- p$data
  numeric_cols <- names(d)[vapply(d, is.numeric, logical(1))]
  if (length(numeric_cols) < 2) return(NULL)
  label_col <- setdiff(names(d), numeric_cols)[1]
  if (is.na(label_col)) return(NULL)

  data.frame(
    cell_group = as.character(d[[label_col]]),
    outgoing   = d[[numeric_cols[1]]],
    incoming   = d[[numeric_cols[2]]],
    group      = group_name,
    stringsAsFactors = FALSE
  )
}

roles <- do.call(rbind, lapply(seq_along(cellchat_list), function(i) {
  extract_signaling_roles(cellchat_list[[i]], names(cellchat_list)[i])
}))

if (!is.null(roles) && nrow(roles) > 0) {
  save_table(roles, "08_cellchat_signaling_roles_table", digits = 4)

  roles$group <- ordered_factor(roles$group, ordered_groups)
  shapes <- rep_len(c(16, 17, 15, 18), length(levels(roles$group)))

  save_figure(
    ggplot2::ggplot(roles, ggplot2::aes(x = .data$outgoing,
                                        y = .data$incoming,
                                        colour = .data$cell_group,
                                        shape = .data$group)) +
      ggplot2::geom_point(size = 3.4, alpha = 0.85) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = paste0(.data$cell_group, " (",
                                    .data$group, ")")),
        size = 2.6, colour = "grey20", max.overlaps = 30,
        seed = cfg$params$seed
      ) +
      ggplot2::scale_colour_manual(values = pal_neuron_subtype,
                                   name = "Cell group") +
      ggplot2::scale_shape_manual(values = shapes,
                                  name = "Condition") +
      ggplot2::labs(x = "Outgoing interaction strength",
                    y = "Incoming interaction strength") +
      theme_pub(),
    "08_cellchat_signaling_role_overlay", width = 6, height = 5
  )
}

# 3h. Interaction table for the supplement (all groups, all pairs).
all_interactions <- do.call(rbind, lapply(names(cellchat_list), function(g) {
  d <- CellChat::subsetCommunication(cellchat_list[[g]])
  if (nrow(d) == 0) return(NULL)
  d$group <- g
  d
}))
if (!is.null(all_interactions)) {
  save_table(all_interactions, "08_cellchat_interactions_all", digits = 5)
}

note("done: CellChat objects and figures written to results/")
