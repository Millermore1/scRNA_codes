# =============================================================================
# R/09_wgcna.R
#
# Co-expression network analysis of the human inhibitory neurons with hdWGCNA.
# The aim is to place PRMT8 in a co-expression module and test how that module
# relates to disease status.
#
# Steps: metacells -> soft-threshold selection -> signed network -> module
# detection -> module connectivity (kME) -> module-trait correlation -> a
# publication-ready hub network around PRMT8.
#
# Outputs
#   results/objects/wgcna_inhibitory.rds
#   results/tables/09_*.csv
#   results/figures/09_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "hdWGCNA", "WGCNA", "ggplot2", "ggraph", "tidygraph", "igraph")

# -----------------------------------------------------------------------------
# 1. Input: inhibitory neurons of the human atlas
# -----------------------------------------------------------------------------

step_header("hdWGCNA on human inhibitory neurons")

inh <- read_seurat(p_processed("human_inhibitory.rds"),
                   "human inhibitory neurons (run 05_inhibitory_neurons.R first)")

parameter <- cfg$params$wgcna

inh <- hdWGCNA::SetupForWGCNA(
  inh,
  gene_select = parameter$gene_select,
  fraction    = parameter$fraction,
  wgcna_name  = "Inhibitory_Network"
)

# Metacells are built within each donor, which keeps the aggregation honest
# with respect to the batch structure.
inh <- hdWGCNA::MetacellsByGroups(
  seurat_obj  = inh,
  group.by    = c("group", "orig.ident"),
  ident.group = "group",
  k           = parameter$k,
  max_shared  = parameter$max_shared,
  min_cells   = 25
)
inh <- hdWGCNA::NormalizeMetacells(inh)

# The expression matrix for the network is taken from the inhibitory neurons.
inh <- hdWGCNA::SetDatExpr(
  inh,
  group_name = "Inhibitory.neurons",
  group.by   = "celltype",
  assay      = "RNA",
  slot       = "data"
)

# -----------------------------------------------------------------------------
# 2. Soft threshold
# -----------------------------------------------------------------------------

step_header("Soft-threshold selection")

inh <- hdWGCNA::TestSoftPowers(inh, networkType = parameter$network_type)
save_graphics("09_wgcna_soft_powers", 6, 4,
              hdWGCNA::PlotSoftPowers(inh))

note("soft power used for the network: ", parameter$soft_power,
     " (inspect 09_wgcna_soft_powers.pdf to confirm the plateau)")

# -----------------------------------------------------------------------------
# 3. Network and modules
# -----------------------------------------------------------------------------

step_header("Network construction")

inh <- hdWGCNA::ConstructNetwork(
  inh,
  soft_power   = parameter$soft_power,
  setDatExpr   = FALSE,
  tom_name     = "Inhib_TOM",
  overwrite_tom = TRUE
)

save_graphics("09_wgcna_module_dendrogram", 8, 4.5,
              hdWGCNA::PlotDendrogram(
                inh, main = "Inhibitory neuron co-expression modules"
              ))

inh <- hdWGCNA::ModuleEigengenes(inh)
inh <- hdWGCNA::ModuleConnectivity(
  inh,
  group.by   = "celltype",
  group_name = "Inhibitory.neurons"
)

modules <- hdWGCNA::GetModules(inh)
save_table(modules, "09_wgcna_modules")

# -----------------------------------------------------------------------------
# 4. Which module contains PRMT8?
# -----------------------------------------------------------------------------

target_gene <- cfg$genes$prmt8
target_module <- as.character(
  modules$module[modules$gene_name == target_gene]
)

if (length(target_module) == 0) {
  warning(target_gene, " was not assigned to any module -- it may have been ",
          "filtered by the fraction criterion.", call. = FALSE)
} else {
  target_module <- target_module[1]
  note(target_gene, " belongs to module: ", target_module)

  hubs <- hdWGCNA::GetHubGenes(inh, n_hubs = parameter$n_hubs)
  save_table(hubs, "09_wgcna_hub_genes")

  module_hubs <- hubs[hubs$module == target_module, , drop = FALSE]
  module_hubs <- module_hubs[order(-module_hubs$kME), , drop = FALSE]
  save_table(module_hubs,
             paste0("09_wgcna_hub_genes_", gsub("[^A-Za-z0-9]", "", target_module)))

  # 4a. hdWGCNA's own circular hub plot. This function writes its own PDF into
  #     `outdir`, so it is called directly rather than through save_graphics().
  hdWGCNA::ModuleNetworkPlot(
    inh,
    outdir   = cfg$dirs$figures,
    mods     = target_module,
    n_inner  = 10,
    n_outer  = 20,
    n_conns  = 300,
    label_center       = FALSE,
    vertex.label.cex   = 1,
    vertex.label.font  = 3,
    vertex.label.color = "grey10",
    vertex.size        = 2,
    edge.alpha         = 0.3
  )

  # 4b. A cleaner ggraph version of the same neighbourhood, which is easier to
  #     edit in Illustrator.
  top_genes <- utils::head(module_hubs$gene_name, 30)
  if (!target_gene %in% top_genes) {
    top_genes <- unique(c(target_gene, top_genes))[1:30]
  }
  top_genes <- intersect(top_genes, rownames(inh))

  tom <- hdWGCNA::GetTOM(inh)
  sub_tom <- as.matrix(tom[top_genes, top_genes])

  graph <- igraph::graph_from_adjacency_matrix(
    sub_tom, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  tg <- tidygraph::as_tbl_graph(graph) |>
    tidygraph::activate(edges) |>
    dplyr::filter(.data$weight > 0.05) |>
    tidygraph::activate(nodes) |>
    dplyr::mutate(
      degree = tidygraph::centrality_degree(weights = .data$weight),
      label  = ifelse(.data$name == target_gene, target_gene, "Other")
    )

  p_hub <- ggraph::ggraph(tg, layout = "stress") +
    ggraph::geom_edge_link(ggplot2::aes(edge_width = .data$weight),
                           edge_colour = "grey80", edge_alpha = 0.5) +
    ggraph::geom_node_point(ggplot2::aes(size = .data$degree,
                                         colour = .data$label),
                            alpha = 0.9) +
    ggraph::geom_node_text(ggplot2::aes(label = .data$name), repel = TRUE,
                           size = 3, max.overlaps = 50,
                           fontface = "italic") +
    ggplot2::scale_size_continuous(range = c(2, 9), guide = "none") +
    ggplot2::scale_colour_manual(
      values = c("PRMT8" = "#E64B35", "Other" = "#4DBBD5"), guide = "none"
    ) +
    ggplot2::ggtitle(paste0(target_gene, " co-expression neighbourhood")) +
    ggraph::theme_graph(base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                      face = "bold"))

  save_figure(p_hub, paste0("09_wgcna_ggraph_", target_gene),
              width = 6, height = 5.5)
}

# -----------------------------------------------------------------------------
# 5. Module-trait correlation
# -----------------------------------------------------------------------------

step_header("Module-trait correlation")

trait_name <- parameter$trait
trait_group <- sub("_Status$", "", trait_name)

if (trait_group %in% unique(as.character(inh$group))) {
  traits <- data.frame(as.integer(as.character(inh$group) == trait_group))
  names(traits) <- trait_name
  rownames(traits) <- colnames(inh)

  inh <- hdWGCNA::AddMetaData(inh, metadata = traits)
  inh <- hdWGCNA::ModuleTraitCorrelation(
    inh,
    traits   = trait_name,
    group.by = "celltype"
  )

  save_graphics("09_wgcna_module_trait", 5, 4,
                hdWGCNA::PlotModuleTraitCorrelation(inh))

  cor_tab <- tryCatch(hdWGCNA::GetModuleTraitCorrelation(inh),
                      error = function(e) NULL)
  if (!is.null(cor_tab)) {
    save_table(as.data.frame(cor_tab), "09_wgcna_module_trait_correlation",
               digits = 4)
  }
} else {
  note("trait group '", trait_group, "' is absent from the object -- skipping")
}

write_object(inh, "wgcna_inhibitory")

note("done: WGCNA module complete")
