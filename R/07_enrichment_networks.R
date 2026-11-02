# =============================================================================
# R/07_enrichment_networks.R
#
# Functional interpretation of the DEG tables:
#   * GO biological-process over-representation (clusterProfiler, BH, simplify)
#     for up- and down-regulated genes of every mouse contrast,
#   * the "GO bar with genes" figure: one bar per term, term name and member
#     genes printed inside the bar,
#   * the term-similarity network: nodes are GO terms sized by gene hits,
#     node pies show the cell-type composition of the hits, and edges join
#     terms that share genes (Jaccard),
#   * the same for KEGG pathways.
#
# Outputs
#   results/tables/07_enrichment_*.csv
#   results/figures/07_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("dplyr", "ggplot2", "tidyr", "stringr", "igraph", "ggraph",
     "scatterpie", "ggforce", "ggnewscale", "patchwork", "clusterProfiler",
     "org.Mm.eg.db", "ggrepel")

# -----------------------------------------------------------------------------
# 1. Figure helpers
#
# The enrichment runners used below -- run_enrich_go(), run_enrich_kegg(),
# enrich_to_bar_df() -- live in R/utils/analysis.R and the bar layout in
# R/utils/plotting.R, so they can be reused and unit-tested. Only the network
# layout, which is specific to this figure, is defined here.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2. Term-similarity network
# -----------------------------------------------------------------------------

#' Term network: nodes = GO terms, edges = shared genes (Jaccard).
#'
#' @param enrich_df Enrichment table.
#' @param deg DEG table used to colour node pies by cell type.
plot_term_network <- function(enrich_df, deg,
                              celltype_cols = pal_celltype_mouse,
                              top_n = cfg$params$enrichment$top_n_terms,
                              edge_jaccard = cfg$params$enrichment$edge_jaccard,
                              seed = cfg$params$seed,
                              node_r = c(0.035, 0.085)) {
  gdata <- build_enrich_graph_data(
    enrich_df, top_n = top_n,
    edge_jaccard = edge_jaccard,
    min_shared = cfg$params$enrichment$min_shared_genes
  )
  terms <- gdata$terms
  ct_names <- names(celltype_cols)

  # Cell-type composition of the genes behind each term.
  gene_ct <- unique(deg[, c("gene", "celltype")])
  nodes <- do.call(rbind, lapply(seq_len(nrow(terms)), function(i) {
    genes <- gdata$term_genes[[i]]
    tab <- table(gene_ct$celltype[gene_ct$gene %in% genes])
    row <- data.frame(ID = terms$ID[i],
                      Term = terms$Description[i],
                      hits = length(unique(genes)),
                      p.adjust = terms$p.adjust[i],
                      stringsAsFactors = FALSE)
    for (ct in ct_names) {
      row[[ct]] <- if (ct %in% names(tab)) as.integer(tab[[ct]]) else 0L
    }
    row
  }))

  # Layout from an igraph object; isolated nodes are placed on a circle.
  g <- igraph::graph_from_data_frame(
    d = gdata$edges, vertices = nodes$ID, directed = FALSE
  )
  set.seed(seed)
  xy <- if (igraph::ecount(g) > 0) {
    igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
  } else {
    igraph::layout_in_circle(g)
  }

  nodes$x <- norm01(xy[, 1])
  nodes$y <- norm01(xy[, 2])
  nodes$r <- norm01(sqrt(nodes$hits), to = node_r)
  nodes$label <- stringr::str_wrap(nodes$Term, width = 20)

  ed <- NULL
  if (!is.null(gdata$edges) && nrow(gdata$edges) > 0) {
    nodes_xy <- nodes[, c("ID", "x", "y")]
    ed <- merge(gdata$edges, nodes_xy, by.x = "from", by.y = "ID")
    ed <- merge(ed, stats::setNames(nodes_xy, c("to", "xend", "yend")),
                by = "to")
  }

  p <- ggplot2::ggplot()
  if (!is.null(ed)) {
    p <- p + ggplot2::geom_segment(
      data = ed,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend),
      colour = "grey65", linewidth = 0.4, alpha = 0.8
    )
  }

  p +
    scatterpie::geom_scatterpie(
      data = nodes, ggplot2::aes(x = .data$x, y = .data$y, r = .data$r),
      cols = ct_names, colour = NA, alpha = 0.95
    ) +
    ggplot2::scale_fill_manual(values = celltype_cols, name = "Cell type") +
    ggrepel::geom_text_repel(
      data = nodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 3, fontface = "bold", lineheight = 0.92,
      segment.colour = "grey50", segment.size = 0.3,
      min.segment.length = 0, max.overlaps = Inf, seed = seed
    ) +
    ggplot2::coord_fixed(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05),
                         expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#F4F6F8",
                                               colour = "grey35",
                                               linewidth = 0.5),
      legend.position = "bottom"
    )
}

# -----------------------------------------------------------------------------
# 3. Run on the mouse DEG tables
# -----------------------------------------------------------------------------

step_header("Enrichment of mouse DEGs")

deg <- read_table_required(p_tab("04_deg_mouse_all.csv"), "mouse DEG table")
deg_sig <- filter_deg(deg)
if (is.null(deg_sig) || nrow(deg_sig) == 0) {
  stop("No significant mouse DEGs -- run 04_differential_expression.R first.",
       call. = FALSE)
}

# Per contrast and direction, within the inhibitory neurons (the paper's focus).
inh_deg <- deg_sig[deg_sig$celltype == "Inhibitory neurons", , drop = FALSE]

for (cmp in unique(inh_deg$comparison)) {
  for (dir in c("Up", "Down")) {
    genes <- unique(inh_deg$gene[inh_deg$comparison == cmp &
                                   inh_deg$direction == dir])
    if (length(genes) < 5) next

    tag <- paste0(cmp, "_", tolower(dir))

    go <- run_enrich_go(genes, "org.Mm.eg.db", "BP")
    if (!is.null(go)) {
      save_table(go, paste0("07_go_", tag), digits = 4)

      bar <- enrich_to_bar_df(go)
      p_bar <- plot_go_bar_with_genes(
        bar, x_label = "Fold enrichment"
      )
      if (!is.null(p_bar)) {
        save_figure(p_bar, paste0("07_go_bar_genes_", tag),
                    width = 5, height = 3.4)
      }

      net <- plot_term_network(go, inh_deg[inh_deg$comparison == cmp &
                                             inh_deg$direction == dir, ])
      save_figure(net, paste0("07_go_network_", tag),
                  width = 6, height = 5.2)
    }

    kegg <- run_enrich_kegg(genes, "mmu")
    if (!is.null(kegg)) save_table(kegg, paste0("07_kegg_", tag), digits = 4)
  }
}

# Signed overview of both directions for the genetic contrast.
up <- run_enrich_go(unique(inh_deg$gene[inh_deg$comparison == "KOF_vs_WT" &
                                          inh_deg$direction == "Up"]),
                    "org.Mm.eg.db", "BP")
dn <- run_enrich_go(unique(inh_deg$gene[inh_deg$comparison == "KOF_vs_WT" &
                                          inh_deg$direction == "Down"]),
                    "org.Mm.eg.db", "BP")

p_signed <- plot_go_signed_bar(up, dn)
if (!is.null(p_signed)) {
  save_figure(p_signed, "07_go_signed_bar_KOF_vs_WT", width = 6, height = 4)
}

# -----------------------------------------------------------------------------
# 4. Cell-type-resolved enrichment (all seven classes)
# -----------------------------------------------------------------------------

step_header("Cell-type-resolved enrichment")

deg <- read_table_required(p_tab("04_deg_mouse_all.csv"), "mouse DEG table")
deg_sig <- filter_deg(deg)

for (cmp in unique(deg_sig$comparison)) {
  for (dir in c("Up", "Down")) {
    genes <- unique(deg_sig$gene[deg_sig$comparison == cmp &
                                   deg_sig$direction == dir])
    go <- run_enrich_go(genes, "org.Mm.eg.db", "BP")
    if (is.null(go)) next

    tag <- paste0(cmp, "_", tolower(dir), "_allcelltypes")
    save_table(go, paste0("07_go_", tag), digits = 4)

    net <- plot_term_network(go, deg_sig[deg_sig$comparison == cmp &
                                           deg_sig$direction == dir, ],
                             celltype_cols = pal_celltype_mouse_alt)
    save_figure(net, paste0("07_go_network_", tag), width = 6, height = 5.2)
  }
}

note("done: enrichment tables and figures written to results/")
