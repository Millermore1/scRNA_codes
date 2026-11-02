# =============================================================================
# R/05_inhibitory_neurons.R
#
# The inhibitory-neuron analysis, the core of the manuscript:
#   * re-cluster the inhibitory compartment and assign PVALB / SST / VIP
#     subtypes,
#   * subtype composition across the four experimental groups,
#   * PRMT8 expression in each subtype and group, with sample-level tests,
#   * subtype-specific differential expression (KOF vs WT, Gq vs vector),
#   * overlap between the two contrasts,
#   * the autism risk-gene module score (HP_AUTISM).
#
# Outputs
#   data/processed/mouse_inhibitory.rds, human_inhibitory.rds
#   results/tables/05_*.csv
#   results/figures/05_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "dplyr", "ggplot2", "patchwork", "ggVennDiagram", "pheatmap",
     "scales")

# =============================================================================
# Part A -- mouse inhibitory neurons
# =============================================================================

step_header("Mouse inhibitory neurons: subclustering")

mouse <- read_seurat(p_processed("mouse_ko_annotated.rds"), "annotated mouse atlas")

inh <- subset(mouse, cells = colnames(mouse)[
  as.character(mouse$celltype) == "Inhibitory neurons"
])
if (ncol(inh) < 50) {
  stop("Only ", ncol(inh), " inhibitory neurons -- check the annotation.",
       call. = FALSE)
}
note("inhibitory neurons: ", ncol(inh), " cells")

inh <- run_seurat_standard(
  inh,
  nfeatures  = cfg$params$inhibitory$nfeatures,
  npcs       = cfg$params$inhibitory$npcs,
  dims       = cfg$params$inhibitory$dims,
  resolution = cfg$params$inhibitory$resolution
)

save_figure(
  plot_umap(inh, "RNA_snn_res.0.1",
            palette = setNames(recycle_colours(30), as.character(0:29)),
            title = "Inhibitory neuron subclusters"),
  "05_umap_inhibitory_clusters_mouse", width = 4.6, height = 4
)

# Marker inspection before annotation: the subtype gene panels from the paper.
inh_panel <- unique(unlist(cfg$genes$inhib_subtypes[c("VIP", "SST", "PV",
                                                      "TH", "LAMP5", "PAX6")]))
inh_panel <- intersect(inh_panel, rownames(inh))

save_figure(
  plot_marker_dotplot(inh, inh_panel, "RNA_snn_res.0.1",
                      title = "Interneuron subtype markers"),
  "05_dotplot_inhibitory_markers_mouse", width = 7, height = 3.4
)

inh_markers <- find_cluster_markers(inh, "RNA_snn_res.0.1", top_n = 100)
save_table(inh_markers, "05_inhibitory_cluster_markers_mouse", digits = 4)

step_header("Annotate PVALB / SST / VIP")

inh <- relabel_clusters(inh, cfg$design$mouse_inhibitory_mapping,
                        new_col = "subtype", cluster_col = "RNA_snn_res.0.1")
inh$inhibit_celltype <- inh$subtype

save_figure(
  plot_umap(inh, "subtype", pal_inhib_subtype,
            title = "Interneuron subtypes"),
  "05_umap_inhibitory_subtypes_mouse", width = 4.4, height = 3.8
)

save_figure(
  plot_marker_dotplot(inh, cfg$genes$mouse_inhib, "subtype",
                      levels = names(pal_inhib_subtype),
                      title = "Subtype marker expression"),
  "05_dotplot_subtype_markers_mouse", width = 5, height = 2.8
)

# -----------------------------------------------------------------------------
# A1. Subtype composition across groups
# -----------------------------------------------------------------------------

comp_sub <- celltype_composition(inh[[]], "group", "subtype",
                                 cfg$design$mouse_group_order)
save_table(comp_sub, "05_subtype_composition_mouse")

save_figure(
  plot_composition(comp_sub, pal_inhib_subtype,
                   levels = cfg$design$mouse_group_order),
  "05_subtype_composition_mouse", width = 4, height = 3.6
)

# -----------------------------------------------------------------------------
# A2. PRMT8 expression by subtype and group
# -----------------------------------------------------------------------------

prmt8_gene <- "Prmt8"
if (!prmt8_gene %in% rownames(inh)) {
  warning("Prmt8 is not in the mouse object; skipping expression panels.",
          call. = FALSE)
} else {
  expr <- Seurat::FetchData(inh, vars = c(prmt8_gene, "group", "subtype",
                                          "orig.ident"))
  names(expr)[1] <- "expr"

  summary_tab <- do.call(rbind, lapply(split(expr, list(expr$subtype,
                                                        expr$group),
                                             drop = TRUE),
    function(d) {
      data.frame(
        subtype        = d$subtype[1],
        group          = d$group[1],
        n_cells        = nrow(d),
        mean_expr      = mean(d$expr),
        median_expr    = stats::median(d$expr),
        pct_expressing = mean(d$expr > 0) * 100,
        stringsAsFactors = FALSE
      )
    }))
  save_table(summary_tab, "05_prmt8_expression_summary_mouse", digits = 4)

  # Cell-level violin (figure) and sample-level test (statistics).
  expr$group   <- ordered_factor(expr$group, cfg$design$mouse_group_order)
  expr$subtype <- ordered_factor(expr$subtype, names(pal_inhib_subtype))

  p_vln <- ggplot2::ggplot(expr, ggplot2::aes(x = .data$group,
                                              y = .data$expr,
                                              fill = .data$group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    ggplot2::geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.1,
                          linewidth = 0.25) +
    ggplot2::facet_wrap(~ subtype, nrow = 1) +
    ggplot2::scale_fill_manual(values = pal_group_mouse) +
    ggplot2::labs(x = NULL, y = "Prmt8 expression (log-normalised)") +
    theme_pub(base_size = 9) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  save_figure(p_vln, "05_prmt8_violin_mouse", width = 6.5, height = 3)

  # Sample-level Wilcoxon: KOF vs WT and Gq vs vector, per subtype.
  sample_means <- stats::aggregate(
    expr ~ orig.ident + group + subtype, data = expr, FUN = mean
  )
  save_table(sample_means, "05_prmt8_per_sample_mouse", digits = 4)

  test_rows <- list()
  for (st in unique(expr$subtype)) {
    d <- sample_means[sample_means$subtype == st, , drop = FALSE]
    for (i in seq_len(nrow(cfg$design$mouse_contrasts))) {
      c1 <- cfg$design$mouse_contrasts$ident1[i]
      c2 <- cfg$design$mouse_contrasts$ident2[i]
      x <- d$expr[d$group == c1]
      y <- d$expr[d$group == c2]
      test_rows[[length(test_rows) + 1]] <- data.frame(
        subtype = st, comparison = paste0(c1, "_vs_", c2),
        n1 = length(x), n2 = length(y),
        mean1 = if (length(x)) mean(x) else NA_real_,
        mean2 = if (length(y)) mean(y) else NA_real_,
        p_value = safe_wilcox(x, y),
        stringsAsFactors = FALSE
      )
    }
  }
  tests <- do.call(rbind, test_rows)
  tests$p_adj <- stats::p.adjust(tests$p_value, method = "BH")
  save_table(tests, "05_prmt8_test_mouse", digits = 4)
}

# -----------------------------------------------------------------------------
# A3. Subtype-specific differential expression
# -----------------------------------------------------------------------------

step_header("Inhibitory subtype DEGs")

subtype_deg <- list()
for (st in names(pal_inhib_subtype)) {
  cells <- colnames(inh)[as.character(inh$subtype) == st]
  if (length(cells) < 50) next

  sub <- subset(inh, cells = cells)
  for (i in seq_len(nrow(cfg$design$mouse_contrasts))) {
    c1 <- cfg$design$mouse_contrasts$ident1[i]
    c2 <- cfg$design$mouse_contrasts$ident2[i]
    res <- run_deg_pair(sub, ident1 = c1, ident2 = c2, group_col = "group")
    if (is.null(res)) next
    res$subtype <- st
    subtype_deg[[length(subtype_deg) + 1]] <- res
  }
}

if (length(subtype_deg) > 0) {
  subtype_deg <- do.call(rbind, subtype_deg)
  save_table(subtype_deg, "05_deg_inhibitory_subtypes_mouse", digits = 4)

  sig <- filter_deg(subtype_deg)
  save_table(sig, "05_deg_inhibitory_subtypes_significant_mouse", digits = 4)

  # One volcano per subtype for the Gq vs vector contrast.
  for (st in unique(sig$subtype)) {
    d <- subtype_deg[subtype_deg$subtype == st &
                       subtype_deg$comparison == "KOF_Gq_CNO_vs_KOF_vector_CNO", ]
    if (nrow(d) == 0) next
    save_figure(
      plot_volcano(d, label_genes = cfg$genes$asd_label_genes,
                   title = paste0(st, ": Gq vs vector")),
      paste0("05_volcano_", st, "_Gq_vs_vector"), width = 4.6, height = 4
    )
  }
}

# -----------------------------------------------------------------------------
# A4. Overlap between the two contrasts
# -----------------------------------------------------------------------------

if (exists("sig") && nrow(sig) > 0) {
  sets <- split(sig$gene, sig$comparison)
  if (length(sets) == 2) {
    venn <- ggVennDiagram::ggVennDiagram(
      sets,
      label_alpha = 0, label = "both",
      category.names = names(sets),
      edge_size = 0.2, set_size = 3.5, label_size = 5
    ) +
      ggplot2::scale_fill_gradient(low = "#FEF1F3", high = "#F9CCE0",
                                   name = "Gene count") +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right")

    save_figure(venn, "05_venn_inhibitory_deg_mouse", width = 4.4, height = 3.6)

    shared <- Reduce(intersect, sets)
    save_table(data.frame(gene = sort(shared), stringsAsFactors = FALSE),
               "05_shared_deg_inhibitory_mouse")
    note(length(shared), " genes shared by both contrasts")
  }
}

# -----------------------------------------------------------------------------
# A5. Autism risk-gene module score
# -----------------------------------------------------------------------------

step_header("HP_AUTISM module score")

hp_genes_file <- p_meta("hp_autism_genes.txt")
if (file.exists(hp_genes_file)) {
  hp_genes <- readLines(hp_genes_file)
  hp_genes <- trimws(hp_genes[nzchar(trimws(hp_genes))])

  scored <- add_module_score_auto(inh, hp_genes, name = "HP_AUTISM")
  inh <- scored$object
  note("module score built from ", scored$n_genes, " genes")

  score_df <- Seurat::FetchData(inh, vars = c(scored$score_col, "group",
                                              "subtype"))
  names(score_df)[1] <- "score"
  score_df$group <- ordered_factor(score_df$group,
                                   cfg$design$mouse_group_order)

  save_figure(
    ggplot2::ggplot(score_df, ggplot2::aes(x = .data$group, y = .data$score,
                                           fill = .data$group)) +
      ggplot2::geom_violin(scale = "width", trim = FALSE, linewidth = 0.2) +
      ggplot2::geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.1,
                            linewidth = 0.25) +
      ggplot2::facet_wrap(~ subtype, nrow = 1) +
      ggplot2::scale_fill_manual(values = pal_group_mouse) +
      ggplot2::labs(x = NULL, y = "HP_AUTISM score") +
      theme_pub(base_size = 9) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      ),
    "05_hp_autism_score_mouse", width = 6.5, height = 3
  )

  # Heat map of mean scores, subtype x group.
  mean_score <- stats::aggregate(score ~ subtype + group, data = score_df,
                                 FUN = mean)
  mat <- stats::xtabs(score ~ subtype + group, data = mean_score)
  mat_df <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(mat_df)[3] <- "mean_score"
  save_table(mat_df, "05_hp_autism_mean_by_subtype_mouse", digits = 4)

  save_graphics("05_hp_autism_heatmap_mouse", 4.2, 3.2,
                pheatmap::pheatmap(mat, cluster_rows = FALSE,
                                   cluster_cols = FALSE,
                                   fontsize = 9,
                                   main = "HP_AUTISM score"))
}

write_object(inh, "mouse_inhibitory")

# =============================================================================
# Part B -- human inhibitory neurons
# =============================================================================

step_header("Human inhibitory neurons")

human <- read_seurat(p_processed("human_atlas_annotated.rds"),
                     "annotated human atlas")

h_inh <- subset(human, cells = colnames(human)[
  as.character(human$celltype) == "Inhibitory.neurons"
])
note("human inhibitory neurons: ", ncol(h_inh), " cells")

h_inh <- run_seurat_standard(
  h_inh,
  nfeatures  = cfg$params$inhibitory$nfeatures,
  npcs       = cfg$params$inhibitory$npcs,
  dims       = cfg$params$inhibitory$dims,
  resolution = cfg$params$inhibitory$resolution
)

save_figure(
  plot_umap(h_inh, "RNA_snn_res.0.1",
            palette = setNames(recycle_colours(30), as.character(0:29)),
            title = "Human inhibitory neuron subclusters"),
  "05_umap_inhibitory_clusters_human", width = 4.6, height = 4
)

save_figure(
  plot_marker_dotplot(h_inh, cfg$genes$human_inhib, "celltype",
                      title = "Interneuron subtype markers (human)"),
  "05_dotplot_inhibitory_markers_human", width = 5, height = 2.6
)

if ("PRMT8" %in% rownames(h_inh)) {
  h_expr <- Seurat::FetchData(h_inh, vars = c("PRMT8", "group", "orig.ident"))
  names(h_expr)[1] <- "expr"
  h_expr$group <- ordered_factor(h_expr$group, cfg$design$human_group_order)

  save_figure(
    ggplot2::ggplot(h_expr, ggplot2::aes(x = .data$group, y = .data$expr,
                                         fill = .data$group)) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
      ggplot2::geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.1,
                            linewidth = 0.25) +
      ggplot2::scale_fill_manual(values = pal_group_human) +
      ggplot2::labs(x = NULL, y = "PRMT8 expression") +
      theme_pub() +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
    "05_prmt8_violin_human_inhibitory", width = 4, height = 3.4
  )

  h_sample <- stats::aggregate(expr ~ orig.ident + group, data = h_expr,
                               FUN = mean)
  save_table(h_sample, "05_prmt8_per_sample_human", digits = 4)
}

write_object(h_inh, "human_inhibitory")

note("done: inhibitory neuron module complete")
