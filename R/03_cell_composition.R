# =============================================================================
# R/03_cell_composition.R
#
# Cell-type composition of both atlases:
#   * stacked proportion bar plots per experimental group,
#   * absolute cell numbers per cell type,
#   * per-sample proportions plus a sample-level Wilcoxon test (this is the
#     conservative test used in the manuscript, because cells within a donor
#     are not independent observations),
#   * a donut chart of cells contributed by each cohort.
#
# Outputs
#   results/tables/03_composition_<dataset>.csv
#   results/tables/03_proportion_test_<dataset>.csv
#   results/figures/03_*.pdf
# =============================================================================

if (!file.exists("R/00_setup.R")) stop("Run this script from the repository root.")
source("R/00_setup.R")

need("Seurat", "dplyr", "ggplot2", "patchwork", "scales")

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------

#' Cell-type proportions computed per sample (donor), not per group.
#'
#' Aggregating per sample first is what makes the downstream test valid: the
#' unit of replication in these cohorts is the donor, not the cell.
proportion_by_sample <- function(meta, group_col = "group",
                                 celltype_col = "celltype",
                                 sample_col = "orig.ident") {
  df <- data.frame(
    sample   = as.character(meta[[sample_col]]),
    group    = as.character(meta[[group_col]]),
    celltype = as.character(meta[[celltype_col]]),
    stringsAsFactors = FALSE
  )
  df <- df[stats::complete.cases(df), , drop = FALSE]

  tab <- as.data.frame(table(sample = df$sample, celltype = df$celltype),
                       stringsAsFactors = FALSE)
  names(tab) <- c("sample", "celltype", "n")
  totals <- stats::aggregate(n ~ sample, data = tab, FUN = sum)
  names(totals)[2] <- "n_total"
  tab <- merge(tab, totals, by = "sample")
  tab$proportion <- tab$n / tab$n_total

  sample_group <- unique(df[, c("sample", "group")])
  merge(tab, sample_group, by = "sample")
}

#' Wilcoxon rank-sum test of per-sample proportions between two groups.
test_proportion_by_group <- function(prop, ident1, ident2,
                                     group_col = "group") {
  celltypes <- sort(unique(prop$celltype))
  out <- do.call(rbind, lapply(celltypes, function(ct) {
    d <- prop[prop$celltype == ct, , drop = FALSE]
    x <- d$proportion[d[[group_col]] == ident1]
    y <- d$proportion[d[[group_col]] == ident2]

    data.frame(
      celltype  = ct,
      n_ident1  = length(x),
      n_ident2  = length(y),
      mean_ident1 = if (length(x)) mean(x) else NA_real_,
      mean_ident2 = if (length(y)) mean(y) else NA_real_,
      p_value   = safe_wilcox(x, y),
      stringsAsFactors = FALSE
    )
  }))
  out$comparison <- paste0(ident1, "_vs_", ident2)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out[out$n_ident1 > 0 & out$n_ident2 > 0, , drop = FALSE]
}

#' Donut chart of the number of cells contributed by each cohort.
plot_cells_donut <- function(meta, group_col, palette, levels = NULL) {
  tab <- as.data.frame(table(group = meta[[group_col]]), stringsAsFactors = FALSE)
  names(tab) <- c("group", "count")
  tab <- tab[tab$count > 0, , drop = FALSE]
  if (!is.null(levels)) tab$group <- ordered_factor(tab$group, levels)

  tab$fraction <- tab$count / sum(tab$count)
  tab$ymax <- cumsum(tab$fraction)
  tab$ymin <- c(0, utils::head(tab$ymax, -1))
  tab$label <- sprintf("%s\n%s cells", tab$group,
                       format(tab$count, big.mark = ","))

  ggplot2::ggplot(tab, ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax,
                                    xmin = 3, xmax = 4, fill = .data$group)) +
    ggplot2::geom_rect(colour = "white", linewidth = 0.6) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(c(2, 4)) +
    ggplot2::scale_fill_manual(values = pal_for(as.character(tab$group),
                                                palette), name = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "right")
}

#' Bar plot of absolute cell numbers per cell type.
plot_celltype_counts <- function(meta, celltype_col, palette, levels = NULL) {
  tab <- as.data.frame(table(celltype = meta[[celltype_col]]),
                       stringsAsFactors = FALSE)
  names(tab) <- c("celltype", "count")
  tab <- tab[tab$count > 0, , drop = FALSE]
  if (!is.null(levels)) tab$celltype <- ordered_factor(tab$celltype, levels)

  ggplot2::ggplot(tab, ggplot2::aes(x = stats::reorder(.data$celltype,
                                                       -.data$count),
                                    y = .data$count, fill = .data$celltype)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = format(.data$count, big.mark = ",")),
                       vjust = -0.4, size = 3) +
    ggplot2::scale_fill_manual(values = pal_for(as.character(tab$celltype),
                                                palette)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12)),
                                labels = scales::comma) +
    ggplot2::labs(x = NULL, y = "Number of cells") +
    theme_pub() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

# -----------------------------------------------------------------------------
# 2. Human atlas
# -----------------------------------------------------------------------------

step_header("Composition: human atlas")

human <- read_seurat(p_processed("human_atlas_annotated.rds"), "annotated human atlas")

comp_human <- celltype_composition(human[[]], "group", "celltype",
                                   cfg$design$human_group_order)
save_table(comp_human, "03_composition_human")

save_figure(
  plot_composition(comp_human, pal_celltype_human,
                   levels = cfg$design$human_group_order),
  "03_composition_human", width = 4.2, height = 4
)

save_figure(
  plot_celltype_counts(human[[]], "celltype", pal_celltype_human,
                       cfg$design$human_celltypes),
  "03_celltype_counts_human", width = 4.6, height = 3.6
)

save_figure(
  plot_cells_donut(human[[]], "group", pal_group_human,
                   cfg$design$human_group_order),
  "03_cells_per_cohort_human", width = 4.6, height = 4
)

# Sample-level test for every disease cohort against controls.
prop_human <- proportion_by_sample(human[[]], "group", "celltype", "orig.ident")
save_table(prop_human, "03_proportion_per_sample_human")

tests_human <- Filter(
  function(d) nrow(d) > 0,
  lapply(cfg$design$human_contrasts$ident1, function(dz) {
    test_proportion_by_group(prop_human, dz, "NC")
  })
)
tests_human <- do.call(rbind, tests_human)
if (!is.null(tests_human)) save_table(tests_human, "03_proportion_test_human", digits = 4)

# -----------------------------------------------------------------------------
# 3. Mouse atlas
# -----------------------------------------------------------------------------

step_header("Composition: mouse atlas")

mouse <- read_seurat(p_processed("mouse_ko_annotated.rds"), "annotated mouse atlas")

comp_mouse <- celltype_composition(mouse[[]], "group", "celltype",
                                   cfg$design$mouse_group_order)
save_table(comp_mouse, "03_composition_mouse")

save_figure(
  plot_composition(comp_mouse, pal_celltype_mouse,
                   levels = cfg$design$mouse_group_order),
  "03_composition_mouse", width = 4.2, height = 4
)

save_figure(
  plot_celltype_counts(mouse[[]], "celltype", pal_celltype_mouse,
                       cfg$design$mouse_celltypes),
  "03_celltype_counts_mouse", width = 4.6, height = 3.6
)

# KOF vs WT and Gq vs vector, at the level of individual samples.
prop_mouse <- proportion_by_sample(mouse[[]], "group", "celltype", "orig.ident")
save_table(prop_mouse, "03_proportion_per_sample_mouse")

tests_mouse <- lapply(seq_len(nrow(cfg$design$mouse_contrasts)), function(i) {
  test_proportion_by_group(prop_mouse,
                           cfg$design$mouse_contrasts$ident1[i],
                           cfg$design$mouse_contrasts$ident2[i])
})
tests_mouse <- do.call(rbind, Filter(function(d) !is.null(d) && nrow(d) > 0,
                                     tests_mouse))
if (!is.null(tests_mouse)) {
  save_table(tests_mouse, "03_proportion_test_mouse", digits = 4)
}

note("done: composition tables and figures written to results/")
